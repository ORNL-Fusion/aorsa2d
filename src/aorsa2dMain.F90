program aorsa2dMain
   
    use constants
    use eqdsk_dlg
    use aorsa2din_mod
    use inv_fourier
    use write_data
    use grid
    use bField
    use interp
    use profiles
    use rotation
    use mat_fill
    use antenna
    use solve
    use parallel
    use timer_mod
    use E_to_lab

    implicit none

#ifdef usepapi
    include "f90papi.h"
#endif
    !include "gptl.inc"

    !   Timer variables

    type ( timer ) :: tFill, tSolve, tTotal

#ifdef usepapi

    !   PAPI variables

    integer :: papi_irc
    real :: papi_rtime, papi_ptime, papi_mflops
    real :: papi_mflops_global
    real :: papi_rtime_zero, papi_ptime_zero
    real :: papi_rtime_fill, papi_rtime_solve
    real :: papi_ptime_fill, papi_ptime_solve
    integer(kind=long) :: papi_flpins

#endif

    !! GPTL vars

    !integer :: stat

    !stat = gptlSetUTR ( GPTLpapitime)
    !write(*,*) 'gptl papiTime stat: ', stat
    !stat = gptlInitialize ()
    !stat = gptlStart ( 'total' )

    !call start_timer ( tTotal )

    integer :: i


!   read namelist input data
!   ------------------------

    if (iAm==0) &
    write(*,*) 'Reading namelist'

    call read_nameList ()


!   initialise the parallel env
!   ---------------------------

    if (iAm==0) &
    write(*,*) 'Initialising the parallel environment'

#ifdef par
    call init_procGrid ()
#endif


!   initialise the spatial grid
!   ---------------------------

    if (iAm==0) &
    write(*,*) 'Creating spatial grid'

    allocate ( allGrids ( nGrid ) )
    nRAll(1) = nPtsX
    nZAll(1) = nPtsY
    rMinAll(1) = rwleft
    rMaxAll(1) = rwright
    zMinAll(1) = yBot
    zMaxAll(1) = yTop

    do i=1,nGrid 

        allGrids(i) = init_GridBlock ( &
            nRAll(i), nZAll(i), rMinAll(i), rMaxAll(i), zMinAll(i), zMaxAll(i) )

    enddo

    !call init_grid ()


!   setup magnetic field 
!   --------------------

    if (iAm==0) &
    write(*,*) 'Reading eqdsk'

    do i=1,nGrid

        if (useEqdsk) then
            call read_geqdsk ( eqdsk, plot = .false. )
            call init_interp ()
            call bFieldEqdsk ( allGrids(i) )
        elseif (useSoloviev) then
            call soloviev ( allGrids(i) )
        elseif (useCircular) then
            call bFieldCircular ( allGrids(i) )
        else
            call bFieldAnalytical ( allGrids(i) )
        endif

    enddo

!   setup profiles
!   --------------

    if (iAm==0) &
    write(*,*) 'Profile setup'
   
    do i=1,nGrid

        call init_profiles ( allGrids(i) )

        if (useFluxProfiles) then
            call flux_profiles ( allGrids(i) )
        elseif (useCircularProfiles) then
            call circular_profiles ( allGrids(i) )
        else
            call flat_profiles ( allGrids(i) )
        endif

    enddo

!   calculate rotation matrix U
!   ---------------------------

    if (iAm==0) &
    write(*,*) 'Building rotation matrix U'

    do i=1,nGrid

        call init_rotation ( allGrids(i) )
        call deriv_rotation ( allGrids(i) )

    enddo


!!   Calculate kx and ky values 
!!   --------------------------
!
!    if (iAM==0) &
!    write(*,*) 'Creating k grid'
!
!    call init_k ()
!
!
!!   Precompute basis functions xx(n,i), yy(m,j)
!!   -------------------------------------------
!
!    call init_basis_functions () 


!   Antenna current
!   ---------------

    if (iAm==0) &
    write(*,*) 'Building antenna current (brhs)'

    call init_brhs ()


!   Write the run input data to disk
!   --------------------------------

    if (iAm==0) &
    write(*,*) 'Writing run input data to file'

    if (iAm==0) &
    call write_runData ( 'runData.nc' )


!   Allocate the matrix
!   -------------------

#ifdef par

    if (iAm == 0) &
    write(*,100), &
        nPtsX*nPtsY*3*nModesX*nModesY*3*2*8.0 / 1024.0**2, &
        nRowLocal*nColLocal*2*8.0 / 1024.0**2
    100 format (' Filling aMat [global size: ',f8.1,' MB, local size: ',f8.1' MB]')

    allocate ( aMat(nRowLocal,nColLocal) )
#else 
    write(*,100), &
        nPtsX*nPtsY*3*nModesX*nModesY*3*2*8.0 / 1024.0**2
    100 format (' Filling aMat [global size: ',f8.1,' MB]')

    allocate ( aMat(nPtsX*nPtsY*3,nModesX*nModesY*3) )

#endif

    aMat = 0

!   Fill matrix 
!   -----------

    call start_timer ( tFill )

#ifdef usepapi
    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_zero = papi_rTime
    papi_pTime_zero = papi_pTime
#endif

    do i=1,nGrid

        call aMat_fill ( allGrids(i) )

    enddo

#ifdef usepapi
    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_fill = papi_rTime - papi_rTime_zero
    papi_pTime_fill = papi_pTime - papi_pTime_zero
#endif

    if (iAm==0) &
    write(*,*) 'Time to fill aMat: ', end_timer ( tFill )

#ifdef usepapi
    !   sum total mflops from all procs

    papi_mflops_global  = papi_mflops

#ifdef par
    call sgSum2D ( iContext, 'All', ' ', 1, 1, &
                papi_mflops_global, 1, -1, -1 )
#endif

    if (iAm==0) then

        write(*,*) 'PAPI data for fill'
        write(*,'(a30,f12.1)') 'real time: ', papi_rTime_fill
        write(*,'(a30,f12.1)') 'proc time: ', papi_pTime_fill
        write(*,'(a30,f12.2)') 'iAm Gflops/s: ', papi_mflops / 1e3
        write(*,'(a30,f12.2)') 'global Gflops/s: ', papi_mflops_global / 1e3
        write(*,'(a30,i1)') 'status: ', papi_irc

    endif

#endif

    !if (iAm==0) &
    !call write_amat ( 'amat.nc' )


!   Solve complex, linear system
!   ----------------------------

    if (iAm==0) &
    write(*,*) 'Solving complex linear system'

#ifdef par
    call init_parallel_aMat ()
#endif

    call start_timer ( tSolve )

#ifdef usepapi
    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_zero = papi_rTime
    papi_pTime_zero = papi_pTime
#endif

#ifdef par
    call solve_lsq_parallel ()
#else
    call solve_lsq ()
#endif

#ifdef usepapi
    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_solve = papi_rTime - papi_rTime_zero
    papi_pTime_solve = papi_pTime - papi_pTime_zero
#endif

    if (iAm==0) &
    write(*,*) 'Time to solve: ', end_timer ( tSolve )

#ifdef usepapi
    !   sum total mflops from all procs

    papi_mflops_global  = papi_mflops

#ifdef par
    call sgSum2D ( iContext, 'All', ' ', 1, 1, &
                papi_mflops_global, 1, -1, -1 )
#endif

    if (iAm==0) then

        write(*,*) 'PAPI data for solve'
        write(*,'(a30,f12.1)') 'real time: ', papi_rTime_solve
        write(*,'(a30,f12.1)') 'proc time: ', papi_pTime_solve
        write(*,'(a30,f12.2)') 'iAm Gflops/s: ', papi_mflops / 1e3
        write(*,'(a30,f12.2)') 'global Gflops/s: ', papi_mflops_global / 1e3
        write(*,'(a30,i1)') 'status: ', papi_irc

    endif
#endif

    call extract_coeffs ()    


!   Inverse Fourier transform the k coefficents
!   to real space
!   -------------------------------------------

    if (iAm==0) &
    write(*,*) 'Constructin E solution from basis set (alpha,beta,b)'
    
    call sftInv2d ( ealphak, f = ealpha )
    call sftInv2d ( ebetak, f = ebeta )
    call sftInv2d ( eBk, f = eB )


!   Rotation E solution to Lab frame 
!   --------------------------------

    if (iAm==0) &
    write(*,*) 'Rotating E solution to lab frame (R,Th,z)'
 
    call rotate_E_to_lab ()


!   Write data to file
!   ------------------

    if (iAm==0) &
    write(*,*) 'Writing solution to file'

    if ( iAm == 0 ) &
    call write_solution ( 'solution.nc' )

#ifdef par
    call release_grid ()
#endif

    !stat = gptlStop ('total')
    !stat = gptlpr (0)
    !stat = gptlFinalize ()

    if (iAm==0) &
    write(*,*) 'Total time: ', end_timer ( tTotal )

end program aorsa2dMain


