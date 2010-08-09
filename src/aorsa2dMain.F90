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
    use power

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

    integer :: i, nPts_tot
    integer :: nModesR_tot, nModesZ_tot


!   read namelist input data
!   ------------------------

    if (iAm==0) &
    write(*,*) 'Reading namelist'

    call read_nameList ()


!   initialise the parallel env
!   ---------------------------

    if (iAm==0) &
    write(*,*) 'Initialising the parallel environment'

    nPts_tot = sum ( nRAll(1:nGrid)*nZAll(1:nGrid) )

#ifdef par
    call init_procGrid ( nPts_tot )
#endif


!   initialise the spatial grid
!   ---------------------------

    if (iAm==0) &
    write(*,*) 'Creating spatial grid'

    allocate ( allGrids ( nGrid ) )

    do i=1,nGrid 

        allGrids(i) = init_GridBlock ( &
            nRAll(i), nZAll(i), &
            rMinAll(i), rMaxAll(i), &
            zMinAll(i), zMaxAll(i) )

        allGrids(i)%gridNumber = i

    enddo


!   setup magnetic field 
!   --------------------

    if (iAm==0) &
    write(*,*) 'Reading eqdsk'

    if(useEqdsk)then
        call read_geqdsk ( eqdsk, plot = .false. )
        call init_interp ()
    endif

    do i=1,nGrid

        if (useEqdsk) then
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
   
    call init_profiles ()

    do i=1,nGrid


        if (useFluxProfiles) then
            call flux_profiles ( allGrids(i) )
        elseif (useCircularProfiles) then
            call circular_profiles ( allGrids(i) )
        else
            call flat_profiles ( allGrids(i), parabolic = parabolic )
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


!   Calculate aMat index for first pt in each grid block
!   ----------------------------------------------------

    do i=1,nGrid

        if(i==1)then
            allGrids(i)%startRow = 1 
            allGrids(i)%startCol = 1
        else
            allGrids(i)%startRow = sum ( allGrids(1:i-1)%nR * allGrids(1:i-1)%nZ * 3 ) + 1
            allGrids(i)%startCol = sum ( allGrids(1:i-1)%nModesR * allGrids(1:i-1)%nModesZ * 3 ) + 1
        endif

    enddo


!   Antenna current
!   ---------------

    if (iAm==0) &
    write(*,*) 'Building antenna current (brhs)'

    call alloc_total_brhs ( nPts_tot )
    do i=1,nGrid
        call init_brhs ( allGrids(i) )
    enddo


!   Write the run input data to disk
!   --------------------------------

    if (iAm==0) &
    write(*,*) 'Writing run input data to file'

    if (iAm==0) then

        do i=1,nGrid
            write(allGrids(i)%fNumber,'(i3.3)') i
            call write_runData ( allGrids(i) )
        enddo

    endif


!   Allocate the matrix
!   -------------------

    call alloc_total_aMat ( nPts_tot )



!   Label each grid blocks boundary points
!   --------------------------------------

    call labelPts ( allGrids, nPts_tot )


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

    call amat_boundaries ( allGrids, nPts_tot )

    !where(abs(aMat)/=0)
    !        aMat = 11
    !endwhere
    !write(*,*) shape(amat)
    
    !write(*,'(3x,48(1x,i3.2))') (/ (i,i=1,48) /)
    !write(*,*) '   --------------------------------------------------'
    !do i=1,48
    !write(*,'(i2.2,x,48(1x,i3.2),4x,i2.2)') i,int((aMat(i,:))),int((brhs(i)))
    !enddo

    !do i=1,nGrid
    !    write(*,*) allGrids(i)%label
    !enddo

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

    do i=1,nGrid
        call extract_coeffs ( allGrids(i) )    
    enddo


!   Inverse Fourier transform the k coefficents
!   to real space
!   -------------------------------------------

    if (iAm==0) &
    write(*,*) 'Constructin E solution from basis set (alpha,beta,b)'
   
    do i=1,nGrid 
        call sftInv2d ( allGrids(i) )
    enddo


!   Calculate the plasma current
!   ----------------------------

    if (iAm==0) &
    write(*,*) 'Calculating plasma current'

    if ( iAm == 0 ) then

        do i=1,nGrid
            call current ( allGrids(i) )
        enddo

    endif


!   Rotation to Lab frame 
!   ---------------------

    if (iAm==0) &
    write(*,*) 'Rotating E solution to lab frame (R,Th,z)'

    do i=1,nGrid 
        call rotate_E_to_lab ( allGrids(i) )
    enddo



!   Calculate the Joule Heating 
!   ---------------------------

    if (iAm==0) &
    write(*,*) 'Calculating Joule heating' 

    if ( iAm == 0 ) then

        do i=1,nGrid
            call jDotE ( allGrids(i) )
        enddo

    endif


!   Write soln to file
!   ------------------

    if (iAm==0) &
    write(*,*) 'Writing solution to file'

    if ( iAm == 0 ) then

        do i=1,nGrid
            call write_solution ( allGrids(i) )
        enddo

    endif


#ifdef par
    call release_grid ()
#endif

    !stat = gptlStop ('total')
    !stat = gptlpr (0)
    !stat = gptlFinalize ()

    if (iAm==0) &
    write(*,*) 'Total time: ', end_timer ( tTotal )

end program aorsa2dMain


