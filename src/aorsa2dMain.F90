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

    implicit none

    include "f90papi.h"

    !   Timer variables

    type ( timer ) :: tFill, tSolve, tTotal

    !   PAPI variables

    integer :: papi_irc
    real :: papi_rtime, papi_ptime, papi_mflops
    real :: papi_mflops_global
    real :: papi_rtime_zero, papi_ptime_zero
    real :: papi_rtime_fill, papi_rtime_solve
    real :: papi_ptime_fill, papi_ptime_solve
    integer(kind=long) :: papi_flpins

    call start_timer ( tTotal )


!   read namelist input data
!   ------------------------

    if (iAm==0) &
    write(*,*) 'Reading namelist'

    call read_nameList ()


!   initialise the parallel env
!   ---------------------------

    if (iAm==0) &
    write(*,*) 'Initialising the parallel environment'

    call init_procGrid ()


!   initialise the spatial grid
!   ---------------------------

    if (iAm==0) &
    write(*,*) 'Creating spatial grid'

    call init_grid ()


!   read g-eqdsk file
!   -----------------

    if (iAm==0) &
    write(*,*) 'Reading eqdsk'

    if (useEqdsk) then
        call read_geqdsk ( eqdsk, plot = .false. )
        call init_interp ()
        call bFieldEqdsk ()
    else
        call bFieldAnalytical ()
    endif

!   setup profiles
!   --------------

    if (iAm==0) &
    write(*,*) 'Profile setup'
    
    call init_profiles ()

    if (useFluxProfiles) then
        call flux_profiles ()
    else
        call flat_profiles ()
    endif


!   calculate rotation matrix U
!   ---------------------------

    if (iAm==0) &
    write(*,*) 'Building rotation matrix U'

    call init_rotation ()
    call deriv_rotation ()


!   Calculate kx and ky values 
!   --------------------------

    if (iAM==0) &
    write(*,*) 'Creating k grid'

    call init_k ()


!   Precompute basis functions xx(n,i), yy(m,j)
!   -------------------------------------------

    call init_basis_functions () 


!   Fill matrix 
!   -----------

    call start_timer ( tFill )

    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_zero = papi_rTime
    papi_pTime_zero = papi_pTime

    call aMat_fill ()

    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_fill = papi_rTime - papi_rTime_zero
    papi_pTime_fill = papi_pTime - papi_pTime_zero

    if (iAm==0) &
    write(*,*) 'Time to fill aMat: ', end_timer ( tFill )

    !   sum total mflops from all procs

    papi_mflops_global  = papi_mflops
    call sgSum2D ( iContext, 'All', ' ', 1, 1, &
                papi_mflops_global, 1, -1, -1 )

    if (iAm==0) then

        write(*,*) 'PAPI data for fill'
        write(*,'(a30,f12.1)') 'real time: ', papi_rTime_fill
        write(*,'(a30,f12.1)') 'proc time: ', papi_pTime_fill
        write(*,'(a30,f12.2)') 'iAm Gflops/s: ', papi_mflops / 1e3
        write(*,'(a30,f12.2)') 'global Gflops/s: ', papi_mflops_global / 1e3
        write(*,'(a30,i1)') 'status: ', papi_irc

    endif

    !call write_amat ( 'amat.nc' )


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


!   Solve complex, linear system
!   ----------------------------

    if (iAm==0) &
    write(*,*) 'Solving complex linear system'

    call init_parallel_aMat ()

    call start_timer ( tSolve )

    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_zero = papi_rTime
    papi_pTime_zero = papi_pTime

    call solve_lsq_parallel ()

    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_solve = papi_rTime - papi_rTime_zero
    papi_pTime_solve = papi_pTime - papi_pTime_zero

    if (iAm==0) &
    write(*,*) 'Time to solve: ', end_timer ( tSolve )

    !   sum total mflops from all procs

    papi_mflops_global  = papi_mflops
    call sgSum2D ( iContext, 'All', ' ', 1, 1, &
                papi_mflops_global, 1, -1, -1 )

    if (iAm==0) then

        write(*,*) 'PAPI data for solve'
        write(*,'(a30,f12.1)') 'real time: ', papi_rTime_solve
        write(*,'(a30,f12.1)') 'proc time: ', papi_pTime_solve
        write(*,'(a30,f12.2)') 'iAm Gflops/s: ', papi_mflops / 1e3
        write(*,'(a30,f12.2)') 'global Gflops/s: ', papi_mflops_global / 1e3
        write(*,'(a30,i1)') 'status: ', papi_irc

    endif


    !call blacs_barrier ( iContext, 'A' )
    call extract_coeffs ()    


!   Inverse Fourier transform the k coefficents
!   to real space
!   -------------------------------------------

    if (iAm==0) &
    write(*,*) 'Inverse Fourier transforming the k coeffs'
    
    call sftInv2d ( ealphak, f = ealpha )
    call sftInv2d ( ebetak, f = ebeta )
    call sftInv2d ( eBk, f = eB )


!   Write data to file
!   ------------------

    if (iAm==0) &
    write(*,*) 'Writing solution to file'

    if ( iAm == 0 ) &
    call write_solution ( 'solution.nc' )

    call release_grid ()

!
!!     ----------------------------------------------
!!     Calculate E in the Lab frame and eplus, eminus
!!     ----------------------------------------------
!      isq2 = SQRT(0.5)
!      do i = 1, nModesX
!         do j = 1, nModesY
!
!            ex(i,j)   = uxx(i,j) * ealpha(i,j) &
!                      + uyx(i,j) * ebeta(i,j) &
!                      + uzx(i,j) * eb(i,j)
!
!            ey(i,j)   = uxy(i,j) * ealpha(i,j) &
!                      + uyy(i,j) * ebeta(i,j) &
!                      + uzy(i,j) * eb(i,j)
!
!            ez(i,j)   = uxz(i,j) * ealpha(i,j) &
!                      + uyz(i,j) * ebeta(i,j) &
!                      + uzz(i,j) * eb(i,j)
!
!            eplus(i,j)  = isq2 * (ealpha(i,j) + zi * ebeta(i,j))
!          eminus(i,j) = isq2 * (ealpha(i,j) - zi * ebeta(i,j))
!
!         end do
!      end do
!

    if (iAm==0) &
    write(*,*) 'Total time: ', end_timer ( tTotal )

end program aorsa2dMain


