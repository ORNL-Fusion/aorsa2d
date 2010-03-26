program aorsa2dMain
    
    use constants
    use eqdsk_dlg
    use aorsasubs_mod
    use aorsa2din_mod
    use fourier_mod
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

    type ( timer ) :: tFill, tSolve, tTotal

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

    !call read_geqdsk ( eqdsk, plot = .false. )
    !call init_interp ()
    !call bFieldEqdsk ()
    call bFieldAnalytical ()

!   setup profiles
!   --------------

    if (iAm==0) &
    write(*,*) 'Profile setup'
    
    call init_profiles ()
    call flat_profiles ()

!   calculate rotation matrix U
!   ---------------------------

    if (iAm==0) &
    write(*,*) 'Building rotation matrix U'

    call init_rotation ()
    call deriv_rotation ()

!   Calculate kx and ky values 
!   --------------------------

    call init_k ()

!   precompute basis functions xx(n,i), yy(m,j)
!   -------------------------------------------

    call init_basis_functions () 



!   Load x, y and z equations for spatial point (i,j) and mode number (n,m)
!   ------------------------------------------------------------------------

    call start_timer ( tFill )

    call aMat_fill ()

    if (iAm==0) &
    write(*,*) 'Time to fill aMat: ', end_timer ( tFill )

    call write_amat ( 'amat.nc' )

!   Antenna current
!   ---------------

    if (iAm==0) &
    write(*,*) 'Building antenna current (brhs)'

    call init_brhs ()

!   Write the run input data to disk
!   --------------------------------

    if (iAm==0) &
    write(*,*) 'Writing run input data to file'

    call write_runData ( 'runData.nc' )


!   Solve complex, linear system
!   ----------------------------

    if (iAm==0) &
    write(*,*) 'Solving complex linear system'

    call init_parallel_aMat ()

    call start_timer ( tSolve )

    call solve_lsq_parallel ()

    if (iAm==0) &
    write(*,*) 'Time to solve: ', end_timer ( tSolve )

    !call blacs_barrier ( iContext, 'A' )
    call release_grid ()
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


