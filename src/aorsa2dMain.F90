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
        
    implicit none

!   Variable list
!   -------------

    integer :: i, j, m, n, s, p, iRow
    integer :: info
    integer, allocatable, dimension(:) :: ipiv
    complex, allocatable, dimension(:,:) :: &
        ealphak, ebetak, eBk
    complex, allocatable, dimension(:,:) :: &
       ealpha, ealphax, ealphay, &
       ebeta, ebetax, ebetay, &
       eB, eBx, eBy


!   read namelist input data
!   ------------------------

    write(*,*) 'Reading namelist'
    call read_nameList ()

!   initialise the spatial grid
!   ---------------------------

    write(*,*) 'Creating spatial grid'
    call init_grid ()

!   read g-eqdsk file
!   -----------------

    write(*,*) 'Reading eqdsk'
    call read_geqdsk ( eqdsk, plot = .false. )
    call init_interp ()
    call bFieldEqdsk ()

!   setup profiles
!   --------------

    write(*,*) 'Profile setup'
    call init_profiles ()
    call flat_profiles ()

!   calculate rotation matrix U
!   ---------------------------

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

    call aMat_fill ()

!   Antenna current
!   ---------------

    write(*,*) 'Building antenna current (brhs)'
    call init_brhs ()

!   Write the run input data to disk
!   --------------------------------

    write(*,*) 'Writing run input data to file'
    call write_runData ( 'runData.nc', &
        capR, y, bxn, byn, bzn, bmod, xjy )


!   Solve complex, linear system
!   ----------------------------

    write(*,*) 'Solving complex linear system'
    
    nRow    = nModesX*nModesY*3
    nCol    = nkx * nky * 3

    allocate ( ipiv ( nRow ) )

    call cgesv ( nRow, 1, aMat, nRow, ipiv, brhs, nRow, info )
            
    write(*,*) '    LAPACK status: ', info


!   Extract the k coefficents from the solution
!   for each field component
!   -------------------------------------------

    allocate ( &
        ealphak(kxL:kxR,kyL:kyR), &
        ebetak(kxL:kxR,kyL:kyR), &
        eBk(kxL:kxR,kyL:kyR) )

    do n=kxL,kxR
        do m=kyL,kyR

            iRow = (m-kyL) * 3 + (n-kxL) * nModesY * 3 + 1
    
            ealphak(n,m)    = brhs(iRow)
            ebetak(n,m)     = brhs(iRow+1)
            eBk(n,m)        = brhs(iRow+2)

        enddo
    enddo


!   Inverse Fourier transform the k coefficents
!   to real space
!   -------------------------------------------

    write(*,*) 'Inverse Fourier transforming the k coeffs'

    call sftinv2d ( xkxsav, xkysav, xx, yy, &
       ealphak, f = ealpha, fx = ealphax, fy = ealphay )
    call sftinv2d ( xkxsav, xkysav, xx, yy, &
       ebetak, f = ebeta, fx = ebetax, fy = ebetay )
    call sftinv2d ( xkxsav, xkysav, xx, yy, &
       eBk, f = eB, fx = eBx, fy = eBy )


!   Write data to file
!   ------------------

    write(*,*) 'Writing solution to file'
    call write_solution ( 'solution.nc', ealpha, ebeta, eB )



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
end program aorsa2dMain


