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
    integer :: info, nRow, nCol
    integer, allocatable, dimension(:) :: ipiv
    complex, allocatable, dimension(:,:) :: &
        ealphak, ebetak, eBk
    complex, allocatable, dimension(:,:) :: &
       ealpha, ealphax, ealphay, &
       ebeta, ebetax, ebetay, &
       eB, eBx, eBy

    !   lapack lsq solve variables

    integer :: M_, N_, NRHS, LDA, LDB, MN_, LWORK, RANK
    real :: RCOND
    real, allocatable :: rWork(:)
    complex, allocatable :: work(:)
    integer, allocatable, dimension(:) :: jpvt


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
    !call read_geqdsk ( eqdsk, plot = .false. )
    !call init_interp ()
    !call bFieldEqdsk ()
    call bFieldAnalytical ()

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
    call write_amat ( 'amat.nc' )


!   Antenna current
!   ---------------

    write(*,*) 'Building antenna current (brhs)'
    call init_brhs ()

!   Write the run input data to disk
!   --------------------------------

    write(*,*) 'Writing run input data to file'
    call write_runData ( 'runData.nc', &
        capR, y, bxn, byn, bzn, bmod, xjy, &
        xkxsav, xkysav )


!   Solve complex, linear system
!   ----------------------------

    write(*,*) 'Solving complex linear system'
    
    nRow    = nPtsX * nPtsY * 3
    nCol    = nModesX * nModesY * 3

    !allocate ( ipiv ( nRow ), jpvt ( nCol ) )

    !call cgesv ( nRow, 1, aMat, nRow, ipiv, brhs, nRow, info )

    M_  = nRow
    N_  = nCol
    NRHS   = 1
    LDA = maxVal ( (/ 1, M_ /) )
    LDB = maxVal ( (/ 1, M_, N_ /) )
    MN_ = minVal ( (/ M_, N_ /) )
    LWORK   = MN_ + maxVal ( (/ 2 * MN_, N_ + 1, MN_ + NRHS /) )
    RCOND   = 1.0E-14

    allocate ( WORK ( maxVal ( (/ 1, LWORK /) ) ) )
    allocate ( RWORK ( 2 * N_ ) )
    allocate ( JPVT ( N_ ) ) 

    JPVT    = 0

    call cgelsy ( M_, N_, NRHS, aMat, LDA, brhs, LDB, JPVT, RCOND, RANK, &
                        WORK, LWORK, RWORK, info )

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
    call write_solution ( 'solution.nc', &
        ealpha, ebeta, eB, &
        ealphak, ebetak, eBk )



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


