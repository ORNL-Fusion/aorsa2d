modUle rotation

implicit none

real, allocatable, dimension(:,:) :: btaU
real :: sqx
real, allocatable, dimension(:,:) :: &
    Uxx, Uxy, Uxz, Uyx, Uyy, Uyz, Uzx, Uzy, Uzz

real, allocatable, dimension(:,:) :: &
    Urr_, Urt_, Urz_, Utr_, Utt_, Utz_, Uzr_, Uzt_, Uzz_

real, allocatable, dimension(:,:,:,:) :: U_xyz, U_cyl

real :: sinTh
real, allocatable, dimension(:,:) :: gradPrlB
!real, allocatable, dimension(:,:) :: &
!    dxUxx, dxxUxx, dxUxy, dxxUxy, dxUxz, dxxUxz, &
!    dxUyx, dxxUyx, dxUyy, dxxUyy, dxUyz, dxxUyz, &
!    dxUzx, dxxUzx, dxUzy, dxxUzy, dxUzz, dxxUzz
!real, allocatable, dimension(:,:) :: &
!    dyUxx, dyyUxx, dyUxy, dyyUxy, dyUxz, dyyUxz, &
!    dyUyx, dyyUyx, dyUyy, dyyUyy, dyUyz, dyyUyz, &
!    dyUzx, dyyUzx, dyUzy, dyyUzy, dyUzz, dyyUzz
!real, allocatable, dimension(:,:) :: &
!        dxyUxx, dxyUxy, dxyUxz, &
!        dxyUyx, dxyUyy, dxyUyz, &
!        dxyUzx, dxyUzy, dxyUzz

! dr first derivatives
real, allocatable, dimension(:,:) :: &
    drUrr, drUrt, drUrz, &
    drUtr, drUtt, drUtz, &
    drUzr, drUzt, drUzz

! dz first derivatives
real, allocatable, dimension(:,:) :: &
    dzUrr, dzUrt, dzUrz, &
    dzUtr, dzUtt, dzUtz, &
    dzUzr, dzUzt, dzUzz

! drr second derivatives
real, allocatable, dimension(:,:) :: &
    drrUrr, drrUrt, drrUrz, &
    drrUtr, drrUtt, drrUtz, &
    drrUzr, drrUzt, drrUzz

! dzz second derivatives
real, allocatable, dimension(:,:) :: &
    dzzUrr, dzzUrt, dzzUrz, &
    dzzUtr, dzzUtt, dzzUtz, &
    dzzUzr, dzzUzt, dzzUzz

!drz derivatives
real, allocatable, dimension(:,:) :: &
    drzUrr, drzUrt, drzUrz, &
    drzUtr, drzUtt, drzUtz, &
    drzUzr, drzUzt, drzUzz

contains

    sUbroUtine init_rotation ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY
        use bField
        use grid, &
        only: capR, y

        implicit none

        integer :: i,j 
        real, allocatable :: sqr(:,:)
        real, allocatable :: det(:,:)

        allocate ( sqr ( nPtsx, nPtsY ) )
        allocate ( btaU ( nPtsX, nPtsY ) )

        allocate ( &
            Uxx( nPtsX, nPtsY ), & 
            Uxy( nPtsX, nPtsY ), &
            Uxz( nPtsX, nPtsY ), &
            Uyx( nPtsX, nPtsY ), & 
            Uyy( nPtsX, nPtsY ), &
            Uyz( nPtsX, nPtsY ), &
            Uzx( nPtsX, nPtsY ), & 
            Uzy( nPtsX, nPtsY ), &
            Uzz( nPtsX, nPtsY ) )

        do i = 1, nPtsX
            do j = 1, nPtsY

                btaU(i,j) = sqrt(bxn(i,j)**2 + byn(i,j)**2)

                sqx = sqrt(1.0 - bxn(i,j)**2)

                Uxx(i, j) =   sqx
                Uxy(i, j) = - bxn(i, j) * byn(i, j) / sqx
                Uxz(i, j) = - bxn(i, j) * bzn(i, j) / sqx
                Uyx(i, j) =   0.0
                Uyy(i, j) =   bzn(i, j) / sqx
                Uyz(i, j) = - byn(i, j) / sqx
                Uzx(i, j) =   bxn(i, j)
                Uzy(i, j) =   byn(i, j)
                Uzz(i, j) =   bzn(i, j)

            enddo
        enddo

        
        ! Create a cylindrical coordinate version of the 
        ! rotation matrix so we have some consistency in the
        ! naming conventions of coordinates and hence preserve
        ! what is left of my sanity

        allocate ( &
            Urr_( nPtsX, nPtsY ), & 
            Urt_( nPtsX, nPtsY ), &
            Urz_( nPtsX, nPtsY ), &
            Utr_( nPtsX, nPtsY ), & 
            Utt_( nPtsX, nPtsY ), &
            Utz_( nPtsX, nPtsY ), &
            Uzr_( nPtsX, nPtsY ), & 
            Uzt_( nPtsX, nPtsY ), &
            Uzz_( nPtsX, nPtsY ) )

        sqr     = sqrt ( 1d0 - brn_**2 )

        Urr_     = sqr
        Urt_    = -brn_ * bthn_ / sqr
        Urz_     = -brn_ * bzn_ / sqr

        Utr_    = 0
        Utt_   = bzn_ / sqr
        Utz_    = -bthn_ / sqr

        Uzr_     = brn_
        Uzt_    = bthn_
        Uzz_     = bzn_


        ! Check the determinant, shoUld = 1
        ! ---------------------------------
    
        allocate(det(nPtsX,nPtsY))

        det = Urr_ * Utt_ * Uzz_ &
            + Urt_ * Utz_ * Uzr_ &
            + Urz_ * Utr_ * Uzt_ &
            - Urr_ * Utz_ * Uzt_ &
            - Urt_ * Utr_ * Uzz_ &
            - Urz_ * Utt_ * Uzr_

        if ( any(1-det>1e-4) ) then

            write(*,*) 'rotation.f90: ERROR det != 1'
            stop

        endif

        deallocate(det)

        allocate ( U_xyz(nPtsX,nPtsY,3,3), U_cyl(nPtsX,nPtsY,3,3) )

        U_xyz(:,:,1,1)  = Uxx
        U_xyz(:,:,2,1)  = Uxy
        U_xyz(:,:,3,1)  = Uxz
        
        U_xyz(:,:,1,2)  = Uyx
        U_xyz(:,:,2,2)  = Uyy
        U_xyz(:,:,3,2)  = Uyz

        U_xyz(:,:,1,3)  = Uzx
        U_xyz(:,:,2,3)  = Uzy
        U_xyz(:,:,3,3)  = Uzz

        U_cyl(:,:,1,1)  = Urr_
        U_cyl(:,:,2,1)  = Urt_
        U_cyl(:,:,3,1)  = Urz_
        
        U_cyl(:,:,1,2)  = Utr_
        U_cyl(:,:,2,2)  = Utt_
        U_cyl(:,:,3,2)  = Utz_

        U_cyl(:,:,1,3)  = Uzr_
        U_cyl(:,:,2,3)  = Uzt_
        U_cyl(:,:,3,3)  = Uzz_

    end sUbroUtine init_rotation


    sUbroUtine deriv_rotation

        Use aorsa2din_mod, &
        only: nPtsX, nPtsY, nZFUn, r0
        Use grid
        Use derivatives 
        Use bField
        Use eqdsk_dlg

        implicit none

        integer :: i, j

        allocate ( gradPrlB (nPtsX,nPtsY) )

        !! XYZ version
        !! -----------

        !allocate ( &  
        !    dxUxx(nPtsX,nPtsY), dxxUxx(nPtsX,nPtsY), &
        !    dxUxy(nPtsX,nPtsY), dxxUxy(nPtsX,nPtsY), &
        !    dxUxz(nPtsX,nPtsY), dxxUxz(nPtsX,nPtsY), &
        !    dxUyx(nPtsX,nPtsY), dxxUyx(nPtsX,nPtsY), &
        !    dxUyy(nPtsX,nPtsY), dxxUyy(nPtsX,nPtsY), &
        !    dxUyz(nPtsX,nPtsY), dxxUyz(nPtsX,nPtsY), &
        !    dxUzx(nPtsX,nPtsY), dxxUzx(nPtsX,nPtsY), &
        !    dxUzy(nPtsX,nPtsY), dxxUzy(nPtsX,nPtsY), &
        !    dxUzz(nPtsX,nPtsY), dxxUzz(nPtsX,nPtsY) )

        !allocate ( &  
        !    dyUxx(nPtsX,nPtsY), dyyUxx(nPtsX,nPtsY), &
        !    dyUxy(nPtsX,nPtsY), dyyUxy(nPtsX,nPtsY), &
        !    dyUxz(nPtsX,nPtsY), dyyUxz(nPtsX,nPtsY), &
        !    dyUyx(nPtsX,nPtsY), dyyUyx(nPtsX,nPtsY), &
        !    dyUyy(nPtsX,nPtsY), dyyUyy(nPtsX,nPtsY), &
        !    dyUyz(nPtsX,nPtsY), dyyUyz(nPtsX,nPtsY), &
        !    dyUzx(nPtsX,nPtsY), dyyUzx(nPtsX,nPtsY), &
        !    dyUzy(nPtsX,nPtsY), dyyUzy(nPtsX,nPtsY), &
        !    dyUzz(nPtsX,nPtsY), dyyUzz(nPtsX,nPtsY) )

        !allocate ( &
        !    dxyUxx(nPtsX,nPtsY), dxyUxy(nPtsX,nPtsY), dxyUxz(nPtsX,nPtsY), &
        !    dxyUyx(nPtsX,nPtsY), dxyUyy(nPtsX,nPtsY), dxyUyz(nPtsX,nPtsY), &
        !    dxyUzx(nPtsX,nPtsY), dxyUzy(nPtsX,nPtsY), dxyUzz(nPtsX,nPtsY) )
        !            
        !do i = 1, nPtsX
        !    do j = 1, nPtsY

        !        !   Brambilla approximation:
        !        !   -----------------------

        !        sinTh = y(j) 
        !        if (abs(y(j))>0) &
        !        sinTh =  y(j) / sqrt ( (capR(i)-r0)**2 + y(j)**2 )

        !        gradprlb(i,j) = bMod(i,j) / capr(i) * abs ( btaU(i,j) * sinTh )

        !        !if (nzfUn == 0) gradPrlB(i,j) = 1.0e-10

        !        call deriv_x(Uxx, i, j, dx, dfdx = dxUxx(i,j), d2fdx2 = dxxUxx(i,j))
        !        call deriv_x(Uxy, i, j, dx, dfdx = dxUxy(i,j), d2fdx2 = dxxUxy(i,j))
        !        call deriv_x(Uxz, i, j, dx, dfdx = dxUxz(i,j), d2fdx2 = dxxUxz(i,j))

        !        call deriv_x(Uyx, i, j, dx, dfdx = dxUyx(i,j), d2fdx2 = dxxUyx(i,j))
        !        call deriv_x(Uyy, i, j, dx, dfdx = dxUyy(i,j), d2fdx2 = dxxUyy(i,j))
        !        call deriv_x(Uyz, i, j, dx, dfdx = dxUyz(i,j), d2fdx2 = dxxUyz(i,j))

        !        call deriv_x(Uzx, i, j, dx, dfdx = dxUzx(i,j), d2fdx2 = dxxUzx(i,j))
        !        call deriv_x(Uzy, i, j, dx, dfdx = dxUzy(i,j), d2fdx2 = dxxUzy(i,j))
        !        call deriv_x(Uzz, i, j, dx, dfdx = dxUzz(i,j), d2fdx2 = dxxUzz(i,j))

        !        call deriv_y(Uxx, i, j, dy, dfdy = dyUxx(i,j), d2fdy2 = dyyUxx(i,j))
        !        call deriv_y(Uxy, i, j, dy, dfdy = dyUxy(i,j), d2fdy2 = dyyUxy(i,j))
        !        call deriv_y(Uxz, i, j, dy, dfdy = dyUxz(i,j), d2fdy2 = dyyUxz(i,j))

        !        call deriv_y(Uyx, i, j, dy, dfdy = dyUyx(i,j), d2fdy2 = dyyUyx(i,j))
        !        call deriv_y(Uyy, i, j, dy, dfdy = dyUyy(i,j), d2fdy2 = dyyUyy(i,j))
        !        call deriv_y(Uyz, i, j, dy, dfdy = dyUyz(i,j), d2fdy2 = dyyUyz(i,j))

        !        call deriv_y(Uzx, i, j, dy, dfdy = dyUzx(i,j), d2fdy2 = dyyUzx(i,j))
        !        call deriv_y(Uzy, i, j, dy, dfdy = dyUzy(i,j), d2fdy2 = dyyUzy(i,j))
        !        call deriv_y(Uzz, i, j, dy, dfdy = dyUzz(i,j), d2fdy2 = dyyUzz(i,j))

        !        call deriv_xy(Uxx, i, j, dx, dy, d2fdxy = dxyUxx(i,j))
        !        call deriv_xy(Uxy, i, j, dx, dy, d2fdxy = dxyUxy(i,j))
        !        call deriv_xy(Uxz, i, j, dx, dy, d2fdxy = dxyUxz(i,j))

        !        call deriv_xy(Uyx, i, j, dx, dy, d2fdxy = dxyUyx(i,j))
        !        call deriv_xy(Uyy, i, j, dx, dy, d2fdxy = dxyUyy(i,j))
        !        call deriv_xy(Uyz, i, j, dx, dy, d2fdxy = dxyUyz(i,j))

        !        call deriv_xy(Uzx, i, j, dx, dy, d2fdxy = dxyUzx(i,j))
        !        call deriv_xy(Uzy, i, j, dx, dy, d2fdxy = dxyUzy(i,j))
        !        call deriv_xy(Uzz, i, j, dx, dy, d2fdxy = dxyUzz(i,j))

        !   enddo
        !enddo

        ! R,Th,z version
        ! --------------

        allocate ( &
            drUrr(nPtsX,nPtsY), drrUrr(nPtsX,nPtsY), &
            drUrt(nPtsX,nPtsY), drrUrt(nPtsX,nPtsY), &
            drUrz(nPtsX,nPtsY), drrUrz(nPtsX,nPtsY), &
            drUtr(nPtsX,nPtsY), drrUtr(nPtsX,nPtsY), &
            drUtt(nPtsX,nPtsY), drrUtt(nPtsX,nPtsY), &
            drUtz(nPtsX,nPtsY), drrUtz(nPtsX,nPtsY), &
            drUzr(nPtsX,nPtsY), drrUzr(nPtsX,nPtsY), &
            drUzt(nPtsX,nPtsY), drrUzt(nPtsX,nPtsY), &
            drUzz(nPtsX,nPtsY), drrUzz(nPtsX,nPtsY) )

        allocate ( &
            dzUrr(nPtsX,nPtsY), dzzUrr(nPtsX,nPtsY), &
            dzUrt(nPtsX,nPtsY), dzzUrt(nPtsX,nPtsY), &
            dzUrz(nPtsX,nPtsY), dzzUrz(nPtsX,nPtsY), &
            dzUtr(nPtsX,nPtsY), dzzUtr(nPtsX,nPtsY), &
            dzUtt(nPtsX,nPtsY), dzzUtt(nPtsX,nPtsY), &
            dzUtz(nPtsX,nPtsY), dzzUtz(nPtsX,nPtsY), &
            dzUzr(nPtsX,nPtsY), dzzUzr(nPtsX,nPtsY), &
            dzUzt(nPtsX,nPtsY), dzzUzt(nPtsX,nPtsY), &
            dzUzz(nPtsX,nPtsY), dzzUzz(nPtsX,nPtsY) )

        allocate ( &
            drzUrr(nPtsX,nPtsY), drzUrt(nPtsX,nPtsY), drzUrz(nPtsX,nPtsY), &
            drzUtr(nPtsX,nPtsY), drzUtt(nPtsX,nPtsY), drzUtz(nPtsX,nPtsY), &
            drzUzr(nPtsX,nPtsY), drzUzt(nPtsX,nPtsY), drzUzz(nPtsX,nPtsY) ) 

        do i = 1, nPtsX
            do j = 1, nPtsY

                call deriv_x(capR, Urr_, i, j, dfdx = drUrr(i,j), d2fdx2 = drrUrr(i,j))
                call deriv_x(capR, Urt_, i, j, dfdx = drUrt(i,j), d2fdx2 = drrUrt(i,j))
                call deriv_x(capR, Urz_, i, j, dfdx = drUrz(i,j), d2fdx2 = drrUrz(i,j))

                call deriv_x(capR, Utr_, i, j, dfdx = drUtr(i,j), d2fdx2 = drrUtr(i,j))
                call deriv_x(capR, Utt_, i, j, dfdx = drUtt(i,j), d2fdx2 = drrUtt(i,j))
                call deriv_x(capR, Utz_, i, j, dfdx = drUtz(i,j), d2fdx2 = drrUtz(i,j))

                call deriv_x(capR, Uzr_, i, j, dfdx = drUzr(i,j), d2fdx2 = drrUzr(i,j))
                call deriv_x(capR, Uzt_, i, j, dfdx = drUzt(i,j), d2fdx2 = drrUzt(i,j))
                call deriv_x(capR, Uzz_, i, j, dfdx = drUzz(i,j), d2fdx2 = drrUzz(i,j))

                call deriv_y(y, Urr_, i, j, dfdy = drUrr(i,j), d2fdy2 = drrUrr(i,j))
                call deriv_y(y, Urt_, i, j, dfdy = drUrt(i,j), d2fdy2 = drrUrt(i,j))
                call deriv_y(y, Urz_, i, j, dfdy = drUrz(i,j), d2fdy2 = drrUrz(i,j))

                call deriv_y(y, Utr_, i, j, dfdy = drUtr(i,j), d2fdy2 = drrUtr(i,j))
                call deriv_y(y, Utt_, i, j, dfdy = drUtt(i,j), d2fdy2 = drrUtt(i,j))
                call deriv_y(y, Utz_, i, j, dfdy = drUtz(i,j), d2fdy2 = drrUtz(i,j))

                call deriv_y(y, Uzr_, i, j, dfdy = drUzr(i,j), d2fdy2 = drrUzr(i,j))
                call deriv_y(y, Uzt_, i, j, dfdy = drUzt(i,j), d2fdy2 = drrUzt(i,j))
                call deriv_y(y, Uzz_, i, j, dfdy = drUzz(i,j), d2fdy2 = drrUzz(i,j))

                call deriv_xy(capR, y, Urr_, i, j, d2fdxy = drzUrr(i,j))
                call deriv_xy(capR, y, Urt_, i, j, d2fdxy = drzUrt(i,j))
                call deriv_xy(capR, y, Urz_, i, j, d2fdxy = drzUrz(i,j))

                call deriv_xy(capR, y, Utr_, i, j, d2fdxy = drzUtr(i,j))
                call deriv_xy(capR, y, Utt_, i, j, d2fdxy = drzUtt(i,j))
                call deriv_xy(capR, y, Utz_, i, j, d2fdxy = drzUtz(i,j))

                call deriv_xy(capR, y, Uzr_, i, j, d2fdxy = drzUzr(i,j))
                call deriv_xy(capR, y, Uzt_, i, j, d2fdxy = drzUzt(i,j))
                call deriv_xy(capR, y, Uzz_, i, j, d2fdxy = drzUzz(i,j))

           enddo
        enddo

        !drUrr =0 
        !drUrt =0 
        !drUrz =0 

        !drUtr =0 
        !drUtt =0 
        !drUtz =0 

        !drUzr =0 
        !drUzt =0 
        !drUzz =0 

        !drUrr =0 
        !drUrt =0 
        !drUrz =0 

        !drUtr =0 
        !drUtt =0 
        !drUtz =0 

        !drUzr =0 
        !drUzt =0 
        !drUzz =0 

        !drrUrr=0 
        !drrUrt=0
        !drrUrz=0

        !drrUtr=0
        !drrUtt=0
        !drrUtz=0

        !drrUzr=0
        !drrUzt=0
        !drrUzz=0

        !drrUrr=0
        !drrUrt=0
        !drrUrz=0

        !drrUtr=0
        !drrUtt=0
        !drrUtz=0

        !drrUzr=0
        !drrUzt=0
        !drrUzz=0

        !drzUrr=0 
        !drzUrt=0
        !drzUrz=0

        !drzUtr=0
        !drzUtt=0
        !drzUtz=0

        !drzUzr=0
        !drzUzt=0
        !drzUzz=0

    end subroutine deriv_rotation

end module rotation
