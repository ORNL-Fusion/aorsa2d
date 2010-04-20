module rotation

implicit none

real, allocatable, dimension(:,:) :: btau
real :: sqx
real, allocatable, dimension(:,:) :: &
    uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz

real, allocatable, dimension(:,:) :: &
    urr_, urth_, urz_, uthr_, uthth_, uthz_, uzr_, uzth_, uzz_

real :: sinTh
real, allocatable, dimension(:,:) :: gradPrlB
real, allocatable, dimension(:,:) :: &
    dxuxx, dxxuxx, dxuxy, dxxuxy, dxuxz, dxxuxz, &
    dxuyx, dxxuyx, dxuyy, dxxuyy, dxuyz, dxxuyz, &
    dxuzx, dxxuzx, dxuzy, dxxuzy, dxuzz, dxxuzz
real, allocatable, dimension(:,:) :: &
    dyuxx, dyyuxx, dyuxy, dyyuxy, dyuxz, dyyuxz, &
    dyuyx, dyyuyx, dyuyy, dyyuyy, dyuyz, dyyuyz, &
    dyuzx, dyyuzx, dyuzy, dyyuzy, dyuzz, dyyuzz
real, allocatable, dimension(:,:) :: &
        dxyuxx, dxyuxy, dxyuxz, &
        dxyuyx, dxyuyy, dxyuyz, &
        dxyuzx, dxyuzy, dxyuzz

! dr first derivatives
real, allocatable, dimension(:,:) :: &
    drUrr, drUrth, drUrz, &
    drUthr, drUthth, drUthz, &
    drUzr, drUzth, drUzz

! dz first derivatives
real, allocatable, dimension(:,:) :: &
    dzUrr, dzUrth, dzUrz, &
    dzUthr, dzUthth, dzUthz, &
    dzUzr, dzUzth, dzUzz

! drr second derivatives
real, allocatable, dimension(:,:) :: &
    drrUrr, drrUrth, drrUrz, &
    drrUthr, drrUthth, drrUthz, &
    drrUzr, drrUzth, drrUzz

! dzz second derivatives
real, allocatable, dimension(:,:) :: &
    dzzUrr, dzzUrth, dzzUrz, &
    dzzUthr, dzzUthth, dzzUthz, &
    dzzUzr, dzzUzth, dzzUzz

!drz derivatives
real, allocatable, dimension(:,:) :: &
    drzUrr, drzUrth, drzUrz, &
    drzUthr, drzUthth, drzUthz, &
    drzUzr, drzUzth, drzUzz

contains

    subroutine init_rotation ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY
        use bField

        implicit none

        integer :: i,j 
        real, allocatable :: sqr(:,:)

        allocate ( sqr ( nPtsx, nPtsY ) )
        allocate ( btau ( nPtsX, nPtsY ) )

        allocate ( &
            uxx( nPtsX, nPtsY ), & 
            uxy( nPtsX, nPtsY ), &
            uxz( nPtsX, nPtsY ), &
            uyx( nPtsX, nPtsY ), & 
            uyy( nPtsX, nPtsY ), &
            uyz( nPtsX, nPtsY ), &
            uzx( nPtsX, nPtsY ), & 
            uzy( nPtsX, nPtsY ), &
            uzz( nPtsX, nPtsY ) )

        do i = 1, nPtsX
            do j = 1, nPtsY

                btau(i,j) = sqrt(bxn(i,j)**2 + byn(i,j)**2)

                sqx = sqrt(1.0 - bxn(i,j)**2)

                uxx(i, j) =   sqx
                uxy(i, j) = - bxn(i, j) * byn(i, j) / sqx
                uxz(i, j) = - bxn(i, j) * bzn(i, j) / sqx
                uyx(i, j) =   0.0
                uyy(i, j) =   bzn(i, j) / sqx
                uyz(i, j) = - byn(i, j) / sqx
                uzx(i, j) =   bxn(i, j)
                uzy(i, j) =   byn(i, j)
                uzz(i, j) =   bzn(i, j)

            enddo
        enddo

        
        ! Create a cylindrical coordinate version of the 
        ! rotation matrix so we have some consistency in the
        ! naming conventions of coordinates and hence preserve
        ! what is left of my sanity

        allocate ( &
            Urr_( nPtsX, nPtsY ), & 
            Urth_( nPtsX, nPtsY ), &
            Urz_( nPtsX, nPtsY ), &
            Uthr_( nPtsX, nPtsY ), & 
            Uthth_( nPtsX, nPtsY ), &
            Uthz_( nPtsX, nPtsY ), &
            Uzr_( nPtsX, nPtsY ), & 
            Uzth_( nPtsX, nPtsY ), &
            Uzz_( nPtsX, nPtsY ) )

        sqr     = sqrt ( 1d0 - brn_**2 )

        Urr_     = sqr
        Urth_    = -brn_ * bthn_ / sqr
        Urz_     = -brn_ * bzn_ / sqr

        Uthr_    = 0
        Uthth_   = bzn_ / sqr
        Uthz_    = -bthn_ / sqr

        Uzr_     = brn_
        Uzth_    = bthn_
        Uzz_     = bzn_

    end subroutine init_rotation


    subroutine deriv_rotation

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nZFun
        use grid
        use derivatives 
        use bField
        use eqdsk_dlg

        implicit none

        integer :: i, j

        allocate ( gradPrlB (nPtsX,nPtsY) )

        ! XYZ version
        ! -----------

        allocate ( &  
            dxuxx(nPtsX,nPtsY), dxxuxx(nPtsX,nPtsY), &
            dxuxy(nPtsX,nPtsY), dxxuxy(nPtsX,nPtsY), &
            dxuxz(nPtsX,nPtsY), dxxuxz(nPtsX,nPtsY), &
            dxuyx(nPtsX,nPtsY), dxxuyx(nPtsX,nPtsY), &
            dxuyy(nPtsX,nPtsY), dxxuyy(nPtsX,nPtsY), &
            dxuyz(nPtsX,nPtsY), dxxuyz(nPtsX,nPtsY), &
            dxuzx(nPtsX,nPtsY), dxxuzx(nPtsX,nPtsY), &
            dxuzy(nPtsX,nPtsY), dxxuzy(nPtsX,nPtsY), &
            dxuzz(nPtsX,nPtsY), dxxuzz(nPtsX,nPtsY) )

        allocate ( &  
            dyuxx(nPtsX,nPtsY), dyyuxx(nPtsX,nPtsY), &
            dyuxy(nPtsX,nPtsY), dyyuxy(nPtsX,nPtsY), &
            dyuxz(nPtsX,nPtsY), dyyuxz(nPtsX,nPtsY), &
            dyuyx(nPtsX,nPtsY), dyyuyx(nPtsX,nPtsY), &
            dyuyy(nPtsX,nPtsY), dyyuyy(nPtsX,nPtsY), &
            dyuyz(nPtsX,nPtsY), dyyuyz(nPtsX,nPtsY), &
            dyuzx(nPtsX,nPtsY), dyyuzx(nPtsX,nPtsY), &
            dyuzy(nPtsX,nPtsY), dyyuzy(nPtsX,nPtsY), &
            dyuzz(nPtsX,nPtsY), dyyuzz(nPtsX,nPtsY) )

        allocate ( &
            dxyuxx(nPtsX,nPtsY), dxyuxy(nPtsX,nPtsY), dxyuxz(nPtsX,nPtsY), &
            dxyuyx(nPtsX,nPtsY), dxyuyy(nPtsX,nPtsY), dxyuyz(nPtsX,nPtsY), &
            dxyuzx(nPtsX,nPtsY), dxyuzy(nPtsX,nPtsY), dxyuzz(nPtsX,nPtsY) )
                    
        do i = 1, nPtsX
            do j = 1, nPtsY

                !   Brambilla approximation:
                !   -----------------------

                sinTh = y(j) / sqrt ( (capR(i)-r0__)**2 + (y(j)-z0__)**2 )
                gradprlb(i,j) = bMod(i,j) / capr(i) * abs ( btau(i,j) * sinTh )

                if (nzfun == 0) gradPrlB(i,j) = 1.0e-10

                call deriv_x(uxx, i, j, dx, dfdx = dxuxx(i,j), d2fdx2 = dxxuxx(i,j))
                call deriv_x(uxy, i, j, dx, dfdx = dxuxy(i,j), d2fdx2 = dxxuxy(i,j))
                call deriv_x(uxz, i, j, dx, dfdx = dxuxz(i,j), d2fdx2 = dxxuxz(i,j))

                call deriv_x(uyx, i, j, dx, dfdx = dxuyx(i,j), d2fdx2 = dxxuyx(i,j))
                call deriv_x(uyy, i, j, dx, dfdx = dxuyy(i,j), d2fdx2 = dxxuyy(i,j))
                call deriv_x(uyz, i, j, dx, dfdx = dxuyz(i,j), d2fdx2 = dxxuyz(i,j))

                call deriv_x(uzx, i, j, dx, dfdx = dxuzx(i,j), d2fdx2 = dxxuzx(i,j))
                call deriv_x(uzy, i, j, dx, dfdx = dxuzy(i,j), d2fdx2 = dxxuzy(i,j))
                call deriv_x(uzz, i, j, dx, dfdx = dxuzz(i,j), d2fdx2 = dxxuzz(i,j))

                call deriv_y(uxx, i, j, dy, dfdy = dyuxx(i,j), d2fdy2 = dyyuxx(i,j))
                call deriv_y(uxy, i, j, dy, dfdy = dyuxy(i,j), d2fdy2 = dyyuxy(i,j))
                call deriv_y(uxz, i, j, dy, dfdy = dyuxz(i,j), d2fdy2 = dyyuxz(i,j))

                call deriv_y(uyx, i, j, dy, dfdy = dyuyx(i,j), d2fdy2 = dyyuyx(i,j))
                call deriv_y(uyy, i, j, dy, dfdy = dyuyy(i,j), d2fdy2 = dyyuyy(i,j))
                call deriv_y(uyz, i, j, dy, dfdy = dyuyz(i,j), d2fdy2 = dyyuyz(i,j))

                call deriv_y(uzx, i, j, dy, dfdy = dyuzx(i,j), d2fdy2 = dyyuzx(i,j))
                call deriv_y(uzy, i, j, dy, dfdy = dyuzy(i,j), d2fdy2 = dyyuzy(i,j))
                call deriv_y(uzz, i, j, dy, dfdy = dyuzz(i,j), d2fdy2 = dyyuzz(i,j))

                call deriv_xy(uxx, i, j, dx, dy, d2fdxy = dxyuxx(i,j))
                call deriv_xy(uxy, i, j, dx, dy, d2fdxy = dxyuxy(i,j))
                call deriv_xy(uxz, i, j, dx, dy, d2fdxy = dxyuxz(i,j))

                call deriv_xy(uyx, i, j, dx, dy, d2fdxy = dxyuyx(i,j))
                call deriv_xy(uyy, i, j, dx, dy, d2fdxy = dxyuyy(i,j))
                call deriv_xy(uyz, i, j, dx, dy, d2fdxy = dxyuyz(i,j))

                call deriv_xy(uzx, i, j, dx, dy, d2fdxy = dxyuzx(i,j))
                call deriv_xy(uzy, i, j, dx, dy, d2fdxy = dxyuzy(i,j))
                call deriv_xy(uzz, i, j, dx, dy, d2fdxy = dxyuzz(i,j))

           enddo
        enddo

        ! R,Th,z version
        ! --------------

        allocate ( &
            drUrr(nPtsX,nPtsY), drrUrr(nPtsX,nPtsY), &
            drUrth(nPtsX,nPtsY), drrUrth(nPtsX,nPtsY), &
            drUrz(nPtsX,nPtsY), drrUrz(nPtsX,nPtsY), &
            drUthr(nPtsX,nPtsY), drrUthr(nPtsX,nPtsY), &
            drUthth(nPtsX,nPtsY), drrUthth(nPtsX,nPtsY), &
            drUthz(nPtsX,nPtsY), drrUthz(nPtsX,nPtsY), &
            drUzr(nPtsX,nPtsY), drrUzr(nPtsX,nPtsY), &
            drUzth(nPtsX,nPtsY), drrUzth(nPtsX,nPtsY), &
            drUzz(nPtsX,nPtsY), drrUzz(nPtsX,nPtsY) )

        allocate ( &
            dzUrr(nPtsX,nPtsY), dzzUrr(nPtsX,nPtsY), &
            dzUrth(nPtsX,nPtsY), dzzUrth(nPtsX,nPtsY), &
            dzUrz(nPtsX,nPtsY), dzzUrz(nPtsX,nPtsY), &
            dzUthr(nPtsX,nPtsY), dzzUthr(nPtsX,nPtsY), &
            dzUthth(nPtsX,nPtsY), dzzUthth(nPtsX,nPtsY), &
            dzUthz(nPtsX,nPtsY), dzzUthz(nPtsX,nPtsY), &
            dzUzr(nPtsX,nPtsY), dzzUzr(nPtsX,nPtsY), &
            dzUzth(nPtsX,nPtsY), dzzUzth(nPtsX,nPtsY), &
            dzUzz(nPtsX,nPtsY), dzzUzz(nPtsX,nPtsY) )

        allocate ( &
            drzUrr(nPtsX,nPtsY), drzUrth(nPtsX,nPtsY), drzUrz(nPtsX,nPtsY), &
            drzUthr(nPtsX,nPtsY), drzUthth(nPtsX,nPtsY), drzUthz(nPtsX,nPtsY), &
            drzUzr(nPtsX,nPtsY), drzUzth(nPtsX,nPtsY), drzUzz(nPtsX,nPtsY) ) 

        do i = 1, nPtsX
            do j = 1, nPtsY

                call deriv_x(Urr_, i, j, dx, dfdx = drUrr(i,j), d2fdx2 = drrUrr(i,j))
                call deriv_x(Urth_, i, j, dx, dfdx = drUrth(i,j), d2fdx2 = drrUrth(i,j))
                call deriv_x(Urz_, i, j, dx, dfdx = drUrz(i,j), d2fdx2 = drrUrz(i,j))

                call deriv_x(Uthr_, i, j, dx, dfdx = drUthr(i,j), d2fdx2 = drrUthr(i,j))
                call deriv_x(Uthth_, i, j, dx, dfdx = drUthth(i,j), d2fdx2 = drrUthth(i,j))
                call deriv_x(Uthz_, i, j, dx, dfdx = drUthz(i,j), d2fdx2 = drrUthz(i,j))

                call deriv_x(Uzr_, i, j, dx, dfdx = drUzr(i,j), d2fdx2 = drrUzr(i,j))
                call deriv_x(Uzth_, i, j, dx, dfdx = drUzth(i,j), d2fdx2 = drrUzth(i,j))
                call deriv_x(Uzz_, i, j, dx, dfdx = drUzz(i,j), d2fdx2 = drrUzz(i,j))

                call deriv_y(Urr_, i, j, dy, dfdy = drUrr(i,j), d2fdy2 = drrUrr(i,j))
                call deriv_y(Urth_, i, j, dy, dfdy = drUrth(i,j), d2fdy2 = drrUrth(i,j))
                call deriv_y(Urz_, i, j, dy, dfdy = drUrz(i,j), d2fdy2 = drrUrz(i,j))

                call deriv_y(Uthr_, i, j, dy, dfdy = drUthr(i,j), d2fdy2 = drrUthr(i,j))
                call deriv_y(Uthth_, i, j, dy, dfdy = drUthth(i,j), d2fdy2 = drrUthth(i,j))
                call deriv_y(Uthz_, i, j, dy, dfdy = drUthz(i,j), d2fdy2 = drrUthz(i,j))

                call deriv_y(Uzr_, i, j, dy, dfdy = drUzr(i,j), d2fdy2 = drrUzr(i,j))
                call deriv_y(Uzth_, i, j, dy, dfdy = drUzth(i,j), d2fdy2 = drrUzth(i,j))
                call deriv_y(Uzz_, i, j, dy, dfdy = drUzz(i,j), d2fdy2 = drrUzz(i,j))

                call deriv_xy(Urr_, i, j, dx, dy, d2fdxy = drzUrr(i,j))
                call deriv_xy(Urth_, i, j, dx, dy, d2fdxy = drzUrth(i,j))
                call deriv_xy(Urz_, i, j, dx, dy, d2fdxy = drzUrz(i,j))

                call deriv_xy(Uthr_, i, j, dx, dy, d2fdxy = drzUthr(i,j))
                call deriv_xy(Uthth_, i, j, dx, dy, d2fdxy = drzUthth(i,j))
                call deriv_xy(Uthz_, i, j, dx, dy, d2fdxy = drzUthz(i,j))

                call deriv_xy(Uzr_, i, j, dx, dy, d2fdxy = drzUzr(i,j))
                call deriv_xy(Uzth_, i, j, dx, dy, d2fdxy = drzUzth(i,j))
                call deriv_xy(Uzz_, i, j, dx, dy, d2fdxy = drzUzz(i,j))

           enddo
        enddo


    end subroutine deriv_rotation

end module rotation
