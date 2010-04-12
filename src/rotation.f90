module rotation

implicit none

real, allocatable, dimension(:,:) :: btau
real :: sqx
real, allocatable, dimension(:,:) :: &
    uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz

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

contains

    subroutine init_rotation ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY
        use bField

        implicit none

        integer :: i,j 

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

    end subroutine deriv_rotation

end module rotation
