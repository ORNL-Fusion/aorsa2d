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
        only: nModesX, nModesY
        use bField

        implicit none

        integer :: i,j 

        allocate ( btau ( nModesX, nModesY ) )

        allocate ( &
            uxx( nModesX, nModesY ), & 
            uxy( nModesX, nModesY ), &
            uxz( nModesX, nModesY ), &
            uyx( nModesX, nModesY ), & 
            uyy( nModesX, nModesY ), &
            uyz( nModesX, nModesY ), &
            uzx( nModesX, nModesY ), & 
            uzy( nModesX, nModesY ), &
            uzz( nModesX, nModesY ) )

        do i = 1, nModesX
            do j = 1, nModesY

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
        only: nModesX, nModesY, nZFun
        use grid
        use aorsasubs_mod, &
        only: deriv_x, deriv_y, deriv_xy
        use bField
        use eqdsk_dlg

        implicit none

        integer :: i, j

        allocate ( gradPrlB (nModesX,nModesY) )

        allocate ( &  
            dxuxx(nModesX,nModesY), dxxuxx(nModesX,nModesY), &
            dxuxy(nModesX,nModesY), dxxuxy(nModesX,nModesY), &
            dxuxz(nModesX,nModesY), dxxuxz(nModesX,nModesY), &
            dxuyx(nModesX,nModesY), dxxuyx(nModesX,nModesY), &
            dxuyy(nModesX,nModesY), dxxuyy(nModesX,nModesY), &
            dxuyz(nModesX,nModesY), dxxuyz(nModesX,nModesY), &
            dxuzx(nModesX,nModesY), dxxuzx(nModesX,nModesY), &
            dxuzy(nModesX,nModesY), dxxuzy(nModesX,nModesY), &
            dxuzz(nModesX,nModesY), dxxuzz(nModesX,nModesY) )

        allocate ( &  
            dyuxx(nModesX,nModesY), dyyuxx(nModesX,nModesY), &
            dyuxy(nModesX,nModesY), dyyuxy(nModesX,nModesY), &
            dyuxz(nModesX,nModesY), dyyuxz(nModesX,nModesY), &
            dyuyx(nModesX,nModesY), dyyuyx(nModesX,nModesY), &
            dyuyy(nModesX,nModesY), dyyuyy(nModesX,nModesY), &
            dyuyz(nModesX,nModesY), dyyuyz(nModesX,nModesY), &
            dyuzx(nModesX,nModesY), dyyuzx(nModesX,nModesY), &
            dyuzy(nModesX,nModesY), dyyuzy(nModesX,nModesY), &
            dyuzz(nModesX,nModesY), dyyuzz(nModesX,nModesY) )

        allocate ( &
            dxyuxx(nModesX,nModesY), dxyuxy(nModesX,nModesY), dxyuxz(nModesX,nModesY), &
            dxyuyx(nModesX,nModesY), dxyuyy(nModesX,nModesY), dxyuyz(nModesX,nModesY), &
            dxyuzx(nModesX,nModesY), dxyuzy(nModesX,nModesY), dxyuzz(nModesX,nModesY) )
                    
        do i = 1, nModesX
            do j = 1, nModesY

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
