module grid

implicit none

!   init_grid
real, allocatable, dimension(:) :: capR, xkphi
real, allocatable, dimension(:) :: y
real :: dx, dy, xRange, yRange

!   init_k
integer :: nkx, nky, nRow, nCol
real, allocatable, dimension(:) :: xkxsav, xkysav
real :: xk_cutOff
integer :: kxL, kxR, kyL, kyR

!   init_basis_functions
complex, allocatable, dimension(:,:) :: xx, xx_inv, yy, yy_inv

contains

    subroutine init_grid ()

        use aorsa2din_mod, &
        only: rwLeft, rwRight, yTop, yBot, &
                nPhi, nModesX, nModesY

        implicit none

        integer :: i, j

        !   define x mesh: capr(i)
        !   --------------------------------------------
        
            allocate ( &
                capR ( nModesX ), &
                xkphi ( nModesX ) )
        
            rwLeft   = 0.2
            rwRight  = 1.7
            xRange  = rwRight - rwLeft
            dx = xRange / nModesX
            do i = 1, nModesX
        
                capr(i) = (i-1) * dx + rwLeft 
                xkphi(i) = nPhi / capr(i)
        
            enddo
        
        
        !   define y mesh: y(j), yprime(j)
        !---------------------------------
        
            allocate ( y ( nModesY ) )
        
            yTop    =  1.1
            yBot    = -1.1
            yRange  = yTop - yBot
            dy = yRange / nModesY
            do j = 1, nModesY
        
                y(j) = (j-1) * dy + yBot
        
            enddo


    end subroutine init_grid


    subroutine init_k ()

        use aorsa2din_mod, &
        only: nModesX, nModesY, xkperp_cutoff
        use constants

        implicit none

        integer :: n, m

        kxL = -nModesX/2+1
        kxR =  nModesX/2
        kyL = -nModesY/2+1
        kyR =  nModesY/2

        allocate ( &
            xkxsav (kxL:kxR), &
            xkysav (kyL:kyR) )

        nkx = size ( xkxsav, 1 )
        nky = size ( xkysav, 1 )

        do n = kxL, kxR 
            xkxsav(n) = 2.0 * pi * n / xRange
        enddo

        do m = kyL, kyR 
            xkysav(m) = 2.0 * pi * m / yRange
        enddo

        xk_cutoff   = sqrt ( xkxsav( kxR )**2 &
                                + xkysav( kyR )**2 ) * xkperp_cutoff

    end subroutine init_k


    subroutine init_basis_functions ()

        use aorsa2din_mod, &
        only: nModesX, nModesY
        use constants

        implicit none

        integer :: i, j, n, m

        allocate ( &
            xx(kxL:kxR,nModesX), xx_inv(kxL:kxR,nModesX), &
            yy(kyL:kyR,nModesY), yy_inv(kyL:kyR,nModesY) )

        do i = 1, nModesX
            do n = kxL, kxR 
                xx(n, i) = exp(zi * xkxsav(n) * capR(i))
                xx_inv(n,i) = 1.0/xx(n,i)
            enddo
        enddo

        do j = 1, nModesY
            do m = kyL, kyR 
                yy(m,j) = exp(zi * xkysav(m) * y(j))
                yy_inv(m,j) = 1.0/yy(m,j)
            enddo
        enddo

    end subroutine init_basis_functions 

end module grid
