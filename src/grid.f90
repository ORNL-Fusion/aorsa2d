module grid

implicit none

!   init_grid
real, allocatable, dimension(:) :: capR, xkphi
real, allocatable, dimension(:) :: y
real :: dx, dy, xRange, yRange

!   init_k
real, allocatable, dimension(:) :: xkxsav, xkysav
real :: xk_cutOff
integer :: kxL, kxR, kyL, kyR

!   init_basis_functions
complex, allocatable, dimension(:,:) :: xx, xx_inv, yy, yy_inv

contains

    subroutine init_grid ()

        use aorsa2din_mod, &
        only: rwLeft, rwRight, yTop, yBot, &
                nPhi, nModesX, nModesY, nPtsX, nPtsY

        implicit none

        integer :: i, j

        !   define x mesh: capr(i)
        !   ----------------------
        
            allocate ( &
                capR ( nPtsX ), &
                xkphi ( nPtsX ) )
        
            xRange  = rwRight - rwLeft
            dx = xRange / (nPtsX-1)
            do i = 1, nPtsX
        
                capr(i) = (i-1) * dx + rwLeft 
                xkphi(i) = nPhi / capr(i)
        
            enddo
        
        !   define y mesh: y(j)
        !----------------------
        
            allocate ( y ( nPtsY ) )
        
            yRange  = yTop - yBot
            dy = yRange / (nPtsY-1)
            do j = 1, nPtsY
        
                y(j) = (j-1) * dy + yBot
        
            enddo

    end subroutine init_grid


    subroutine init_k ()

        use aorsa2din_mod, &
        only: nModesX, nModesY, xkperp_cutoff
        use constants
        use parallel

        implicit none

        integer :: n, m

        kxL = -nModesX/2
        kxR =  nModesX/2
        kyL = -nModesY/2
        kyR =  nModesY/2


        ! Catch for even number of modes
        ! for producing a square matrix
        ! ------------------------------

        if (mod(nModesX,2)==0) kXL = kXL+1
        if (mod(nModesY,2)==0) kYL = kYL+1

        if (iAm==0) then
            write(*,*) '    kx: ', kxL, kxR
            write(*,*) '    ky: ', kyL, kyR
        endif

        allocate ( &
            xkxsav (kxL:kxR), &
            xkysav (kyL:kyR) )

        do n = kxL, kxR 
            xkxsav(n) = 2.0 * pi * n / xRange 
        enddo

        do m = kyL, kyR 
            xkysav(m) = 2.0 * pi * m / yRange 
        enddo

        !write(*,*) 'WARNING: Running with /4 on the mode resolution'

        xk_cutoff   = sqrt ( xkxsav( kxR )**2 &
                                + xkysav( kyR )**2 ) * xkperp_cutoff

    end subroutine init_k


    subroutine init_basis_functions ()

        use aorsa2din_mod, &
        only: nModesX, nModesY, nPtsX, nPtsY
        use constants

        implicit none

        integer :: i, j, n, m

        allocate ( &
            xx(kxL:kxR,nPtsX), xx_inv(kxL:kxR,nPtsX), &
            yy(kyL:kyR,nPtsY), yy_inv(kyL:kyR,nPtsY) )

        do i = 1, nPtsX
            do n = kxL, kxR 
                xx(n, i) = exp(zi * xkxsav(n) * (capR(i)-minVal(capR)))
                xx_inv(n,i) = 1.0/xx(n,i)
            enddo
        enddo

        do j = 1, nPtsY
            do m = kyL, kyR 
                yy(m,j) = exp(zi * xkysav(m) * (y(j)-minVal(y)))
                yy_inv(m,j) = 1.0/yy(m,j)
            enddo
        enddo

    end subroutine init_basis_functions 

end module grid
