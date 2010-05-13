module grid

implicit none

!   init_grid
real, allocatable, dimension(:) :: capR, xkphi
real, allocatable, dimension(:) :: y
real :: dx, dy, xRange, yRange

!   init_k
real, allocatable, dimension(:) :: kxsav, kysav
real :: xk_cutOff
integer :: kxL, kxR, kyL, kyR

!   init_basis_functions
complex, allocatable, dimension(:,:) :: xx, yy

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
            kxsav (kxL:kxR), &
            kysav (kyL:kyR) )

        do n = kxL, kxR 
            kxsav(n) = 2.0 * pi * n / (xRange+dx)
        enddo

        do m = kyL, kyR 
            kysav(m) = 2.0 * pi * m / (yRange+dy) 
        enddo

        xk_cutOff   = sqrt ( maxVal ( abs(kxsav) )**2 &
                                + maxVal ( abs(kysav) )**2 ) * xkPerp_cutOff

    end subroutine init_k


    subroutine init_basis_functions ()

        use aorsa2din_mod, &
        only: nModesX, nModesY, nPtsX, nPtsY, &
                rwLeft, yBot
        use constants

        implicit none

        integer :: i, j, n, m

        allocate ( &
            xx(kxL:kxR,nPtsX), &
            yy(kyL:kyR,nPtsY) )

        do i = 1, nPtsX
            do n = kxL, kxR 
                xx(n, i) = exp(zi * kxsav(n) * ( capR(i)-rwLeft+dx/2 ) )
            enddo
        enddo

        do j = 1, nPtsY
            do m = kyL, kyR 
                yy(m,j) = exp(zi * kysav(m) * ( y(j)-yBot+dy/2 ) )
            enddo
        enddo

    end subroutine init_basis_functions 

end module grid
