module grid

use constants
use chebyshev_mod

implicit none

!   init_grid
real, allocatable, dimension(:) :: capR, xkphi
real, allocatable, dimension(:) :: y, xGrid_basis, yGrid_basis
real :: dx, dy, xRange, yRange

!   init_k
!real, allocatable, dimension(:) :: kxsav, kysav
!real :: k_cutOff, kx_cutOff, ky_cutOff
integer :: nMin, nMax, mMin, mMax

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
                xkphi ( nPtsX ), &
                xGrid_basis(nPtsX) )
      
            if(nPtsX>1) then 
                xRange  = rwRight - rwLeft
                dx = xRange / (nPtsX-1)
                do i = 1, nPtsX
      
                    !xGrid_basis(i)  = -cos(pi*(i-1)/(nPtsX-1))
                    xGrid_basis(i)  = (i-1) * 2 * pi / ( nPtsX - 1 ) - pi

                    capR(i) = (i-1) * dx + rwLeft 
                    xkphi(i) = nPhi / capR(i)
        
                enddo
            else

                capR(1) = rwLeft
                xkphi(1) = nPhi / capR(1)
                
            endif

       
        !   define y mesh: y(j)
        !----------------------
        
            allocate ( y ( nPtsY ), &
                yGrid_basis(nPtsY) )
       
            if(nPtsY>1) then 

                yRange  = yTop - yBot
                dy = yRange / (nPtsY-1)
                do j = 1, nPtsY

                    !yGrid_basis(j)  = -cos(pi*(j-1)/(nPtsY-1))
                    yGrid_basis(j)  = (j-1) * 2 * pi / ( nPtsY - 1 ) - pi

                    y(j) = (j-1) * dy + yBot 
                    !y(j) = (yGrid_basis(j)+1)/2 * yRange + yBot
        
                enddo
            else
                y(1)    = 0
                yGrid_basis(1) = 0
            endif

    end subroutine init_grid


    subroutine init_k ()

        use aorsa2din_mod, &
        only: nModesX, nModesY, nPtsX, nPtsY, &
            xkperp_cutOff, xkx_cutOff, xky_cutOff
        use constants
        use parallel

        implicit none

        integer :: n, m

        nMin = -nModesX/2
        nMax =  nModesX/2
        mMin = -nModesY/2
        mMax =  nModesY/2

        ! Catch for even number of modes
        ! for producing a square matrix
        ! ------------------------------

        if (mod(nModesX,2)==0) nMin = nMin+1
        if (mod(nModesY,2)==0) mMin = mMin+1

        !nMin    = 0
        !nMax    = nModesX-1
        !mMin    = 0
        !mMax    = nModesY-1

        if (iAm==0) then
            write(*,*) '    n: ', nMin, nMax
            write(*,*) '    m: ', mMin, mMax
        endif

        !allocate ( &
        !    kxSav (nMin:nMax), &
        !    kySav (mMin:mMax) )

        !if(nPtsX/=1)then

        !    do n = nMin, nMax 
        !        !kxSav(n) = 2.0 * pi * n / (xRange+dx)
        !        kxSav(n) = 2.0 * pi * n / xRange
        !    enddo
        !else
        !        kxSav(0) = 0
        !endif

        !if(nPtsY/=1)then
        !    do m = mMin, mMax 
        !        !kySav(m) = 2.0 * pi * m / (yRange+dy) 
        !        kySav(m) = 2.0 * pi * m / yRange 
        !    enddo
        !else
        !        kySav(0) = 0
        !endif

        !k_cutOff   = sqrt ( maxVal ( abs(kxSav) )**2 &
        !                        + maxVal ( abs(kySav) )**2 ) * xkPerp_cutOff
        !kx_cutOff   = maxVal ( abs(kxSav) ) * xkx_cutOff
        !ky_cutOff   = maxVal ( abs(kySav) ) * xky_cutOff


    end subroutine init_k


    subroutine init_basis_functions ()

        use aorsa2din_mod, &
        only: nModesX, nModesY, nPtsX, nPtsY, &
                rwLeft, yBot
        use constants

        implicit none

        integer :: i, j, n, m

        allocate ( &
            xx(nMin:nMax,nPtsX), &
            yy(mMin:mMax,nPtsY) )

        do i = 1, nPtsX
            do n = nMin, nMax 
                !xx(n,i) = exp(zi * kxsav(n) * ( capR(i)-rwLeft-xRange/2 ) )
                !xx(n,i) = exp(zi * n * xGrid_basis(i) )
                xx(n,i) = xBasis(n,xGrid_basis(i))
                !xx(n,i) = exp(zi * kxsav(n) * ( capR(i)-rwLeft ) )
                !xx(n,i) = chebT(n+nMin,xGrid_basis(i))
            enddo
        enddo

        do j = 1, nPtsY
            do m = mMin, mMax 
                !yy(m,j) = exp(zi * kysav(m) * ( y(j)-yBot-yRange/2 ) )
                !yy(m,j) = exp(zi * m * yGrid_basis(j) )
                yy(m,j) = yBasis(m,yGrid_basis(j))
                !yy(m,j) = chebT(m+mMin,yGrid_basis(j))
            enddo
        enddo

    end subroutine init_basis_functions 


    ! Basis function derivative routines
    ! for Fourier & Chebyshev
    ! ----------------------------------

    function drBfn_bfn(i,j,n,m)

        implicit none

        integer, intent(in) :: i, j, n, m
        complex :: drBfn_bfn
        real :: kx

        kx = 2 * pi * n / xRange
        drBfn_bfn = zi * kx
        !drBfn_bfn = n * chebT(-1+n,xGrid_basis(i)) &
        !    / chebT(n,xGrid_basis(i))

    end function drBfn_bfn


    function dzBfn_bfn(i,j,n,m)

        implicit none

        integer, intent(in) :: i, j, n, m
        complex :: dzBfn_bfn
        real :: ky

        ky = 2 * pi * m / yRange
        dzBfn_bfn = zi * ky
        !dzBfn_bfn = m * chebU(-1+m,yGrid_basis(j)) &
        !    / chebT(m,yGrid_basis(j))

    end function dzBfn_bfn


    function drrBfn_bfn(i,j,n,m)

        implicit none

        integer, intent(in) :: i, j, n, m
        complex :: drrBfn_bfn
        real :: kx

        kx = 2 * pi * n / xRange
        drrBfn_bfn = - kx**2
        !drrBfn_bfn = n * &
        !    ( -n * chebU(-2+n,xGrid_basis(i)) &
        !        + (-1+n)*xGrid_basis(i) * chebU(-1+n,xGrid_basis(i) )) &
        !    / ((-1+xGrid_basis(i)**2)*chebT(n,xGrid_basis(i)))

    end function drrBfn_bfn


    function dzzBfn_bfn(i,j,n,m)

        implicit none

        integer, intent(in) :: i, j, n, m
        complex :: dzzBfn_bfn
        real :: ky

        ky = 2 * pi * m / yRange
        dzzBfn_bfn = - ky**2
        !dzzBfn_bfn = m * &
        !    ( -m*chebU(-2+m,yGrid_basis(j)) &
        !        + (-1+m)*yGrid_basis(j)*chebU(-1+m,yGrid_basis(j))) &
        !    / ( (-1+yGrid_basis(j)**2)*chebT(m,yGrid_basis(j)))

    end function dzzBfn_bfn


    function drzBfn_bfn(i,j,n,m)

        implicit none

        integer, intent(in) :: i, j, n, m
        complex :: drzBfn_bfn
        real :: kx, ky

        kx = 2 * pi * n / xRange
        ky = 2 * pi * m / yRange

        drzBfn_bfn = - kx * ky
        !drzBfn_bfn = ( m * n * &
        !    chebU(-1+m,yGrid_basis(j)) * chebU(-1+n,xGrid_basis(i)) ) &
        !    / ( chebT(m,yGrid_basis(j))*chebT(n,xGrid_basis(i)))

    end function drzBfn_bfn


    function xBasis (n,x)

        implicit none
        integer, intent(in) :: n
        real, intent(in) :: x
        complex :: xBasis

        xBasis = exp ( zi * n * x )
        !xBasis = chebT ( n, x )

    end function xBasis


    function yBasis (m,y)

        implicit none
        integer, intent(in) :: m
        real, intent(in) :: y
        complex :: yBasis

        yBasis = exp ( zi * m * y )
        !yBasis = chebT ( m, y )

    end function yBasis


end module grid
