module grid

use constants
use chebyshev_mod
use aorsa2din_mod, &
only: chebyshev

implicit none

!   init_grid
real, allocatable, dimension(:) :: capR, xkphi
real, allocatable, dimension(:) :: y, xGrid_basis, yGrid_basis
real :: xRange, yRange, normFacX, normFacY

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

        ! define x mesh: capr(i)
        ! ----------------------
        
            allocate ( &
                capR ( nPtsX ), &
                xkphi ( nPtsX ), &
                xGrid_basis(nPtsX) )
      
            if(nPtsX>1) then 
                xRange  = rwRight - rwLeft
                !dx = xRange / (nPtsX-1)
                do i = 1, nPtsX
     
                    if(chebyshev) then
                        ! Gauss-Lobatto grid [-1,1] (Chebyshev basis)
                        xGrid_basis(i)  = -cos(pi*(i-1)/(nPtsX-1))
                    else
                        ! Uniform grid [0,2pi] (Fourier basis)
                        xGrid_basis(i)  = (i-1) * 2 * pi / ( nPtsX - 1 )
                    endif
       
                enddo
            else

                capR(1) = rwLeft
                xGrid_basis(1) = 0
                
            endif

            if(nPtsX>1) &
            capR = (xGrid_basis-xGrid_basis(1)) &
                / (xGrid_basis(nPtsX)-xGrid_basis(1)) * xRange + rwLeft

            xkphi = nPhi / capR

            if(chebyshev)then
                normFacX = 2 / xRange
            else
                normFacX = 2 * pi / xRange
            endif

       
        ! define y mesh: y(j)
        !--------------------
        
            allocate ( y ( nPtsY ), &
                yGrid_basis(nPtsY) )
       
            if(nPtsY>1) then 

                yRange  = yTop - yBot
                !dy = yRange / (nPtsY-1)
                do j = 1, nPtsY

                    if(chebyshev) then
                        ! Gauss-Lobatto grid [-1,1] (Chebyshev basis)
                        yGrid_basis(j)  = -cos(pi*(j-1)/(nPtsY-1))
                    else
                        ! Uniform grid [0,2pi] (Fourier basis)
                        yGrid_basis(j)  = (j-1) * 2 * pi / ( nPtsY - 1 )
                    endif

                enddo
            else
                y(1)    = 0
                yGrid_basis(1) = 0
            endif

            if(nPtsY>1) &
            y = (yGrid_basis-yGrid_basis(1)) &
                / (yGrid_basis(nPtsY)-yGrid_basis(1)) * yRange + yBot

            if(chebyshev) then
                normFacY = 2 / yRange
            else 
                normFacY = 2 * pi / yRange
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

        if(chebyshev) then

            nMin    = 0
            nMax    = nModesX-1
            mMin    = 0
            mMax    = nModesY-1

        else

            nMin = -nModesX/2
            nMax =  nModesX/2
            mMin = -nModesY/2
            mMax =  nModesY/2

            ! Catch for even number of modes
            ! for producing a square matrix
            ! ------------------------------

            if (mod(nModesX,2)==0) nMin = nMin+1
            if (mod(nModesY,2)==0) mMin = mMin+1

        endif

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
                xx(n,i) = xBasis(n,xGrid_basis(i))
            enddo
        enddo

        do j = 1, nPtsY
            do m = mMin, mMax 
                yy(m,j) = yBasis(m,yGrid_basis(j))
            enddo
        enddo

    end subroutine init_basis_functions 


    ! Basis function and their derivative routines
    ! for Fourier & Chebyshev
    ! --------------------------------------------

    function xBasis (n,x)

        implicit none
        integer, intent(in) :: n
        real, intent(in) :: x
        complex :: xBasis

        if(chebyshev) then
            xBasis = chebT ( n, x )
        else
            xBasis = exp ( zi * n * x )
        endif

    end function xBasis


    function yBasis (m,y)

        implicit none
        integer, intent(in) :: m
        real, intent(in) :: y
        complex :: yBasis

        if(chebyshev) then
            yBasis = chebT ( m, y )
        else
            yBasis = exp ( zi * m * y )
        endif

    end function yBasis


    function drBfn_bfn(i,j,n,m)

        implicit none

        integer, intent(in) :: i, j, n, m
        complex :: drBfn_bfn

        if(chebyshev) then
            drBfn_bfn = n * chebT(-1+n,xGrid_basis(i)) &
                / chebT(n,xGrid_basis(i))
        else
            drBfn_bfn = zi * n
        endif

        drBfn_bfn = drBfn_bfn * normFacX

    end function drBfn_bfn


    function dzBfn_bfn(i,j,n,m)

        implicit none

        integer, intent(in) :: i, j, n, m
        complex :: dzBfn_bfn
        
        if(chebyshev) then
            dzBfn_bfn = m * chebU(-1+m,yGrid_basis(j)) &
                / chebT(m,yGrid_basis(j))
        else
            dzBfn_bfn = zi * m
        endif

        dzBfn_bfn = dzBfn_bfn * normFacY

    end function dzBfn_bfn


    function drrBfn_bfn(i,j,n,m)

        implicit none

        integer, intent(in) :: i, j, n, m
        complex :: drrBfn_bfn

        if(chebyshev) then

            drrBfn_bfn = n * &
                ( -n * chebU(-2+n,xGrid_basis(i)) &
                    + (-1+n)*xGrid_basis(i) * chebU(-1+n,xGrid_basis(i) )) &
                / ((-1+xGrid_basis(i)**2)*chebT(n,xGrid_basis(i)))
        else
            drrBfn_bfn = - n**2
        endif

        drrBfn_bfn = drrBfn_bfn * normFacX**2

    end function drrBfn_bfn


    function dzzBfn_bfn(i,j,n,m)

        implicit none

        integer, intent(in) :: i, j, n, m
        complex :: dzzBfn_bfn

        if(chebyshev) then 
            dzzBfn_bfn = m * &
                ( -m*chebU(-2+m,yGrid_basis(j)) &
                    + (-1+m)*yGrid_basis(j)*chebU(-1+m,yGrid_basis(j))) &
                / ( (-1+yGrid_basis(j)**2)*chebT(m,yGrid_basis(j)))
        else
            dzzBfn_bfn = - m**2
        endif

        dzzBfn_bfn = dzzBfn_bfn * normFacY**2

    end function dzzBfn_bfn


    function drzBfn_bfn(i,j,n,m)

        implicit none

        integer, intent(in) :: i, j, n, m
        complex :: drzBfn_bfn

        if(chebyshev) then
            drzBfn_bfn = ( m * n * &
                chebU(-1+m,yGrid_basis(j)) * chebU(-1+n,xGrid_basis(i)) ) &
                / ( chebT(m,yGrid_basis(j))*chebT(n,xGrid_basis(i)))
        else
            drzBfn_bfn = - n * m 
        endif

        drzBfn_bfn = drzBfn_bfn * ( normFacX * normFacY )

    end function drzBfn_bfn

end module grid
