module grid

use constants
use chebyshev_mod
use aorsa2din_mod, &
only: chebyshevX, chebyshevY

implicit none

type :: gridBlock

    integer :: nR, nZ
    real, allocatable, dimension(:) :: rNorm, zNorm, R, z, kPhi
    real :: rMin, rMax, zMin, zMax, rRange, zRange
    real :: normFacR, normFacZ
    integer :: nMin, nMax, mMin, mMax
    real :: k_cutOff
    complex, allocatable, dimension(:,:) :: xx, yy    

end type gridBlock

type(gridBlock), allocatable :: allGrids(:)

!   init_grid
real, allocatable, dimension(:) :: capR, kPhi
real, allocatable, dimension(:) :: y, xGrid_basis, yGrid_basis
real :: xRange, yRange, normFacX, normFacY

!   init_k
real :: k_cutOff!, kx_cutOff, ky_cutOff
integer :: nMin, nMax, mMin, mMax

!   init_basis_functions
complex, allocatable, dimension(:,:) :: xx, yy

contains

    function init_gridBlock ( nR, nZ, rMin, rMax, zMin, zMax ) result ( grid )

        use aorsa2din_mod, &
        only : nPhi, xkPerp_cutOff
        use parallel

        implicit none

        type(gridBlock) :: grid

        integer, intent(in) :: nR, nZ
        real, intent(in) :: rMin, rMax, zMin, zMax

        integer :: i, j, n, m

        allocate ( grid%rNorm(nR), &
                    grid%zNorm(nZ), &
                    grid%R(nR), &
                    grid%z(nZ), &
                    grid%kPhi(nR) )

            ! Populate grid block parameters
            ! ------------------------------

            grid%rRange  = rMax - rMin    

            grid%rMin = rMin
            grid%rMax = rMax
            grid%zMin = zMin
            grid%zMax = zMax
            grid%nR = nR
            grid%nZ = nZ

            
            ! Create the grid block grids
            ! ---------------------------

            if(nR>1) then 

                do i = 1, nR
     
                    if(chebyshevX) then
                        ! Gauss-Lobatto grid [-1,1] (Chebyshev basis)
                        grid%rNorm(i)  = -cos(pi*(i-1)/(nR-1))
                    else
                        ! Uniform grid [0,2pi] (Fourier basis)
                        grid%rNorm(i)  = (i-1) * 2 * pi / ( nR - 1 )
                    endif
       
                enddo
            else

                grid%R(1) = rMin
                grid%rNorm(1) = 0
                
            endif

            if(nR>1) &
            grid%R = (grid%rNorm-grid%rNorm(1)) &
                / (grid%rNorm(nR)-grid%rNorm(1)) * grid%rRange + rMin


            if(chebyshevX)then
                grid%normFacR = 2 / grid%rRange
            else
                grid%normFacR = 2 * pi / grid%rRange
            endif

       
            grid%zRange  = zMax - zMin
            if(nZ>1) then 

                do j = 1, nZ

                    if(chebyshevY) then
                        ! Gauss-Lobatto grid [-1,1] (Chebyshev basis)
                        grid%zNorm(j)  = -cos(pi*(j-1)/(nZ-1))
                    else
                        ! Uniform grid [0,2pi] (Fourier basis)
                        grid%zNorm(j)  = (j-1) * 2 * pi / ( nZ - 1 )
                    endif

                enddo
            else
                grid%z(1)    = 0
                grid%zNorm(1) = 0
            endif

            if(nZ>1) &
            grid%z = (grid%zNorm-grid%zNorm(1)) &
                / (grid%zNorm(nZ)-grid%zNorm(1)) * grid%zRange + zMin

            if(chebyshevY) then
                grid%normFacZ = 2 / grid%zRange
            else 
                grid%normFacZ = 2 * pi / grid%zRange
            endif


            ! Set the grid block kPhi
            ! -----------------------

            grid%kPhi = nPhi / grid%R


            ! Set the grid block n,m ranges
            ! -----------------------------

            if(chebyshevX) then

                grid%nMin    = 0
                grid%nMax    = nR-1

            else

                grid%nMin = -nR/2
                grid%nMax =  nR/2

                ! Catch for even number of modes
                ! for producing a square matrix
                ! ------------------------------

                if (mod(nR,2)==0) grid%nMin = grid%nMin+1

            endif


            if(chebyshevY) then

                grid%mMin    = 0
                grid%mMax    = nZ-1

            else

                grid%mMin = -nZ/2
                grid%mMax =  nZ/2
                if (mod(nZ,2)==0) grid%mMin = grid%mMin+1

            endif

            if (iAm==0) then
                write(*,*) '    n: ', grid%nMin, grid%nMax
                write(*,*) '    m: ', grid%mMin, grid%mMax
            endif


            ! Set the k above which will be damped
            ! ------------------------------------

            grid%k_cutOff = xkPerp_cutOff * sqrt(&
               (grid%nMax * grid%normFacR)**2+(grid%mMax*grid%normFacZ)**2)


            ! Fill the grid block basis functions
            ! -----------------------------------

            allocate ( &
                grid%xx(grid%nMin:nMax,nR), &
                grid%yy(grid%mMin:mMax,nZ) )

            do i = 1, nR
                do n = grid%nMin, grid%nMax 
                    grid%xx(n,i) = xBasis(n,grid%rNorm(i))
                enddo
            enddo

            do j = 1, nZ
                do m = grid%mMin, grid%mMax 
                    grid%yy(m,j) = yBasis(m,grid%zNorm(j))
                enddo
            enddo


    end function init_gridBlock


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
                kPhi ( nPtsX ), &
                xGrid_basis(nPtsX) )

            xRange  = rwRight - rwLeft     
            if(nPtsX>1) then 

                do i = 1, nPtsX
     
                    if(chebyshevX) then
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

            kPhi = nPhi / capR

            if(chebyshevX)then
                normFacX = 2 / xRange
            else
                normFacX = 2 * pi / xRange
            endif

       
        ! define y mesh: y(j)
        ! -------------------
        
            allocate ( y ( nPtsY ), &
                yGrid_basis(nPtsY) )
       
            yRange  = yTop - yBot
            if(nPtsY>1) then 

                do j = 1, nPtsY

                    if(chebyshevY) then
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

            if(chebyshevY) then
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

        if(chebyshevX) then

            nMin    = 0
            nMax    = nModesX-1

        else

            nMin = -nModesX/2
            nMax =  nModesX/2

            ! Catch for even number of modes
            ! for producing a square matrix
            ! ------------------------------

            if (mod(nModesX,2)==0) nMin = nMin+1

        endif


        if(chebyshevY) then

            mMin    = 0
            mMax    = nModesY-1

        else

            mMin = -nModesY/2
            mMax =  nModesY/2
            if (mod(nModesY,2)==0) mMin = mMin+1

        endif


        if (iAm==0) then
            write(*,*) '    n: ', nMin, nMax
            write(*,*) '    m: ', mMin, mMax
        endif

        k_cutOff = xkPerp_cutOff * sqrt(&
           (nMax * normFacX)**2+(mMax*normFacY)**2)

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

        if(chebyshevX) then
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

        if(chebyshevY) then
            yBasis = chebT ( m, y )
        else
            yBasis = exp ( zi * m * y )
        endif

    end function yBasis


    function drBfn_bfn(i,j,n,m)

        implicit none

        integer, intent(in) :: i, j, n, m
        complex :: drBfn_bfn

        if(chebyshevX) then
            drBfn_bfn = n * chebU(-1+n,xGrid_basis(i)) &
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
        
        if(chebyshevY) then
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

        if(chebyshevX) then

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

        if(chebyshevY) then 
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
        complex :: drzBfn_bfn, drBfn_bfn, dzBfn_bfn

        if(chebyshevX) then
            drBfn_bfn = n * chebU(-1+n,xGrid_basis(i)) &
                / chebT(n,xGrid_basis(i))
        else
            drBfn_bfn = zi * n
        endif

        if(chebyshevY) then
            dzBfn_bfn = m * chebU(-1+m,yGrid_basis(j)) &
                / chebT(m,yGrid_basis(j))
        else
            dzBfn_bfn = zi * m
        endif

        drzBfn_bfn = drBfn_bfn * dzBfn_bfn * ( normFacX * normFacY )

    end function drzBfn_bfn

end module grid
