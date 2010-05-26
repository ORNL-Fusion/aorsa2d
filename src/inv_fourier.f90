module inv_fourier

implicit none

contains

    subroutine sftinv2d( a, f )

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, chebyshevX, chebyshevY
        use grid, &
        only: nMin, nMax, mMin, mMax, &
            xGrid_basis, yGrid_basis, &
            xBasis, yBasis
 
        implicit none
        
        complex, intent(in) :: a(:,:)
        complex, intent(inout), optional, allocatable :: &
            f(:,:)

        complex :: bFn
        integer :: i, j, n, m
        integer :: nS, nF, mS, mF

        if (.not. allocated ( f ) ) allocate ( f(nPtsX,nPtsY) )

        f = 0

        if(chebyshevX)then
            nS = nMin
            nF = nMax*2/3
        else
            nS = nMin*2/3
            nF = nMax*2/3
        endif

        if(chebyshevY)then
            mS = mMin
            mF = mMax*2/3
        else
            mS = mMin*2/3
            mF = mMax*2/3
        endif


        do i = 1, nPtsX
            do j = 1, nPtsY
        
                do n = nS, nF
                    do m = mS, mF

                      bFn = xBasis(n,xGrid_basis(i)) * yBasis(m,yGrid_basis(j))

                      f(i,j) = f(i,j) + a(n-nMin+1,m-mMin+1) * bFn

                    enddo
                enddo

            enddo
        enddo
    
    
    end subroutine sftinv2d

end module inv_fourier

