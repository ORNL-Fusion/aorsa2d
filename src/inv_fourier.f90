module inv_fourier

implicit none

contains

    subroutine sftinv2d( a, f )

        use aorsa2din_mod, &
        only: nPtsX, nPtsY
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

        if (.not. allocated ( f ) ) allocate ( f(nPtsX,nPtsY) )

        f = 0

        do i = 1, nPtsX
            do j = 1, nPtsY
        
                do n = nMin, nMax
                    do m = mMin, mMax

                      bFn = xBasis(n,xGrid_basis(i)) * yBasis(m,yGrid_basis(j))

                      f(i,j) = f(i,j) + a(n-nMin+1,m-mMin+1) * bFn

                    enddo
                enddo

            enddo
        enddo
    
    
    end subroutine sftinv2d

end module inv_fourier

