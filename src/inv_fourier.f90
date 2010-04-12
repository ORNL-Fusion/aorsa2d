module inv_fourier

implicit none

contains

    subroutine sftinv2d( a, f, fx, fy )

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nModesX, nModesY
        use grid, &
        only: kxL, kxR, kyL, kyR, xkxsav, xkysav, xx, yy
 
        implicit none
        
        complex, intent(in) :: a(:,:)
        complex, intent(inout), optional, allocatable :: &
            f(:,:), fx(:,:), fy(:,:)

        complex :: cexpkxky
        integer :: i, j, n, m

        if (.not. allocated ( f ) ) allocate ( f(nPtsX,nPtsY) )
        if ( present ( fx ) ) then
            if (.not. allocated ( fx ) ) allocate ( fx(nPtsX,nPtsY) )
        endif
        if ( present ( fy ) ) then 
            if (.not. allocated ( fy ) ) allocate ( fy(nPtsX,nPtsY) )
        endif

        f = 0
        if (present(fx)) fx = 0
        if (present(fy)) fy = 0

        do i = 1, nPtsX
            do j = 1, nPtsY
        
                do n = 1, nModesX 
                    do m = 1, nModesY 

                      !----------------------------------------------------
                      !cexpkxky = exp(zi * (xkx(n) * x(i) + xky(m) * y(j)))
                      !----------------------------------------------------

                      cexpkxky = xx(kxL+n-1, i) * yy(kyL+m-1, j)

                      f(i,j) = f(i,j) + a(n,m) * cexpkxky
                      if (present(fx)) fx(i,j) = f(i,j) + xkxsav(kxL+n-1) * a(n,m) * cexpkxky
                      if (present(fy)) fy(i,j) = f(i,j) + xkysav(kyL+m-1) * a(n,m) * cexpkxky

                    enddo
                enddo

            enddo
        enddo
    
    
    end subroutine sftinv2d

end module inv_fourier

