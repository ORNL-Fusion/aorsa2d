module inv_fourier

implicit none

contains

    subroutine sftinv2d( g, rhs )

        use aorsaNamelist, &
            only: chebyshevX, chebyshevY, cosX, cosY, fracOfModesInSolution
        use grid
        use parallel, only : iAm, ICTXT
 
        implicit none

        type(gridBlock), intent(inout) :: g
        integer, intent(in) :: rhs
        
        complex :: bFn
        integer :: i, j, n, m, w
        integer :: nS, nF, mS, mF

        if (.not. allocated ( g%eAlpha ) ) allocate ( g%eAlpha(g%nR,g%nZ) )
        if (.not. allocated ( g%eBeta ) ) allocate ( g%eBeta(g%nR,g%nZ) )
        if (.not. allocated ( g%eB ) ) allocate ( g%eB(g%nR,g%nZ) )

        g%eAlpha = 0
        g%eBeta = 0
        g%eB = 0

        if(chebyshevX)then
            g%nS = g%nMin
            g%nF = g%nMax
        elseif(cosX)then
            g%nS = g%nMin
            g%nF = g%nMax
        else
            g%nS = g%nMin*fracOfModesInSolution
            g%nF = g%nMax*fracOfModesInSolution
        endif

        if(chebyshevY)then
            g%mS = g%mMin
            g%mF = g%mMax
        elseif(cosY)then
            g%mS = g%mMin
            g%mF = g%mMax
        else
            g%mS = g%mMin*fracOfModesInSolution
            g%mF = g%mMax*fracOfModesInSolution
            if(g%nZ==1)then
                g%mS = g%mMin
                g%mF = g%mMax
            endif
        endif


        !! NOTE: There may be an indexing error here after changing to the g
        !! type routine. If so it means the g%eAlphak etc have retained their
        !! neg through pos indexing that was not retained in the previous
        !! approach.

        do w=1,size(g%wl)

            n = g%wl(w)%n
            m = g%wl(w)%m
            i = g%wl(w)%i
            j = g%wl(w)%j

            twoThirdsRule: &
            if(m >= g%mMin*fracOfModesInSolution .and. m <= g%mMax*fracOfModesInSolution &
                .and. n >= g%nMin*fracOfModesInSolution .and. n <= g%nMax*fracOfModesInSolution ) then

                bFn = g%xx(n,i) * g%yy(m,j)
                g%eAlpha(i,j) = g%eAlpha(i,j) + g%eAlphak(n,m,rhs) * bFn
                g%eBeta(i,j) = g%eBeta(i,j) + g%eBetak(n,m,rhs) * bFn
                g%eB(i,j) = g%eB(i,j) + g%eBk(n,m,rhs) * bFn

                write(*,*) 'Ea: ', m, n, i, j, g%eAlphak(n,m,rhs)
                write(*,*) 'Eb: ', m, n, i, j, g%eBk(n,m,rhs)

            endif twoThirdsRule

        enddo

#ifdef par
        call blacs_barrier ( ICTXT, 'All' ) 

        call cGSUM2D ( ICTXT, 'All', ' ', g%nR, g%nZ, g%eAlpha(:,:), g%nR, -1, -1 )
        call cGSUM2D ( ICTXT, 'All', ' ', g%nR, g%nZ, g%eBeta(:,:), g%nR, -1, -1 )
        call cGSUM2D ( ICTXT, 'All', ' ', g%nR, g%nZ, g%eB(:,:), g%nR, -1, -1 )

        call blacs_barrier ( ICTXT, 'All' ) 
#endif
    end subroutine sftinv2d

end module inv_fourier

