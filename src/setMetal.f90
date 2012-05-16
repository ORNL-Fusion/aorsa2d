module setMetal

contains

    subroutine setMetalRegions ( g )

        use grid
        use aorsaNamelist, only : &
            metalLeft, metalRight, metalTop, metalBot, &
            limiter_boundary, useEqdsk, UseAR2Input
        use eqdsk_dlg, only: is_inside_lim, rLim__,zLim
        use IsInside, only: IsInsideOf
        use ar2Input, only: ar2_rLim=>rLim, ar2_zLim=>zLim
        use parallel, only: iAm

        implicit none

        type(gridBlock), intent(inout) :: g

        integer :: w, i, j

        !allocate (g%isMetal(g%nR,g%nZ))
        allocate (g%isMetal(size(g%pt)))

        g%isMetal = .false.

        !if(useEqdsk)then
        !where(grid%rho>=0.99)
        !        isMetal = .true.
        !endwhere
        !endif

        if(limiter_boundary .and. useEqdsk)then
    
            !do i=1,g%nR
            !    do j=1,g%nZ
            !       g%isMetal(i,j) = .not. IsInsideOf ( g%R(i), g%z(j), rLim__, zLim )
            !    enddo
            !enddo
            stop

        elseif(limiter_boundary .and. UseAR2Input)then

             do w=1,size(g%pt)
                i = g%pt(w)%i
                j = g%pt(w)%j
             !do i=1,g%nR
             !   do j=1,g%nZ
                    g%isMetal(w) = .not. IsInsideOf ( g%R(i), g%z(j), ar2_rLim, ar2_zLim )
             !   enddo
             !enddo
             enddo

        else ! square box defined by metalLeft, metalRight, metalTop, metalBot

            do w=1,size(g%pt)
                i = g%pt(w)%i
                j = g%pt(w)%j
            !do i=1,g%nR
            !    do j=1,g%nZ
                    if ( g%R(i) < metalLeft .or. g%R(i) > metalRight &
                            .or. g%Z(j) > metalTop .or. g%Z(j) < metalBot ) &
                        g%isMetal(w) = .true.
            !    enddo
            !enddo
            enddo

        endif

        if(iAM==0)write(*,*) '    No. metal points: ', count(g%isMetal)

    end subroutine setMetalRegions

end module setMetal
