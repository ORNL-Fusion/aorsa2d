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

        integer :: i, j

        allocate (g%isMetal(g%nR,g%nZ))
        g%isMetal = .false.

        !if(useEqdsk)then
        !where(grid%rho>=0.99)
        !        isMetal = .true.
        !endwhere
        !endif

        if(limiter_boundary .and. useEqdsk)then
    
            do i=1,g%nR
                do j=1,g%nZ
                    g%isMetal(i,j) = .not. IsInsideOf ( g%R(i), g%z(j), rLim__, zLim )
                enddo
            enddo

        elseif(limiter_boundary .and. UseAR2Input)then
            
             do i=1,g%nR
                do j=1,g%nZ
                    
                    g%isMetal(i,j) = .not. IsInsideOf ( g%R(i), g%z(j), ar2_rLim, ar2_zLim )
                enddo
            enddo

        else ! square box defined by metalLeft, metalRight, metalTop, metalBot

            do i=1,g%nR
                do j=1,g%nZ
                    if ( g%R(i) < metalLeft .or. g%R(i) > metalRight &
                            .or. g%Z(j) > metalTop .or. g%Z(j) < metalBot ) &
                        g%isMetal(i,j) = .true.
                enddo
            enddo

        endif

        if(iAM==0)write(*,*) 'No. metal points: ', count(g%isMetal)

    end subroutine setMetalRegions

end module setMetal
