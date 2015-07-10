module setMetal

contains

    subroutine setMetalRegions ( g )

        use grid
        use aorsaNamelist, only : &
            metalLeft, metalRight, metalTop, metalBot, &
            limiter_boundary, useEqdsk, UseAR2Input, lcfs_boundary, useMetalBox
        use ar2Input, only: ar2_BbbMask=>BbbMask, ar2_LimMask=>LimMask, &
           ar2_nR=>nR, ar2_nZ=>nZ, ar2_r=>r, ar2_z=>z, ar2_LCFSMask=>bbbMask
        use parallel, only: iAm

        implicit none

        type(gridBlock), intent(inout) :: g

        integer :: w, i, j
        real :: iTmp,jTmp

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

                iTmp = (g%r(i)-ar2_r(1))/(ar2_r(ar2_nR)-ar2_r(1))*(ar2_nR-1)+1.0
                jTmp = (g%z(j)-ar2_z(1))/(ar2_z(ar2_nZ)-ar2_z(1))*(ar2_nZ-1)+1.0

#if __DebugSetMetal__ > 0 

                if(iTmp>ar2_nR)then 
                        write(*,*) 'iTmp error: ', iTmp, jTmp, i, j, g%nR, g%nZ, ar2_nR, ar2_nZ
                        stop 
                endif
                if(iTmp<1)then
                        write(*,*) 'iTmp error: ',  iTmp, jTmp, i, j, g%nR, g%nZ, ar2_nR, ar2_nZ
                        stop 
                endif
 
                if(jTmp>ar2_nZ)then
                        write(*,*) 'jTmp error: ', iTmp, jTmp, i, j, g%nR, g%nZ, ar2_nR, ar2_nZ
                        stop 
                endif
 
                if(jTmp<1)then
                        write(*,*) 'jTmp error: ', iTmp, jTmp, i, j, g%nR, g%nZ, ar2_nR, ar2_nZ
                        stop 
                endif
#endif 

                if(any( (/  ar2_LimMask(floor(iTmp),floor(jTmp)),&
                            ar2_LimMask(ceiling(iTmp),floor(jTmp)),&
                            ar2_LimMask(ceiling(iTmp),ceiling(jTmp)),&
                            ar2_LimMask(floor(iTmp),ceiling(jTmp)) /)==0 ) ) g%isMetal(w) = .true.

             enddo

        elseif(lcfs_boundary .and. UseAR2Input)then

             do w=1,size(g%pt)
                i = g%pt(w)%i
                j = g%pt(w)%j

                iTmp = (g%r(i)-ar2_r(1))/(ar2_r(ar2_nR)-ar2_r(1))*(ar2_nR-1)+1.0
                jTmp = (g%z(j)-ar2_z(1))/(ar2_z(ar2_nZ)-ar2_z(1))*(ar2_nZ-1)+1.0

                    if(i>1.and.i<g%nR.and.j>1.and.j<g%nZ)then
                        if(any( (/  ar2_LCFSMask(floor(iTmp),floor(jTmp)),&
                                    ar2_LCFSMask(ceiling(iTmp),floor(jTmp)),&
                                    ar2_LCFSMask(ceiling(iTmp),ceiling(jTmp)),&
                                    ar2_LCFSMask(floor(iTmp),ceiling(jTmp)) /)==0 ) ) g%isMetal(w) = .true.
                    else
                            g%isMetal(w) = .true.
                    endif
             enddo

        endif ! square box defined by metalLeft, metalRight, metalTop, metalBot

        if(useMetalBox)then 
            do w=1,size(g%pt)
                i = g%pt(w)%i
                j = g%pt(w)%j
                    if ( g%R(i) < metalLeft .or. g%R(i) > metalRight &
                            .or. g%Z(j) > metalTop .or. g%Z(j) < metalBot ) &
                        g%isMetal(w) = .true.
            enddo

        endif

        if(iAM==0)write(*,*) '    No. metal points: ', count(g%isMetal)

    end subroutine setMetalRegions

end module setMetal
