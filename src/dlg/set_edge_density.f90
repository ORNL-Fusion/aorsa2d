module set_edge_density
use eqdsk_dlg

implicit none

contains

    real function density_by_gradient ( R, z, offset, gradient, limit )

        implicit none
        real, intent(in) :: R, z, offset, gradient, limit

        real :: rDist, zDist, dist
        real, allocatable :: newR(:), newZ(:)
        integer :: nInterp, nbbbs, i, j, cnt
        real :: m, b, d, dStep

        ! create an interpolated bbbs boundary

        nInterp = 10
        nbbbs    = size ( rbbbs )
        cnt = 1

        allocate ( newR(nInterp*nbbbs+nbbbs), &
            newZ(nInterp*nbbbs+nbbbs) )
    
        newR(1:nbbbs)    = rbbbs
        newZ(1:nbbbs)    = zbbbs

        do i=1,nbbbs-1

            !   get slope

            m   = ( zbbbs(i+1)-zbbbs(i) ) &
                    / ( rbbbs(i+1) - rbbbs(i) )
            b   = zbbbs(i) - m * rbbbs(i)

            !   distance

            d   = sqrt ( ( rbbbs(i+1) - rbbbs(i) )**2 &
                    + ( zbbbs(i+1) - zbbbs(i) )**2 )

            dStep   = (rbbbs(i+1) - rbbbs(i)) / nInterp

            do j = 1, nInterp

                if (abs(dStep) .gt. 0 ) then 
                    newR(nbbbs+cnt)  = rbbbs(i) + dStep*j
                    newZ(nbbbs+cnt)  = m * (rbbbs(i) + dStep*j) + b 
                    cnt = cnt + 1
                endif

            enddo

        enddo

 
        ! get distance from lcfs to point

        rDist   = minVal ( abs ( newR - R ) )
        zDist   = minVal ( abs ( newZ - z ) )

        dist    = minVal ( sqrt ( (newR-R)**2 + (newZ-z)**2 ) )
        deallocate ( newR, newZ )

        density_by_gradient = gradient * dist + offset
        
        if (density_by_gradient .lt. limit) &
            density_by_gradient = limit

        write(*,*) '***********', dist, density_by_gradient, gradient, offset
    
    end function density_by_gradient

end module set_edge_density
