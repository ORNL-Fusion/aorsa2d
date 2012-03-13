module isInside

implicit none

contains

function IsInsideOf ( rIn, zIn, rArr, zArr )

    implicit none

    logical :: IsInsideOf 
    real, intent(in) :: rIn, zIn
    real, intent(in) :: rArr(:), zArr(:)

    integer :: q1, q2, q3, q4
    real, allocatable :: newR(:), newZ(:)
    integer :: nInterp, nLim, i, j, cnt
    real :: m, b, d, dStep

    ! create an interpolated limiter boundary

    nInterp = 100
    nLim    = size ( rArr )
    cnt = 1

    allocate ( newR(nInterp*nLim+nLim), &
        newZ(nInterp*nLim+nLim) )

    newR(1:nLim)    = rArr
    newZ(1:nLim)    = zArr

    do i=1,nLim-1

        !   get slope

        dStep   = (rArr(i+1) - rArr(i)) / nInterp

        if (abs(dStep) .gt. 0 ) then 

            m   = ( zArr(i+1)-zArr(i) ) &
                    / ( rArr(i+1) - rArr(i) )
            b   = zArr(i) - m * rArr(i)

            !   distance

            d   = sqrt ( ( rArr(i+1) - rArr(i) )**2 &
                    + ( zArr(i+1) - zArr(i) )**2 )

            do j = 1, nInterp

                    newR(nLim+cnt)  = rArr(i) + dStep*j
                    newZ(nLim+cnt)  = m * (rArr(i) + dStep*j) + b 
                    cnt = cnt + 1

            enddo

        endif

    enddo

    !   interpolated eqdsk limite boundary

    q1  = count ( rIn - newR(1:nLim+cnt-1) .ge. 0 .and. zIn - newZ(1:nLim+cnt-1) .ge. 0 )
    q2  = count ( rIn - newR(1:nLim+cnt-1) .ge. 0 .and. zIn - newZ(1:nLim+cnt-1) .le. 0 )
    q3  = count ( rIn - newR(1:nLim+cnt-1) .le. 0 .and. zIn - newZ(1:nLim+cnt-1) .ge. 0 )
    q4  = count ( rIn - newR(1:nLim+cnt-1) .le. 0 .and. zIn - newZ(1:nLim+cnt-1) .le. 0 )

    !!   coarse eqdsk limiter boundary

    !q1  = count ( rIn - rArr .ge. 0 .and. zIn - zArr .ge. 0 )
    !q2  = count ( rIn - rArr .ge. 0 .and. zIn - zArr .le. 0 )
    !q3  = count ( rIn - rArr .le. 0 .and. zIn - zArr .ge. 0 )
    !q4  = count ( rIn - rArr .le. 0 .and. zIn - zArr .le. 0 )

    if ( q1 > 0 .and. q2 > 0 .and. q3 > 0 .and. q4 > 0 ) then

       IsInsideOf    = .true. 

    else

       IsInsideOf   = .false.

    endif
  
    deallocate ( newR, newZ )

    return

end function IsInsideOf 
 

end module isInside
