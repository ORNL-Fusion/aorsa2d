module antenna

implicit none

real :: antSigX, antSigY
complex, allocatable :: brhs(:)
complex, allocatable, dimension(:,:) :: &
    xjx, xjy, xjz

contains

    subroutine init_brhs ()

        use aorsa2din_mod, &
        only: nModesX, nModesY, rAnt
        use grid
        use constants
        use profiles

        implicit none

        integer :: i, j, iRow

        allocate ( brhs(nModesX*nModesY*3) )
        allocate ( &
            xjx(nModesX,nModesY), &
            xjy(nModesX,nModesY), &
            xjz(nModesX,nModesY) )

        !   note curden is in Amps per meter of toroidal length (2.*pi*rt).

        antSigX = 0.1
        antSigY = 0.2

        do i = 1, nModesX
            do j = 1, nModesY

                xjx(i,j) = 0.0
                xjy(i,j) = 1.0 / dx &
                    * exp ( &
                    -( (capR(i)-rant)**2/antSigX**2 + (y(j)-0.0)**2/antSigY**2 ) &
                          )
                xjz(i,j) = 0.0

                !   boundary conditions
                !   -------------------

                if ( i==1 .or. i==nModesX &
                        .or. j==1 .or. j==nModesY ) then
                    xjx(i,j)    = 0
                    xjy(i,j)    = 0
                    xjz(i,j)    = 0
                endif

           enddo
        enddo

        xjx = -zi / omgrf / eps0 * xjx
        xjy = -zi / omgrf / eps0 * xjy
        xjz = -zi / omgrf / eps0 * xjz

        do i = 1, nModesX
            do j = 1, nModesY

                iRow = (j-1) * 3 + (i-1) * nModesY * 3 + 1

                brhs(iRow)      = xjx(i,j)
                brhs(iRow+1)    = xjy(i,j)
                brhs(iRow+2)    = xjz(i,j)

            enddo
        enddo

    end subroutine

end module antenna
