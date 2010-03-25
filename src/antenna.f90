module antenna

implicit none

real :: antSigX, antSigY
complex, allocatable :: brhs(:), brhs_global(:)
complex, allocatable, dimension(:,:) :: &
    xjx, xjy, xjz

contains

    subroutine init_brhs ()

        use aorsa2din_mod, &
        only: rAnt, nPtsX, nPtsY, npRow, npCol
        use grid
        use constants
        use profiles
        use parallel

        implicit none

        integer :: i, j, iRow, iCol

        ! scalapack index variables

        integer :: pr_sp, pc_sp, l_sp, m_sp, x_sp, y_sp
        integer :: localRow, localCol, ii, jj

        allocate ( brhs_global(nPtsX*nPtsY*3) )
        allocate ( brhs(nRowLocal) )

        brhs        = 0
        brhs_global = 0

        allocate ( &
            xjx(nPtsX,nPtsY), &
            xjy(nPtsX,nPtsY), &
            xjz(nPtsX,nPtsY) )

        !   note curden is in Amps per meter of toroidal length (2.*pi*rt).

        antSigX = 0.1
        antSigY = 0.9

        do i = 1, nPtsX
            do j = 1, nPtsY

                xjx(i,j) = 0.0
                xjy(i,j) = 1.0 / dx &
                    * exp ( &
                    -( (capR(i)-rant)**2/antSigX**2 + (y(j)-0.0)**2/antSigY**2 ) &
                          )
                xjz(i,j) = 0.0

                !!   boundary conditions
                !!   -------------------

                !if ( i==1 .or. i==nPtsX &
                !        .or. j==1 .or. j==nPtsY ) then
                !    xjx(i,j)    = 0
                !    xjy(i,j)    = 0
                !    xjz(i,j)    = 0
                !endif

           enddo
        enddo

        xjy = 0
        xjy(nPtsX/2,nPtsY/2)    = 1

        xjx = -zi / omgrf / eps0 * xjx
        xjy = -zi / omgrf / eps0 * xjy
        xjz = -zi / omgrf / eps0 * xjz

        do i = 1, nPtsX
            do j = 1, nPtsY

                iRow = (j-1) * 3 + (i-1) * nPtsY * 3 + 1

                do ii = 0, 2

                    !   2D (with only 1 col) block cyclic storage, see:
                    !   http://www.netlib.org/scalapack/slug/node76.html
                    !       and
                    !   http://acts.nersc.gov/scalapack/hands-on/example4.html
                    !
                    !   note that brhs is distributed only in col 0 of the 
                    !   process grid. All other columns possess an empty
                    !   local portion of brhs

                    iCol    = 1 
                    jj      = 0
                    l_sp    = ( iRow-1+ii ) / ( npRow * rowBlockSize )
                    m_sp    = ( iCol-1+jj ) / ( npCol * colBlockSize )

                    pr_sp   = mod ( rowStartProc + (iRow-1+ii)/rowBlockSize, npRow )
                    pc_sp   = mod ( colStartProc + (iCol-1+jj)/colBlockSize, npCol )

                    x_sp    = mod ( iRow-1+ii, rowBlockSize ) + 1
                    y_sp    = mod ( iCol-1+jj, colBlockSize ) + 1

                    localRow    = l_sp*rowBlockSize+x_sp
                    localCol    = m_sp*colBlockSize+y_sp

                    if (myRow==pr_sp .and. myCol==pc_sp) then

                        if (ii==0) brhs(localRow)    = xjx(i,j)
                        if (ii==1) brhs(localRow)    = xjy(i,j)
                        if (ii==2) brhs(localRow)    = xjz(i,j)

                    endif

                enddo

            enddo
        enddo
        
    end subroutine

end module antenna
