module antenna

use constants

implicit none

real :: antSigX, antSigY
#ifndef dblprec
    complex, allocatable :: brhs(:)
#else
    complex(kind=dbl), allocatable :: brhs(:)
#endif
complex, allocatable :: brhs_global(:)

contains

    subroutine alloc_total_brhs ( nPts_tot )

        use parallel

        implicit none

        integer, intent(in) :: nPts_tot

        allocate ( brhs_global(nPts_tot*3) )
#ifdef par
        allocate ( brhs(nRowLocal) )
#else
        allocate ( brhs(nPts_tot*3) )
#endif
        brhs        = 0
        brhs_global = 0

    end subroutine alloc_total_brhs


    subroutine init_brhs ( g )

        use aorsa2din_mod, &
        only: rAnt, zAnt, npRow, npCol, &
            antSigX, antSigY, &
            metalLeft, metalRight, metalTop, metalBot, &
            useEqdsk, r0, rhoAnt, antSigRho
        use grid
        use parallel
        use profiles, &
        only: omgrf

        implicit none

        type(gridBlock), intent(inout) :: g

        integer :: i, j, iRow, iCol

        ! scalapack index variables

        integer :: pr_sp, pc_sp, l_sp, m_sp, x_sp, y_sp
        integer :: localRow, localCol, ii, jj


        allocate ( &
            g%jR(g%nR,g%nZ), &
            g%jT(g%nR,g%nZ), &
            g%jZ(g%nR,g%nZ) )

        !   note curden is in Amps per meter of toroidal length (2.*pi*rt).

        do i = 1, g%nR
            do j = 1, g%nZ

                g%jR(i,j) = 0.0
                g%jT(i,j) = 0.0
                if(useEqdsk) then
                    if(g%R(i)>r0) &
                    g%jZ(i,j) = exp ( &
                    -( (g%rho(i,j)-rhoAnt)**2/antSigRho**2 + (g%Z(j)-zAnt)**2/antSigY**2 ) &
                          )
                else
                    g%jZ(i,j) = exp ( &
                    -( (g%R(i)-rAnt)**2/antSigX**2 + (g%Z(j)-zAnt)**2/antSigY**2 ) &
                          )
                endif

                !   boundary conditions
                !   -------------------

                if(g%nR/=1) then
                if ( i==1 .or. i==g%nR ) then
                    g%jR(i,j)    = 0
                    g%jT(i,j)    = 0
                    g%jZ(i,j)    = 0
                endif
                endif

                if(g%nZ/=1) then
                if ( j==1 .or. j==g%nZ ) then
                    g%jR(i,j)    = 0
                    g%jT(i,j)    = 0
                    g%jZ(i,j)    = 0
                endif
                endif


                !if ( capR(i) < metalLeft .or. capR(i) > metalRight &
                !        .or. y(j) > metalTop .or. y(j) < metalBot ) then 
                !    jR(i,j)    = 0
                !    jT(i,j)    = 0
                !    jZ(i,j)    = 0
                !endif


           enddo
        enddo

        !jT = 0
        !jT(nPtsX/2,nPtsY/2)    = 1

        g%jR = -zi / omgrf / eps0 * g%jR 
        g%jT = -zi / omgrf / eps0 * g%jT
        g%jZ = -zi / omgrf / eps0 * g%jZ

        do i = 1, g%nR
            do j = 1, g%nZ

                iRow = (j-1) * 3 + (i-1) * g%nZ * 3 + 1
                iRow = iRow + ( g%startRow-1 )

                do ii = 0, 2

#ifdef par
                    !   2D (with only 1 col) block cyclic storage, see:
                    !   http://www.netlib.org/scalapack/slug/node76.html
                    !       and
                    !   http://acts.nersc.gov/scalapack/hands-on/example4.html
                    !
                    !   note that brhs is distributed only in col 0 of the 
                    !   process grid. All other columns possess an empty
                    !   local portion of brhs

                    !iCol    = 1 
                    !jj      = 0
                    l_sp    = ( iRow-1+ii ) / ( npRow * rowBlockSize )
                    !m_sp    = ( iCol-1+jj ) / ( npCol * colBlockSize )

                    pr_sp   = mod ( rowStartProc + (iRow-1+ii)/rowBlockSize, npRow )
                    !pc_sp   = mod ( colStartProc + (iCol-1+jj)/colBlockSize, npCol )

                    x_sp    = mod ( iRow-1+ii, rowBlockSize ) + 1
                    !y_sp    = mod ( iCol-1+jj, colBlockSize ) + 1

                    localRow    = l_sp*rowBlockSize+x_sp
                    !localCol    = m_sp*colBlockSize+y_sp

                    if (myRow==pr_sp .and. myCol==pc_sp) then

                        if (ii==0) brhs(localRow)    = g%jR(i,j)
                        if (ii==1) brhs(localRow)    = g%jT(i,j)
                        if (ii==2) brhs(localRow)    = g%jZ(i,j)

                    endif
#else
                    if (ii==0) brhs(iRow+0)    = g%jR(i,j)
                    if (ii==1) brhs(iRow+1)    = g%jT(i,j)
                    if (ii==2) brhs(iRow+2)    = g%jZ(i,j)
#endif
                enddo

            enddo
        enddo
        
    end subroutine

end module antenna
