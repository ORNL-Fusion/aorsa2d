module antenna

use constants

implicit none

real :: antSigX, antSigY
integer :: NRHS = 3

#ifndef dblprec
    complex, allocatable :: brhs(:,:)
#else
    complex(kind=dbl), allocatable :: brhs(:,:)
#endif
complex, allocatable :: brhs_global(:,:)

contains

    subroutine alloc_total_brhs ( nPts_tot )

        use parallel
        use AR2SourceLocationsInput

        implicit none

        integer(kind=long), intent(in) :: nPts_tot

        if(NRHS<1)then
                write(*,*) 'ERROR: NRHS<1'
                stop
        endif

        allocate ( brhs_global(nPts_tot*3,NRHS) )
#ifdef par
        allocate ( brhs(nRowLocal,NRHS) )
#else
        allocate ( brhs(nPts_tot*3,NRHS) )
#endif
        brhs        = 0
        brhs_global = 0

    end subroutine alloc_total_brhs


    subroutine init_brhs ( g, rhs )

        use aorsaNamelist, &
        only: rAnt, zAnt, npRow, npCol, &
            antSigX, antSigY, &
            metalLeft, metalRight, metalTop, metalBot, &
            useEqdsk, r0, rhoAnt, antSigRho, &
            AntennaJ_R, AntennaJ_T, AntennaJ_Z, &
            antSigUnit, useAR2SourceLocationsFile, &
            useAllRHSsSource
        use grid
        use parallel
        use profiles, only: omgrf
        use constants
        use scalapack_mod, only: IndxL2G
        use ar2SourceLocationsInput, &
            only: CurrentSource_r, CurrentSource_z, &
            CurrentSource_ComponentID

        implicit none

        type(gridBlock), intent(inout) :: g
        integer, intent(in) :: rhs

        integer :: i, j, iRow, iCol
        complex :: TmpAntJ

        ! scalapack index variables

        integer :: pr_sp, pc_sp, l_sp, m_sp, x_sp, y_sp
        integer :: localRow, localCol, ii, jj


        if(.not.allocated(g%jR))allocate ( &
            g%jR(g%nR,g%nZ), &
            g%jT(g%nR,g%nZ), &
            g%jZ(g%nR,g%nZ) )

        g%jR = 0
        g%jT = 0
        g%jZ = 0
 
        !   note curden is in Amps per meter of toroidal length (2.*pi*rt).

        do i = 1, g%nR
            do j = 1, g%nZ

                if(useAR2SourceLocationsFile)then

                    TmpAntJ = 1*exp ( -( &
                        (g%R(i)-CurrentSource_r(rhs))**2/antSigUnit**2 &
                            + (g%Z(j)-CurrentSource_z(rhs))**2/antSigUnit**2 ) )

                    ! ID mapping to components:
                    !   0 -> r_re
                    !   1 -> t_re
                    !   2 -> z_re
                    !   3 -> r_im
                    !   4 -> t_im
                    !   5 -> z_im

                    if(CurrentSource_ComponentID(rhs).eq.0) g%jR(i,j) = TmpAntJ+zi*TmpAntJ
                    if(CurrentSource_ComponentID(rhs).eq.1) g%jT(i,j) = TmpAntJ+zi*TmpAntJ
                    if(CurrentSource_ComponentID(rhs).eq.2) g%jZ(i,j) = TmpAntJ+zi*TmpAntJ

                    !if(CurrentSource_ComponentID(rhs).eq.3) g%jR(i,j) = zi*TmpAntJ
                    !if(CurrentSource_ComponentID(rhs).eq.4) g%jT(i,j) = zi*TmpAntJ
                    !if(CurrentSource_ComponentID(rhs).eq.5) g%jZ(i,j) = zi*TmpAntJ

                endif 

                if(useAllRHSsSource)then

                    TmpAntJ = 50*exp ( -( (g%R(i)-rAnt)**2/antSigX**2 + (g%Z(j)-zAnt)**2/antSigY**2 ) )

                    if (AntennaJ_R) g%jR(i,j) = TmpAntJ
                    if (AntennaJ_T) g%jT(i,j) = TmpAntJ
                    if (AntennaJ_Z) g%jZ(i,j) = TmpAntJ

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


        do i = 1, g%nR
            do j = 1, g%nZ

                iRow = (j-1) * 3 + (i-1) * g%nZ * 3 + 1
                iRow = iRow + ( g%startRow-1 )

                do ii = 0, 2

#ifdef par
                    ! now below in it's own section  
#else
                    if (ii==0) brhs(iRow+0,rhs)    = -zi*omgrf*mu0*g%jR(i,j)
                    if (ii==1) brhs(iRow+1,rhs)    = -zi*omgrf*mu0*g%jT(i,j)
                    if (ii==2) brhs(iRow+2,rhs)    = -zi*omgrf*mu0*g%jZ(i,j)
#endif
                enddo
            enddo
        enddo
       
#ifdef par
        if(MyCol==0)then
            do ii=1,nRowLocal,3
       
                iRow = IndxL2G ( ii, RowBlockSize, MyRow, 0, NpRow )

                i = (iRow-1)/(3*g%nZ)+1
                j = (mod(iRow-1,3*g%nZ)+1-1)/3+1
          
                brhs(ii+0,rhs) =  -zi*omgrf*mu0*g%jR(i,j)
                brhs(ii+1,rhs) =  -zi*omgrf*mu0*g%jT(i,j)
                brhs(ii+2,rhs) =  -zi*omgrf*mu0*g%jZ(i,j)

            enddo
        endif
#endif

    end subroutine

end module antenna
