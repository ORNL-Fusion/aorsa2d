module antenna

use constants

implicit none

real :: antSigX, antSigY

#ifndef dblprec
    complex, allocatable :: brhs(:,:)
#else
    complex(kind=dbl), allocatable :: brhs(:,:)
#endif
complex, allocatable :: brhs_global(:,:)

contains

    subroutine alloc_total_brhs ( nPts_tot )

        use parallel
        use AR2SourceLocationsInput, NRHS_FromInputFile=>NRHS

        implicit none

        integer, intent(in) :: nPts_tot

        if(NRHS<1)then
                write(*,*) 'ERROR: NRHS<1'
                stop
        endif

#ifdef par
        allocate ( brhs(LM_B,LN_B) )
        allocate ( brhs_global(M_B,N_B) )
#else
        allocate ( brhs(nPts_tot*3,NRHS) )
        allocate ( brhs_global(nPts_tot*3,NRHS) )
#endif
        brhs        = 0
        brhs_global = 0

    end subroutine alloc_total_brhs


    function get_jA (g, i, j, rhs)

        use aorsaNamelist, &
        only: rAnt, zAnt, npRow, npCol, &
            antSigX, antSigY, &
            metalLeft, metalRight, metalTop, metalBot, &
            useEqdsk, r0, rhoAnt, antSigRho, &
            AntennaJ_R, AntennaJ_T, AntennaJ_Z, &
            antSigUnit, useAR2SourceLocationsFile, &
            useAllRHSsSource, ZeroJp, &
            ZeroJp_rMin, ZeroJp_rMax, ZeroJp_zMin, ZeroJp_zMax, &
            useAntennaFromAR2Input
        use grid
        use profiles, only: omgrf
        use constants
        use ar2SourceLocationsInput, &
            only: CurrentSource_r, CurrentSource_z, &
            CurrentSource_ComponentID, &
            CS_RealFac, CS_ImagFac

        use ar2Input, only: &
            ar2_ja=>jant, &
            ar2_ja_r=>jant_r, &
            ar2_ja_t=>jant_t, &
            ar2_ja_z=>jant_z

        implicit none

        type(gridBlock), intent(in) :: g
        integer, intent(in) :: rhs, i, j
        complex, dimension(3) :: get_jA
        real :: TmpAntJ
        logical :: DoThisPoint

        get_jA = 0

        if(useAR2SourceLocationsFile)then
        
            TmpAntJ = 1*exp ( -( &
                (g%R(i)-CurrentSource_r(rhs))**2/antSigX**2 &
                    + (g%Z(j)-CurrentSource_z(rhs))**2/antSigY**2 ) )

            ! ID mapping to components:
            !   0 -> r
            !   1 -> t
            !   2 -> z

            DoThisPoint = .true.
            !if(ZeroJp)then
            !        DoThisPoint = .false.
            !        if(g%R(i)>=ZeroJp_rMin &
            !                .and.g%R(i)<=ZeroJp_rMax &
            !                .and.g%z(j)>=ZeroJp_zMin &
            !                .and.g%z(j)<=ZeroJp_zMax) DoThisPoint = .true.
            !endif
      
            if(DoThisPoint)then 
                if(CurrentSource_ComponentID(rhs).eq.0) &
                        get_jA(1) = CS_RealFac(rhs)*TmpAntJ + CS_ImagFac(rhs)*zi*TmpAntJ
                if(CurrentSource_ComponentID(rhs).eq.1) &
                        get_jA(2) = CS_RealFac(rhs)*TmpAntJ + CS_ImagFac(rhs)*zi*TmpAntJ
                if(CurrentSource_ComponentID(rhs).eq.2) &
                        get_jA(3) = CS_RealFac(rhs)*TmpAntJ + CS_ImagFac(rhs)*zi*TmpAntJ
            endif
        
        endif 
        
        if(useAllRHSsSource)then ! This will ADD to the multiple RHS sources if enabled
        
            TmpAntJ = 1*exp ( -( (g%R(i)-rAnt)**2/antSigX**2 + (g%Z(j)-zAnt)**2/antSigY**2 ) )
        
            if (AntennaJ_R) get_jA(1) = get_jA(1) + cmplx(TmpAntJ,TmpAntJ*0)
            if (AntennaJ_T) get_jA(2) = get_jA(2) + cmplx(TmpAntJ,TmpAntJ*0)
            if (AntennaJ_Z) get_jA(3) = get_jA(3) + cmplx(TmpAntJ,TmpAntJ*0)
        
        endif

        if(useAntennaFromAR2Input)then 
                ! nothing to do since these data are read in with the profile read :) 
        endif
        
        !   boundary conditions
        !   -------------------
        
        if(g%nR/=1) then
        if ( i==1 .or. i==g%nR ) then
            get_jA  = 0
        endif
        endif
        
        if(g%nZ/=1) then
        if ( j==1 .or. j==g%nZ ) then
            get_jA = 0
        endif
        endif
        
        !if ( capR(i) < metalLeft .or. capR(i) > metalRight &
        !        .or. y(j) > metalTop .or. y(j) < metalBot ) then 
        !    jR(i,j)    = 0
        !    jT(i,j)    = 0
        !    jZ(i,j)    = 0
        !endif

    end function get_jA

    subroutine fill_jA ( g, rhs )

        use grid
        use aorsaNamelist, only: useAntennaFromAR2Input

        implicit none

        type(gridBlock), intent(inout) :: g
        integer, intent(in) :: rhs
        integer :: i, j
        complex :: This_jA(3)

        if(.not.allocated(g%jR))allocate ( &
            g%jR(g%nR,g%nZ), &
            g%jT(g%nR,g%nZ), &
            g%jZ(g%nR,g%nZ) )

        ! overwrite jA from ar2Input.nc if required 
        if(.not.useAntennaFromAR2Input)then 

            g%jR = 0
            g%jT = 0
            g%jZ = 0

            do i = 1,g%nR
                do j = 1,g%nZ

                    This_jA = get_jA( g, i, j, rhs )

                    g%jR(i,j) = This_jA(1)
                    g%jT(i,j) = This_jA(2)
                    g%jZ(i,j) = This_jA(3)

                enddo
            enddo
        endif
       
    end subroutine fill_jA


    subroutine init_brhs ( g, nPts_tot )

        use aorsaNamelist, &
        only: rAnt, zAnt, npRow, npCol, &
            antSigX, antSigY, &
            metalLeft, metalRight, metalTop, metalBot, &
            useEqdsk, r0, rhoAnt, antSigRho, &
            AntennaJ_R, AntennaJ_T, AntennaJ_Z, &
            antSigUnit, useAR2SourceLocationsFile, &
            useAllRHSsSource, useAntennaFromAR2Input
        use grid
        use parallel
        use profiles, only: omgrf
        use constants
        use scalapack_mod
        use ar2SourceLocationsInput, &
            only: CurrentSource_r, CurrentSource_z, &
            CurrentSource_ComponentID

        implicit none

        type(gridBlock), intent(inout) :: g

        integer :: i, j, iRow, iCol, rhs
        integer, intent(in) :: nPts_tot
        complex :: This_jA(3)

        ! scalapack index variables

        integer :: pr_sp, pc_sp, l_sp, m_sp, x_sp, y_sp
        integer :: localRow, localCol, ii, jj

        if(.not.allocated(g%jR))allocate ( &
            g%jR(g%nR,g%nZ), &
            g%jT(g%nR,g%nZ), &
            g%jZ(g%nR,g%nZ) )

        if(.not.useAntennaFromAR2Input)then
            g%jR = 0
            g%jT = 0
            g%jZ = 0
        endif
 
#ifndef par
        do ii = 1,nPts_tot*3,3
            do rhs = 1,NRHS 

                !iRow = (j-1) * 3 + (i-1) * g%nZ * 3 + 1
                !iRow = iRow + ( g%startRow-1 )

                i = (ii-1)/(3*g%nZ)+1
                j = (mod(ii-1,3*g%nZ)+1-1)/3+1

                if(.not.useAntennaFromAR2Input)then  
                    This_jA = get_jA( g, i, j, rhs )
                else
                    This_jA = (/g%jR(i,j),g%jT(i,j),g%jZ(i,j)/)
                endif

                brhs(ii+0,rhs) = -zi*omgrf*mu0*This_jA(1)
                brhs(ii+1,rhs) = -zi*omgrf*mu0*This_jA(2)
                brhs(ii+2,rhs) = -zi*omgrf*mu0*This_jA(3)

            enddo
        enddo
#endif
       
#ifdef par
        do ii=1,LM_B,3
            do jj=1,LN_B
        
                iRow = IndxL2G ( ii, desc_B(MB_), myRow, 0, NpRow )
                rhs = IndxL2G ( jj, desc_B(NB_), myCol, 0, NpCol )

                i = (iRow-1)/(3*g%nZ)+1
                j = (mod(iRow-1,3*g%nZ)+1-1)/3+1

                This_jA = get_jA( g, i, j, rhs )

                brhs(ii+0,jj) =  -zi*omgrf*mu0*This_jA(1)
                brhs(ii+1,jj) =  -zi*omgrf*mu0*This_jA(2)
                brhs(ii+2,jj) =  -zi*omgrf*mu0*This_jA(3)

!#ifdef __GFORTRAN__
!                if(ISNAN(REAL(brhs(ii+0,jj)))) stop '"re brhs(ii+0,jj)" is a NaN'
!                if(ISNAN(REAL(brhs(ii+1,jj)))) stop '"re brhs(ii+1,jj)" is a NaN'
!                if(ISNAN(REAL(brhs(ii+2,jj)))) stop '"re brhs(ii+2,jj)" is a NaN'
!
!                if(ISNAN(AIMAG(brhs(ii+0,jj)))) stop '"im brhs(ii+0,jj)" is a NaN'
!                if(ISNAN(AIMAG(brhs(ii+1,jj)))) stop '"im brhs(ii+1,jj)" is a NaN'
!                if(ISNAN(AIMAG(brhs(ii+2,jj)))) stop '"im brhs(ii+2,jj)" is a NaN'
!#endif

            enddo
        enddo
#endif

    end subroutine

end module antenna
