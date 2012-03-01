module mat_fill

use constants

implicit none

complex :: bFn, dRBFn, dZBFn, d2RbFn, d2ZbFn, &
    bFnHere, d3RbFn, d4RbFn
complex(kind=dbl) :: &
    sigAlpAlp, sigAlpBet, sigAlpPrl, &
    sigBetAlp, sigBetBet, sigBetPrl, &
    sigPrlAlp, sigPrlBet, sigPrlPrl
complex(kind=dbl) :: &
    sigAlpAlpTmp, sigAlpBetTmp, sigAlpPrlTmp, &
    sigBetAlpTmp, sigBetBetTmp, sigBetPrlTmp, &
    sigPrlAlpTmp, sigPrlBetTmp, sigPrlPrlTmp


complex :: &
    kAlpAlp, kAlpBet, kAlpPrl, &
    kBetAlp, kBetBet, kBetPrl, &
    kPrlAlp, kPrlBet, kPrlPrl
real :: rnx, rny, rnPhi 
complex :: &
    fdk, fek, ffk, &
    fgk, fak, fpk, &
    frk, fqk, fsk
complex, allocatable, dimension(:) :: &
    sss, ttt, qqq

#ifndef dblprec
    complex, allocatable :: aMat(:,:)
#else
    complex(kind=dbl), allocatable :: aMat(:,:)
#endif

type :: workListEntry
        integer :: i,j,m,n
end type workListEntry


contains

    subroutine alloc_total_aMat ( nPts_tot )

        use parallel

        implicit none

        integer, intent(in) :: nPts_tot

        write(*,*) '    nPts_tot: ', nPts_tot

#ifdef par

        if (iAm == 0) &
        write(*,100), &
            nPts_tot*3*nPts_tot*3*2*8.0 / 1024.0**2, &
            nRowLocal*nColLocal*2*8.0 / 1024.0**2
        100 format (' Filling aMat [global size: ',f8.1,' MB, local size: ',f8.1' MB]')

        allocate ( aMat(nRowLocal,nColLocal) )
#else 
        write(*,100), &
            nPts_tot*3*nPts_tot*3*2*8.0 / 1024.0**2
        100 format (' Filling aMat [global size: ',f8.1,' MB]')

        allocate ( aMat(nPts_tot*3,nPts_tot*3) )
        write(*,*) '    ',nPts_tot*3,nPts_tot*3  

#endif

        aMat = 0

    end subroutine alloc_total_aMat


    subroutine amat_boundaries ( gAll, nPts_tot )

        use grid
        use antenna
        use aorsa2din_mod, &
        only: overlap
        use spline_dlg

        implicit none

        type(gridBlock), intent(in) :: gAll(:)
        integer, intent(in) :: nPts_tot

        type(gridBlock) :: me, nbr
        integer :: i, j, n, m, ii, iRow, iCol 
        complex(kind=dbl) :: aMatBlock(3,3)
        real :: r, kt
        complex :: bFn_iL, bFn_iR, bFn_iRR
        complex :: bFn_i0, bFn_i1, bFn_i2, bFn_i3, bFn_i4, bFn_i5
        complex :: bFn_R, bFn_Z, bFn_0
        real :: rNorm, zNorm


        real :: h1, h2, alpha, coeffL, coeffR, splines(3,4)
        integer :: iiL, iiR, iiRR, iiArr(1)

        do ii=1,nPts_tot

            if(bndryBlockID(3,ii)/=0)then

                write(*,*) 'BndryType: ', bndryType(ii)
                write(*,*) 'BndryBlockID: ', bndryBlockID(:,ii)

                me = gAll(bndryBlockID(3,ii))

                do n=me%nMin,me%nMax
                    do m=me%mMin,me%mMax

                        aMatBlock = 0

                        if(bndryType(ii)==-999) then

                            ! create E=E @ overlap point

                            i = me%nR-overlap
                            j = bndryBlockID(2,ii)

                            bFn = me%xx(n, i) * me%yy(m, j)

                            aMatBlock(1,1) = -bFn
                            aMatBlock(2,2) = -bFn
                            aMatBlock(3,3) = -bFn

                            !aMatBlock = -99

                        elseif(bndryType(ii)<=-1 .and. bndryType(ii)>=-20) then

                            ! Dirichlet 
                            ! ---------

                            nbr = gAll(bndryBlockID(4,ii))
                            i = bndryBlockID(1,ii)
                            j = bndryBlockID(2,ii)

                            ! linear interpolation
                            ! --------------------

                            ! determine left index for linear interpolation
                            iiArr = maxLoc ( me%R-nbr%R(i), (me%R-nbr%R(i))<0)
                            iiL = iiArr(1)
                            iiR = iiL + 1

                            coeffL = 1 - nbr%R(i) / ( -me%R(iiL)+me%R(iiR) ) &
                                + me%R(iiL) / ( -me%R(iiL)+me%R(iiR) )
                            coeffR = nbr%R(i) / ( -me%R(iiL)+me%R(iiR) ) &
                                - me%R(iiL) / ( -me%R(iiL)+me%R(iiR) )

                            bFn_iL =me%xx(n, iiL) * me%yy(m, j)
                            bFn_iR =me%xx(n, iiR) * me%yy(m, j)

                            aMatBlock(1,1) = -(coeffL*bFn_iL + coeffR*bFn_iR)
                            aMatBlock(2,2) = -(coeffL*bFn_iL + coeffR*bFn_iR)
                            aMatBlock(3,3) = -(coeffL*bFn_iL + coeffR*bFn_iR)

                            !!write(*,*) n, aMatBlock(1,1)
                            !! clamped cubic spline interpolation
                            !! ----------------------------------

                            !! need 3 pts 
                            !iiRR = iiR + 1
                            !bFn_iRR =me%xx(n, iiRR) * me%yy(m, j)
                            !splines = spline ( me%R(iiL:iiRR), &
                            !    realPart((/ bFn_iL, bFn_iR, bFn_iRR /)), &
                            !    realPart(me%dRbfn_bFn(n,iiL)*bFn_iL), &
                            !    realPart(me%dRbFn_bFn(n,iiRR)*bFn_iRR) )

                            !aMatBlock(1,1) = -( &
                            !    splines(1,1) &
                            !    + splines(1,2) * ( nbr%R(i) - me%R(iiL) ) &
                            !    + splines(1,3) * ( nbr%R(i) - me%R(iiL) )**2 &
                            !    + splines(1,4) * ( nbr%R(i) - me%R(iiL) )**3 )

                            !aMatBlock(2,2) = aMatBlock(1,1)
                            !aMatBlock(3,3) = aMatBlock(1,1)
                            !!write(*,*) n, aMatBlock(1,1)


                            !! direct fn evaluation
                            !! --------------------

                            !if(chebyshevX)then
                            !    rNorm = (nbr%R(i)-me%rMin)/me%rRange*2-1
                            !    if(nbr%nR==1) rNorm = nbR%rNorm(1)
                            !    bFn_R = xBasis(n,rNorm)
                            !else
                            !    rNorm = (nbr%R(i)-me%rMin)/me%rRange*2*pi
                            !    if(nbr%nR==1) rNorm = nbR%rNorm(1)
                            !    bFn_R = xBasis(n,rNorm)
                            !endif

                            !if(chebyshevY)then
                            !    zNorm = (nbr%z(j)-me%zMin)/me%zRange*2-1
                            !    if(nbr%nZ==1) zNorm = nbr%zNorm(1)
                            !    bFn_Z = yBasis(m,zNorm)
                            !else
                            !    zNorm = (nbr%z(j)-me%zMin)/me%zRange*2*pi
                            !    if(nbr%nZ==1) zNorm = nbr%zNorm(1)
                            !    bFn_Z = yBasis(m,zNorm)
                            !endif

                            !bFn_0 = bFn_R * bFn_Z 

                            !aMatBlock(1,1) = -(bFn_0)
                            !aMatBlock(2,2) = -(bFn_0)
                            !aMatBlock(3,3) = -(bFn_0)

                            !!write(*,*) n, aMatBlock(1,1)
                            !!write(*,*)

                            !! Match single (or subset) basis functions
                            !! ----------------------------------------

                            !if(mod(abs(bndryType(ii)),2)==1) then 
                            !    i = me%nR-overlap
                            !else
                            !    i = 1+overlap
                            !endif
                            !j = bndryBlockID(2,ii)

                            !if(n<=me%nR/(abs(bndryType(ii))+1))then
                            !!if(n==abs(bndryType(ii))+1)then
                            !    bFn = me%xx(n, i) * me%yy(m, j)
                            !else
                            !    bFn = 0
                            !endif

                            !aMatBlock(1,1) = -(bFn)
                            !aMatBlock(2,2) = -(bFn)
                            !aMatBlock(3,3) = -(bFn)

                        endif

                        ! but couple with the neighbour block

                        i = bndryBlockID(1,ii)
                        j = bndryBlockID(2,ii)

                        nbr = gAll(bndryBlockID(4,ii))

                        iRow = (i-1) * 3 * nbr%nZ + (j-1) * 3 + 1
                        iRow = iRow + ( nbr%startRow-1 )

                        iCol = (n-me%nMin) * 3 * me%nModesZ + (m-me%mMin) * 3 + 1
                        iCol = iCol + ( me%startCol-1 )

                        aMat(iRow:iRow+2,iCol:iCol+2) = aMatBlock
                        brhs(iRow:iRow+2) = 0


                    enddo
                enddo

            endif

        enddo

    end subroutine amat_boundaries



    subroutine amat_fill ( g )

        use aorsa2din_mod, &
        only: &
            delta0, nSpec, &
            iSigma, npRow, npCol, &
            metalLeft, metalRight, &
            metalTop, metalBot, nPhi, square, lsWeightFac, &
            useEqdsk, overlap
        use grid
        use sigma_mod
        use rotation
        use profiles
        use bField
        use parallel
        use chebyshev_mod
        use write_data

        implicit none

        type(gridBlock), intent(in) :: g

        type(dBfnArg) :: d
        integer :: iRow, iCol, i, j, n, m, p, s, ii, jj, iOL, jOL
        integer :: localRow, localCol
        complex :: metal
        logical, allocatable :: isMetal(:,:)
        real :: kr, kt, kz, r, z, kVec_stix(3)

        complex :: &
            mat_r_alp, mat_r_bet, mat_r_prl, &
            mat_th_alp, mat_th_bet, mat_th_prl, &
            mat_z_alp, mat_z_bet, mat_z_prl

        complex :: bFn_iL, bFn_iR, bFn_i0, bFn_i1, bFn_i2, bFn_i3, bFn_i4, bFn_i5 
        real :: h1, h2, alpha

        complex(kind=dbl) :: sigma_tmp(3,3), sigma_tmp_neg(3,3)
        complex, allocatable :: sigma_write(:,:,:,:,:,:,:)

        integer :: nr, nz, iStat
        integer :: sigma_nc_id, sigma_re_id, sigma_im_id

        real :: &
            Urr, Urt, Urz, &
            Utr, Utt, Utz, &
            Uzr, Uzt, Uzz

        real :: &
            drUrr, drUrt, drUrz, &
            drUtr, drUtt, drUtz, &
            drUzr, drUzt, drUzz
        
        real :: &
            dzUrr, dzUrt, dzUrz, &
            dzUtr, dzUtt, dzUtz, &
            dzUzr, dzUzt, dzUzz
        
        real :: &
            drrUrr, drrUrt, drrUrz, &
            drrUtr, drrUtt, drrUtz, &
            drrUzr, drrUzt, drrUzz
        
        real :: &
            dzzUrr, dzzUrt, dzzUrz, &
            dzzUtr, dzzUtt, dzzUtz, &
            dzzUzr, dzzUzt, dzzUzz
        
        real :: &
            drzUrr, drzUrt, drzUrz, &
            drzUtr, drzUtt, drzUtz, &
            drzUzr, drzUzt, drzUzz

#ifdef par
        !   scalapack indicies
        !   see http://www.netlib.org/scalapack/slug/node76.html

        integer :: l_sp, m_sp, pr_sp, pc_sp, x_sp, y_sp
        integer :: pr_sp_thisPt(3), pc_sp_thisPt(3)
#endif

        integer :: workListPosition
        type(workListEntry) :: thisWorkList
        type(workListEntry), allocatable :: workList(:), workListTooLong(:)

        metal   = ( 1e8,1e8 )


        ! Set the metal regions
        ! ---------------------

        allocate (isMetal(g%nR,g%nZ))
        isMetal = .false.

        !if(useEqdsk)then
        !where(g%rho>=0.99)
        !        isMetal = .true.
        !endwhere
        !endif

        do i=1,g%nR
            do j=1,g%nZ
                if ( g%R(i) < metalLeft .or. g%R(i) > metalRight &
                        .or. g%Z(j) > metalTop .or. g%Z(j) < metalBot ) &
                    isMetal(i,j) = .true.
            enddo
        enddo

        if(square) lsWeightFac = 1


        ! Initialize grid sigma file
        ! --------------------------

        if(iAm==0) &
        call init_sigma_file ( g, 'sigma'//g%fNumber//'.nc', &
            sigma_nc_id, sigma_re_id, sigma_im_id, nSpec )

        allocate ( sigma_write(g%nR,g%nZ,g%nMin:g%nMax,g%mMin:g%mMax,3,3,nSpec), stat = iStat )
        if(iStat/=0)then
                write(*,*) 'ERROR src/mat_fill.f90 - allocation failed :('
                stop
        endif

        sigma_write = 0

        ! Calculate sigma seperately outside main loop
        ! --------------------------------------------

        do s=1,nSpec

        if(iAm==0) &
        write(*,*) 'Calculating sigma for species ', s, ' of', nSpec

            do m=g%mMin,g%mMax
                do n=g%nMin,g%nMax
#ifndef par
                !   progress indicator
                !   ------------------
                do p=1,7 
                    write(*,'(a)',advance='no') char(8)
                enddo
                write(*,'(1x,f5.1,a)',advance='no') &
                    real((m-g%mMin)*g%nR+(n-g%nMin))/(g%nR*g%nZ)*100, '%'
#endif
                    do j=1,g%nZ
                        do i=1,g%nR

                            if(chebyshevX) then
                                if(n>1) then
                                    kr = n / sqrt ( sin ( pi * (g%rNorm(i)+1)/2  ) ) * g%normFacR 
                                else
                                    kr = n * g%normFacR
                                endif
                            else
                                kr = n * g%normFacR
                            endif

                            if(chebyshevY) then
                                if(m>1) then
                                    kz = m / sqrt ( sin ( pi * (g%zNorm(j)+1)/2 ) ) * g%normFacZ 
                                else
                                    kz = m * g%normFacZ
                                endif
                            else
                                kz = m * g%normFacZ
                            endif

  
                            hotPlasma:& 
                            if (iSigma==1 .and. (.not. isMetal(i,j)) ) then        

                                kVec_stix = matMul( g%U_RTZ_to_ABb(i,j,:,:), (/ kr, g%kPhi(i), kz /) ) 

                                sigma_tmp = sigmaHot_maxwellian&
                                    ( mSpec(s), &
                                    g%ktSpec(i,j,s), g%omgc(i,j,s), g%omgp2(i,j,s), &
                                    kVec_stix, g%R(i), &
                                    omgrf, k0, &
                                    g%k_cutoff, s, &
                                    g%sinTh(i,j), g%bPol(i,j), g%bMag(i,j), g%gradPrlB(i,j), &
                                    g%nuOmg(i,j) )

                                !if(cosX)then
                                !kVec_stix = matMul( g%U_RTZ_to_ABb(i,j,:,:), (/ -kr, g%kPhi(i), kz /) ) 

                                !sigma_tmp_neg = sigmaHot_maxwellian&
                                !    ( mSpec(s), &
                                !    g%ktSpec(i,j,s), g%omgc(i,j,s), g%omgp2(i,j,s), &
                                !    kVec_stix, g%R(i), &
                                !    omgrf, k0, &
                                !    g%k_cutoff, s, &
                                !    g%sinTh(i,j), g%bPol(i,j), g%bMag(i,j), g%gradPrlB(i,j), &
                                !    g%nuOmg(i,j) )

                                !sigma_tmp = ( sigma_tmp + sigma_tmp_neg ) / 2
                                !endif

                            endif hotPlasma


                            coldPlasma: &
                            if (iSigma==0 .and. (.not. isMetal(i,j)) ) then 

                                sigma_tmp = sigmaCold_stix &
                                    ( g%omgc(i,j,s), g%omgp2(i,j,s), omgrf, &
                                    g%nuOmg(i,j) )

                            endif coldPlasma

                            sigma_write(i,j,n,m,:,:,s) = sigma_tmp

                        enddo
                    enddo

                enddo
            enddo

        enddo
#ifndef par
        write(*,*)
#endif

        ! Write sigma
        ! -----------
        if(iAm==0) &
        write(*,*) 'Writing sigma ...'
        if(iAm==0) &
        call write_sigma_pt ( sigma_write, &
            sigma_nc_id, sigma_re_id, sigma_im_id )
        if(iAm==0) &
        write(*,*) 'DONE'

        ! Close sigma file
        ! ----------------

        if(iAm==0) &
        call close_sigma_file ( sigma_nc_id )


        ! Create a work list for this processor
        ! -------------------------------------

        allocate(workListTooLong(nRowLocal*nColLocal))
        workListPosition = 0

        i_workList: &
        do i=1,g%nR
            j_workList: &
            do j=1,g%nZ

                iRow = (i-1) * 3 * g%nZ + (j-1) * 3 + 1
                iRow = iRow + ( g%startRow-1 )

                n_workList: &
                do n=g%nMin,g%nMax
                    m_workList: &
                    do m=g%mMin,g%mMax

                        iCol = (n-g%nMin) * 3 * g%nModesZ + (m-g%mMin) * 3 + 1
                        iCol = iCol + ( g%startCol-1 )
#ifdef par
                        pr_sp_thisPt   = mod ( rowStartProc + (iRow-1+(/0,1,2/))/rowBlockSize, npRow )
                        pc_sp_thisPt   = mod ( colStartProc + (iCol-1+(/0,1,2/))/colBlockSize, npCol )

                        addToWorkList: &
                        if ( any(pr_sp_thisPt==myRow) .and. any(pc_sp_thisPt==myCol) ) then
#endif
                            workListPosition = workListPosition + 1
                            if(workListPosition > size(workListTooLong)) &
                            write(*,*) iAm, workListPosition, size(workListTooLong)
                            thisWorkList = workListEntry(i,j,m,n)
                            workListTooLong(workListPosition) = thisWorkList
#ifdef par
                        endif addToWorkList
#endif
                    enddo m_workList
                enddo n_workList

            enddo j_workList
        enddo i_workList

        allocate(workList(workListPosition))
        workList = workListTooLong(1:workListPosition)
        deallocate(workListTooLong)



        ! Begin loop
        ! ----------

        i_loop: &
        do i=1,g%nR

#ifndef par
            !   progress indicator
            !   ------------------
            do p=1,7 
                write(*,'(a)',advance='no') char(8)
            enddo
            write(*,'(1x,f5.1,a)',advance='no') real(i)/g%nR*100, '%'
#endif
            j_loop: &
            do j=1,g%nZ

                iRow = (i-1) * 3 * g%nZ + (j-1) * 3 + 1
                iRow = iRow + ( g%startRow-1 )

                n_loop: &
                do n=g%nMin,g%nMax
                    m_loop: &
                    do m=g%mMin,g%mMax

                        iCol = (n-g%nMin) * 3 * g%nModesZ + (m-g%mMin) * 3 + 1
                        iCol = iCol + ( g%startCol-1 )
#ifdef par
                        pr_sp_thisPt   = mod ( rowStartProc + (iRow-1+(/0,1,2/))/rowBlockSize, npRow )
                        pc_sp_thisPt   = mod ( colStartProc + (iCol-1+(/0,1,2/))/colBlockSize, npCol )

                        doWork: &
                        if ( any(pr_sp_thisPt==myRow) .and. any(pc_sp_thisPt==myCol) ) then
#endif
                            sigAlpAlp = 0.0
                            sigAlpBet = 0.0
                            sigAlpPrl = 0.0
                            sigBetAlp = 0.0
                            sigBetBet = 0.0
                            sigBetPrl = 0.0
                            sigPrlAlp = 0.0
                            sigPrlBet = 0.0
                            sigPrlPrl = 0.0

                            z   = g%z(j)
                            r   = g%R(i)
                            kt  = nPhi!g%kPhi(i)

                            bFn = g%xx(n, i) * g%yy(m, j)

        interior: &
        if(g%label(i,j)==0)then


                            !   interior plasma region:
                            !   ----------------------

                            ! The chebyshev k will be infinite at the
                            ! boundaries. However, these should never be
                            ! calculated anyway. I am still not sure where
                            ! the sqrt comes from in the below expressions
                            ! for k? If you figure it out let me know. Also
                            ! not sure about the n,m == 1 values for k. For
                            ! n,m == 0 we have simply k = 0 but n,m == 1 the
                            ! chebT is a straight line. For now leaving it
                            ! as the Fourier equiv.

                            if(chebyshevX) then
                                if(n>1) then
                                    kr = n / sqrt ( sin ( pi * (g%rNorm(i)+1)/2  ) ) * g%normFacR 
                                else
                                    kr = n * g%normFacR
                                endif
                            else
                                kr = n * g%normFacR
                            endif

                            if(chebyshevY) then
                                if(m>1) then
                                    kz = m / sqrt ( sin ( pi * (g%zNorm(j)+1)/2 ) ) * g%normFacZ 
                                else
                                    kz = m * g%normFacZ
                                endif
                            else
                                kz = m * g%normFacZ
                            endif


                            !! Short of k's above the xkPerp_cutOff 
                            !! ------------------------------------

                            !if ( abs(kxSav(n))> kx_cutOff &
                            !    .or. abs(kySav(m))> ky_cutOff ) then 
                            !      
                            !    sigAlpAlp = metal 
                            !    sigAlpBet = 0
                            !    sigAlpPrl = 0 
                            !            
                            !    sigBetAlp = 0 
                            !    sigBetBet = metal 
                            !    sigBetPrl = 0 
                            !            
                            !    sigPrlAlp = 0 
                            !    sigPrlBet = 0 
                            !    sigPrlPrl = metal 

                            !endif

                            ! Sum sigma over species
                            ! ----------------------
                    
                            sigAlpAlp = sum ( sigma_write(i,j,n,m,1,1,:) )
                            sigAlpBet = sum ( sigma_write(i,j,n,m,1,2,:) )
                            sigAlpPrl = sum ( sigma_write(i,j,n,m,1,3,:) )

                            sigBetAlp = sum ( sigma_write(i,j,n,m,2,1,:) )
                            sigBetBet = sum ( sigma_write(i,j,n,m,2,2,:) )
                            sigBetPrl = sum ( sigma_write(i,j,n,m,2,3,:) )

                            sigPrlAlp = sum ( sigma_write(i,j,n,m,3,1,:) )
                            sigPrlBet = sum ( sigma_write(i,j,n,m,3,2,:) )
                            sigPrlPrl = sum ( sigma_write(i,j,n,m,3,3,:) )
 
                            ! Metal
                            ! -----

                            if (isMetal(i,j)) then 

                                sigAlpAlp = metal 
                                sigAlpBet = 0
                                sigAlpPrl = 0 
                                        
                                sigBetAlp = 0 
                                sigBetBet = metal 
                                sigBetPrl = 0 
                                        
                                sigPrlAlp = 0 
                                sigPrlBet = 0 
                                sigPrlPrl = metal 

                            endif

                            kAlpAlp = 1.0 + zi / (eps0 * omgrf) * sigAlpAlp
                            kAlpBet =       zi / (eps0 * omgrf) * sigAlpBet
                            kAlpPrl =       zi / (eps0 * omgrf) * sigAlpPrl

                            kBetAlp =       zi / (eps0 * omgrf) * sigBetAlp
                            kBetBet = 1.0 + zi / (eps0 * omgrf) * sigBetBet
                            kBetPrl =       zi / (eps0 * omgrf) * sigBetPrl

                            kPrlAlp =       zi / (eps0 * omgrf) * sigPrlAlp
                            kPrlBet =       zi / (eps0 * omgrf) * sigPrlBet
                            kPrlPrl = 1.0 + zi / (eps0 * omgrf) * sigPrlPrl

                            d%n = n
                            d%m = m
                            d%xNorm = g%rNorm(i)
                            d%yNorm = g%zNorm(j)
                            d%normFacX = g%normFacR
                            d%normFacY = g%normFacZ

                            Urr = g%Urr(i,j)
                            Urt = g%Urt(i,j)
                            Urz = g%Urz(i,j)

                            Utr = g%Utr(i,j)
                            Utt = g%Utt(i,j)
                            Utz = g%Utz(i,j)

                            Uzr = g%Uzr(i,j)
                            Uzt = g%Uzt(i,j)
                            Uzz = g%Uzz(i,j)


                            drUrr = g%drUrr(i,j)
                            drUrt = g%drUrt(i,j)
                            drUrz = g%drUrz(i,j)

                            drUtr = g%drUtr(i,j)
                            drUtt = g%drUtt(i,j)
                            drUtz = g%drUtz(i,j)

                            drUzr = g%drUzr(i,j)
                            drUzt = g%drUzt(i,j)
                            drUzz = g%drUzz(i,j)

                            dzUrr = g%dzUrr(i,j)
                            dzUrt = g%dzUrt(i,j)
                            dzUrz = g%dzUrz(i,j)

                            dzUtr = g%dzUtr(i,j)
                            dzUtt = g%dzUtt(i,j)
                            dzUtz = g%dzUtz(i,j)

                            dzUzr = g%dzUzr(i,j)
                            dzUzt = g%dzUzt(i,j)
                            dzUzz = g%dzUzz(i,j)


                            drrUrr = g%drrUrr(i,j)
                            drrUrt = g%drrUrt(i,j)
                            drrUrz = g%drrUrz(i,j)
                            
                            drrUtr = g%drrUtr(i,j)
                            drrUtt = g%drrUtt(i,j)
                            drrUtz = g%drrUtz(i,j)

                            drrUzr = g%drrUzr(i,j)
                            drrUzt = g%drrUzt(i,j)
                            drrUzz = g%drrUzz(i,j)

                            dzzUrr = g%dzzUrr(i,j)
                            dzzUrt = g%dzzUrt(i,j)
                            dzzUrz = g%dzzUrz(i,j)
                            
                            dzzUtr = g%dzzUtr(i,j)
                            dzzUtt = g%dzzUtt(i,j)
                            dzzUtz = g%dzzUtz(i,j)

                            dzzUzr = g%dzzUzr(i,j)
                            dzzUzt = g%dzzUzt(i,j)
                            dzzUzz = g%dzzUzz(i,j)

                            drzUrr = g%drzUrr(i,j)
                            drzUrt = g%drzUrt(i,j)
                            drzUrz = g%drzUrz(i,j)
                            
                            drzUtr = g%drzUtr(i,j)
                            drzUtt = g%drzUtt(i,j)
                            drzUtz = g%drzUtz(i,j)

                            drzUzr = g%drzUzr(i,j)
                            drzUzt = g%drzUzt(i,j)
                            drzUzz = g%drzUzz(i,j)


        ! Matrix elements. See mathematica worksheet for calculation of 
        ! these. They are for a general basis set.
        ! This section will fail if you try and fill the boundary pts
        ! using the Chebyshev basis. So don't do it. The boundary pts
        ! are filled later. See below. Runs fine for Fourier basis as
        ! the boundary conds are periodic.
        ! -------------------------------------------------------------

       mat_r_alp = (-((kt**2*Urr)/r**2) + k0**2*KAlpAlp*Urr - (zi*kt*Urt)/r**2 + &
          k0**2*KBetAlp*Utr + k0**2*KPrlAlp*Uzr + &
          (2*dzBfn_bfn(d)*dzUrr) + &
          (Urr*dzzBfn_bfn(d)) + dzzUrr - &
          (zi*kt*Urt*drBfn_bfn(d))/r - &
          (dzUrz*drBfn_bfn(d)) - &
          (zi*kt*drUrt)/r - &
          (dzBfn_bfn(d)*drUrz) - &
          (Urz*drzBfn_bfn(d)) - drzUrz)

       mat_r_bet = (k0**2*KAlpBet*Urr - (kt**2*Utr)/r**2 + k0**2*KBetBet*Utr - &
          (zi*kt*Utt)/r**2 + k0**2*KPrlBet*Uzr + &
          (2*dzBfn_bfn(d)*dzUtr) + &
          (Utr*dzzBfn_bfn(d)) + dzzUtr - &
          (zi*kt*Utt*drBfn_bfn(d))/r - &
          (dzUtz*drBfn_bfn(d)) - &
          (zi*kt*drUtt)/r - &
          (dzBfn_bfn(d)*drUtz) - &
          (Utz*drzBfn_bfn(d)) - drzUtz)

       mat_r_prl = (k0**2*KAlpPrl*Urr + k0**2*KBetPrl*Utr - (kt**2*Uzr)/r**2 + &
          k0**2*KPrlPrl*Uzr - (zi*kt*Uzt)/r**2 + &
          (2*dzBfn_bfn(d)*dzUzr) + &
          (Uzr*dzzBfn_bfn(d)) + dzzUzr - &
          (zi*kt*Uzt*drBfn_bfn(d))/r - &
          (dzUzz*drBfn_bfn(d)) - &
          (zi*kt*drUzt)/r - &
          (dzBfn_bfn(d)*drUzz) - &
          (Uzz*drzBfn_bfn(d)) - drzUzz)

       mat_th_alp = ((zi*kt*Urr)/r**2 - Urt/r**2 + k0**2*KAlpAlp*Urt + &
          k0**2*KBetAlp*Utt + k0**2*KPrlAlp*Uzt - &
          (zi*kt*Urz*dzBfn_bfn(d))/r + &
          (2*dzBfn_bfn(d)*dzUrt) - &
          (zi*kt*dzUrz)/r + (Urt*dzzBfn_bfn(d)) + &
          dzzUrt - (zi*kt*Urr*drBfn_bfn(d))/r + &
          (Urt*drBfn_bfn(d))/r - (zi*kt*drUrr)/r + &
          drUrt/r + (2*drBfn_bfn(d)*drUrt) + &
          (Urt*drrBfn_bfn(d)) + drrUrt)

       mat_th_bet = (k0**2*KAlpBet*Urt + (zi*kt*Utr)/r**2 - Utt/r**2 + &
          k0**2*KBetBet*Utt + k0**2*KPrlBet*Uzt - &
          (zi*kt*Utz*dzBfn_bfn(d))/r + &
          (2*dzBfn_bfn(d)*dzUtt) - &
          (zi*kt*dzUtz)/r + (Utt*dzzBfn_bfn(d)) + &
          dzzUtt - (zi*kt*Utr*drBfn_bfn(d))/r + &
          (Utt*drBfn_bfn(d))/r - (zi*kt*drUtr)/r + &
          drUtt/r + (2*drBfn_bfn(d)*drUtt) + &
          (Utt*drrBfn_bfn(d)) + drrUtt)

       mat_th_prl = (k0**2*KAlpPrl*Urt + k0**2*KBetPrl*Utt + (zi*kt*Uzr)/r**2 - &
          Uzt/r**2 + k0**2*KPrlPrl*Uzt - &
          (zi*kt*Uzz*dzBfn_bfn(d))/r + &
          (2*dzBfn_bfn(d)*dzUzt) - &
          (zi*kt*dzUzz)/r + (Uzt*dzzBfn_bfn(d)) + &
          dzzUzt - (zi*kt*Uzr*drBfn_bfn(d))/r + &
          (Uzt*drBfn_bfn(d))/r - (zi*kt*drUzr)/r + &
          drUzt/r + (2*drBfn_bfn(d)*drUzt) + &
          (Uzt*drrBfn_bfn(d)) + drrUzt)

       mat_z_alp = (-((kt**2*Urz)/r**2) + k0**2*KAlpAlp*Urz + k0**2*KBetAlp*Utz + &
          k0**2*KPrlAlp*Uzz - (Urr*dzBfn_bfn(d))/r - &
          (zi*kt*Urt*dzBfn_bfn(d))/r - dzUrr/r - &
          (zi*kt*dzUrt)/r + (Urz*drBfn_bfn(d))/r - &
          (dzUrr*drBfn_bfn(d)) - &
          (dzBfn_bfn(d)*drUrr) + drUrz/r + &
          (2*drBfn_bfn(d)*drUrz) - &
          (Urr*drzBfn_bfn(d)) - drzUrr + &
          (Urz*drrBfn_bfn(d)) + drrUrz)

       mat_z_bet = (k0**2*KAlpBet*Urz - (kt**2*Utz)/r**2 + k0**2*KBetBet*Utz + &
          k0**2*KPrlBet*Uzz - (Utr*dzBfn_bfn(d))/r - &
          (zi*kt*Utt*dzBfn_bfn(d))/r - dzUtr/r - &
          (zi*kt*dzUtt)/r + (Utz*drBfn_bfn(d))/r - &
          (dzUtr*drBfn_bfn(d)) - &
          (dzBfn_bfn(d)*drUtr) + drUtz/r + &
          (2*drBfn_bfn(d)*drUtz) - &
          (Utr*drzBfn_bfn(d)) - drzUtr + &
          (Utz*drrBfn_bfn(d)) + drrUtz) 

       mat_z_prl = (k0**2*KAlpPrl*Urz + k0**2*KBetPrl*Utz - (kt**2*Uzz)/r**2 + &
          k0**2*KPrlPrl*Uzz - (Uzr*dzBfn_bfn(d))/r - &
          (zi*kt*Uzt*dzBfn_bfn(d))/r - dzUzr/r - &
          (zi*kt*dzUzt)/r + (Uzz*drBfn_bfn(d))/r - &
          (dzUzr*drBfn_bfn(d)) - &
          (dzBfn_bfn(d)*drUzr) + drUzz/r + &
          (2*drBfn_bfn(d)*drUzz) - &
          (Uzr*drzBfn_bfn(d)) - drzUzr + &
          (Uzz*drrBfn_bfn(d)) + drrUzz)

            endif interior

                            ii_loop: &
                            do ii=0,2
                                jj_loop: &
                                do jj=0,2
#ifdef par
                                    pr_sp   = mod ( rowStartProc + (iRow-1+ii)/rowBlockSize, npRow )
                                    pc_sp   = mod ( colStartProc + (iCol-1+jj)/colBlockSize, npCol )

                                    ! scalapack indicies for 2D block cyclic data format
                                    ! see http://www.netlib.org/scalapack/slug/node76.html
                                    !       and
                                    ! http://acts.nersc.gov/scalapack/hands-on/example4.html
 
                                    myProc: &
                                    if ( myRow==pr_sp .and. myCol==pc_sp ) then

                                        l_sp    = ( iRow-1+ii ) / ( npRow * rowBlockSize )
                                        m_sp    = ( iCol-1+jj ) / ( npCol * colBlockSize )

                                        x_sp    = mod ( iRow-1+ii, rowBlockSize ) + 1
                                        y_sp    = mod ( iCol-1+jj, colBlockSize ) + 1

                                        localRow    = l_sp*rowBlockSize+x_sp
                                        localCol    = m_sp*colBlockSize+y_sp
#else
                                        localRow    = iRow+ii
                                        localCol    = iCol+jj
#endif
                                    
                                        ! Interior points
                                        ! ---------------

                                        if(g%label(i,j)==0) then 

                                            if (ii==0 .and. jj==0) aMat(localRow,localCol) = mat_r_alp * bFn  
                                            if (ii==0 .and. jj==1) aMat(localRow,localCol) = mat_r_bet * bFn   
                                            if (ii==0 .and. jj==2) aMat(localRow,localCol) = mat_r_prl * bFn   
                                                                                                  
                                            if (ii==1 .and. jj==0) aMat(localRow,localCol) = mat_th_alp * bFn   
                                            if (ii==1 .and. jj==1) aMat(localRow,localCol) = mat_th_bet * bFn   
                                            if (ii==1 .and. jj==2) aMat(localRow,localCol) = mat_th_prl * bFn   
                                                                                                  
                                            if (ii==2 .and. jj==0) aMat(localRow,localCol) = mat_z_alp * bFn   
                                            if (ii==2 .and. jj==1) aMat(localRow,localCol) = mat_z_bet * bFn   
                                            if (ii==2 .and. jj==2) aMat(localRow,localCol) = mat_z_prl * bFn   

                                        endif


                                        ! Outer boundary points
                                        ! ---------------------

                                        if (g%label(i,j)==888) then
                                            
                                            if (ii==0 .and. jj==0) aMat(localRow,localCol) = bFn * lsWeightFac 
                                            if (ii==0 .and. jj==1) aMat(localRow,localCol) = 0
                                            if (ii==0 .and. jj==2) aMat(localRow,localCol) = 0
                                   
                                            if (ii==1 .and. jj==0) aMat(localRow,localCol) = 0
                                            if (ii==1 .and. jj==1) aMat(localRow,localCol) = bFn * lsWeightFac  
                                            if (ii==1 .and. jj==2) aMat(localRow,localCol) = 0
                                   
                                            if (ii==2 .and. jj==0) aMat(localRow,localCol) = 0
                                            if (ii==2 .and. jj==1) aMat(localRow,localCol) = 0
                                            if (ii==2 .and. jj==2) aMat(localRow,localCol) = bFn * lsWeightFac 

                                        endif

                                        if (g%label(i,j)==999) then

                                            !iOL = 1+overlap
                                            !jOL = j
                                       
                                            !bFnHere = g%xx(n, iOL) * g%yy(m, jOL)

                                            if (ii==0 .and. jj==0) aMat(localRow,localCol) = bFn!Here
                                            if (ii==0 .and. jj==1) aMat(localRow,localCol) = 0  
                                            if (ii==0 .and. jj==2) aMat(localRow,localCol) = 0  
                                   
                                            if (ii==1 .and. jj==0) aMat(localRow,localCol) = 0  
                                            if (ii==1 .and. jj==1) aMat(localRow,localCol) = bFn!Here
                                            if (ii==1 .and. jj==2) aMat(localRow,localCol) = 0   
                                   
                                            if (ii==2 .and. jj==0) aMat(localRow,localCol) = 0   
                                            if (ii==2 .and. jj==1) aMat(localRow,localCol) = 0   
                                            if (ii==2 .and. jj==2) aMat(localRow,localCol) = bFn!Here

                                        endif


                                        ! Mesh-Mesh boundary for R bFn 
                                        ! ----------------------------

                                        if (g%label(i,j)>=1 .and. g%label(i,j)<=20) then

                                            !if(mod(g%label(i,j),2)==1) then 
                                            !    iOL = 1+overlap
                                            !else
                                            !    iOL = g%nR-overlap
                                            !endif
 
                                            !jOL = j
                                      
                                            !if(n<=g%nR/(g%label(i,j)+1))then 
                                            !!if(n==(g%label(i,j)+1))then 
                                            !    bFnHere = g%xx(n, iOL) * g%yy(m, jOL)
                                            !else
                                            !    bFnHere = 0
                                            !endif
                                         
                                            if (ii==0 .and. jj==0) aMat(localRow,localCol) = bFn!Here 
                                            if (ii==0 .and. jj==1) aMat(localRow,localCol) = 0  
                                            if (ii==0 .and. jj==2) aMat(localRow,localCol) = 0  
                                   
                                            if (ii==1 .and. jj==0) aMat(localRow,localCol) = 0  
                                            if (ii==1 .and. jj==1) aMat(localRow,localCol) = bFn!Here  
                                            if (ii==1 .and. jj==2) aMat(localRow,localCol) = 0   
                                   
                                            if (ii==2 .and. jj==0) aMat(localRow,localCol) = 0   
                                            if (ii==2 .and. jj==1) aMat(localRow,localCol) = 0   
                                            if (ii==2 .and. jj==2) aMat(localRow,localCol) = bFn!Here 

                                        endif

#ifdef par
                                    endif myProc
#endif
                                enddo jj_loop
                            enddo ii_loop
#ifdef par
                        endif doWork
#endif
                    enddo m_loop
                enddo n_loop 

            enddo j_loop
        enddo i_loop 


        deallocate ( sigma_write )

        !if(iAm==0) then
        !write(*,*) 
        !do i=1,nRowLocal
        !    write(*,*) real ( aMat(i,:) )
        !enddo
        !endif

        !aMat_   = aMat
    end subroutine amat_fill

end module mat_fill
