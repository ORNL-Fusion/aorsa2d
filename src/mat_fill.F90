module mat_fill

use constants

implicit none

complex :: bFn
complex :: &
    sigAlpAlp, sigAlpBet, sigAlpPrl, &
    sigBetAlp, sigBetBet, sigBetPrl, &
    sigPrlAlp, sigPrlBet, sigPrlPrl
complex :: &
    sigAlpAlpTmp, sigAlpBetTmp, sigAlpPrlTmp, &
    sigBetAlpTmp, sigBetBetTmp, sigBetPrlTmp, &
    sigPrlAlpTmp, sigPrlBetTmp, sigPrlPrlTmp

complex :: sigma_tmp(3,3)

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
    complex(kind=dbl), allocatable :: aMat_(:,:)
#endif

contains

    subroutine amat_fill

        use aorsa2din_mod, &
        only: nModesX, nModesY, &
            delta0, nSpec, &
            iSigma, nPtsX, nPtsY, npRow, npCol, &
            metalLeft, metalRight, &
            metalTop, metalBot, nPhi, square, lsWeightFac, &
            useEqdsk
        use grid
        use sigma_mod
        use rotation
        use profiles
        use bField
        use parallel
        use chebyshev_mod

        implicit none

        integer :: iRow, iCol, i, j, n, m, p, s, ii, jj
        integer :: localRow, localCol
        complex :: metal
        logical, allocatable :: isMetal(:,:)
        real :: kr, kt, kz, r, z
        complex :: &
            dxx, dxy, dxz, &
            dyx, dyy, dyz, &
            dzx, dzy, dzz 

        complex :: &
            mat_r_alp, mat_r_bet, mat_r_prl, &
            mat_th_alp, mat_th_bet, mat_th_prl, &
            mat_z_alp, mat_z_bet, mat_z_prl

        integer :: nr, nz

#ifdef par

        !   scalapack indicies
        !   see http://www.netlib.org/scalapack/slug/node76.html

        integer :: l_sp, m_sp, pr_sp, pc_sp, x_sp, y_sp
        integer :: pr_sp_thisPt(3), pc_sp_thisPt(3)


        if (iAm == 0) &
        write(*,100), &
            nPtsX*nPtsY*3*nModesX*nModesY*3*2*8.0 / 1024.0**2, &
            nRowLocal*nColLocal*2*8.0 / 1024.0**2
        100 format (' Filling aMat [global size: ',f8.1,' MB, local size: ',f8.1' MB]')

        allocate ( aMat(nRowLocal,nColLocal) )
#else 
        write(*,100), &
            nPtsX*nPtsY*3*nModesX*nModesY*3*2*8.0 / 1024.0**2
        100 format (' Filling aMat [global size: ',f8.1,' MB]')

        allocate ( aMat(nPtsX*nPtsY*3,nModesX*nModesY*3) )

#endif

        aMat    = 0
        metal   = ( 1e8,1e8 )


        ! Set the metal regions
        ! ---------------------

        allocate (isMetal(nPtsX,nPtsY))
        isMetal = .false.

        if(useEqdsk)then
        where(rho>=0.99)
                isMetal = .true.
        endwhere
        endif

        do i=1,nPtsX
            do j=1,nPtsY
                if ( capR(i) < metalLeft .or. capR(i) > metalRight &
                        .or. y(j) > metalTop .or. y(j) < metalBot ) &
                    isMetal(i,j) = .true.
            enddo
        enddo

        if(square) lsWeightFac = 1


        ! Begin loop
        ! ----------

        i_loop: &
        do i=1,nPtsX

#ifndef par
            !   progress indicator
            !   ------------------
            do p=1,7 
                write(*,'(a)',advance='no') char(8)
            enddo
            write(*,'(1x,f5.1,a)',advance='no') real(i)/nPtsX*100, '%'
#endif
            j_loop: &
            do j=1,nPtsY

                iRow = (i-1) * 3 * nPtsY + (j-1) * 3 + 1

                n_loop: &
                do n=nMin,nMax
                    m_loop: &
                    do m=mMin,mMax

                        iCol = (n-nMin) * 3 * nModesY + (m-mMin) * 3 + 1
#ifdef par
                        pr_sp_thisPt   = mod ( rowStartProc + (iRow-1+(/0,1,2/))/rowBlockSize, npRow )
                        pc_sp_thisPt   = mod ( colStartProc + (iCol-1+(/0,1,2/))/colBlockSize, npCol )

                        doWork: &
                        if ( any(pr_sp_thisPt==myRow) .and.  any(pc_sp_thisPt==myCol) ) then
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

                            bFn = xx(n, i) * yy(m, j)

        if(  (i/=1 .and. i/=nPtsX .and. j/=1 .and. j/=nPtsY) &
            .or. (nPtsY==1 .and. i/=1 .and. i/=nPtsX) &
            .or. (nPtsX==1 .and. j/=1 .and. j/=nPtsY) ) then


                            !   interior plasma region:
                            !   ----------------------
                            
                            species: &
                            do s=1,nSpec

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
                                        kr = n / sqrt ( sin ( pi * (xGrid_basis(i)+1)/2  ) ) * normFacX 
                                    else
                                        kr = n * normFacX
                                    endif
                                else
                                    kr = n * normFacX
                                endif

                                if(chebyshevY) then
                                    if(m>1) then
                                        kz = m / sqrt ( sin ( pi * (yGrid_basis(j)+1)/2 ) ) * normFacY 
                                    else
                                        kz = m * normFacY
                                    endif
                                else
                                    kz = m * normFacY
                                endif

                                if (iSigma==1 .and. (.not. isMetal(i,j)) ) & ! hot plasma        
                                sigma_tmp = sigmaHot_maxwellian&
                                    ( mSpec(s), &
                                    ktSpec(i,j,s), omgc(i,j,s), omgp2(i,j,s), &
                                    kr, kz, kPhi(i), capr(i), &
                                    omgrf, k0, &
                                    k_cutoff, s, &
                                    bMod(i,j), gradPrlB(i,j), &
                                    U_xyz(i,j,:,:), U_cyl(i,j,:,:), nuOmg2D(i,j) )
                              
                                if (iSigma==0) & ! cold plasma 
                                sigma_tmp = sigmaCold_stix &
                                    ( i, j, omgc(i,j,s), omgp2(i,j,s), omgrf, &
                                    nuOmg2D(i,j) )

                                sigAlpAlp = sigAlpAlp + sigma_tmp(1,1) 
                                sigAlpBet = sigAlpBet + sigma_tmp(1,2) 
                                sigAlpPrl = sigAlpPrl + sigma_tmp(1,3) 
                                               
                                sigBetAlp = sigBetAlp + sigma_tmp(2,1) 
                                sigBetBet = sigBetBet + sigma_tmp(2,2) 
                                sigBetPrl = sigBetPrl + sigma_tmp(2,3) 
                                               
                                sigPrlAlp = sigPrlAlp + sigma_tmp(3,1) 
                                sigPrlBet = sigPrlBet + sigma_tmp(3,2) 
                                sigPrlPrl = sigPrlPrl + sigma_tmp(3,3) 

                            enddo species

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

                            kAlpAlp = 1.0 + zi / (eps0 * omgrf) * sigAlpAlp
                            kAlpBet =       zi / (eps0 * omgrf) * sigAlpBet
                            kAlpPrl =       zi / (eps0 * omgrf) * sigAlpPrl

                            kBetAlp =       zi / (eps0 * omgrf) * sigBetAlp
                            kBetBet = 1.0 + zi / (eps0 * omgrf) * sigBetBet
                            kBetPrl =       zi / (eps0 * omgrf) * sigBetPrl

                            kPrlAlp =       zi / (eps0 * omgrf) * sigPrlAlp
                            kPrlBet =       zi / (eps0 * omgrf) * sigPrlBet
                            kPrlPrl = 1.0 + zi / (eps0 * omgrf) * sigPrlPrl

                            z   = y(j)
                            r   = capR(i)
                            kt  = kPhi(i)

        ! Matrix elements. See mathematica worksheet for calculation of 
        ! these. They are for a general basis set.
        ! This section will fail if you try and fill the boundary pts
        ! using the Chebyshev basis. So don't do it. The boundary pts
        ! are filled later. See below. Runs fine for Fourier basis as
        ! the boundary conds are periodic.
        ! -------------------------------------------------------------

       mat_r_alp = (-((kt**2*Urr_(i,j))/r**2) + k0**2*KAlpAlp*Urr_(i,j) - (zi*kt*Urt_(i,j))/r**2 + &
          k0**2*KBetAlp*Utr_(i,j) + k0**2*KPrlAlp*Uzr_(i,j) + &
          (2*dzBfn_bfn(i,j,n,m)*dzUrr(i,j)) + &
          (Urr_(i,j)*dzzBfn_bfn(i,j,n,m)) + dzzUrr(i,j) - &
          (zi*kt*Urt_(i,j)*drBfn_bfn(i,j,n,m))/r - &
          (dzUrz(i,j)*drBfn_bfn(i,j,n,m)) - &
          (zi*kt*drUrt(i,j))/r - &
          (dzBfn_bfn(i,j,n,m)*drUrz(i,j)) - &
          (Urz_(i,j)*drzBfn_bfn(i,j,n,m)) - drzUrz(i,j))

       mat_r_bet = (k0**2*KAlpBet*Urr_(i,j) - (kt**2*Utr_(i,j))/r**2 + k0**2*KBetBet*Utr_(i,j) - &
          (zi*kt*Utt_(i,j))/r**2 + k0**2*KPrlBet*Uzr_(i,j) + &
          (2*dzBfn_bfn(i,j,n,m)*dzUtr(i,j)) + &
          (Utr_(i,j)*dzzBfn_bfn(i,j,n,m)) + dzzUtr(i,j) - &
          (zi*kt*Utt_(i,j)*drBfn_bfn(i,j,n,m))/r - &
          (dzUtz(i,j)*drBfn_bfn(i,j,n,m)) - &
          (zi*kt*drUtt(i,j))/r - &
          (dzBfn_bfn(i,j,n,m)*drUtz(i,j)) - &
          (Utz_(i,j)*drzBfn_bfn(i,j,n,m)) - drzUtz(i,j))

       mat_r_prl = (k0**2*KAlpPrl*Urr_(i,j) + k0**2*KBetPrl*Utr_(i,j) - (kt**2*Uzr_(i,j))/r**2 + &
          k0**2*KPrlPrl*Uzr_(i,j) - (zi*kt*Uzt_(i,j))/r**2 + &
          (2*dzBfn_bfn(i,j,n,m)*dzUzr(i,j)) + &
          (Uzr_(i,j)*dzzBfn_bfn(i,j,n,m)) + dzzUzr(i,j) - &
          (zi*kt*Uzt_(i,j)*drBfn_bfn(i,j,n,m))/r - &
          (dzUzz(i,j)*drBfn_bfn(i,j,n,m)) - &
          (zi*kt*drUzt(i,j))/r - &
          (dzBfn_bfn(i,j,n,m)*drUzz(i,j)) - &
          (Uzz_(i,j)*drzBfn_bfn(i,j,n,m)) - drzUzz(i,j))

       mat_th_alp = ((zi*kt*Urr_(i,j))/r**2 - Urt_(i,j)/r**2 + k0**2*KAlpAlp*Urt_(i,j) + &
          k0**2*KBetAlp*Utt_(i,j) + k0**2*KPrlAlp*Uzt_(i,j) - &
          (zi*kt*Urz_(i,j)*dzBfn_bfn(i,j,n,m))/r + &
          (2*dzBfn_bfn(i,j,n,m)*dzUrt(i,j)) - &
          (zi*kt*dzUrz(i,j))/r + (Urt_(i,j)*dzzBfn_bfn(i,j,n,m)) + &
          dzzUrt(i,j) - (zi*kt*Urr_(i,j)*drBfn_bfn(i,j,n,m))/r + &
          (Urt_(i,j)*drBfn_bfn(i,j,n,m))/r - (zi*kt*drUrr(i,j))/r + &
          drUrt(i,j)/r + (2*drBfn_bfn(i,j,n,m)*drUrt(i,j)) + &
          (Urt_(i,j)*drrBfn_bfn(i,j,n,m)) + drrUrt(i,j))

       mat_th_bet = (k0**2*KAlpBet*Urt_(i,j) + (zi*kt*Utr_(i,j))/r**2 - Utt_(i,j)/r**2 + &
          k0**2*KBetBet*Utt_(i,j) + k0**2*KPrlBet*Uzt_(i,j) - &
          (zi*kt*Utz_(i,j)*dzBfn_bfn(i,j,n,m))/r + &
          (2*dzBfn_bfn(i,j,n,m)*dzUtt(i,j)) - &
          (zi*kt*dzUtz(i,j))/r + (Utt_(i,j)*dzzBfn_bfn(i,j,n,m)) + &
          dzzUtt(i,j) - (zi*kt*Utr_(i,j)*drBfn_bfn(i,j,n,m))/r + &
          (Utt_(i,j)*drBfn_bfn(i,j,n,m))/r - (zi*kt*drUtr(i,j))/r + &
          drUtt(i,j)/r + (2*drBfn_bfn(i,j,n,m)*drUtt(i,j)) + &
          (Utt_(i,j)*drrBfn_bfn(i,j,n,m)) + drrUtt(i,j))

       mat_th_prl = (k0**2*KAlpPrl*Urt_(i,j) + k0**2*KBetPrl*Utt_(i,j) + (zi*kt*Uzr_(i,j))/r**2 - &
          Uzt_(i,j)/r**2 + k0**2*KPrlPrl*Uzt_(i,j) - &
          (zi*kt*Uzz_(i,j)*dzBfn_bfn(i,j,n,m))/r + &
          (2*dzBfn_bfn(i,j,n,m)*dzUzt(i,j)) - &
          (zi*kt*dzUzz(i,j))/r + (Uzt_(i,j)*dzzBfn_bfn(i,j,n,m)) + &
          dzzUzt(i,j) - (zi*kt*Uzr_(i,j)*drBfn_bfn(i,j,n,m))/r + &
          (Uzt_(i,j)*drBfn_bfn(i,j,n,m))/r - (zi*kt*drUzr(i,j))/r + &
          drUzt(i,j)/r + (2*drBfn_bfn(i,j,n,m)*drUzt(i,j)) + &
          (Uzt_(i,j)*drrBfn_bfn(i,j,n,m)) + drrUzt(i,j))

       mat_z_alp = (-((kt**2*Urz_(i,j))/r**2) + k0**2*KAlpAlp*Urz_(i,j) + k0**2*KBetAlp*Utz_(i,j) + &
          k0**2*KPrlAlp*Uzz_(i,j) - (Urr_(i,j)*dzBfn_bfn(i,j,n,m))/r - &
          (zi*kt*Urt_(i,j)*dzBfn_bfn(i,j,n,m))/r - dzUrr(i,j)/r - &
          (zi*kt*dzUrt(i,j))/r + (Urz_(i,j)*drBfn_bfn(i,j,n,m))/r - &
          (dzUrr(i,j)*drBfn_bfn(i,j,n,m)) - &
          (dzBfn_bfn(i,j,n,m)*drUrr(i,j)) + drUrz(i,j)/r + &
          (2*drBfn_bfn(i,j,n,m)*drUrz(i,j)) - &
          (Urr_(i,j)*drzBfn_bfn(i,j,n,m)) - drzUrr(i,j) + &
          (Urz_(i,j)*drrBfn_bfn(i,j,n,m)) + drrUrz(i,j))

       mat_z_bet = (k0**2*KAlpBet*Urz_(i,j) - (kt**2*Utz_(i,j))/r**2 + k0**2*KBetBet*Utz_(i,j) + &
          k0**2*KPrlBet*Uzz_(i,j) - (Utr_(i,j)*dzBfn_bfn(i,j,n,m))/r - &
          (zi*kt*Utt_(i,j)*dzBfn_bfn(i,j,n,m))/r - dzUtr(i,j)/r - &
          (zi*kt*dzUtt(i,j))/r + (Utz_(i,j)*drBfn_bfn(i,j,n,m))/r - &
          (dzUtr(i,j)*drBfn_bfn(i,j,n,m)) - &
          (dzBfn_bfn(i,j,n,m)*drUtr(i,j)) + drUtz(i,j)/r + &
          (2*drBfn_bfn(i,j,n,m)*drUtz(i,j)) - &
          (Utr_(i,j)*drzBfn_bfn(i,j,n,m)) - drzUtr(i,j) + &
          (Utz_(i,j)*drrBfn_bfn(i,j,n,m)) + drrUtz(i,j)) 

       mat_z_prl = (k0**2*KAlpPrl*Urz_(i,j) + k0**2*KBetPrl*Utz_(i,j) - (kt**2*Uzz_(i,j))/r**2 + &
          k0**2*KPrlPrl*Uzz_(i,j) - (Uzr_(i,j)*dzBfn_bfn(i,j,n,m))/r - &
          (zi*kt*Uzt_(i,j)*dzBfn_bfn(i,j,n,m))/r - dzUzr(i,j)/r - &
          (zi*kt*dzUzt(i,j))/r + (Uzz_(i,j)*drBfn_bfn(i,j,n,m))/r - &
          (dzUzr(i,j)*drBfn_bfn(i,j,n,m)) - &
          (dzBfn_bfn(i,j,n,m)*drUzr(i,j)) + drUzz(i,j)/r + &
          (2*drBfn_bfn(i,j,n,m)*drUzz(i,j)) - &
          (Uzr_(i,j)*drzBfn_bfn(i,j,n,m)) - drzUzr(i,j) + &
          (Uzz_(i,j)*drrBfn_bfn(i,j,n,m)) + drrUzz(i,j))

                            dxx = mat_r_alp / k0**2
                            dxy = mat_r_bet / k0**2
                            dxz = mat_r_prl / k0**2

                            dyx = mat_th_alp / k0**2
                            dyy = mat_th_bet / k0**2
                            dyz = mat_th_prl / k0**2

                            dzx = mat_z_alp / k0**2
                            dzy = mat_z_bet / k0**2
                            dzz = mat_z_prl / k0**2

           endif

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
                                        if (ii==0 .and. jj==0) aMat(localRow,localCol) = dxx * bFn  
                                        if (ii==1 .and. jj==0) aMat(localRow,localCol) = dxy * bFn   
                                        if (ii==2 .and. jj==0) aMat(localRow,localCol) = dxz * bFn   
                                                                                                        
                                        if (ii==0 .and. jj==1) aMat(localRow,localCol) = dyx * bFn   
                                        if (ii==1 .and. jj==1) aMat(localRow,localCol) = dyy * bFn   
                                        if (ii==2 .and. jj==1) aMat(localRow,localCol) = dyz * bFn   
                                                                                                        
                                        if (ii==0 .and. jj==2) aMat(localRow,localCol) = dzx * bFn   
                                        if (ii==1 .and. jj==2) aMat(localRow,localCol) = dzy * bFn   
                                        if (ii==2 .and. jj==2) aMat(localRow,localCol) = dzz * bFn   


                                        ! Boundary Conditions
                                        ! The seperate X and Y sections are
                                        ! there to catch the scenario where the
                                        ! 2D code is run as 1D and there is only
                                        ! a single Y or X grid pt. 
                                        ! --------------------------------------

                                        if (nPtsX/=1) then
                                        if ( i==1 .or. i==nPtsX ) then 

                                            if (ii==0 .and. jj==0) aMat(localRow,localCol) = bFn * lsWeightFac 
                                            if (ii==1 .and. jj==0) aMat(localRow,localCol) = 0  
                                            if (ii==2 .and. jj==0) aMat(localRow,localCol) = 0  
                                   
                                            if (ii==0 .and. jj==1) aMat(localRow,localCol) = 0  
                                            if (ii==1 .and. jj==1) aMat(localRow,localCol) = bFn * lsWeightFac  
                                            if (ii==2 .and. jj==1) aMat(localRow,localCol) = 0   
                                   
                                            if (ii==0 .and. jj==2) aMat(localRow,localCol) = 0   
                                            if (ii==1 .and. jj==2) aMat(localRow,localCol) = 0   
                                            if (ii==2 .and. jj==2) aMat(localRow,localCol) = bFn * lsWeightFac 

                                        endif
                                        endif

                                        if (nPtsY/=1) then
                                        if ( j==1 .or. j==nPtsY ) then 

                                            if (ii==0 .and. jj==0) aMat(localRow,localCol) = bFn * lsWeightFac 
                                            if (ii==1 .and. jj==0) aMat(localRow,localCol) = 0  
                                            if (ii==2 .and. jj==0) aMat(localRow,localCol) = 0  
                                   
                                            if (ii==0 .and. jj==1) aMat(localRow,localCol) = 0  
                                            if (ii==1 .and. jj==1) aMat(localRow,localCol) = bFn * lsWeightFac  
                                            if (ii==2 .and. jj==1) aMat(localRow,localCol) = 0   
                                   
                                            if (ii==0 .and. jj==2) aMat(localRow,localCol) = 0   
                                            if (ii==1 .and. jj==2) aMat(localRow,localCol) = 0   
                                            if (ii==2 .and. jj==2) aMat(localRow,localCol) = bFn * lsWeightFac 

                                        endif
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

        !if(iAm==0) then
        !write(*,*) 
        !do i=1,nRowLocal
        !    write(*,*) real ( aMat(i,:) )
        !enddo
        !endif

        !aMat_   = aMat
    end subroutine amat_fill

end module mat_fill
