module mat_fill

use constants

implicit none

complex :: cexpkxky
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
#endif

contains

    subroutine amat_fill

        use aorsa2din_mod, &
        only: nModesX, nModesY, &
            delta0, nSpec, &
            iSigma, nPtsX, nPtsY, npRow, npCol, &
            metalLeft, metalRight, &
            metalTop, metalBot, nPhi, square, lsWeightFac
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

        if(square) lsWeightFac = 1

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

                        !write(*,*) i, j, m, n, nMin, nMax, mMin, mMax

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


                            !   interior plasma region:
                            !   ----------------------
                            
                            species: &
                            do s=1,nSpec

                                !if (iSigma==1) & ! hot plasma        
                                !sigma_tmp = sigmaHot_maxwellian&
                                !    ( mSpec(s), &
                                !    ktSpec(i,j,s), omgc(i,j,s), omgp2(i,j,s), &
                                !    kxsav(n), kysav(m), capr(i), &
                                !    omgrf, k0, &
                                !    k_cutoff, s, &
                                !    bMod(i,j), gradPrlB(i,j), &
                                !    U_xyz(i,j,:,:), U_cyl(i,j,:,:), nuOmg2D(i,j) )
                              
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

                            if ( capR(i) < metalLeft .or. capR(i) > metalRight &
                                .or. y(j) > metalTop .or. y(j) < metalBot ) then 

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


                            cexpkxky = xx(n, i) * yy(m, j)
                            !if(n==nMin .or. n==nMax) cexpkxky = cexpkxky / 2
                            !if(m==mMin .or. m==mMax) cexpkxky = cexpkxky / 2

                            kAlpAlp = 1.0 + zi / (eps0 * omgrf) * sigAlpAlp
                            kAlpBet =       zi / (eps0 * omgrf) * sigAlpBet
                            kAlpPrl =       zi / (eps0 * omgrf) * sigAlpPrl

                            kBetAlp =       zi / (eps0 * omgrf) * sigBetAlp
                            kBetBet = 1.0 + zi / (eps0 * omgrf) * sigBetBet
                            kBetPrl =       zi / (eps0 * omgrf) * sigBetPrl

                            kPrlAlp =       zi / (eps0 * omgrf) * sigPrlAlp
                            kPrlBet =       zi / (eps0 * omgrf) * sigPrlBet
                            kPrlPrl = 1.0 + zi / (eps0 * omgrf) * sigPrlPrl

                            !rnx = kxsav(n) / k0
                            !rny = kysav(m) / k0
                            !rnPhi = xkphi(i) / k0

                            !kz  = kysav(m)  
                            !kr  = kxsav(n) 

                            z   = y(j)
                            r   = capR(i)
                            kt  = nPhi / r

        !if(i/=1 .and. i/=nPtsX .and. j/=1 .and. j/=nPtsY) then

        mat_r_alp = (-((kt**2*Urr_(i,j))/r**2) + k0**2*KAlpAlp*Urr_(i,j) - (zi*kt*Urt_(i,j))/r**2 + &
          k0**2*KBetAlp*Utr_(i,j) + k0**2*KPrlAlp*Uzr_(i,j) + &
          (2*dzbfn_bfn(i,j,n,m)*dzUrr(i,j)) + &
          (Urr_(i,j)*dzzbfn_bfn(i,j,n,m)) + dzzUrr(i,j) - &
          (zi*kt*Urt_(i,j)*drbfn_bfn(i,j,n,m))/(r) - &
          (dzUrz(i,j)*drbfn_bfn(i,j,n,m)) - &
          (zi*kt*drUrt(i,j))/r - &
          (dzbfn_bfn(i,j,n,m)*drUrz(i,j)) - &
          (Urz_(i,j)*drzbfn_bfn(i,j,n,m)) - drzUrz(i,j)) 

       mat_r_bet = (k0**2*KAlpBet*Urr_(i,j) - (kt**2*Utr_(i,j))/r**2 + k0**2*KBetBet*Utr_(i,j) - &
          (zi*kt*Utt_(i,j))/r**2 + k0**2*KPrlBet*Uzr_(i,j) + &
          (2*dzbfn_bfn(i,j,n,m)*dzUtr(i,j)) + &
          (Utr_(i,j)*dzzbfn_bfn(i,j,n,m)) + dzzUtr(i,j) - &
          (zi*kt*Utt_(i,j)*drbfn_bfn(i,j,n,m))/(r) - &
          (dzUtz(i,j)*drbfn_bfn(i,j,n,m)) - &
          (zi*kt*drUtt(i,j))/r - &
          (dzbfn_bfn(i,j,n,m)*drUtz(i,j)) - &
          (Utz_(i,j)*drzbfn_bfn(i,j,n,m)) - drzUtz(i,j)) 

       mat_r_prl = (k0**2*KAlpPrl*Urr_(i,j) + k0**2*KBetPrl*Utr_(i,j) - (kt**2*Uzr_(i,j))/r**2 + &
          k0**2*KPrlPrl*Uzr_(i,j) - (zi*kt*Uzt_(i,j))/r**2 + &
          (2*dzbfn_bfn(i,j,n,m)*dzUzr(i,j)) + &
          (Uzr_(i,j)*dzzbfn_bfn(i,j,n,m)) + dzzUzr(i,j) - &
          (zi*kt*Uzt_(i,j)*drbfn_bfn(i,j,n,m))/(r) - &
          (dzUzz(i,j)*drbfn_bfn(i,j,n,m)) - &
          (zi*kt*drUzt(i,j))/r - &
          (dzbfn_bfn(i,j,n,m)*drUzz(i,j)) - &
          (Uzz_(i,j)*drzbfn_bfn(i,j,n,m)) - drzUzz(i,j))

       mat_th_alp = ((zi*kt*Urr_(i,j))/r**2 - Urt_(i,j)/r**2 + k0**2*KAlpAlp*Urt_(i,j) + k0**2*KBetAlp*Utt_(i,j) + &
          k0**2*KPrlAlp*Uzt_(i,j) - (zi*kt*Urz_(i,j)*dzbfn_bfn(i,j,n,m))/(r) + &
          (2*dzbfn_bfn(i,j,n,m)*dzUrt(i,j)) - &
          (zi*kt*dzUrz(i,j))/r + (Urt_(i,j)*dzzbfn_bfn(i,j,n,m)) + &
          dzzUrt(i,j) - (zi*kt*Urr_(i,j)*drbfn_bfn(i,j,n,m))/(r) + &
          (Urt_(i,j)*drbfn_bfn(i,j,n,m))/(r) - (zi*kt*drUrr(i,j))/r + &
          drUrt(i,j)/r + (2*drbfn_bfn(i,j,n,m)*drUrt(i,j)) + &
          (Urt_(i,j)*drrbfn_bfn(i,j,n,m)) + drrUrt(i,j))

       mat_th_bet = (k0**2*KAlpBet*Urt_(i,j) + (zi*kt*Utr_(i,j))/r**2 - Utt_(i,j)/r**2 + k0**2*KBetBet*Utt_(i,j) + &
          k0**2*KPrlBet*Uzt_(i,j) - (zi*kt*Utz_(i,j)*dzbfn_bfn(i,j,n,m))/(r) + &
          (2*dzbfn_bfn(i,j,n,m)*dzUtt(i,j)) - &
          (zi*kt*dzUtz(i,j))/r + (Utt_(i,j)*dzzbfn_bfn(i,j,n,m)) + &
          dzzUtt(i,j) - (zi*kt*Utr_(i,j)*drbfn_bfn(i,j,n,m))/(r) + &
          (Utt_(i,j)*drbfn_bfn(i,j,n,m))/(r) - (zi*kt*drUtr(i,j))/r + &
          drUtt(i,j)/r + (2*drbfn_bfn(i,j,n,m)*drUtt(i,j)) + &
          (Utt_(i,j)*drrbfn_bfn(i,j,n,m)) + drrUtt(i,j))

       mat_th_prl = (k0**2*KAlpPrl*Urt_(i,j) + k0**2*KBetPrl*Utt_(i,j) + (zi*kt*Uzr_(i,j))/r**2 - Uzt_(i,j)/r**2 + &
          k0**2*KPrlPrl*Uzt_(i,j) - (zi*kt*Uzz_(i,j)*dzbfn_bfn(i,j,n,m))/(r) + &
          (2*dzbfn_bfn(i,j,n,m)*dzUzt(i,j)) - &
          (zi*kt*dzUzz(i,j))/r + (Uzt_(i,j)*dzzbfn_bfn(i,j,n,m)) + &
          dzzUzt(i,j) - (zi*kt*Uzr_(i,j)*drbfn_bfn(i,j,n,m))/(r) + &
          (Uzt_(i,j)*drbfn_bfn(i,j,n,m))/(r) - (zi*kt*drUzr(i,j))/r + &
          drUzt(i,j)/r + (2*drbfn_bfn(i,j,n,m)*drUzt(i,j)) + &
          (Uzt_(i,j)*drrbfn_bfn(i,j,n,m)) + drrUzt(i,j)) 


        mat_z_alp = (-((kt**2*Urz_(i,j))/r**2) + k0**2*KAlpAlp*Urz_(i,j) + k0**2*KBetAlp*Utz_(i,j) + &
          k0**2*KPrlAlp*Uzz_(i,j) - (Urr_(i,j)*dzbfn_bfn(i,j,n,m))/(r) - &
          (zi*kt*Urt_(i,j)*dzbfn_bfn(i,j,n,m))/(r) - dzUrr(i,j)/r - &
          (zi*kt*dzUrt(i,j))/r + (Urz_(i,j)*drbfn_bfn(i,j,n,m))/(r) - &
          (dzUrr(i,j)*drbfn_bfn(i,j,n,m)) - &
          (dzbfn_bfn(i,j,n,m)*drUrr(i,j)) + drUrz(i,j)/r + &
          (2*drbfn_bfn(i,j,n,m)*drUrz(i,j)) - &
          (Urr_(i,j)*drzbfn_bfn(i,j,n,m)) - drzUrr(i,j) + &
          (Urz_(i,j)*drrbfn_bfn(i,j,n,m)) + drrUrz(i,j))

       mat_z_bet = (k0**2*KAlpBet*Urz_(i,j) - (kt**2*Utz_(i,j))/r**2 + k0**2*KBetBet*Utz_(i,j) + &
          k0**2*KPrlBet*Uzz_(i,j) - (Utr_(i,j)*dzbfn_bfn(i,j,n,m))/(r) - &
          (zi*kt*Utt_(i,j)*dzbfn_bfn(i,j,n,m))/(r) - dzUtr(i,j)/r - &
          (zi*kt*dzUtt(i,j))/r + (Utz_(i,j)*drbfn_bfn(i,j,n,m))/(r) - &
          (dzUtr(i,j)*drbfn_bfn(i,j,n,m)) - &
          (dzbfn_bfn(i,j,n,m)*drUtr(i,j)) + drUtz(i,j)/r + &
          (2*drbfn_bfn(i,j,n,m)*drUtz(i,j)) - &
          (Utr_(i,j)*drzbfn_bfn(i,j,n,m)) - drzUtr(i,j) + &
          (Utz_(i,j)*drrbfn_bfn(i,j,n,m)) + drrUtz(i,j))

       mat_z_prl = (k0**2*KAlpPrl*Urz_(i,j) + k0**2*KBetPrl*Utz_(i,j) - (kt**2*Uzz_(i,j))/r**2 + &
          k0**2*KPrlPrl*Uzz_(i,j) - (Uzr_(i,j)*drbfn_bfn(i,j,n,m))/(r) - &
          (zi*kt*Uzt_(i,j)*dzbfn_bfn(i,j,n,m))/(r) - dzUzr(i,j)/r - &
          (zi*kt*dzUzt(i,j))/r + (Uzz_(i,j)*drbfn_bfn(i,j,n,m))/(r) - &
          (dzUzr(i,j)*drbfn_bfn(i,j,n,m)) - &
          (dzbfn_bfn(i,j,n,m)*drUzr(i,j)) + drUzz(i,j)/r + &
          (2*drbfn_bfn(i,j,n,m)*drUzz(i,j)) - &
          (Uzr_(i,j)*drzbfn_bfn(i,j,n,m)) - drzUzr(i,j) + &
          (Uzz_(i,j)*drrbfn_bfn(i,j,n,m)) + drrUzz(i,j))

                            dxx = mat_r_alp / k0**2
                            dxy = mat_r_bet / k0**2
                            dxz = mat_r_prl / k0**2

                            dyx = mat_th_alp / k0**2
                            dyy = mat_th_bet / k0**2
                            dyz = mat_th_prl / k0**2

                            dzx = mat_z_alp / k0**2
                            dzy = mat_z_bet / k0**2
                            dzz = mat_z_prl / k0**2

           !endif

                            !dxx = (kAlpAlp - rny**2 - rnphi**2) * uxx(i,j) &
                            !    +  kBetAlp * uyx(i,j) &
                            !    +  kPrlAlp * uzx(i,j) &
                            !    + rnx * (rny * uxy(i,j) + rnphi * uxz(i,j)) &
                            !    - zi * rnphi / k0 * &
                            !             (uxz(i,j) / capr(i) + dxuxz(i,j)) &
                            !    - zi * rny / k0 * (dxuxy(i,j) - 2. * dyuxx(i,j)) &
                            !    - zi * rnx / k0 * dyuxy(i,j) &
                            !    + 1. / k0**2 * (dyyuxx(i,j) - dxyuxy(i,j))

                            !dxy =  kAlpBet * uxx(i,j) &
                            !    + (kBetBet - rny**2 - rnphi**2) * uyx(i,j) &
                            !    +  kPrlBet * uzx(i,j) &
                            !    + rnx * (rny * uyy(i,j) + rnphi * uyz(i,j)) &
                            !    - zi * rnphi / k0 * &
                            !             (uyz(i,j) / capr(i) + dxuyz(i,j)) &
                            !    - zi * rny / k0 * (dxuyy(i,j) - 2. * dyuyx(i,j)) &
                            !    - zi * rnx / k0 * dyuyy(i,j) &
                            !    + 1. / k0**2 * (dyyuyx(i,j) - dxyuyy(i,j))
                            !dxz =  kAlpPrl * uxx(i,j) &
                            !    +  kBetPrl * uyx(i,j) &
                            !    + (kPrlPrl - rny**2 - rnphi**2) * uzx(i,j) &
                            !    + rnx * (rny * uzy(i,j) + rnphi * uzz(i,j)) &
                            !    - zi * rnphi / k0 * &
                            !             (uzz(i,j) / capr(i) + dxuzz(i,j)) &
                            !    - zi * rny / k0 * (dxuzy(i,j) - 2. * dyuzx(i,j)) &
                            !    - zi * rnx / k0 * dyuzy(i,j) &
                            !    + 1. / k0**2 * (dyyuzx(i,j) - dxyuzy(i,j))
                            !dyx = (kAlpAlp - rnx**2 - rnphi**2) * uxy(i,j) &
                            !    +  kBetAlp * uyy(i,j) &
                            !    +  kPrlAlp * uzy(i,j) &
                            !    + rny * (rnx * uxx(i,j) + rnphi * uxz(i,j)) &
                            !    - zi * rny / k0 * &
                            !             (dxuxx(i,j) + uxx(i,j) / capr(i)) &
                            !    - zi * rnphi / k0 * dyuxz(i,j) &
                            !    - zi * rnx / k0 * &
                            !    (dyuxx(i,j) - uxy(i,j) / capr(i) - 2.* dxuxy(i,j)) &
                            !    + 1. / k0**2 * (dxxuxy(i,j) - dxyuxx(i,j) &
                            !    - dyuxx(i,j)/ capr(i) + dxuxy(i,j) / capr(i))
                            !dyy =  kAlpBet * uxy(i,j) &
                            !    + (kBetBet - rnx**2 - rnphi**2) * uyy(i,j) &
                            !    +  kPrlBet * uzy(i,j) &
                            !    + rny * (rnx * uyx(i,j) + rnphi * uyz(i,j)) &
                            !    - zi * rny / k0 * &
                            !             (dxuyx(i,j) + uyx(i,j) / capr(i)) &
                            !    - zi * rnphi / k0 * dyuyz(i,j) &
                            !    - zi * rnx / k0 * &
                            !    (dyuyx(i,j) - uyy(i,j) / capr(i) - 2.* dxuyy(i,j)) &
                            !    + 1. / k0**2 * (dxxuyy(i,j) - dxyuyx(i,j) &
                            !    - dyuyx(i,j)/ capr(i) + dxuyy(i,j) / capr(i))
                            !dyz =  kAlpPrl * uxy(i,j) &
                            !    +  kBetPrl * uyy(i,j) &
                            !    + (kPrlPrl - rnx**2 - rnphi**2) * uzy(i,j) &
                            !    + rny * (rnx * uzx(i,j) + rnphi * uzz(i,j)) &
                            !    - zi * rny / k0 * &
                            !             (dxuzx(i,j) + uzx(i,j) / capr(i)) &
                            !    - zi * rnphi / k0 * dyuzz(i,j) &
                            !    - zi * rnx / k0 * &
                            !    (dyuzx(i,j) - uzy(i,j) / capr(i) - 2.* dxuzy(i,j)) &
                            !    + 1. / k0**2 * (dxxuzy(i,j) - dxyuzx(i,j) &
                            !    - dyuzx(i,j)/ capr(i) + dxuzy(i,j) / capr(i))

                            !dzx = (kAlpAlp - rnx**2 - rny**2) * uxz(i,j) &
                            !    +  kBetAlp * uyz(i,j) &
                            !    +  kPrlAlp * uzz(i,j) &
                            !    + rnphi * (rny * uxy(i,j) + rnx * uxx(i,j)) &
                            !    + zi * rny / k0 * 2. * dyuxz(i,j) &
                            !    - zi * rnphi / k0 * &
                            !       (dyuxy(i,j) + dxuxx(i,j) - uxx(i,j) / capr(i)) &
                            !    + zi * rnx / k0 * &
                            !       (uxz(i,j) / capr(i) + 2.* dxuxz(i,j)) &
                            !    - 1. / (k0**2 * capr(i)) * &
                            !       (uxz(i,j)/ capr(i) - dxuxz(i,j)) &
                            !    + 1. / k0**2  * (dxxuxz(i,j) + dyyuxz(i,j))

                            !dzy =  kAlpBet * uxz(i,j) &
                            !    + (kBetBet - rnx**2 - rny**2) * uyz(i,j) &
                            !    +  kPrlBet * uzz(i,j) &
                            !    + rnphi * (rny * uyy(i,j) + rnx * uyx(i,j)) &
                            !    + zi * rny / k0 * 2. * dyuyz(i,j) &
                            !    - zi * rnphi / k0 * &
                            !       (dyuyy(i,j) + dxuyx(i,j) - uyx(i,j) / capr(i)) &
                            !    + zi * rnx / k0 * &
                            !       (uyz(i,j) / capr(i) + 2.* dxuyz(i,j)) &
                            !    - 1. / (k0**2 * capr(i)) * &
                            !       (uyz(i,j)/ capr(i) - dxuyz(i,j)) &
                            !    + 1. / k0**2  * (dxxuyz(i,j) + dyyuyz(i,j))

                            !dzz =  kAlpPrl * uxz(i,j) &
                            !    +  kBetPrl * uyz(i,j) &
                            !    + (kPrlPrl - rnx**2 - rny**2) * uzz(i,j) &
                            !    + rnphi * (rny * uzy(i,j) + rnx * uzx(i,j)) &
                            !    + zi * rny / k0 * 2. * dyuzz(i,j) &
                            !    - zi * rnphi / k0 * &
                            !       (dyuzy(i,j) + dxuzx(i,j) - uzx(i,j) / capr(i)) &
                            !    + zi * rnx / k0 * &
                            !       (uzz(i,j) / capr(i) + 2.* dxuzz(i,j)) &
                            !    - 1. / (k0**2 * capr(i)) * &
                            !       (uzz(i,j)/ capr(i) - dxuzz(i,j)) &
                            !    + 1. / k0**2  * (dxxuzz(i,j) + dyyuzz(i,j))


                            !if(iAm==0) then
                            !   write(*,*) dxx, dxy, dxz
                            !   write(*,*) dyx, dyy, dyz
                            !   write(*,*) dzx, dzy, dzz

                            !   write(*,*)

                            !   write(*,*) mat_r_alp/k0**2, mat_r_bet/k0**2, mat_r_prl/k0**2
                            !   write(*,*) mat_th_alp/k0**2, mat_th_bet/k0**2, mat_th_prl/k0**2
                            !   write(*,*) mat_z_alp/k0**2, mat_z_bet/k0**2, mat_z_prl/k0**2
                            !endif
                            !stop

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
                                        !write(*,*) iRow+ii, iCol+jj, l_sp, m_sp, &
                                        !    pr_sp, pc_sp, x_sp, y_sp, localRow, localCol

                                        if (ii==0 .and. jj==0) aMat(localRow,localCol) = cexpkxky * dxx 
                                        if (ii==1 .and. jj==0) aMat(localRow,localCol) = cexpkxky * dxy  
                                        if (ii==2 .and. jj==0) aMat(localRow,localCol) = cexpkxky * dxz  
                                   
                                        if (ii==0 .and. jj==1) aMat(localRow,localCol) = cexpkxky * dyx  
                                        if (ii==1 .and. jj==1) aMat(localRow,localCol) = cexpkxky * dyy  
                                        if (ii==2 .and. jj==1) aMat(localRow,localCol) = cexpkxky * dyz  
                                   
                                        if (ii==0 .and. jj==2) aMat(localRow,localCol) = cexpkxky * dzx  
                                        if (ii==1 .and. jj==2) aMat(localRow,localCol) = cexpkxky * dzy  
                                        if (ii==2 .and. jj==2) aMat(localRow,localCol) = cexpkxky * dzz  


                                        !!   boundary conditions
                                        !!   -------------------

                                        !if (nPtsX/=1) then
                                        !if ( i==1 .or. i==nPtsX &
                                        !        .or. (.not. mask(i,j)) ) then 

                                        !    if (ii==0 .and. jj==0) aMat(localRow,localCol) = cexpkxky * lsWeightFac 
                                        !    if (ii==1 .and. jj==0) aMat(localRow,localCol) = 0  
                                        !    if (ii==2 .and. jj==0) aMat(localRow,localCol) = 0  
                                   
                                        !    if (ii==0 .and. jj==1) aMat(localRow,localCol) = 0  
                                        !    if (ii==1 .and. jj==1) aMat(localRow,localCol) = cexpkxky * lsWeightFac  
                                        !    if (ii==2 .and. jj==1) aMat(localRow,localCol) = 0   
                                   
                                        !    if (ii==0 .and. jj==2) aMat(localRow,localCol) = 0   
                                        !    if (ii==1 .and. jj==2) aMat(localRow,localCol) = 0   
                                        !    if (ii==2 .and. jj==2) aMat(localRow,localCol) = cexpkxky * lsWeightFac 

                                        !endif
                                        !endif

                                        !if (nPtsY/=1) then
                                        !if ( j==1 .or. j==nPtsY &
                                        !        .or. (.not. mask(i,j)) ) then 

                                        !    if (ii==0 .and. jj==0) aMat(localRow,localCol) = cexpkxky * lsWeightFac 
                                        !    if (ii==1 .and. jj==0) aMat(localRow,localCol) = 0  
                                        !    if (ii==2 .and. jj==0) aMat(localRow,localCol) = 0  
                                   
                                        !    if (ii==0 .and. jj==1) aMat(localRow,localCol) = 0  
                                        !    if (ii==1 .and. jj==1) aMat(localRow,localCol) = cexpkxky * lsWeightFac  
                                        !    if (ii==2 .and. jj==1) aMat(localRow,localCol) = 0   
                                   
                                        !    if (ii==0 .and. jj==2) aMat(localRow,localCol) = 0   
                                        !    if (ii==1 .and. jj==2) aMat(localRow,localCol) = 0   
                                        !    if (ii==2 .and. jj==2) aMat(localRow,localCol) = cexpkxky * lsWeightFac 

                                        !endif
                                        !endif


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

    end subroutine amat_fill

end module mat_fill
