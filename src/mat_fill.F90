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
            metalTop, metalBot, nPhi
        use grid
        use sigma_mod
        use rotation
        use profiles
        use bField
        use parallel

        implicit none

        integer :: iRow, iCol, i, j, n, m, p, s, ii, jj
        integer :: localRow, localCol
        complex :: metal
        real :: kr, kz, z
        complex :: &
            dxx, dxy, dxz, &
            dyx, dyy, dyz, &
            dzx, dzy, dzz 

        complex :: &
            mat_r_alp, mat_r_bet, mat_r_prl, &
            mat_th_alp, mat_th_bet, mat_th_prl, &
            mat_z_alp, mat_z_bet, mat_z_prl

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
        metal   = ( 0,1e4 )

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
                do n=kxL,kxR
                    m_loop: &
                    do m=kyL,kyR

                        !write(*,*) i, j, m, n, kxL, kxR, kyL, kyR

                        iCol = (n-kxL) * 3 * nModesY + (m-kyL) * 3 + 1
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

                                if (iSigma==1) & ! hot plasma        
                                sigma_tmp = sigmaHot_maxwellian&
                                    ( mSpec(s), &
                                    ktSpec(i,j,s), omgc(i,j,s), omgp2(i,j,s), &
                                    kxsav(n), kysav(m), capr(i), &
                                    omgrf, k0, &
                                    xk_cutoff, s, &
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

                            cexpkxky = xx(n, i) * yy(m, j)

                            kAlpAlp = 1.0 + zi / (eps0 * omgrf) * sigAlpAlp
                            kAlpBet =       zi / (eps0 * omgrf) * sigAlpBet
                            kAlpPrl =       zi / (eps0 * omgrf) * sigAlpPrl

                            kBetAlp =       zi / (eps0 * omgrf) * sigBetAlp
                            kBetBet = 1.0 + zi / (eps0 * omgrf) * sigBetBet
                            kBetPrl =       zi / (eps0 * omgrf) * sigBetPrl

                            kPrlAlp =       zi / (eps0 * omgrf) * sigPrlAlp
                            kPrlBet =       zi / (eps0 * omgrf) * sigPrlBet
                            kPrlPrl = 1.0 + zi / (eps0 * omgrf) * sigPrlPrl

                            rnx = kxsav(n) / k0
                            rny = kysav(m) / k0
                            rnPhi = xkphi(i) / k0

                            kz  = kysav(m)  
                            kr  = kxsav(n) 
                            z   = y(j)
#ifdef cylProper
                            mat_r_alp = -((nPhi**2*Urr_(i,j))/capR(i)**2) &
                                - kz**2*Urr_(i,j) + &
                                k0**2*KAlpAlp*Urr_(i,j) + &
                                ((-zi + capR(i)*kr)*nPhi*Urth_(i,j))/capR(i)**2 &
                                + kr*kz*Urz_(i,j) + &
                                k0**2*KBetAlp*Uthr_(i,j) + &
                                k0**2*KPrlAlp*Uzr_(i,j) + &
                                2*zi*kz*dzUrr(i,j) - &
                                zi*kr*dzUrz(i,j) + &
                                dzzUrr(i,j) - &
                                (zi*nPhi*drUrth(i,j))/capR(i) - &
                                zi*kz*drUrz(i,j) - drzUrz(i,j)

                            mat_r_bet = k0**2*KAlpBet*Urr_(i,j) &
                                - (nPhi**2*Uthr_(i,j))/capR(i)**2 - &
                                kz**2*Uthr_(i,j) + k0**2*KBetBet*Uthr_(i,j) - &
                                (zi*nPhi*Uthth_(i,j))/capR(i)**2 &
                                + (kr*nPhi*Uthth_(i,j))/capR(i) + &
                                kr*kz*Uthz_(i,j) + k0**2*KPrlBet*Uzr_(i,j) + &
                                2*zi*kz*dzUthr(i,j) &
                                - zi*kr*dzUthz(i,j) + &
                                dzzUthr(i,j) &
                                - (zi*nPhi*drUthth(i,j))/capR(i) - &
                                zi*kz*drUthz(i,j) &
                                - drzUthz(i,j)

                            mat_r_prl = k0**2*KAlpPrl*Urr_(i,j) + &
                                k0**2*KBetPrl*Uthr_(i,j) - &
                                (nPhi**2*Uzr_(i,j))/capR(i)**2 - kz**2*Uzr_(i,j) + &
                                k0**2*KPrlPrl*Uzr_(i,j) &
                                - (zi*nPhi*Uzth_(i,j))/capR(i)**2 + &
                                (kr*nPhi*Uzth_(i,j))/capR(i) + kr*kz*Uzz_(i,j) + &
                                2*zi*kz*dzUzr(i,j) &
                                - zi*kr*dzUzz(i,j) + &
                                dzzUzr(i,j) &
                                - (zi*nPhi*drUzth(i,j))/capR(i) - &
                                zi*kz*drUzz(i,j) &
                                - drzUzz(i,j)

                            mat_th_alp = (((zi + capR(i)*kr)*nPhi*Urr_(i,j) &
                                - Urth_(i,j) + zi*capR(i)*kr*Urth_(i,j) - &
                                capR(i)**2*kr**2*Urth_(i,j) - capR(i)**2*kz**2*Urth_(i,j) + &
                                capR(i)**2*k0**2*KAlpAlp*Urth_(i,j) &
                                + capR(i)*kz*nPhi*Urz_(i,j) + &
                                capR(i)**2*k0**2*KBetAlp*Uthth_(i,j) + &
                                capR(i)**2*k0**2*KPrlAlp*Uzth_(i,j) + &
                                2*zi*capR(i)**2*kz*dzUrth(i,j) - &
                                zi*capR(i)*nPhi*dzUrz(i,j) + &
                                capR(i)**2*dzzUrth(i,j) - &
                                zi*capR(i)*nPhi*drUrr(i,j) + &
                                capR(i)*drUrth(i,j) + &
                                2*zi*capR(i)**2*kr*drUrth(i,j) + &
                                capR(i)**2*drrUrth(i,j)))/capR(i)**2 

                            mat_th_bet = ((capR(i)**2*k0**2*KAlpBet*Urth_(i,j) &
                                + zi*nPhi*Uthr_(i,j) + &
                                capR(i)*kr*nPhi*Uthr_(i,j) &
                                - Uthth_(i,j) + zi*capR(i)*kr*Uthth_(i,j) - &
                                capR(i)**2*kr**2*Uthth_(i,j) &
                                - capR(i)**2*kz**2*Uthth_(i,j) + &
                                capR(i)**2*k0**2*KBetBet*Uthth_(i,j) &
                                + capR(i)*kz*nPhi*Uthz_(i,j) + &
                                capR(i)**2*k0**2*KPrlBet*Uzth_(i,j) + &
                                2*zi*capR(i)**2*kz*dzUthth(i,j) - &
                                zi*capR(i)*nPhi*dzUthz(i,j) + &
                                capR(i)**2*dzzUthth(i,j) - &
                                zi*capR(i)*nPhi*drUthr(i,j) + &
                                capR(i)*drUthth(i,j) + &
                                2*zi*capR(i)**2*kr*drUthth(i,j) + &
                                capR(i)**2*drrUthth(i,j)))/capR(i)**2 

                            mat_th_prl = ((capR(i)**2*k0**2*KAlpPrl*Urth_(i,j) + &
                                capR(i)**2*k0**2*KBetPrl*Uthth_(i,j) &
                                + zi*nPhi*Uzr_(i,j) + &
                                capR(i)*kr*nPhi*Uzr_(i,j) - Uzth_(i,j) &
                                + zi*capR(i)*kr*Uzth_(i,j) - &
                                capR(i)**2*kr**2*Uzth_(i,j) - capR(i)**2*kz**2*Uzth_(i,j) + &
                                capR(i)**2*k0**2*KPrlPrl*Uzth_(i,j) &
                                + capR(i)*kz*nPhi*Uzz_(i,j) + &
                                2*zi*capR(i)**2*kz*dzUzth(i,j) - &
                                zi*capR(i)*nPhi*dzUzz(i,j) + &
                                capR(i)**2*dzzUzth(i,j) - &
                                zi*capR(i)*nPhi*drUzr(i,j) + &
                                capR(i)*drUzth(i,j) + &
                                2*zi*capR(i)**2*kr*drUzth(i,j) + &
                                capR(i)**2*drrUzth(i,j)))/capR(i)**2 

                            mat_z_alp = ((-zi + capR(i)*kr)*kz*Urr_(i,j))/capR(i) &
                                + (kz*nPhi*Urth_(i,j))/capR(i) + &
                                (zi*kr*Urz_(i,j))/capR(i) - kr**2*Urz_(i,j) &
                                - (nPhi**2*Urz_(i,j))/capR(i)**2 + &
                                k0**2*KAlpAlp*Urz_(i,j) &
                                + k0**2*KBetAlp*Uthz_(i,j) + &
                                k0**2*KPrlAlp*Uzz_(i,j) &
                                - dzUrr(i,j)/capR(i) - &
                                zi*kr*dzUrr(i,j) - &
                                (zi*nPhi*dzUrth(i,j))/capR(i) - &
                                zi*kz*drUrr(i,j) &
                                + drUrz(i,j)/capR(i) + &
                                2*zi*kr*drUrz(i,j) &
                                - drzUrr(i,j) + &
                                drrUrz(i,j)

                            mat_z_bet = k0**2*KAlpBet*Urz_(i,j) &
                                - (zi*kz*Uthr_(i,j))/capR(i) + &
                                kr*kz*Uthr_(i,j) + (kz*nPhi*Uthth_(i,j))/capR(i) &
                                + (zi*kr*Uthz_(i,j))/capR(i) - &
                                kr**2*Uthz_(i,j) - (nPhi**2*Uthz_(i,j))/capR(i)**2 + &
                                k0**2*KBetBet*Uthz_(i,j) &
                                + k0**2*KPrlBet*Uzz_(i,j) - &
                                dzUthr(i,j)/capR(i) &
                                - zi*kr*dzUthr(i,j) - &
                                (zi*nPhi*dzUthth(i,j))/capR(i) - &
                                zi*kz*drUthr(i,j) &
                                + drUthz(i,j)/capR(i) + &
                                2*zi*kr*drUthz(i,j) &
                                - drzUthr(i,j) + &
                                drrUthz(i,j)

                            mat_z_prl =k0**2*KAlpPrl*Urz_(i,j) &
                                + k0**2*KBetPrl*Uthz_(i,j) - &
                                (zi*kz*Uzr_(i,j))/capR(i) + kr*kz*Uzr_(i,j) &
                                + (kz*nPhi*Uzth_(i,j))/capR(i) + &
                                (zi*kr*Uzz_(i,j))/capR(i) - kr**2*Uzz_(i,j) &
                                - (nPhi**2*Uzz_(i,j))/capR(i)**2 + &
                                k0**2*KPrlPrl*Uzz_(i,j) &
                                - dzUzr(i,j)/capR(i) - &
                                zi*kr*dzUzr(i,j) - &
                                (zi*nPhi*dzUzth(i,j))/capR(i) - &
                                zi*kz*drUzr(i,j) &
                                + drUzz(i,j)/capR(i) + &
                                2*zi*kr*drUzz(i,j) &
                                - drzUzr(i,j) + &
                                drrUzz(i,j) 
                                
                            dxx = mat_r_alp / k0**2
                            dxy = mat_r_bet / k0**2
                            dxz = mat_r_prl / k0**2

                            dyx = mat_th_alp / k0**2
                            dyy = mat_th_bet / k0**2
                            dyz = mat_th_prl / k0**2

                            dzx = mat_z_alp / k0**2
                            dzy = mat_z_bet / k0**2
                            dzz = mat_z_prl / k0**2

#else
                            dxx = (kAlpAlp - rny**2 - rnphi**2) * uxx(i,j) &
                                +  kBetAlp * uyx(i,j) &
                                +  kPrlAlp * uzx(i,j) &
                                + rnx * (rny * uxy(i,j) + rnphi * uxz(i,j)) &
                                - zi * rnphi / k0 * &
                                         (uxz(i,j) / capr(i) + dxuxz(i,j)) &
                                - zi * rny / k0 * (dxuxy(i,j) - 2. * dyuxx(i,j)) &
                                - zi * rnx / k0 * dyuxy(i,j) &
                                + 1. / k0**2 * (dyyuxx(i,j) - dxyuxy(i,j))

                            dxy =  kAlpBet * uxx(i,j) &
                                + (kBetBet - rny**2 - rnphi**2) * uyx(i,j) &
                                +  kPrlBet * uzx(i,j) &
                                + rnx * (rny * uyy(i,j) + rnphi * uyz(i,j)) &
                                - zi * rnphi / k0 * &
                                         (uyz(i,j) / capr(i) + dxuyz(i,j)) &
                                - zi * rny / k0 * (dxuyy(i,j) - 2. * dyuyx(i,j)) &
                                - zi * rnx / k0 * dyuyy(i,j) &
                                + 1. / k0**2 * (dyyuyx(i,j) - dxyuyy(i,j))
                            dxz =  kAlpPrl * uxx(i,j) &
                                +  kBetPrl * uyx(i,j) &
                                + (kPrlPrl - rny**2 - rnphi**2) * uzx(i,j) &
                                + rnx * (rny * uzy(i,j) + rnphi * uzz(i,j)) &
                                - zi * rnphi / k0 * &
                                         (uzz(i,j) / capr(i) + dxuzz(i,j)) &
                                - zi * rny / k0 * (dxuzy(i,j) - 2. * dyuzx(i,j)) &
                                - zi * rnx / k0 * dyuzy(i,j) &
                                + 1. / k0**2 * (dyyuzx(i,j) - dxyuzy(i,j))
                           dyx = (kAlpAlp - rnx**2 - rnphi**2) * uxy(i,j) &
                                +  kBetAlp * uyy(i,j) &
                                +  kPrlAlp * uzy(i,j) &
                                + rny * (rnx * uxx(i,j) + rnphi * uxz(i,j)) &
                                - zi * rny / k0 * &
                                         (dxuxx(i,j) + uxx(i,j) / capr(i)) &
                                - zi * rnphi / k0 * dyuxz(i,j) &
                                - zi * rnx / k0 * &
                                (dyuxx(i,j) - uxy(i,j) / capr(i) - 2.* dxuxy(i,j)) &
                                + 1. / k0**2 * (dxxuxy(i,j) - dxyuxx(i,j) &
                                - dyuxx(i,j)/ capr(i) + dxuxy(i,j) / capr(i))
                            dyy =  kAlpBet * uxy(i,j) &
                                + (kBetBet - rnx**2 - rnphi**2) * uyy(i,j) &
                                +  kPrlBet * uzy(i,j) &
                                + rny * (rnx * uyx(i,j) + rnphi * uyz(i,j)) &
                                - zi * rny / k0 * &
                                         (dxuyx(i,j) + uyx(i,j) / capr(i)) &
                                - zi * rnphi / k0 * dyuyz(i,j) &
                                - zi * rnx / k0 * &
                                (dyuyx(i,j) - uyy(i,j) / capr(i) - 2.* dxuyy(i,j)) &
                                + 1. / k0**2 * (dxxuyy(i,j) - dxyuyx(i,j) &
                                - dyuyx(i,j)/ capr(i) + dxuyy(i,j) / capr(i))
                            dyz =  kAlpPrl * uxy(i,j) &
                                +  kBetPrl * uyy(i,j) &
                                + (kPrlPrl - rnx**2 - rnphi**2) * uzy(i,j) &
                                + rny * (rnx * uzx(i,j) + rnphi * uzz(i,j)) &
                                - zi * rny / k0 * &
                                         (dxuzx(i,j) + uzx(i,j) / capr(i)) &
                                - zi * rnphi / k0 * dyuzz(i,j) &
                                - zi * rnx / k0 * &
                                (dyuzx(i,j) - uzy(i,j) / capr(i) - 2.* dxuzy(i,j)) &
                                + 1. / k0**2 * (dxxuzy(i,j) - dxyuzx(i,j) &
                                - dyuzx(i,j)/ capr(i) + dxuzy(i,j) / capr(i))

                            dzx = (kAlpAlp - rnx**2 - rny**2) * uxz(i,j) &
                                +  kBetAlp * uyz(i,j) &
                                +  kPrlAlp * uzz(i,j) &
                                + rnphi * (rny * uxy(i,j) + rnx * uxx(i,j)) &
                                + zi * rny / k0 * 2. * dyuxz(i,j) &
                                - zi * rnphi / k0 * &
                                   (dyuxy(i,j) + dxuxx(i,j) - uxx(i,j) / capr(i)) &
                                + zi * rnx / k0 * &
                                   (uxz(i,j) / capr(i) + 2.* dxuxz(i,j)) &
                                - 1. / (k0**2 * capr(i)) * &
                                   (uxz(i,j)/ capr(i) - dxuxz(i,j)) &
                                + 1. / k0**2  * (dxxuxz(i,j) + dyyuxz(i,j))

                            dzy =  kAlpBet * uxz(i,j) &
                                + (kBetBet - rnx**2 - rny**2) * uyz(i,j) &
                                +  kPrlBet * uzz(i,j) &
                                + rnphi * (rny * uyy(i,j) + rnx * uyx(i,j)) &
                                + zi * rny / k0 * 2. * dyuyz(i,j) &
                                - zi * rnphi / k0 * &
                                   (dyuyy(i,j) + dxuyx(i,j) - uyx(i,j) / capr(i)) &
                                + zi * rnx / k0 * &
                                   (uyz(i,j) / capr(i) + 2.* dxuyz(i,j)) &
                                - 1. / (k0**2 * capr(i)) * &
                                   (uyz(i,j)/ capr(i) - dxuyz(i,j)) &
                                + 1. / k0**2  * (dxxuyz(i,j) + dyyuyz(i,j))

                            dzz =  kAlpPrl * uxz(i,j) &
                                +  kBetPrl * uyz(i,j) &
                                + (kPrlPrl - rnx**2 - rny**2) * uzz(i,j) &
                                + rnphi * (rny * uzy(i,j) + rnx * uzx(i,j)) &
                                + zi * rny / k0 * 2. * dyuzz(i,j) &
                                - zi * rnphi / k0 * &
                                   (dyuzy(i,j) + dxuzx(i,j) - uzx(i,j) / capr(i)) &
                                + zi * rnx / k0 * &
                                   (uzz(i,j) / capr(i) + 2.* dxuzz(i,j)) &
                                - 1. / (k0**2 * capr(i)) * &
                                   (uzz(i,j)/ capr(i) - dxuzz(i,j)) &
                                + 1. / k0**2  * (dxxuzz(i,j) + dyyuzz(i,j))
#endif

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


                                        !   boundary conditions
                                        !   -------------------

                                        if ( i==1 .or. i==nPtsX &
                                                .or. j==1 .or. j==nPtsY &
                                                .or. (.not. mask(i,j)) ) then 

                                            if (ii==0 .and. jj==0) aMat(localRow,localCol) = cexpkxky 
                                            if (ii==1 .and. jj==0) aMat(localRow,localCol) = 0  
                                            if (ii==2 .and. jj==0) aMat(localRow,localCol) = 0  
                                   
                                            if (ii==0 .and. jj==1) aMat(localRow,localCol) = 0  
                                            if (ii==1 .and. jj==1) aMat(localRow,localCol) = cexpkxky   
                                            if (ii==2 .and. jj==1) aMat(localRow,localCol) = 0   
                                   
                                            if (ii==0 .and. jj==2) aMat(localRow,localCol) = 0   
                                            if (ii==1 .and. jj==2) aMat(localRow,localCol) = 0   
                                            if (ii==2 .and. jj==2) aMat(localRow,localCol) = cexpkxky   

                                        endif

                                        if ( n<kxL/2 .or. n>kxR/2 &
                                            .or. m<kyL/2 .or. m>kyR/2 ) then 
                                            
                                            if (ii==0 .and. jj==0) aMat(localRow,localCol) = cexpkxky 
                                            if (ii==1 .and. jj==0) aMat(localRow,localCol) = 0  
                                            if (ii==2 .and. jj==0) aMat(localRow,localCol) = 0  
                                   
                                            if (ii==0 .and. jj==1) aMat(localRow,localCol) = 0  
                                            if (ii==1 .and. jj==1) aMat(localRow,localCol) = cexpkxky   
                                            if (ii==2 .and. jj==1) aMat(localRow,localCol) = 0   
                                   
                                            if (ii==0 .and. jj==2) aMat(localRow,localCol) = 0   
                                            if (ii==1 .and. jj==2) aMat(localRow,localCol) = 0   
                                            if (ii==2 .and. jj==2) aMat(localRow,localCol) = cexpkxky   

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

    end subroutine amat_fill

end module mat_fill
