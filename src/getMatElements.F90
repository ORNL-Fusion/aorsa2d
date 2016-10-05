module getMatElements

contains

function get3by3Block( g, w)!, r, z)

    use grid
    use aorsaNamelist, only: &
       chebyshevX, chebyshevY, iSigma, nSpec, nPhi, ZeroJp, &
       ZeroJp_rMin, ZeroJp_rMax, ZeroJp_zMin, ZeroJp_zMax, &
       useJpFromFile, coldIons

    use constants
    use profiles, only: omgrf, mSpec
    use sigma
    use generic_biLinearInterp

    implicit none

    type(gridBlock), intent(in) :: g
    integer, intent(in) :: w

    logical :: HotSpecies

    complex :: get3by3Block(3,3)

    complex(kind=dbl) :: &
        sigAlpAlp, sigAlpBet, sigAlpPrl, &
        sigBetAlp, sigBetBet, sigBetPrl, &
        sigPrlAlp, sigPrlBet, sigPrlPrl

    type(dBfnArg) :: d
    complex(kind=dbl) :: sigma_tmp(3,3), sigma_tmp_neg(3,3), sigmaHere(3,3), &
            sigma_tmp2(3,3),sigma_tmp3(3,3)

    real :: r, z, kr, kt, kz, kVec_stix(3)

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

    complex :: &
        kAlpAlp, kAlpBet, kAlpPrl, &
        kBetAlp, kBetBet, kBetPrl, &
        kPrlAlp, kPrlBet, kPrlPrl

    complex :: &
        mat_r_alp, mat_r_bet, mat_r_prl, &
        mat_t_alp, mat_t_bet, mat_t_prl, &
        mat_z_alp, mat_z_bet, mat_z_prl

    integer :: s

    type(spatialSigmaInput_cold) :: sigmaIn_cold
    real :: R_(3,3)
    complex(kind=dbl) :: k0

    z   = g%z(g%wl(w)%j)
    r   = g%R(g%wl(w)%i)
    k0  = g%k0(g%wl(w)%iPt)
    kt  = nPhi!g%kPhi(i)

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
            if(g%wl(w)%n>1) then
                kr = g%wl(w)%n / sqrt ( sin ( real(pi) * (g%rNorm(g%wl(w)%i)+1)/2  ) ) * g%normFacR 
            else
                kr = g%wl(w)%n * g%normFacR
            endif
        else
            kr = g%wl(w)%n * g%normFacR
        endif
        
        if(chebyshevY) then
            if(g%wl(w)%m>1) then
                kz = g%wl(w)%m / sqrt ( sin ( real(pi) * (g%zNorm(g%wl(w)%j)+1)/2 ) ) * g%normFacZ 
            else
                kz = g%wl(w)%m * g%normFacZ
            endif
        else
            kz = g%wl(w)%m * g%normFacZ
        endif
        
#if __sigma__ == 2

        sigmaHere = 0

        do s=1,nSpec

            if(coldIons.and.s>1) then
                HotSpecies = .false. 
            else
                HotSpecies = iSigma
            endif

            hotPlasma:& 
            if (HotSpecies .and. (.not. g%isMetal(g%wl(w)%iPt)) ) then        

                kVec_stix = matMul( g%U_RTZ_to_ABb(g%wl(w)%iPt,:,:), &
                    (/ kr, g%kPhi(g%wl(w)%i), kz /) ) 
                sigma_tmp = sigmaHot_maxwellian &
                    ( mSpec(s), &
                    g%ktSpec(g%wl(w)%iPt,s), &
                    g%omgc(g%wl(w)%iPt,s), &
                    g%omgp2(g%wl(w)%iPt,s), &
                    kVec_stix, g%R(g%wl(w)%i), &
                    omgrf, k0, &
                    g%k_cutoff, s, &
                    g%sinTh(g%wl(w)%iPt), &
                    g%bPol(g%wl(w)%iPt), g%bMag(g%wl(w)%iPt), &
                    g%gradPrlB(g%wl(w)%iPt), &
                    g%nuOmg(g%wl(w)%iPt,s) )

            endif hotPlasma

            coldPlasma: &
            if ((.not.HotSpecies) .and. (.not. g%isMetal(g%wl(w)%iPt)) ) then 

                sigmaIn_cold = spatialSigmaInput_cold( &
                    g%omgc(g%wl(w)%iPt,s), &
                    g%omgp2(g%wl(w)%iPt,s), &
                    omgrf, &
                    g%nuOmg(g%wl(w)%iPt,s) )

                sigma_tmp = sigmaCold_stix ( sigmaIn_cold )

#if PRINT_SIGMA==1
                if(g%wl(w)%i==5)then 

                    R_ = transpose(g%U_RTZ_to_ABb(g%wl(w)%iPt,:,:))
                    sigma_tmp2 = sigma_tmp
                    sigma_tmp3 = matmul(R_,matmul(sigma_tmp,transpose(R_)))
                    write(*,*) 'Species: ', s
                    write(*,*) 'r: ', g%R(g%wl(w)%i)

                    write(*,*) 'abp :'
                    write(*,'(3(e11.5,2x,e11.5,5x))') sigma_tmp2(:,1)
                    write(*,'(3(e11.5,2x,e11.5,5x))') sigma_tmp2(:,2)
                    write(*,'(3(e11.5,2x,e11.5,5x))') sigma_tmp2(:,3)

                    write(*,*) 'rtz :'
                    write(*,'(3(e11.5,2x,e11.5,5x))') sigma_tmp3(:,1)
                    write(*,'(3(e11.5,2x,e11.5,5x))') sigma_tmp3(:,2)
                    write(*,'(3(e11.5,2x,e11.5,5x))') sigma_tmp3(:,3)
                    write(*,*)

                endif
#endif
            endif coldPlasma

            ! Sum sigma over species
            sigmaHere = sigmaHere + sigma_tmp

        enddo

 
#if __noU__==1
        ! Rotate sigma from alp,bet,prl to r,t,z
        R_ = transpose(g%U_RTZ_to_ABb(g%wl(w)%iPt,:,:))
        sigmaHere = matmul(R_,matmul(sigmaHere,transpose(R_)))
        write(*,*) sigmaHere
#endif

        sigAlpAlp = sigmaHere(1,1)
        sigAlpBet = sigmaHere(1,2)
        sigAlpPrl = sigmaHere(1,3)
                         
        sigBetAlp = sigmaHere(2,1)
        sigBetBet = sigmaHere(2,2)
        sigBetPrl = sigmaHere(2,3)
                         
        sigPrlAlp = sigmaHere(3,1)
        sigPrlBet = sigmaHere(3,2)
        sigPrlPrl = sigmaHere(3,3)

#else 
! __sigma__ == 2

        ! Sum sigma over species
        ! ----------------------
        
        sigAlpAlp = sum ( sigma_write(g%wl(w)%i,g%wl(w)%j,g%wl(w)%n,g%wl(w)%m,1,1,:) )
        sigAlpBet = sum ( sigma_write(g%wl(w)%i,g%wl(w)%j,g%wl(w)%n,g%wl(w)%m,1,2,:) )
        sigAlpPrl = sum ( sigma_write(g%wl(w)%i,g%wl(w)%j,g%wl(w)%n,g%wl(w)%m,1,3,:) )
                                                                    
        sigBetAlp = sum ( sigma_write(g%wl(w)%i,g%wl(w)%j,g%wl(w)%n,g%wl(w)%m,2,1,:) )
        sigBetBet = sum ( sigma_write(g%wl(w)%i,g%wl(w)%j,g%wl(w)%n,g%wl(w)%m,2,2,:) )
        sigBetPrl = sum ( sigma_write(g%wl(w)%i,g%wl(w)%j,g%wl(w)%n,g%wl(w)%m,2,3,:) )
                                                                    
        sigPrlAlp = sum ( sigma_write(g%wl(w)%i,g%wl(w)%j,g%wl(w)%n,g%wl(w)%m,3,1,:) )
        sigPrlBet = sum ( sigma_write(g%wl(w)%i,g%wl(w)%j,g%wl(w)%n,g%wl(w)%m,3,2,:) )
        sigPrlPrl = sum ( sigma_write(g%wl(w)%i,g%wl(w)%j,g%wl(w)%n,g%wl(w)%m,3,3,:) )
 
#endif 
! __sigma__ == 2

        ! Metal
        ! -----

        if (g%isMetal(g%wl(w)%iPt)) then 

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

        ! This will zero the plasma current in regions
        ! where we want to specify it manually in the
        ! RHS.

        if(ZeroJp)then 
            if(R>=ZeroJp_rMin &
                    .and.R<=ZeroJp_rMax &
                    .and.z>=ZeroJp_zMin &
                    .and.z<=ZeroJp_zMax) then

                sigAlpAlp = 0
                sigAlpBet = 0
                sigAlpPrl = 0 
                        
                sigBetAlp = 0 
                sigBetBet = 0
                sigBetPrl = 0 
                        
                sigPrlAlp = 0 
                sigPrlBet = 0 
                sigPrlPrl = 0

            endif
        endif

        if(useJpFromFile)then

            sigAlpAlp = 0
            sigAlpBet = 0
            sigAlpPrl = 0
                     
            sigBetAlp = 0
            sigBetBet = 0
            sigBetPrl = 0
                     
            sigPrlAlp = 0
            sigPrlBet = 0
            sigPrlPrl = 0

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

        ! These ARE used below bitch.
        d%n = g%wl(w)%n
        d%m = g%wl(w)%m
        ! FIX: Also, these may require some interpolation.
        ! But ONLY for cos and chebyshec basis functions.
        ! i.e., OK, for NOW.
        d%xNorm = g%rNorm(g%wl(w)%i)
        d%yNorm = g%zNorm(g%wl(w)%j)
        d%normFacX = g%normFacR
        d%normFacY = g%normFacZ

#if __noU__ == 1
        Urr = 1
        Urt = 0
        Urz = 0
        
        Utr = 0
        Utt = 1
        Utz = 0
              
        Uzr = 0
        Uzt = 0
        Uzz = 1
#else
        Urr = g%Urr(g%wl(w)%iPt)
        Urt = g%Urt(g%wl(w)%iPt)
        Urz = g%Urz(g%wl(w)%iPt)
        
        Utr = g%Utr(g%wl(w)%iPt)
        Utt = g%Utt(g%wl(w)%iPt)
        Utz = g%Utz(g%wl(w)%iPt)
                    
        Uzr = g%Uzr(g%wl(w)%iPt)
        Uzt = g%Uzt(g%wl(w)%iPt)
        Uzz = g%Uzz(g%wl(w)%iPt)
#endif
       
#if __noU__==1
        drUrr = 0 
        drUrt = 0 
        drUrz = 0 

        drUtr = 0 
        drUtt = 0 
        drUtz = 0 

        drUzr = 0 
        drUzt = 0 
        drUzz = 0 
#else
        drUrr = g%drUrr(g%wl(w)%iPt)
        drUrt = g%drUrt(g%wl(w)%iPt)
        drUrz = g%drUrz(g%wl(w)%iPt)

        drUtr = g%drUtr(g%wl(w)%iPt)
        drUtt = g%drUtt(g%wl(w)%iPt)
        drUtz = g%drUtz(g%wl(w)%iPt)

        drUzr = g%drUzr(g%wl(w)%iPt)
        drUzt = g%drUzt(g%wl(w)%iPt)
        drUzz = g%drUzz(g%wl(w)%iPt)
#endif
        !drUrr = biLinearInterp(r,z,g,g%drUrr)
        !drUrt = biLinearInterp(r,z,g,g%drUrt)
        !drUrz = biLinearInterp(r,z,g,g%drUrz)

        !drUtr = biLinearInterp(r,z,g,g%drUtr)
        !drUtt = biLinearInterp(r,z,g,g%drUtt)
        !drUtz = biLinearInterp(r,z,g,g%drUtz)

        !drUzr = biLinearInterp(r,z,g,g%drUzr)
        !drUzt = biLinearInterp(r,z,g,g%drUzt)
        !drUzz = biLinearInterp(r,z,g,g%drUzz)
#if __noU__==1
        dzUrr = 0 
        dzUrt = 0 
        dzUrz = 0 

        dzUtr = 0 
        dzUtt = 0 
        dzUtz = 0 

        dzUzr = 0 
        dzUzt = 0 
        dzUzz = 0 
#else
        dzUrr = g%dzUrr(g%wl(w)%iPt)
        dzUrt = g%dzUrt(g%wl(w)%iPt)
        dzUrz = g%dzUrz(g%wl(w)%iPt)

        dzUtr = g%dzUtr(g%wl(w)%iPt)
        dzUtt = g%dzUtt(g%wl(w)%iPt)
        dzUtz = g%dzUtz(g%wl(w)%iPt)

        dzUzr = g%dzUzr(g%wl(w)%iPt)
        dzUzt = g%dzUzt(g%wl(w)%iPt)
        dzUzz = g%dzUzz(g%wl(w)%iPt)
#endif

        !dzUrr = biLinearInterp(r,z,g,g%dzUrr)
        !dzUrt = biLinearInterp(r,z,g,g%dzUrt)
        !dzUrz = biLinearInterp(r,z,g,g%dzUrz)

        !dzUtr = biLinearInterp(r,z,g,g%dzUtr)
        !dzUtt = biLinearInterp(r,z,g,g%dzUtt)
        !dzUtz = biLinearInterp(r,z,g,g%dzUtz)

        !dzUzr = biLinearInterp(r,z,g,g%dzUzr)
        !dzUzt = biLinearInterp(r,z,g,g%dzUzt)
        !dzUzz = biLinearInterp(r,z,g,g%dzUzz)
#if __noU__==1
        drrUrr = 0 
        drrUrt = 0 
        drrUrz = 0 

        drrUtr = 0 
        drrUtt = 0 
        drrUtz = 0 

        drrUzr = 0 
        drrUzt = 0 
        drrUzz = 0 
#else
        drrUrr = g%drrUrr(g%wl(w)%iPt)
        drrUrt = g%drrUrt(g%wl(w)%iPt)
        drrUrz = g%drrUrz(g%wl(w)%iPt)

        drrUtr = g%drrUtr(g%wl(w)%iPt)
        drrUtt = g%drrUtt(g%wl(w)%iPt)
        drrUtz = g%drrUtz(g%wl(w)%iPt)

        drrUzr = g%drrUzr(g%wl(w)%iPt)
        drrUzt = g%drrUzt(g%wl(w)%iPt)
        drrUzz = g%drrUzz(g%wl(w)%iPt)
#endif

        !drrUrr = biLinearInterp(r,z,g,g%drrUrr)
        !drrUrt = biLinearInterp(r,z,g,g%drrUrt)
        !drrUrz = biLinearInterp(r,z,g,g%drrUrz)

        !drrUtr = biLinearInterp(r,z,g,g%drrUtr)
        !drrUtt = biLinearInterp(r,z,g,g%drrUtt)
        !drrUtz = biLinearInterp(r,z,g,g%drrUtz)
        !
        !drrUzr = biLinearInterp(r,z,g,g%drrUzr)
        !drrUzt = biLinearInterp(r,z,g,g%drrUzt)
        !drrUzz = biLinearInterp(r,z,g,g%drrUzz)
#if __noU__==1 
        dzzUrr = 0 
        dzzUrt = 0 
        dzzUrz = 0 
        
        dzzUtr = 0 
        dzzUtt = 0 
        dzzUtz = 0 

        dzzUzr = 0 
        dzzUzt = 0 
        dzzUzz = 0 
#else
        dzzUrr = g%dzzUrr(g%wl(w)%iPt)
        dzzUrt = g%dzzUrt(g%wl(w)%iPt)
        dzzUrz = g%dzzUrz(g%wl(w)%iPt)
        
        dzzUtr = g%dzzUtr(g%wl(w)%iPt)
        dzzUtt = g%dzzUtt(g%wl(w)%iPt)
        dzzUtz = g%dzzUtz(g%wl(w)%iPt)

        dzzUzr = g%dzzUzr(g%wl(w)%iPt)
        dzzUzt = g%dzzUzt(g%wl(w)%iPt)
        dzzUzz = g%dzzUzz(g%wl(w)%iPt)
#endif
        !dzzUrr = biLinearInterp(r,z,g,g%dzzUrr)
        !dzzUrt = biLinearInterp(r,z,g,g%dzzUrt)
        !dzzUrz = biLinearInterp(r,z,g,g%dzzUrz)

        !dzzUtr = biLinearInterp(r,z,g,g%dzzUtr)
        !dzzUtt = biLinearInterp(r,z,g,g%dzzUtt)
        !dzzUtz = biLinearInterp(r,z,g,g%dzzUtz)
        !
        !dzzUzr = biLinearInterp(r,z,g,g%dzzUzr)
        !dzzUzt = biLinearInterp(r,z,g,g%dzzUzt)
        !dzzUzz = biLinearInterp(r,z,g,g%dzzUzz)
#if __noU__==1 
        drzUrr = 0 
        drzUrt = 0 
        drzUrz = 0 
        
        drzUtr = 0 
        drzUtt = 0 
        drzUtz = 0 

        drzUzr = 0 
        drzUzt = 0 
        drzUzz = 0 
#else
        drzUrr = g%drzUrr(g%wl(w)%iPt)
        drzUrt = g%drzUrt(g%wl(w)%iPt)
        drzUrz = g%drzUrz(g%wl(w)%iPt)
        
        drzUtr = g%drzUtr(g%wl(w)%iPt)
        drzUtt = g%drzUtt(g%wl(w)%iPt)
        drzUtz = g%drzUtz(g%wl(w)%iPt)

        drzUzr = g%drzUzr(g%wl(w)%iPt)
        drzUzt = g%drzUzt(g%wl(w)%iPt)
        drzUzz = g%drzUzz(g%wl(w)%iPt)
#endif

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

        mat_t_alp = ((zi*kt*Urr)/r**2 - Urt/r**2 + k0**2*KAlpAlp*Urt + &
           k0**2*KBetAlp*Utt + k0**2*KPrlAlp*Uzt - &
           (zi*kt*Urz*dzBfn_bfn(d))/r + &
           (2*dzBfn_bfn(d)*dzUrt) - &
           (zi*kt*dzUrz)/r + (Urt*dzzBfn_bfn(d)) + &
           dzzUrt - (zi*kt*Urr*drBfn_bfn(d))/r + &
           (Urt*drBfn_bfn(d))/r - (zi*kt*drUrr)/r + &
           drUrt/r + (2*drBfn_bfn(d)*drUrt) + &
           (Urt*drrBfn_bfn(d)) + drrUrt)

        mat_t_bet = (k0**2*KAlpBet*Urt + (zi*kt*Utr)/r**2 - Utt/r**2 + &
           k0**2*KBetBet*Utt + k0**2*KPrlBet*Uzt - &
           (zi*kt*Utz*dzBfn_bfn(d))/r + &
           (2*dzBfn_bfn(d)*dzUtt) - &
           (zi*kt*dzUtz)/r + (Utt*dzzBfn_bfn(d)) + &
           dzzUtt - (zi*kt*Utr*drBfn_bfn(d))/r + &
           (Utt*drBfn_bfn(d))/r - (zi*kt*drUtr)/r + &
           drUtt/r + (2*drBfn_bfn(d)*drUtt) + &
           (Utt*drrBfn_bfn(d)) + drrUtt)

        mat_t_prl = (k0**2*KAlpPrl*Urt + k0**2*KBetPrl*Utt + (zi*kt*Uzr)/r**2 - &
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

    !endif interior

    get3by3Block(1,1) = mat_r_alp
    get3by3Block(1,2) = mat_r_bet
    get3by3Block(1,3) = mat_r_prl

    get3by3Block(2,1) = mat_t_alp
    get3by3Block(2,2) = mat_t_bet
    get3by3Block(2,3) = mat_t_prl

    get3by3Block(3,1) = mat_z_alp
    get3by3Block(3,2) = mat_z_bet
    get3by3Block(3,3) = mat_z_prl

end function get3by3Block

end module getMatElements
