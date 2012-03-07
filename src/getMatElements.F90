module getMatElements

contains

function get3by3Block( g, w)

    use grid
    use aorsa2din_mod, only: &
       chebyshevX, chebyshevY, iSigma, nSpec, nPhi
    use constants
    use profiles, only: k0, omgrf, mSpec
    use sigma_mod

    implicit none

    type(gridBlock), intent(in) :: g
    integer, intent(in) :: w
    
    complex :: get3by3Block(3,3)

    complex(kind=dbl) :: &
        sigAlpAlp, sigAlpBet, sigAlpPrl, &
        sigBetAlp, sigBetBet, sigBetPrl, &
        sigPrlAlp, sigPrlBet, sigPrlPrl

    type(dBfnArg) :: d
    complex(kind=dbl) :: sigma_tmp(3,3), sigma_tmp_neg(3,3), sigmaHere(3,3)

    real :: kr, kt, kz, r, z, kVec_stix(3)

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

    z   = g%z(g%wl(w)%j)
    r   = g%R(g%wl(w)%i)
    kt  = nPhi!g%kPhi(i)

    interior: &
    if(g%label(g%wl(w)%i,g%wl(w)%j)==0)then
    
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
                kr = g%wl(w)%n / sqrt ( sin ( pi * (g%rNorm(g%wl(w)%i)+1)/2  ) ) * g%normFacR 
            else
                kr = g%wl(w)%n * g%normFacR
            endif
        else
            kr = g%wl(w)%n * g%normFacR
        endif
        
        if(chebyshevY) then
            if(g%wl(w)%m>1) then
                kz = g%wl(w)%m / sqrt ( sin ( pi * (g%zNorm(g%wl(w)%j)+1)/2 ) ) * g%normFacZ 
            else
                kz = g%wl(w)%m * g%normFacZ
            endif
        else
            kz = g%wl(w)%m * g%normFacZ
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

#if __sigma__ == 2

        sigmaHere = 0

        do s=1,nSpec

            hotPlasma:& 
            if (iSigma==1 .and. (.not. g%isMetal(g%wl(w)%i,g%wl(w)%j)) ) then        

                kVec_stix = matMul( g%U_RTZ_to_ABb(g%wl(w)%i,g%wl(w)%j,:,:), &
                    (/ kr, g%kPhi(g%wl(w)%i), kz /) ) 

                sigma_tmp = sigmaHot_maxwellian &
                    ( mSpec(s), &
                    g%ktSpec(g%wl(w)%i,g%wl(w)%j,s), &
                    g%omgc(g%wl(w)%i,g%wl(w)%j,s), &
                    g%omgp2(g%wl(w)%i,g%wl(w)%j,s), &
                    kVec_stix, g%R(g%wl(w)%i), &
                    omgrf, k0, &
                    g%k_cutoff, s, &
                    g%sinTh(g%wl(w)%i,g%wl(w)%j), &
                    g%bPol(g%wl(w)%i,g%wl(w)%j), g%bMag(g%wl(w)%i,g%wl(w)%j), &
                    g%gradPrlB(g%wl(w)%i,g%wl(w)%j), &
                    g%nuOmg(g%wl(w)%i,g%wl(w)%j) )

            endif hotPlasma

            coldPlasma: &
            if (iSigma==0 .and. (.not. g%isMetal(g%wl(w)%i,g%wl(w)%j)) ) then 

                sigmaIn_cold = spatialSigmaInput_cold( &
                    g%omgc(g%wl(w)%i,g%wl(w)%j,s), &
                    g%omgp2(g%wl(w)%i,g%wl(w)%j,s), &
                    omgrf, &
                    g%nuOmg(g%wl(w)%i,g%wl(w)%j) )

                !sigma_tmp = sigmaCold_stix &
                !    ( g%omgc(g%wl(w)%i,g%wl(w)%j,s), &
                !    g%omgp2(g%wl(w)%i,g%wl(w)%j,s), omgrf, &
                !    g%nuOmg(g%wl(w)%i,g%wl(w)%j) )
                sigma_tmp = sigmaCold_stix ( sigmaIn_cold )

            endif coldPlasma

            ! Sum sigma over species
            sigmaHere = sigmaHere + sigma_tmp

        enddo

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

        if (g%isMetal(g%wl(w)%i,g%wl(w)%j)) then 

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

        d%n = g%wl(w)%n
        d%m = g%wl(w)%m
        d%xNorm = g%rNorm(g%wl(w)%i)
        d%yNorm = g%zNorm(g%wl(w)%j)
        d%normFacX = g%normFacR
        d%normFacY = g%normFacZ

        Urr = g%Urr(g%wl(w)%i,g%wl(w)%j)
        Urt = g%Urt(g%wl(w)%i,g%wl(w)%j)
        Urz = g%Urz(g%wl(w)%i,g%wl(w)%j)
                                  
        Utr = g%Utr(g%wl(w)%i,g%wl(w)%j)
        Utt = g%Utt(g%wl(w)%i,g%wl(w)%j)
        Utz = g%Utz(g%wl(w)%i,g%wl(w)%j)
                                  
        Uzr = g%Uzr(g%wl(w)%i,g%wl(w)%j)
        Uzt = g%Uzt(g%wl(w)%i,g%wl(w)%j)
        Uzz = g%Uzz(g%wl(w)%i,g%wl(w)%j)


        drUrr = g%drUrr(g%wl(w)%i,g%wl(w)%j)
        drUrt = g%drUrt(g%wl(w)%i,g%wl(w)%j)
        drUrz = g%drUrz(g%wl(w)%i,g%wl(w)%j)

        drUtr = g%drUtr(g%wl(w)%i,g%wl(w)%j)
        drUtt = g%drUtt(g%wl(w)%i,g%wl(w)%j)
        drUtz = g%drUtz(g%wl(w)%i,g%wl(w)%j)

        drUzr = g%drUzr(g%wl(w)%i,g%wl(w)%j)
        drUzt = g%drUzt(g%wl(w)%i,g%wl(w)%j)
        drUzz = g%drUzz(g%wl(w)%i,g%wl(w)%j)

        dzUrr = g%dzUrr(g%wl(w)%i,g%wl(w)%j)
        dzUrt = g%dzUrt(g%wl(w)%i,g%wl(w)%j)
        dzUrz = g%dzUrz(g%wl(w)%i,g%wl(w)%j)

        dzUtr = g%dzUtr(g%wl(w)%i,g%wl(w)%j)
        dzUtt = g%dzUtt(g%wl(w)%i,g%wl(w)%j)
        dzUtz = g%dzUtz(g%wl(w)%i,g%wl(w)%j)

        dzUzr = g%dzUzr(g%wl(w)%i,g%wl(w)%j)
        dzUzt = g%dzUzt(g%wl(w)%i,g%wl(w)%j)
        dzUzz = g%dzUzz(g%wl(w)%i,g%wl(w)%j)


        drrUrr = g%drrUrr(g%wl(w)%i,g%wl(w)%j)
        drrUrt = g%drrUrt(g%wl(w)%i,g%wl(w)%j)
        drrUrz = g%drrUrz(g%wl(w)%i,g%wl(w)%j)
        
        drrUtr = g%drrUtr(g%wl(w)%i,g%wl(w)%j)
        drrUtt = g%drrUtt(g%wl(w)%i,g%wl(w)%j)
        drrUtz = g%drrUtz(g%wl(w)%i,g%wl(w)%j)

        drrUzr = g%drrUzr(g%wl(w)%i,g%wl(w)%j)
        drrUzt = g%drrUzt(g%wl(w)%i,g%wl(w)%j)
        drrUzz = g%drrUzz(g%wl(w)%i,g%wl(w)%j)

        dzzUrr = g%dzzUrr(g%wl(w)%i,g%wl(w)%j)
        dzzUrt = g%dzzUrt(g%wl(w)%i,g%wl(w)%j)
        dzzUrz = g%dzzUrz(g%wl(w)%i,g%wl(w)%j)
        
        dzzUtr = g%dzzUtr(g%wl(w)%i,g%wl(w)%j)
        dzzUtt = g%dzzUtt(g%wl(w)%i,g%wl(w)%j)
        dzzUtz = g%dzzUtz(g%wl(w)%i,g%wl(w)%j)

        dzzUzr = g%dzzUzr(g%wl(w)%i,g%wl(w)%j)
        dzzUzt = g%dzzUzt(g%wl(w)%i,g%wl(w)%j)
        dzzUzz = g%dzzUzz(g%wl(w)%i,g%wl(w)%j)

        drzUrr = g%drzUrr(g%wl(w)%i,g%wl(w)%j)
        drzUrt = g%drzUrt(g%wl(w)%i,g%wl(w)%j)
        drzUrz = g%drzUrz(g%wl(w)%i,g%wl(w)%j)
        
        drzUtr = g%drzUtr(g%wl(w)%i,g%wl(w)%j)
        drzUtt = g%drzUtt(g%wl(w)%i,g%wl(w)%j)
        drzUtz = g%drzUtz(g%wl(w)%i,g%wl(w)%j)

        drzUzr = g%drzUzr(g%wl(w)%i,g%wl(w)%j)
        drzUzt = g%drzUzt(g%wl(w)%i,g%wl(w)%j)
        drzUzz = g%drzUzz(g%wl(w)%i,g%wl(w)%j)


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

    endif interior

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
