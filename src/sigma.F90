module sigma_mod

contains

    function sigmaHot_maxwellian( &
        xm, kt, omgc, omgp2, &
        kVec, capr, &
        omgrf, k0, &
        k_cutoff, specNo, &
        sinTh, bPol, bMod, gradPrlB, nuOmg )

        use constants
        use zfunction_mod
        use aorsa2din_mod, &
        only: upshift, nzfun, damping, lmax, &
            nPhi, delta0, toroidalBroadening, &
            kPrlEffLimit
        use bessel_mod

        ! Calculate the hot Maxwellian sigma at a single
        ! point in space in Stix frame including toriodal
        ! broadening using either Smithe or Brambilla
        ! methods. Stix frame is (alp,bet,prl).
        !------------------------------------------------

        implicit none

        !include "gptl.inc"

        real, intent(in) :: sinTh, bPol, bMod, gradPrlB
        complex(kind=dbl) :: sigmaHot_maxwellian(3,3)
        integer, intent(in) :: specNo
        real, intent(in) :: k_cutOff, nuOmg
        real, intent(in) :: kVec(:), capr

        real(kind=dbl), intent(in) :: omgc, omgp2, omgrf
        real(kind=dbl), intent(in) :: kt, k0

        integer :: labs, i, j
        integer(kind=long) :: l
        real :: kPerp, kPrl, xm, kPrlTmp
        real :: kPrlEff, fgam, y0, y
        real(kind=dbl) :: sgn_kPrl
        real :: descrim
        real :: nu_coll
        real :: vTh
        real(kind=dbl) :: alpha_smithe(-lmax:lmax), gamma_smithe(-lmax:lmax)
        real(kind=dbl) :: gamma_brambilla(-lmax:lmax)
        real :: rhol
        real :: kAlp, kBet, rgamma
        real(kind=dbl) :: kr, step
        complex(kind=dbl) :: omgrfc
        complex(kind=dbl), allocatable :: Z(:), ZPrime(:), zetaZPrime(:)
        complex(kind=dbl) :: sig0, sig1, sig2, sig3, sig4, sig5
        complex(kind=dbl), dimension(-lmax:lmax) :: &
            sig0l, sig1l, sig2l, sig3l, sig4l, sig5l
        integer, parameter :: lmaxdim = 99
        complex(kind=dbl), allocatable, dimension(:) :: &
          expBesselI, expBesselIPrime, expBesselIOverGam
        complex(kind=dbl) :: zetal(-lmax:lmax)
        complex :: zieps0 
        complex(kind=dbl) :: P, PPrime 
        complex(kind=dbl) ::gamma_
        real :: dR, Lpar
        integer :: stat
        real :: sinPh
        real :: cosPsi, sinPsi

        zetal = 0

        zieps0 = zi * eps0
        vTh = sqrt(2. * kt / xm)
        rhol = vTh / omgc
        omgrfc = omgrf * (1d0 + zi * nuOmg)

        ! k in Stix frame
        ! ---------------

        kAlp = kVec(1)
        kBet = kVec(2)
        kPrl = kVec(3)

        kPerp = sqrt(kAlp**2 + kBet**2)

        if (kPrl  <= 1d-2) kPrl  = 1d-2
        if (kPerp == 0.0) kPerp = 1.0e-03

        sgn_kPrl    = sign ( 1.0, kPrl )


        ! Calculate Z function argument
        ! -----------------------------

        nu_coll =  .01 * omgrf

        do l = -lmax, lmax


            ! Two toroidal broadening terms from 
            ! Smithe et al., PRL Vol 60, No. 9, Feb 1988
            ! ------------------------------------------

            !alpha_smithe(l) = l * omgc / ( kPrl * abs( kPrl ) * gradPrlB * vTh )
            !gamma_smithe(l) = nu_coll / (abs(kPrl) * vTh)


            ! Brambilla has a single approximate broadening term
            ! Phys. Lett. A, 188, 376-383, 1994
            ! ---------------------------------
    
            sinPh = bPol / bMod
            gamma_brambilla(l) = omgrf / ( 2.0 * kPrl**2 * capR * vTh ) * abs ( sinTh * sinPh ) 


            ! ions only 
            ! ---------

            if(specNo==1) then 

                alpha_smithe(l) = 0
                gamma_brambilla(l) = 0

            endif


            ! Create and effective kPrl that includes toroidal broadening
            ! ----------------------------------------------------------

            if(gamma_brambilla(l)>0 .and. toroidalBroadening)  then

                kPrlEff = kPrl * ( sqrt ( 1 + 4 * gamma_brambilla(l) ) - 1 ) &
                                / ( 2 * gamma_brambilla(l))
            else

                kPrlEff = kPrl

            endif

            if (kPrlEff  <= kPrlEffLimit) kPrlEff  = kPrlEffLimit

            zetal(l) = (omgrfc - l * omgc) / ( abs( kPrlEff ) * vTh) 

        enddo


        ! Exp(-Gamma) * Il etc 
        ! --------------------

        gamma_ = 0.5 * kPerp**2 * rhol**2
        
        allocate ( &
            expBesselI(0:lMax), &
            expBesselIPrime(0:lMax), &
            expBesselIOverGam(0:lMax) )

        if ( real (gamma_) >= 1.0e-08 ) then
            call besIExp ( gamma_, lMax, expBesselI, expBesselIPrime, expBesselIOverGam )
        else
            call bes_expand ( gamma_, lMax, expBesselI, expBesselIPrime, expBesselIOverGam )
        endif

        
        ! Calculate Z function
        ! --------------------

        allocate ( Z(-lMax:lMax), ZPrime(-lMax:lMax), zetaZPrime(-lMax:lMax) )

        call z_approx_dlg ( zetal, Z, ZPrime )
        !call z_approx(sgn_kprl, zetal, gamma_brambilla*0, z, zPrime, zetaZPrime)


        ! Calculate sigmas according to Swanson
        ! -------------------------------------

        do l = -lmax, lmax

            labs = abs(l)

            P = omgp2 / (abs(kPrlEff) * vTh) * Z(l)
            PPrime = omgp2 / (abs(kPrlEff) * vTh) * ZPrime(l)

            sig0l(l) = - zieps0 * rhol**2 * (expBesselI(labs) - expBesselIPrime(labs)) * P 
            sig1l(l) = - zieps0 * l**2 * expBesselIOverGam(labs) * P
            sig2l(l) = - eps0 * l * (expBesselI(labs) - expBesselIPrime(labs)) * P
            sig3l(l) =   zieps0 * expBesselI(labs) * zetal(l) * PPrime
            sig4l(l) =   0.5d0 * zieps0 * rhol * l * expBesselIOverGam(labs) * PPrime
            sig5l(l) =   0.5d0 * eps0 * rhol * (expBesselI(labs) - expBesselIPrime(labs)) * PPrime

        enddo

        sig0    = sum ( sig0l )
        sig1    = sum ( sig1l )
        sig2    = sum ( sig2l )
        sig3    = sum ( sig3l )
        sig4    = sum ( sig4l )
        sig5    = sum ( sig5l )

        deallocate ( expBesselI, expBesselIPrime, expBesselIOverGam )


        ! Add numerical damping for Bernstein wave (delta0 ~ 1e-4)
        ! --------------------------------------------------------

        sig1 = sig1 + delta0 * eps0 * omgrfc * kPerp**2 / k0**2
        sig3 = sig3 + delta0 * eps0 * omgrfc * kPerp**2 / k0**2


        ! electron damping
        ! ----------------

        if (specNo==1) then

            kr = kPerp / k_cutOff
            step = damping * kr**16 / (1. + kr**16)
            sig3 = sig3 * (1.0 + step)

        endif

        
        ! Swanson's sigma tensor in the 
        ! (alp,bet,prl) system where alp and bet
        ! are NOT related to kPerp as in the Stix
        ! case where Psi=0. See pg 176 of Swanson.
        ! ----------------------------------------

        cosPsi = kAlp / kPerp
        sinPsi = kBet / kPerp

        sigmaHot_maxwellian(1,1) = sig1 + sig0 * sinPsi**2
        sigmaHot_maxwellian(1,2) = sig2 - sig0 * cosPsi * sinPsi
        sigmaHot_maxwellian(1,3) = sig4 * cosPsi + sig5 * sinPsi

        sigmaHot_maxwellian(2,1) = - sig2 - sig0 * cosPsi * sinPsi
        sigmaHot_maxwellian(2,2) =   sig1 + sig0 * cosPsi**2
        sigmaHot_maxwellian(2,3) =   sig4 * sinPsi - sig5 * cosPsi

        sigmaHot_maxwellian(3,1) = sig4 * cosPsi - sig5 * sinPsi
        sigmaHot_maxwellian(3,2) = sig4 * sinPsi + sig5 * cosPsi
        sigmaHot_maxwellian(3,3) = sig3

        return

    end function sigmaHot_maxwellian


    function sigmaCold_stix &
        ( omgC, omgP2, omgRF, nuOmg )

        use constants 

        ! This routine calculates sigma_cold in the Stix frame
        ! (alp,bet,prl)
        ! ----------------------------------------------------

        implicit none

        real(kind=dbl), intent(in) :: omgc, omgp2, omgrf
        real, intent(in) :: nuOmg

        complex :: omgrfc
        complex :: sig1, sig2, sig3
        complex :: sigmaCold_stix(3,3)
        complex :: zieps0

        zieps0 = zi * eps0
        omgRFc = omgRF * (1.0 + zi * nuOmg)

        sig1 = zieps0 * omgRFc * omgP2 / (omgRFc**2 - omgC**2)
        sig2 = eps0 * omgC   * omgP2 / (omgC**2 - omgRFc**2)
        sig3 = zieps0 * omgp2 / omgRFc

        sigmaCold_stix(1,1) = sig1 
        sigmaCold_stix(1,2) = sig2 
        sigmaCold_stix(1,3) = 0 

        sigmaCold_stix(2,1) = -sig2 
        sigmaCold_stix(2,2) = sig1 
        sigmaCold_stix(2,3) = 0 

        sigmaCold_stix(3,1) = 0
        sigmaCold_stix(3,2) = 0 
        sigmaCold_stix(3,3) = sig3

        return

    end function sigmaCold_stix

end module sigma_mod
