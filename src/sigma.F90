#define ERROR(string) write(*,'(/,a,a,/,a,a,a,i5)') \
    "Error: ", string, "In ",__FILE__,":", __LINE__ ;stop

module sigma

use constants

type :: spatialSigmaInput_cold

        real(kind=dbl) :: omgC
        real(kind=dbl) :: omgP2
        real :: omgrf
        real :: nuOmg

end type spatialSigmaInput_cold

contains

    function sigmaHot_maxwellian( &
        xm, kt, omgc, omgp2, &
        kVec, capr, &
        omgrf, k0, &
        k_cutoff, specNo, &
        sinTh, bPol, bMod, gradPrlB, nuOmg )

        use constants
        use zfunction_mod
        use aorsaNamelist, &
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
        integer :: l
        real :: kPerp, kPrl, xm, kPrlTmp
        real :: kPrlEff, fgam, y0, y
        real(kind=dbl) :: sgn_kPrl
        real :: descrim
        !real :: nu_coll
        real :: vTh
        real(kind=dbl) :: alpha_smithe(-lmax:lmax), gamma_smithe(-lmax:lmax)
        real(kind=dbl) :: gamma_brambilla(-lmax:lmax)
        real :: rhol
        real :: kAlp, kBet, rgamma
        real(kind=dbl) :: kr, step
        complex(kind=dbl) :: omgrfc
        complex(kind=dbl), allocatable :: Z(:), zPrime(:), zeta_zPrime(:)
        complex(kind=dbl), allocatable :: Z_swan(:), ZPrime_swan(:), ZetaZPrime_swan(:)
        complex(kind=dbl) :: sig0, sig1, sig2, sig3, sig4, sig5
        complex(kind=dbl), dimension(-lmax:lmax) :: &
            sig0l, sig1l, sig2l, sig3l, sig4l, sig5l
        integer, parameter :: lmaxdim = 99
        complex(kind=dbl), allocatable, dimension(:) :: &
          expBesselI, expBesselIPrime, expBesselIOverGam
        complex(kind=dbl) :: zetal(-lmax:lmax),Zeta_l_swan(-lmax:lmax)
        complex :: zieps0 
        complex(kind=dbl) :: P, PPrime 
        complex(kind=dbl) ::gamma_
        complex(kind=dbl) ::lambda_swan
        real :: dR, Lpar
        integer :: stat
        real :: sinPh
        real :: cosPsi, sinPsi

        !real :: e_swan
        !complex(kind=DBL) :: factor1,factor2
        !complex(kind=DBL) :: K0_swan,K1_swan,K2_swan,K3_swan,K4_swan,K5_swan
        !complex(kind=dbl) :: SigmaHotMaxwellian_swan(3,3),EpsilonHotMaxwellian_swan(3,3)
        !integer :: IdentMat(3,3)
        !type(SpatialSigmaInput_cold) :: scIn
        !complex :: sc(3,3)

        !zetal = 0
        !Zeta_l_swan = 0

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

        if (kPrl  == 0.0) kPrl  = 1.0e-08
        if (kPerp == 0.0) kPerp = 1.0e-08

        sgn_kPrl    = sign ( 1.0, kPrl )
        

        ! Calculate Z function argument
        ! -----------------------------

        !nu_coll =  .01 * omgrf

        do l = -lmax, lmax


            ! Two toroidal broadening terms from 
            ! Smithe et al., PRL Vol 60, No. 9, Feb 1988
            ! ------------------------------------------

            !alpha_smithe(l) = l * omgc / ( kPrl * abs( kPrl ) * gradPrlB * vTh )
            !gamma_smithe(l) = nu_coll / (abs(kPrl) * vTh)


            ! Effective kPrl for ions 
            ! -----------------------

            if(toroidalBroadening .and. abs(sinTh)>0 .and. specNo>1)  then

                ! Brambilla has a single approximate broadening term
                ! Phys. Lett. A, 188, 376-383, 1994
                ! ---------------------------------
    
                sinPh = bPol / bMod
                gamma_brambilla(l) = omgrf / ( 2.0 * kPrl**2 * capR * vTh ) * abs ( sinTh * sinPh ) 

                if(gamma_brambilla(l)<=0)then
                    write(*,*) kPrl,capR,vTh,gamma_brambilla(l),sinTh,sinPh
                    ERROR("gamma_brambilla<0")
                endif

                if(kPrl==0.0)then
                    kPrlEff = sqrt ( omgrf / ( 2*capR*vTh) * sinPh )
                else
                    kPrlEff = kPrl*( sqrt(1 + 4*gamma_brambilla(l)) - 1 ) / (2*gamma_brambilla(l))
                endif

            else
                kPrlEff = kPrl
            endif

            if(abs(kPrlEff)<kPrlEffLimit)then
                    kPrlEff = sign(kPrlEffLimit,kPrlEff)
            endif

            zetal(l)        = (omgrfc-l*omgc)/(kPrlEff*vTh) 
            !Zeta_l_swan(l)  = (omgrfc+l*omgc)/(kPrlEff*vTh)

        enddo



        ! Exp(-Gamma) * Il etc 
        ! --------------------

        gamma_ = 0.5 * kPerp**2 * rhol**2

        allocate ( &
            expBesselI(0:lMax), &
            expBesselIPrime(0:lMax), &
            expBesselIOverGam(0:lMax) )

        if ( abs(gamma_) >= 1.0e-08 ) then
            call besIExp ( gamma_, lMax, expBesselI, expBesselIPrime, expBesselIOverGam )
        else
            call bes_expand ( gamma_, lMax, expBesselI, expBesselIPrime, expBesselIOverGam )
        endif

        ! Calculate Z function
        ! --------------------

        allocate ( Z(-lMax:lMax), zPrime(-lMax:lMax), zeta_zPrime(-lMax:lMax) )
        !allocate ( Z_swan(-lMax:lMax), ZPrime_swan(-lMax:lMax), ZetaZPrime_swan(-lMax:lMax) )

        call z_approx_dlg ( sgn_kprl, zetal, z, zPrime, zeta_zPrime )
        !call z_approx_dlg ( sgn_kprl, Zeta_l_swan, Z_swan, ZPrime_swan, ZetaZPrime_swan )

        ! Calculate sigmas according to Swanson
        ! -------------------------------------

        do l = -lmax, lmax

            labs = abs(l)

            P = omgp2 / (abs(kPrlEff) * vTh) * Z(l)
            PPrime = omgp2 / (abs(kPrlEff) * vTh) * zPrime(l)

            sig0l(l) = - zieps0 * rhol**2 * (expBesselI(labs) - expBesselIPrime(labs)) * P 
            sig1l(l) = - zieps0 * l**2 * expBesselIOverGam(labs) * P
            sig2l(l) = - eps0 * l * (expBesselI(labs) - expBesselIPrime(labs)) * P
            sig3l(l) =   zieps0 * expBesselI(labs) * zetal(l) * PPrime
            sig4l(l) =   0.5d0 * zieps0 * rhol * l * expBesselIOverGam(labs) * PPrime
            sig5l(l) =   0.5d0 * eps0 * rhol * (expBesselI(labs) - expBesselIPrime(labs)) * PPrime

        enddo

        !! Try using the Swanson approach

        !factor1 = omgp2/(omgRFc*abs(kPrlEff)*vTh)
        !factor2 = kPerp*omgP2/(2*abs(kPrlEff)*omgRFc*omgC)
        !e_swan = sign(1d0,omgc) ! this is q/|q|

        !K0_swan=(0,0) 
        !K1_swan=(0,0) 
        !K2_swan=(0,0) 
        !K3_swan=(0,0) 
        !K4_swan=(0,0) 
        !K5_swan=(0,0)

        !do l=-lmax,lmax

        !    labs = abs(l)

        !    K0_swan = K0_swan + gamma_*(expBesselI(labs)-expBesselIPrime(labs))*Z_swan(l)
        !    K1_swan = K1_swan + l**2*expBesselIOverGam(labs)*Z_swan(l)
        !    K2_swan = K2_swan + l*(expBesselI(labs)-expBesselIPrime(labs))*Z_swan(l)
        !    K3_swan = K3_swan + expBesselI(labs)*Zeta_l_swan(l)*ZPrime_swan(l)
        !    K4_swan = K4_swan + l*expBesselIOverGam(labs)*ZPrime_swan(l)
        !    K5_swan = K5_swan + (expBesselI(labs)-expBesselIPrime(labs))*ZPrime_swan(l)

        !enddo

        !K0_swan = 2*factor1*K0_swan
        !K1_swan = 1+factor1*K1_swan
        !K2_swan = zi*e_swan*factor1*K2_swan
        !K3_swan = 1-factor1*K3_swan
        !K4_swan = factor2*K4_swan
        !K5_swan = zi*e_swan*factor2*K5_swan

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


        !EpsilonHotMaxwellian_swan(1,1) = K1_swan+sinPsi**2*K0_swan
        !EpsilonHotMaxwellian_swan(1,2) = K2_swan-cosPsi*sinPsi*K0_swan
        !EpsilonHotMaxwellian_swan(1,3) = cosPsi*K4_swan+sinPsi*K5_swan

        !EpsilonHotMaxwellian_swan(2,1) = -K2_swan-cosPsi*sinPsi*K0_swan
        !EpsilonHotMaxwellian_swan(2,2) = K1_swan+cosPsi**2*K0_swan 
        !EpsilonHotMaxwellian_swan(2,3) = sinPsi*K4_swan-cosPsi*K5_swan

        !EpsilonHotMaxwellian_swan(3,1) = cosPsi*K4_swan-sinPsi*K5_swan
        !EpsilonHotMaxwellian_swan(3,2) = sinPsi*K4_swan+cosPsi*K5_swan 
        !EpsilonHotMaxwellian_swan(3,3) = K3_swan 

        !IdentMat = 0
        !IdentMat(1,1) = 1
        !IdentMat(2,2) = 1
        !IdentMat(3,3) = 1

        !SigmaHotMaxwellian_swan = -(EpsilonHotMaxwellian_swan-IdentMat)*omgrfc*eps0*zi

        !scIn = SpatialSigmaInput_cold(omgc,omgp2,omgrf,nuOmg)
        !sc = SigmaCold_stix(scIn)

        !write(*,*) 'kVec: ', kVec
        !write(*,*) '1,1  ', sigmaHot_maxwellian(1,1), SigmaHotMaxwellian_swan(1,1), sc(1,1)
        !!write(*,*) '1,2  ', sigmaHot_maxwellian(1,2), SigmaHotMaxwellian_swan(1,2), sc(1,2)
        !!write(*,*) '1,3  ', sigmaHot_maxwellian(1,3), SigmaHotMaxwellian_swan(1,3), sc(1,3)

        !!write(*,*) '2,1  ', sigmaHot_maxwellian(2,1), SigmaHotMaxwellian_swan(2,1), sc(2,1)
        !write(*,*) '2,2  ', sigmaHot_maxwellian(2,2), SigmaHotMaxwellian_swan(2,2), sc(2,2)
        !!write(*,*) '2,3  ', sigmaHot_maxwellian(2,3), SigmaHotMaxwellian_swan(2,3), sc(2,3)

        !!write(*,*) '3,1  ', sigmaHot_maxwellian(3,1), SigmaHotMaxwellian_swan(3,1), sc(3,1)
        !!write(*,*) '3,2  ', sigmaHot_maxwellian(3,2), SigmaHotMaxwellian_swan(3,2), sc(3,2)
        !write(*,*) '3,3  ', sigmaHot_maxwellian(3,3), SigmaHotMaxwellian_swan(3,3), sc(3,3)

        !sigmaHot_maxwellian = sc

        return

    end function sigmaHot_maxwellian


!    function sigmaCold_stix &
!        ( omgC, omgP2, omgRF, nuOmg )
    function sigmaCold_stix ( a )

        use constants 

        ! This routine calculates sigma_cold in the Stix frame
        ! (alp,bet,prl)
        ! ----------------------------------------------------

        implicit none

        !real(kind=dbl), intent(in) :: omgc, omgp2, omgrf
        !real, intent(in) :: nuOmg
        type(spatialSigmaInput_cold) :: a 

        complex :: omgrfc
        complex :: sig1, sig2, sig3
        complex :: sigmaCold_stix(3,3)
        complex :: zieps0

        complex :: K1_swan,K2_swan,K3_swan
        real(kind=DBL) :: e_swan
        complex :: sigmaCold_swan(3,3),epsilonCold_swan(3,3)
        integer :: IdentMat(3,3)

        zieps0 = zi * eps0
        omgRFc = a%omgRF * (1.0 + zi * a%nuOmg)

        sig1 = zieps0 * omgRFc * a%omgP2 / (omgRFc**2 - a%omgC**2)
        sig2 = eps0 * a%omgC   * a%omgP2 / (a%omgC**2 - omgRFc**2)
        sig3 = zieps0 * a%omgp2 / omgRFc

        sigmaCold_stix(1,1) = sig1 
        sigmaCold_stix(1,2) = sig2 
        sigmaCold_stix(1,3) = 0 

        sigmaCold_stix(2,1) = -sig2 
        sigmaCold_stix(2,2) = sig1 
        sigmaCold_stix(2,3) = 0 

        sigmaCold_stix(3,1) = 0
        sigmaCold_stix(3,2) = 0 
        sigmaCold_stix(3,3) = sig3


        ! Swanson version

        e_swan = sign(1d0,a%omgc) ! this is q/|q|
        K1_swan = 1 - a%omgP2 / (omgRFc**2-a%omgC**2)
        K2_swan = e_swan*a%omgC*a%omgP2/(omgRFc*(omgRFc**2-a%omgC**2))/zi
        K3_swan = 1-a%omgP2/omgRFc**2

        EpsilonCold_swan = (0,0)
        EpsilonCold_swan(1,1) = +K1_swan
        EpsilonCold_swan(1,2) = +K2_swan
        EpsilonCold_swan(2,1) = -K2_swan
        EpsilonCold_swan(2,2) = +K1_swan
        EpsilonCold_swan(3,3) = +K3_swan

        IdentMat = 0
        IdentMat(1,1) = 1
        IdentMat(2,2) = 1
        IdentMat(3,3) = 1

        SigmaCold_swan = -(EpsilonCold_swan-IdentMat)*omgrfc*eps0*zi

        write(*,*) '1,1  ', SigmaCold_stix(1,1), SigmaCold_swan(1,1)
        write(*,*) '2,2  ', SigmaCold_stix(2,2), SigmaCold_swan(2,2)
        write(*,*) '3,3  ', SigmaCold_stix(3,3), SigmaCold_swan(3,3)

        return

    end function sigmaCold_stix

end module sigma
