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

        implicit none

        real, intent(in) :: sinTh, bPol, bMod, gradPrlB
        complex(kind=dbl) :: sigmaHot_maxwellian(3,3)
        integer, intent(in) :: specNo
        real, intent(in) :: k_cutOff, nuOmg
        real, intent(in) :: kVec(:), capr

        real(kind=dbl), intent(in) :: omgc, omgp2
        real(kind=dbl), intent(in) :: omgrf
        real(kind=dbl), intent(in) :: kt
        complex(kind=dbl), intent(in) :: k0

        integer :: labs, i, j
        integer :: l, lMaxTmp
        real(kind=dbl) :: kPer, kPrl, xm, kPrlTmp
        real(kind=dbl) :: kPrlEff, fgam, y0, y
        real(kind=dbl) :: sgn_kPrl
        real :: descrim
        real(kind=dbl) :: vTh
        real(kind=dbl) :: alpha_smithe(-lmax:lmax), gamma_smithe(-lmax:lmax)
        real(kind=dbl) :: gamma_brambilla(-lmax:lmax)
        real(kind=dbl) :: rhol
        real(kind=dbl) :: kAlp, kBet, rgamma
        real(kind=dbl) :: kr, step
        complex(kind=dbl) :: omgrfc
        complex(kind=dbl), allocatable :: Z_bram(:), Zp_bram(:), zeta_Zp_bram(:)
        complex(kind=dbl), allocatable :: Z_swan(:), Zp_swan(:), zeta_Zp_swan(:)
        complex(kind=dbl) :: sig0, sig1, sig2, sig3, sig4, sig5
        complex(kind=dbl), dimension(-lmax:lmax) :: &
            sig0l, sig1l, sig2l, sig3l, sig4l, sig5l
        integer, parameter :: lmaxdim = 99
        complex(kind=dbl), allocatable, dimension(:) :: &
          expBesselI, expBesselIPrime, expBesselIOverGam
        complex(kind=dbl) :: zeta_l_bram(-lmax:lmax),zeta_l_swan(-lmax:lmax)
        complex :: zieps0 
        complex(kind=dbl) :: P, PPrime 
        complex(kind=dbl) ::gamma_
        complex(kind=dbl) ::lambda
        real :: dR, Lpar
        integer :: stat
        real :: sinPh
        real :: cosPsi, sinPsi

        real :: e_swan
        real(kind=dbl) :: eps_swan
        complex(kind=DBL) :: factor1,factor2
        complex(kind=DBL) :: K0_HarmSum,K1_HarmSum,K2_HarmSum,K3_HarmSum,K4_HarmSum,K5_HarmSum
        complex(kind=DBL) :: K0_swan,K1_swan,K2_swan,K3_swan,K4_swan,K5_swan
        complex(kind=dbl) :: SigmaHotMaxwellian_swan(3,3),EpsilonHotMaxwellian_swan(3,3)
        integer :: IdentMat(3,3)
        type(SpatialSigmaInput_cold) :: scIn
        complex :: sc(3,3)

        complex(kind=dbl), allocatable :: In_(:), Inp(:)
        integer :: ier
        real :: omgC_swan
        complex(kind=dbl) :: g1_swan, g2_swan

        zeta_l_bram = 0
        Zeta_l_swan = 0

        zieps0 = zi * eps0
        vTh = sqrt(2. * kt / xm)
        rhol = vTh / omgc
        omgrfc = omgrf * (1d0 + zi * nuOmg)
        omgC_swan = abs(omgC)

        ! k in Stix frame
        ! ---------------

        kAlp = kVec(1)
        kBet = kVec(2)
        kPrl = kVec(3)

        kPer = sqrt(kAlp**2 + kBet**2)

        if (kPrl  == 0.0) kPrl  = 1.0e-05
        if (kPer == 0.0) kPer = 1.0e-05

        sgn_kPrl    = sign ( 1d0, kPrl )
        

        ! Calculate Zeta (Z function argument + broadening terms)
        ! -------------------------------------------------------

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
                    stop
                endif

                if(kPrl==0.0)then
                    kPrlEff = sqrt ( omgrf / ( 2*capR*vTh) * sinPh )
                else
                    kPrlEff = kPrl*( sqrt(1 + 4*gamma_brambilla(l)) - 1 ) / (2*gamma_brambilla(l))
                endif

            else
                kPrlEff = kPrl
            endif

            ! I'm seeing problems (i.e., non-convergence with cold plasma
            ! sigma), for abs(kPrl)<~1e-3, so setting kPrlEffLimit~1e-1.
            ! The problems originate in there being a 1/kPrl in the eta
            ! argument to the Z function. And in the 3,3 sigma element 
            ! there is a multiplication by this factor. I think there must
            ! be lack of precision that occurs since there are many numbers
            ! requiring large precision and cancellation resulting in an
            ! error when kPrl gets too small. I would not think these terems
            ! are all that important anyway, since kPrl=1e-1 is a parallel
            ! wavelength of ~60 m, which is more than a trip around the 
            ! machine.
            ! -----------------------------------------------------------

            if(abs(kPrlEff)<kPrlEffLimit)then
                    kPrlEff = sign(kPrlEffLimit,kPrlEff)
            endif

            zeta_l_bram(l)  = (omgrfc-l*omgC)/(kPrlEff*vTh) 
            zeta_l_swan(l)  = (omgrfc+l*omgC_swan)/(kPrlEff*vTh)

        enddo


        ! Bessel functions ... exp(-Gamma) * I_l etc 
        ! ------------------------------------------

        gamma_ = 0.5 * kPer**2 * rhol**2
        lambda = gamma_

        allocate ( &
            expBesselI(0:lMax), &
            expBesselIPrime(0:lMax), &
            expBesselIOverGam(0:lMax), &
            In_(0:lMax), Inp(0:lMax) )

        ! Somewhere in the bessel function routines the lMax argument
        ! is being modified, hence the lMaxTmp temp variable.

        if ( abs(gamma_) >= 1.0e-08 ) then
            lMaxTmp = lMax
            call besIExp    ( gamma_, lMaxTmp, expBesselI, expBesselIPrime, expBesselIOverGam )
            lMaxTmp = lMax
            call besic ( gamma_, lMaxTmp, In_, ier )
        else
            lMaxTmp = lMax
            call bes_expand ( gamma_, lMaxTmp, expBesselI, expBesselIPrime, expBesselIOverGam )
            lMaxTmp = lMax
            call besic ( gamma_, lMaxTmp, In_, ier )
        endif

        do l=0,lMax

            if (l==0) then 
                Inp(l) = In_(1) ! recall the I(x)_n = I(x)_-n symmetry for integer n.
            else
                Inp(l) = In_(l-1) - l / gamma_ * In_(l)
            endif 

        enddo

        !write(*,*)
        !write(*,*) expBesselI
        !write(*,*) exp(-lambda)*In_
        !write(*,*) expBesselIPrime
        !write(*,*) exp(-lambda)*Inp


        ! Z function
        ! ----------

        allocate ( Z_bram(-lMax:lMax), Zp_bram(-lMax:lMax), zeta_Zp_bram(-lMax:lMax) )
        allocate ( Z_swan(-lMax:lMax), Zp_swan(-lMax:lMax), zeta_Zp_swan(-lMax:lMax) )

        do l=-lMax,lMax
            call z_from_table ( real(zeta_l_bram(l)), z_bram(l), Zp_bram(l), zeta_Zp_bram(l) )
            call z_from_table ( real(zeta_l_swan(l)), z_swan(l), Zp_swan(l), zeta_Zp_swan(l) )
        enddo

        write(*,*)
        write(*,*) zeta_l_bram
        write(*,*) zeta_l_swan

        write(*,*)
        write(*,*) z_bram
        write(*,*) z_swan

        write(*,*)
        write(*,*) zp_bram
        write(*,*) zp_swan

 
        ! Brambilla
        ! ---------

        do l = -lmax, lmax

            labs = abs(l) ! For integer order, I(x)_l = I(x)_-l 

            P = omgp2 / (kPrlEff * vTh) * Z_bram(l)
            PPrime = omgp2 / (kPrlEff * vTh) * Zp_bram(l)

            sig0l(l) = - zieps0 * rhol**2 * (expBesselI(labs) - expBesselIPrime(labs)) * P 
            sig1l(l) = - zieps0 * l**2 * expBesselIOverGam(labs) * P
            sig2l(l) = - eps0 * l * (expBesselI(labs) - expBesselIPrime(labs)) * P
            sig3l(l) =   zieps0 * expBesselI(labs) * zeta_l_bram(l) * PPrime
            sig4l(l) =   0.5d0 * zieps0 * rhol * l * expBesselIOverGam(labs) * PPrime
            sig5l(l) =   0.5d0 * eps0 * rhol * (expBesselI(labs) - expBesselIPrime(labs)) * PPrime

        enddo

        sig0    = sum ( sig0l )
        sig1    = sum ( sig1l )
        sig2    = sum ( sig2l )
        sig3    = sum ( sig3l )
        sig4    = sum ( sig4l )
        sig5    = sum ( sig5l )

        ! Add numerical damping for Bernstein wave (delta0 ~ 1e-4)
        ! --------------------------------------------------------

        sig1 = sig1 + delta0 * eps0 * omgrfc * kPer**2 / k0**2
        sig3 = sig3 + delta0 * eps0 * omgrfc * kPer**2 / k0**2


        ! electron damping
        ! ----------------

        if (specNo==1) then

            kr = kPer / k_cutOff
            step = damping * kr**16 / (1. + kr**16)
            sig3 = sig3 * (1.0 + step)

        endif

        cosPsi = kAlp / kPer
        sinPsi = kBet / kPer

        sigmaHot_maxwellian(1,1) = sig1 + sig0 * sinPsi**2
        sigmaHot_maxwellian(1,2) = sig2 - sig0 * cosPsi * sinPsi
        sigmaHot_maxwellian(1,3) = sig4 * cosPsi + sig5 * sinPsi

        sigmaHot_maxwellian(2,1) = - sig2 - sig0 * cosPsi * sinPsi
        sigmaHot_maxwellian(2,2) =   sig1 + sig0 * cosPsi**2
        sigmaHot_maxwellian(2,3) =   sig4 * sinPsi - sig5 * cosPsi

        sigmaHot_maxwellian(3,1) = sig4 * cosPsi - sig5 * sinPsi
        sigmaHot_maxwellian(3,2) = sig4 * sinPsi + sig5 * cosPsi
        sigmaHot_maxwellian(3,3) = sig3


        ! Swanson Dielectric 
        ! ------------------

        K0_HarmSum=(0,0) 
        K1_HarmSum=(0,0) 
        K2_HarmSum=(0,0) 
        K3_HarmSum=(0,0) 
        K4_HarmSum=(0,0) 
        K5_HarmSum=(0,0)

        do l=-lmax,lmax
        
            labs = abs(l) ! For integer order, I(x)_l = I(x)_-l

            K0_HarmSum = K0_HarmSum + lambda * ( In_(labs) - Inp(labs)  ) * Z_swan(l)
            K1_HarmSum = K1_HarmSum + l**2 * In_(labs)  / lambda * Z_swan(l)
            K2_HarmSum = K2_HarmSum + l * ( In_(labs)  - Inp(labs)  ) * Z_swan(l)
            K3_HarmSum = K3_HarmSum + In_(labs)  * zeta_l_swan(l) * Zp_swan(l) 
            K4_HarmSum = K4_HarmSum + l * In_(labs)  / lambda * Zp_swan(l)  
            K5_HarmSum = K5_HarmSum + ( In_(labs)  - Inp(labs)  ) * Zp_swan(l)

        enddo

        eps_swan = omgC / abs(omgC) 
        g1_swan = omgP2 * exp(-lambda) / ( omgRFc * kPrlEff * vTh )
        g2_swan = kPer * omgP2 * exp(-lambda) / ( 2 * kPrlEff * omgRFc * omgC_swan )

        K0_swan = 2 * g1_swan * K0_HarmSum 
        K1_swan = 1 + g1_swan * K1_HarmSum
        K2_swan = zi * eps_swan * g1_swan * K2_HarmSum
        K3_swan = 1 - g1_swan * K3_HarmSum
        K4_swan = g2_swan * K4_HarmSum
        K5_swan = zi * eps_swan * g2_swan * K5_HarmSum

        deallocate ( expBesselI, expBesselIPrime, expBesselIOverGam )

        EpsilonHotMaxwellian_swan(1,1) = +K1_swan
        EpsilonHotMaxwellian_swan(1,2) = +K2_swan
        EpsilonHotMaxwellian_swan(1,3) = +K4_swan

        EpsilonHotMaxwellian_swan(2,1) = -K2_swan
        EpsilonHotMaxwellian_swan(2,2) = +K1_swan+K0_swan 
        EpsilonHotMaxwellian_swan(2,3) = -K5_swan

        EpsilonHotMaxwellian_swan(3,1) = +K4_Swan
        EpsilonHotMaxwellian_swan(3,2) = +K5_Swan
        EpsilonHotMaxwellian_swan(3,3) = +K3_Swan 


        !EpsilonHotMaxwellian_swan = transpose(EpsilonHotMaxwellian_swan)

        ! Convert dielectrics to conductivities
        ! -------------------------------------

        IdentMat = 0
        IdentMat(1,1) = 1
        IdentMat(2,2) = 1
        IdentMat(3,3) = 1

        SigmaHotMaxwellian_swan = -(EpsilonHotMaxwellian_swan-IdentMat)*omgrfc*eps0*zi

        !scIn = SpatialSigmaInput_cold(omgc,omgp2,omgrf,nuOmg)
        !sc = SigmaCold_stix(scIn)

        !sigmahot_maxwellian = sc
        
        write(*,*)
        write(*,*) 'bram 1,1', sigmaHot_Maxwellian(1,1)
        write(*,*) 'swan 1,1', sigmaHotMaxwellian_swan(1,1)

        write(*,*)
        write(*,*) 'bram 2,2', sigmaHot_Maxwellian(2,2)
        write(*,*) 'swan 2,2', sigmaHotMaxwellian_swan(2,2)

        write(*,*)
        write(*,*) 'bram 3,1', sigmaHot_Maxwellian(3,1)
        write(*,*) 'swan 3,1', sigmaHotMaxwellian_swan(3,1)

        write(*,*)
        write(*,*) 'bram 3,2', sigmaHot_Maxwellian(3,2)
        write(*,*) 'swan 3,2', sigmaHotMaxwellian_swan(3,2)

        write(*,*)
        write(*,*) 'bram 3,3', sigmaHot_Maxwellian(3,3)
        write(*,*) 'swan 3,3', sigmaHotMaxwellian_swan(3,3)

        SigmaHot_Maxwellian = SigmaHotMaxwellian_swan 

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

        type(spatialSigmaInput_cold) :: a 

        complex(kind=DBL) :: omgrfc
        complex :: sig1, sig2, sig3
        complex :: sigmaCold_stix(3,3)
        complex :: zieps0

        complex :: K1_HarmSum,K2_HarmSum,K3_HarmSum
        real(kind=DBL) :: e_swan
        complex :: sigmaCold_swan(3,3),epsilonCold_swan(3,3)
        integer :: IdentMat(3,3)
        real(kind=DBL) :: nu

        zieps0 = zi * eps0
        nu = a%nuOmg * a%omgRF

        sig1 = zieps0 * omgRFc * a%omgP2 / (omgRFc**2 - a%omgC**2) ! Stix_S
        sig2 = eps0 * a%omgC   * a%omgP2 / (a%omgC**2 - omgRFc**2) ! Stix_D
        sig3 = zieps0 * a%omgp2 / omgRFc ! Stix_P

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

        e_swan = sign(1d0,a%omgC) ! this is q/|q|
        K1_HarmSum = 1d0 - a%omgP2*(a%omgRF+zi*nu) / ( a%omgRF * ( a%omgRF * (a%omgRF+zi*nu)**2-a%omgC**2) )
        K2_HarmSum = e_swan*a%omgC*a%omgP2 / ( a%omgRF * ( (a%omgRF+zi*nu)**2 - a%omgC**2 ) ) / zi
        K3_HarmSum = 1d0 - a%omgP2/(a%omgRF*(a%omgRF+zi*nu))

        EpsilonCold_swan = (0,0)
        EpsilonCold_swan(1,1) = +K1_HarmSum
        EpsilonCold_swan(1,2) = +K2_HarmSum
        EpsilonCold_swan(2,1) = -K2_HarmSum
        EpsilonCold_swan(2,2) = +K1_HarmSum
        EpsilonCold_swan(3,3) = +K3_HarmSum

        IdentMat = 0
        IdentMat(1,1) = 1
        IdentMat(2,2) = 1
        IdentMat(3,3) = 1

        SigmaCold_swan = -(EpsilonCold_swan-IdentMat)*a%omgRF*eps0*zi
        SigmaCold_stix = SigmaCold_swan

!        write(*,*) '1,1  ', SigmaCold_stix(1,1), SigmaCold_swan(1,1)
!        write(*,*) '2,2  ', SigmaCold_stix(2,2), SigmaCold_swan(2,2)
!        write(*,*) '3,3  ', SigmaCold_stix(3,3), SigmaCold_swan(3,3)
!
        return

    end function sigmaCold_stix

end module sigma
