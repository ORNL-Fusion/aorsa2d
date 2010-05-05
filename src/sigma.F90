module sigma_mod

contains

    function sigmaHot_maxwellian( &
        xm, kt, omgc, omgp2, &
        kxsav, kysav, capr, &
        omgrf, k0, &
        k_cutoff, specNo, &
        bMod, gradPrlB, U_xyz, U_cyl )

        use constants
        use zfunction_mod
        use aorsa2din_mod, &
        only: upshift, nzfun, damping, lmax, xnuomg, &
            nPhi, delta0
        use bessel_mod

        ! Calculate the hot Maxwellian sigma at a single
        ! point in space
        !-----------------------------------------------

        ! This routine uses the modified Z functions Z0, Z1, Z2
        ! with the appropriate sign changes for k_parallel < 0.0
        ! No rotation is made.  Result is in the Stix frame,
        ! (alp,bet,prl)
        ! ---------------------------------------------------------

        implicit none

        real, intent(in) :: bMod, gradPrlB, U_xyz(:,:), U_cyl(:,:)
        complex :: sigmaHot_maxwellian(3,3)
        integer, intent(in) :: specNo
        real, intent(in) :: k_cutOff

        integer :: l, labs, i, j
        real :: kPerp, kPrl, xm, kt, omgc, omgp2, kPrlTmp
        real :: kPrl_eff, fgam, y0, y
        real(kind=dbl) :: sgn_kPrl
        real :: descrim
        real :: nu_coll
        real :: akPrl,  alpha, omgrf
        real(kind=dbl) :: gammaBroaden(-lmax:lmax), gamma_coll(-lmax:lmax)
        real :: rhol
        real :: kxsav, kysav, capr
        real :: kPhi
        real :: kAlp, kBet, k0, rgamma, kr, step
        complex :: omgrfc
        complex(kind=dbl) :: z0, z1, z2
        complex :: sig0, sig1, sig2, sig3, sig4, sig5
        complex :: sig0l, sig1l, sig2l, sig3l, sig4l, sig5l
        integer, parameter :: lmaxdim = 99
        complex(kind=dbl), allocatable, dimension(:) :: &
          expBesselI, expBesselIPrime, expBesselIOverGam
        complex(kind=dbl) :: zetal(-lmax:lmax)
        complex :: zieps0, al, bl, cl 
        complex(kind=dbl) ::gamma_
        real :: dR

        !nu_coll =  .01 * omgrf
        zieps0 = zi * eps0
        alpha = sqrt(2. * kt / xm)
        rhol = alpha / omgc
        kPhi = nphi / capr
        omgrfc = omgrf * (1. + zi * xnuomg)


        ! k in Stix frame
        ! ---------------

#ifdef cylProper 

        !kAlp = Urr_(i,j) * kxsav + Urth_(i,j) * kPhi + Urz_(i,j) * kysav
        !kBet = Uthr_(i,j) * kxsav + Uthth_(i,j) * kPhi + Uthz_(i,j) * kysav
        !kPrl = Uzr_(i,j) * kxsav + Uzth_(i,j) * kPhi + Uzz_(i,j) * kysav

        kAlp = U_cyl(1,1) * kxsav + U_cyl(1,2) * kPhi + U_cyl(1,3) * kysav
        kBet = U_cyl(2,1) * kxsav + U_cyl(2,2) * kPhi + U_cyl(2,3) * kysav
        kPrl = U_cyl(3,1) * kxsav + U_cyl(3,2) * kPhi + U_cyl(3,3) * kysav

        !if ( upShift == 0 ) kPrl = Uzth_(i,j) * kPhi
        if ( upShift == 0 ) kPrl = U_cyl(3,2) * kPhi
#else 
        !kAlp = uxx(i,j) * kxsav + uxy(i,j) * kysav + uxz(i,j) * kPhi
        !kBet = uyx(i,j) * kxsav + uyy(i,j) * kysav + uyz(i,j) * kPhi
        !kPrl = uzx(i,j) * kxsav + uzy(i,j) * kysav + uzz(i,j) * kPhi

        kAlp = U_xyz(1,1) * kxsav + U_xyz(1,2) * kPhi + U_xyz(1,3) * kysav
        kBet = U_xyz(2,1) * kxsav + U_xyz(2,2) * kPhi + U_xyz(2,3) * kysav
        kPrl = U_xyz(3,1) * kxsav + U_xyz(3,2) * kPhi + U_xyz(3,3) * kysav

        !if (upshift == 0)  kPrl = uzz(i,j) * kPhi
        if (upshift == 0)  kPrl = U_xyz(3,3) * kPhi

#endif

        kPerp = sqrt(kAlp**2 + kBet**2)


        if (kPrl  == 0.0) kPrl  = 1.0e-08
        if (kPerp == 0.0) kPerp = 1.0e-08

        sgn_kPrl    = sign ( 1.0, kPrl )
        akPrl       = abs ( kPrl )


        ! Calculate zetal(l) and gammaBroaden(l)
        ! ---------------------------------

        do l = -lmax, lmax

            kPrlTmp = kPrl
            dR      = real(omgc) / ( real(kPrl) * alpha )
            !if((abs(dR)/10.0)<abs(dx)) 
            !write(*,*) dx, dR, kPrl, alpha
            !if(abs(kPrlTmp)>25) kPrlTmp=sgn_kPrl * 25 
            zetal(l) = (omgrfc - l * omgc) / (kPrlTmp * alpha) 

            gammaBroaden(l) = abs(l * omgc / (2.0 * alpha * kPrlTmp**2) &
                                                    * gradPrlB / bMod)
            !gamma_coll(l) = nu_coll / (akPrl * alpha)

            ! electrons
            if(specNo==1) &
            gammaBroaden(l) = 0.0
           
            ! not sure why this is nescessary yet ... ask EFJ 
            if ( abs(gammaBroaden(l) ) < .01 ) &
            gammaBroaden(l) = .01

        enddo


        ! Maxwellian distribution
        ! -----------------------

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

        sig0 = 0.0
        sig1 = 0.0
        sig2 = 0.0
        sig3 = 0.0
        sig4 = 0.0
        sig5 = 0.0

        do l = -lmax, lmax

            labs = abs(l)

            !if(nzfun .eq. 0) call z_approx(sgn_kPrl, zetal(l), 0.0, z0, z1, z2)
            if(nzfun .eq. 1) call z_approx ( sgn_kPrl, zetal(l), gammaBroaden(l), z0, z1, z2)
            !if(nzfun .eq. 2) call z_smithe(sgn_kPrl,zetal(l),gammaBroaden(l), z0, z1, z2)
            !if(nzfun .eq. 3) call z_table(sgn_kPrl,zetal(l),gammaBroaden(l), gamma_coll(l), z0, z1, z2)

            al = 1.0 / (kPrl * alpha) * z0
            bl = 1.0 / (kPrl * alpha) * z1
            cl = 1.0 / (kPrl * alpha) * z2

            sig0l = - zieps0 * omgp2 * rhol**2 * (expBesselI(labs) - expBesselIPrime(labs)) * al
            sig1l = - zieps0 * omgp2 * l**2 * expBesselIOverGam(labs) * al
            sig2l = - eps0 * omgp2 * l * (expBesselI(labs) - expBesselIPrime(labs)) * al
            sig3l = - zieps0 * omgp2 * 2.0 * expBesselI(labs) * cl
            sig4l = - zieps0 * omgp2 * rhol * l * expBesselIOverGam(labs) * bl
            sig5l = - eps0 * omgp2 * rhol * (expBesselI(labs) - expBesselIPrime(labs)) * bl

            sig0 = sig0 + sig0l
            sig1 = sig1 + sig1l
            sig2 = sig2 + sig2l
            sig3 = sig3 + sig3l
            sig4 = sig4 + sig4l
            sig5 = sig5 + sig5l

        enddo
           
        deallocate ( expBesselI, expBesselIPrime, expBesselIOverGam )


        ! Add numerical damping for Bernstein wave (delta0 ~ 1e-4)
        ! --------------------------------------------------------

        sig1 = sig1 + delta0 * eps0 * omgrfc * kPerp**2 / k0**2
        sig3 = sig3 + delta0 * eps0 * omgrfc * kPerp**2 / k0**2


        ! electron damping
        ! ----------------

        if (specNo==1) then

            kr = kPerp / k_cutoff
            step = damping * kr**16 / (1. + kr**16)
            sig3 = sig3 * (1.0 + step)

        end if

        
        ! Swanson's rotation (original), to
        ! (alp,bet,prl)
        ! -----------------------------

        sigmaHot_maxwellian(1,1) = sig1 + sig0 * kBet**2
        sigmaHot_maxwellian(1,2) = sig2 - sig0 * kBet * kAlp
        sigmaHot_maxwellian(1,3) = sig4 * kAlp + sig5 * kBet

        sigmaHot_maxwellian(2,1) = - sig2 - sig0 * kBet * kAlp
        sigmaHot_maxwellian(2,2) =   sig1 + sig0 * kAlp**2
        sigmaHot_maxwellian(2,3) =   sig4 * kBet - sig5 * kAlp

        sigmaHot_maxwellian(3,1) = sig4 * kAlp - sig5 * kBet
        sigmaHot_maxwellian(3,2) = sig4 * kBet + sig5 * kAlp
        sigmaHot_maxwellian(3,3) = sig3

        return

    end function sigmaHot_maxwellian


    function sigmaCold_stix &
        ( i, j, omgC, omgP2, omgRF )

        use constants 
        use aorsa2din_mod, &
        only: xnuomg

        ! This routine calculates sigma_cold in the Stix frame
        ! (alp,bet,prl)
        ! ----------------------------------------------------

        implicit none

        integer, intent(in) :: i, j
        real, intent(in) :: omgc, omgp2
        real, intent(in) :: omgrf
        complex :: omgrfc
        complex :: sig1, sig2, sig3
        complex :: sigmaCold_stix(3,3)
        complex :: zieps0

        zieps0 = zi * eps0
        omgRFc = omgRF * (1.0 + zi * xNuomg)

        sig1 = zieps0 * omgRFc * omgP2 / (omgRFc**2 - omgC**2)
        sig2 = - eps0 * omgC   * omgP2 / (omgRFc**2 - omgC**2)
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
