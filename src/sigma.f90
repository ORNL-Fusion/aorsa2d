module sigma_mod

contains

    function sigmaHot_maxwellian(i, j, &
        xm, &
        kt, omgc, omgp2, &
        kxsav, kysav, capr, &
        delta0, &
        omgrf, k0, &
        k_cutoff, specNo )

        use constants
        use zfunction_mod
        use rotation
        use bField
        use aorsa2din_mod, &
        only: upshift, nzfun, damping, lmax, xnuomg, nPhi

        ! This routine uses the modified Z functions Z0, Z1, Z2
        ! with the appropriate sign changes for k_parallel < 0.0
        ! No rotation is made.  Result is in the Stix frame,
        ! (alp,bet,prl)
        ! ---------------------------------------------------------

        implicit none

        complex :: sigmaHot_maxwellian(3,3)
        integer, intent(in) :: specNo
        real, intent(in) :: k_cutOff
        integer :: l, labs, i, j
        real :: kperp, kprl, xm, kt, omgc, omgp2
        real :: kprl_eff, fgam, y0, y, sgn_kprl
        real :: descrim
        real :: nu_coll
        real :: akprl,  alpha, omgrf
        real :: gammab(-lmax:lmax), gamma_coll(-lmax:lmax)
        real :: delta0, rhol
        real :: kxsav, kysav, capr
        real :: kphi
        real :: kalp, kbet, k0, rgamma, kr, step
        complex :: omgrfc
        complex :: z0, z1, z2
        complex :: sig0, sig1, sig2, sig3, sig4, sig5
        complex :: sig0l, sig1l, sig2l, sig3l, sig4l, sig5l
        integer, parameter :: lmaxdim = 99
        complex, allocatable, dimension(:) :: &
          exil, exilp, exilovergam
        complex :: zetal(-lmax:lmax)
        complex :: zieps0, al, bl, cl, gamma_

        nu_coll =  .01 * omgrf
        zieps0 = zi * eps0
        alpha = sqrt(2. * kt / xm)
        rhol = alpha / omgc
        kphi = nphi / capr
        omgrfc = omgrf * (1. + zi * xnuomg)


        ! k in Stix frame
        ! ---------------

        kalp = uxx(i,j) * kxsav + uxy(i,j) * kysav + uxz(i,j) * kphi
        kbet = uyx(i,j) * kxsav + uyy(i,j) * kysav + uyz(i,j) * kphi
        kprl = uzx(i,j) * kxsav + uzy(i,j) * kysav + uzz(i,j) * kphi
        kperp = sqrt(kalp**2 + kbet**2)


        ! Optional: leave out upshift in kprl
        ! --------------------------------- --

        if (upshift .eq. 0)  kprl = uzz(i,j) * kphi
        if (kprl  .eq. 0.0) kprl  = 1.0e-08
        if (kperp .eq. 0.0) kperp = 1.0e-08

        sgn_kprl = sign(1.0, kprl)
        akprl = abs(kprl)


        ! Calculate zetal(l) and gammab(l)
        ! ---------------------------------

        do l = -lmax, lmax

            labs = abs(l)
            zetal(l) = (omgrfc - l * omgc) / (kprl * alpha)
            gammab(l) = abs(l * omgc / (2.0 * alpha * kprl**2) &
                                                    * gradprlb(i,j) / bmod(i,j))
            gamma_coll(l) = nu_coll / (akprl * alpha)

            !if(xm .eq. xme)gammab(l) = 0.0

            ! electrons
            ! ---------
            if(specNo==1) gammab(l) = 0.0

            if(abs(gammab(l)) .lt. .01)gammab(l) = .01

        enddo


        ! Calculate Brambillas krpl_eff using l = 1 only
        ! ------------------------------------------------

        y0 = 1.5
        y = y0

        if(sgn_kprl .ge. 0.0)then

            fgam = 1.0

            if(gammab(1) .gt. 1.0e-05)then
                y = y0
                fgam = (sqrt(1. +  4. * gammab(1) * y) - 1.) &
                  / (2. * gammab(1) * y)
            endif

            kprl_eff = kprl / fgam

        endif


        if(sgn_kprl .lt. 0.0)then

            fgam = 1.0

            if(gammab(1) .gt. 1.0e-05)then

                descrim = 1. - 4. * gammab(1) * y0
                if (descrim .ge. 0.0) y =   y0
                if (descrim .lt. 0.0) y = - y0
                fgam = (1. - sqrt(1. -  4. * gammab(1) * y) ) &
                   / (2. * gammab(1) * y)

            endif

            kprl_eff = kprl / fgam

        endif


        ! Maxwellian distribution
        ! -----------------------

        gamma_ = 0.5 * kperp**2 * rhol**2
        rgamma = real(gamma_)
        
        if (.not. allocated(exil) ) &
        allocate ( exil(0:lMax), exilp(0:lMax), exilOverGam(0:lMax) )

        if(rgamma .ge. 1.0e-08) &
        call besiexp(gamma_, lmax, exil, exilp, lmaxdim, exilovergam)

        if(rgamma .lt. 1.0e-08) &
        call bes_expand(gamma_, lmax, exil, exilp, lmaxdim, exilovergam)

        sig0 = 0.0
        sig1 = 0.0
        sig2 = 0.0
        sig3 = 0.0
        sig4 = 0.0
        sig5 = 0.0

        do l = -lmax, lmax

            labs = abs(l)

            !if(nzfun .eq. 0) call z_approx(sgn_kprl, zetal(l), 0.0, z0, z1, z2)
            if(nzfun .eq. 1) call z_approx(sgn_kprl,zetal(l),gammab(l), z0, z1, z2)
            !if(nzfun .eq. 2) call z_smithe(sgn_kprl,zetal(l),gammab(l), z0, z1, z2)
            !if(nzfun .eq. 3) call z_table(sgn_kprl,zetal(l),gammab(l), gamma_coll(l), z0, z1, z2)

            al = 1.0 / (kprl * alpha) * z0
            bl = 1.0 / (kprl * alpha) * z1
            cl = 1.0 / (kprl * alpha) * z2

            sig0l = - zieps0 * omgp2 * rhol**2 * (exil(labs) - exilp(labs)) * al
            sig1l = - zieps0 * omgp2 * l**2 * exilovergam(labs) * al
            sig2l = - eps0 * omgp2 * l * (exil(labs) - exilp(labs)) * al
            sig3l = - zieps0 * omgp2 * 2.0 * exil(labs) * cl
            sig4l = - zieps0 * omgp2 * rhol * l * exilovergam(labs) * bl
            sig5l = - eps0 * omgp2 * rhol * (exil(labs) - exilp(labs)) * bl

            sig0 = sig0 + sig0l
            sig1 = sig1 + sig1l
            sig2 = sig2 + sig2l
            sig3 = sig3 + sig3l
            sig4 = sig4 + sig4l
            sig5 = sig5 + sig5l

        enddo
            
        sig1 = sig1 + delta0 * eps0 * omgrfc * kperp**2 / k0**2
        sig3 = sig3 + delta0 * eps0 * omgrfc * kperp**2 / k0**2

        !if (xm .eq. xme) then

        ! electrons
        ! ---------
        if (specNo==1) then

            kr = kperp / k_cutoff
            step = damping * kr**16 / (1. + kr**16)
            sig3 = sig3 * (1.0 + step)

        end if

        
        ! Swanson's rotation (original), to
        ! (alp,bet,prl)
        ! -----------------------------

        sigmaHot_maxwellian(1,1) = sig1 + sig0 * kbet**2
        sigmaHot_maxwellian(1,2) = sig2 - sig0 * kbet * kalp
        sigmaHot_maxwellian(1,3) = sig4 * kalp + sig5 * kbet

        sigmaHot_maxwellian(2,1) = - sig2 - sig0 * kbet * kalp
        sigmaHot_maxwellian(2,2) =   sig1 + sig0 * kalp**2
        sigmaHot_maxwellian(2,3) =   sig4 * kbet - sig5 * kalp

        sigmaHot_maxwellian(3,1) = sig4 * kalp - sig5 * kbet
        sigmaHot_maxwellian(3,2) = sig4 * kbet + sig5 * kalp
        sigmaHot_maxwellian(3,3) = sig3

        return

    end function sigmaHot_maxwellian


    function sigmaCold_stix &
        ( i, j, omgC, omgP2, omgRF )

        use constants 
        use rotation
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

!
!***************************************************************************
!


      subroutine besiexp(gamma_, lmax, expbes, expbesp, lmaxdim, &
         expbesovergam)

        use bessel_mod
        use constants

!     ----------------------------------------------------------
!     Calculates exp(-gamma) times the modified bessel functions
!     (and derivatives) of order up to lmax
!     ----------------------------------------------------------


      implicit none

      integer, intent(in) :: lmax
      integer :: nmax, ier, l, lmaxdim
      complex, intent(in) :: gamma_
      complex, intent(out) :: expbes(:), expbesp(:), expbesovergam(:)
      complex(kind=dbl), allocatable, dimension(:) :: xil, xilp, b
      complex(kind=dbl) :: exgam
      real :: gammod
      integer :: status
      complex(kind=dbl) :: gamma_dbl

      gamma_dbl = gamma_

      allocate ( xil(0:lMax), xilp(0:lMax), b(lMax+1), stat = status )
      b = 0

      exgam = exp(-gamma_dbl)
      gammod = cabs(gamma_)

      if(gammod .le. 700.)then
         nmax = lmax
         call besic(gamma_dbl, nmax, b, ier)
         if(ier .ne. 0)write(6,100) ier

         do l = 0, lMax
            xil(l) = b(l+1)
         end do

         do l = 0, lmax
           if(l .eq. 0) xilp(0) = xil(1)
           if(l .ne. 0) xilp(l) = xil(l-1) - l / gamma_dbl * xil(l)
           expbes(l+1) = exgam * xil(l)
           expbesp(l+1) = exgam * xilp(l)
         end do
      end if

      if(gammod .gt. 700.)then
         do l = 0, lmax
            call bes_asym(gamma_, l, expbes(l+1), expbesp(l+1))
         end do
      end if

      do l = 0, lmax
         expbesovergam(l+1) = expbes(l+1) / gamma_
      end do

  100 format('ier = ', i5, 'besic failed')

      deallocate ( b )
      deallocate ( xil )
      deallocate ( xilp )
      end subroutine besiexp


!
!***************************************************************************
!


      subroutine bes_expand(gamma, lmax, expbes, expbesp, lmaxdim, &
         expbesovergam)

!-------------------------------------------------------------------
!     Calculates exp(-gamma) times the modified bessel functions
!     (and derivatives) of order up to lmax using second order
!     expansion for small argument
!-------------------------------------------------------------------

        use bessel_mod

      implicit none

      integer lmax, l, lmaxdim
      real :: factl
      complex gamma, expbes(0:lmaxdim), expbesp(0:lmaxdim), &
         expbesovergam(0:lmaxdim), xilovergam, &
         xil(0:lmaxdim), xilp(0:lmaxdim), exgam

      exgam = 1.0 - gamma + gamma**2 / 2.0

      do l = 0, lmax
         factl = factrl(l)
         xil(l) = gamma**l / (2**l * factl) * &
                                       ( 1. + gamma**2 / (4. * (l+1)))
         xilp(l) = gamma**(l-1) / (2**l * factl) * &
                           (l + (l+2) * gamma**2 / (4. * (l+1)))

         xilovergam = gamma**(l-1) / (2**l * factl) * &
                                       ( 1. + gamma**2 / (4. * (l+1)))

         expbes(l) = exgam * xil(l)
         expbesp(l) = exgam * xilp(l)
         expbesovergam(l) = exgam * xilovergam

!         write(6, 100)l, expbes(l), expbesp(l)
      end do


      return
      end subroutine bes_expand


!
!***************************************************************************
!

      subroutine bes_asym(z, n, exil, exilp)

      implicit none

      integer mu, n
      real pi
      complex z, exil, exilp
      data pi/3.141592654/

      mu = 4 * n**2
      exil =  1.0 / csqrt(2.0 * pi * z) &
         * (1.0 &
         - (mu - 1)/(8.0 * z) &
         + (mu - 1) * (mu - 9) / (2.0 * (8.0 * z)**2) &
         - (mu - 1) * (mu - 9) * (mu - 25) / (6.0 * (8.0 * z)**3)  )
      exilp = 1.0 / csqrt(2.0 * pi * z) &
         * (1.0 &
         - (mu + 3)/(8.0 * z) &
         + (mu - 1) * (mu + 15) / (2.0 * (8.0 * z)**2) &
         - (mu - 1) * (mu - 9) * (mu + 35) / (6.0 * (8.0 * z)**3)  )

      return
      end subroutine bes_asym

!
!***************************************************************************
!


      function fzeta (arg)
      complex a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3, arg, aux0 &
         , aux1, term, z, zz, fzeta
      data d1r/0.0/
      data d1i/1.77245385090551/
      data d2r/0.0/
      data d2i/3.54490770181103/
      data d3r/0.0/
      data d3i/7.08981540362206/
      data d4/0.33333333333333/
      data eps/1.0E-07/

!     data d4/0.33333333333333/

!ray  code analysis
!ray  optimize
!
      i = 0
      z = arg
      zz = z*z
      x = real(z)
      y = aimag(z)
      d1 = cmplx(d1r,d1i)
      d2 = cmplx(d2r,d2i)
      d3 = cmplx(d3r,d3i)
      ymag = abs(y)
      if (ymag - 1.0 .ge. 0.) then
!
!     continued fraction method: abs(y).ge.1.0
!
         y0 = y
         y = ymag
         aux1 = 1.5 - z*z
         aux2 = 0.0
         del = 1.5
         a1 = 0.0
         a2 = -1.0
         b1 = 1.0
         b2 = aux1
         c1 = a2/b2
!
  100    continue
         aux1 = aux1 + 2.0
         aux2 = aux2 - del
         del = del + 2.0
         a3 = aux1*a2 + aux2*a1
         b3 = aux1*b2 + aux2*b1
         c2 = a3/b3
         c3 = c2 - c1
         c3r = real(c3)
         c3i = aimag(c3)
         if (abs(c3r) + abs(c3i) .lt. eps) go to 110
         a1 = a2
         a2 = a3
         b1 = b2
         b2 = b3
         c1 = c2
         go to 100
  110    continue
         if (y0 .lt. 0.) then
            y = y0
            c2 = conjg(c2) - d3*z*exp(-zz)
         endif
         aux0 = -(0.5*c2 + 1.0)/z
      else
!
!     asymptotic series method: abs(x).ge.4.0 and abs(y).lt.1.0
!
         xmag = abs(x)
         if (xmag - 4.0 .lt. 0.) go to 130
         term = 1.0/z
         aux0 = -term
         aux1 = 0.5*term**2
         p = 1.0
         if (y .le. 0.) then
            if (y .ne. 0.) then
               aux0 = aux0 + d2*exp(-zz)
            else
               aux0 = aux0 + d1*exp(-zz)
            endif
         endif
  120    continue
         term = aux1*term*p
         aux0 = aux0 - term
         p = p + 2.0
         termr = real(term)
         termi = aimag(term)
!     if(abs(termr)+abs(termi).lt.eps)30,18
         if (abs(termr) + abs(termi) .lt. eps) go to 160
         go to 120
!
!     power series method: abs(x).lt.4.0 and abs(y).lt.1.0
!
  130    continue
         aux0 = 1.0
         aux1 = -(zz + zz)
         aux2 = eps/(eps + xmag + ymag)
         term = d4*aux1
         p = 3.0
  140    continue
         aux0 = aux0 + term
         termr = real(term)
         termi = aimag(term)
!     if(abs(termr)+abs(termi).lt.aux2)26,24
         if (abs(termr) + abs(termi) .lt. aux2) go to 150
         p = p + 2.0
         term = aux1*term/p
         go to 140
  150    continue
         aux0 = d1*exp(-zz) - 2.0*z*aux0
      endif
  160 continue
      fzeta = aux0
      if (i .le. 0) return
      fzeta = -2.0*(1.0 + arg*aux0)
      return
      end function fzeta
!
!***********************************************************************
!


end module sigma_mod
