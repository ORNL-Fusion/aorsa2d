module sigma_mod

contains

      subroutine sigma_maxwellian(i, j, n, m, &
         gradprlb, bmod, bmod0, &
         xm, xn, dlg_xnuomg, &
         xkt, omgc, omgp2, &
         lmin, lmax, nzfun, &
         xkxsav, xkysav, nphi, capr, &
         bx, by, bz, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         sigxx, sigxy, sigxz, &
         sigyx, sigyy, sigyz, &
         sigzx, sigzy, sigzz, &
         delta0, ndist, &
         omgrf, xk0, &
         upshift, damping, xk_cutoff )

      use constants
      use zfunction_mod

!     ---------------------------------------------------------
!     This routine uses the modified Z functions Z0, Z1, Z2
!     with the appropriate sign changes for k_parallel < 0.0
!     No rotation is made.  Result is in the Stix frame.
!     ---------------------------------------------------------

      implicit none

      real, dimension(:,:), allocatable :: DFDUPER0, DFDUPAR0

      integer lmin, lmax, nzfun, lmaxdim, l, labs, &
          i, j, n, m, nphi, ndist, myid, iflag, nproc, ni, mi, &
          ni0, mi0, upshift

      real u0, u2, fnorm, f_cql, eNormIN

      real xkperp, xkprl, xm, xn, xkt, omgc, omgp2
      real xkprl_eff, fgam, y0, y, sgn_kprl, reson, duperp, dupara
      real xkprl_eff0
      real dzetal(lmin:lmax), descrim
      real dakbdkb, gradprlb, bmod, bmod0, nu_coll
      real, intent(in) :: dlg_xnuomg
      real akprl,  alpha, omgrf, emax
      real gammab(lmin:lmax), gamma_coll(lmin:lmax)
      real a, b, xnurf, delta0, rhol
      real bx, by, bz, bratio, denom
      real dfdth, dfdupar_check, dfduper_check, dfdth_check
      real uperp0_grid, upara0_grid, zeta, eta, ai, bi, ci, di
      real dfduper0_intplt, dfdupar0_intplt

      real xkxsav, xkysav, capr, argd
      real xkphi
      real xkalp, xkbet, xk0, rgamma, xk_cutoff, damping, kr, step

      complex zfunct, fzeta, omgrfc

      complex zfunct0, zeta0, sig3cold, z0, z1, z2, dz0, dz1, dz2

      complex sig0, sig1, sig2, sig3, sig4, sig5
      complex sig0_a, sig1_a, sig2_a, sig3_a, sig4_a, sig5_a
      complex sig0_h, sig1_h, sig2_h, sig3_h, sig4_h, sig5_h

      complex sig0l, sig1l, sig2l, sig3l, sig4l, sig5l
      logical :: l_interp   !new
      logical :: l_first  !new


      complex sigxx, sigxy, sigxz, &
              sigyx, sigyy, sigyz, &
              sigzx, sigzy, sigzz

      real uxx, uxy, uxz, &
           uyx, uyy, uyz, &
           uzx, uzy, uzz

      parameter (lmaxdim = 99)

      complex xil(0: lmaxdim), xilp(0: lmaxdim)
      complex exil(0: lmaxdim), exilp(0: lmaxdim), &
                          exilovergam(0: lmaxdim)

      complex zetal(lmin:lmax)
      complex  zieps0, arg, &
         al, bl, cl, &
         gamma, zeta_eff

      integer NBESSJ

      real :: UPERP0, UPARA0

      real :: W,K1(3),K2(3),KPER1,KPER2,ENORM,ZSPEC,ASPEC,BMAG,DensSPEC
      integer :: NSBESSJ,IFAIL
      COMPLEX WSPEC(3,3)
      complex :: factor

!DLG:   Define variables
        real :: genunf
      character(len=100) :: ncFileName
      integer :: nc_id, n_uper_id, n_upar_id, dfduper_id, scalar_id, &
       dfdupar_id, i_id, j_id, R_id, z_id, uper_id, upar_id, &
       fvv_id

      integer :: i_uperp, i_upara, i_psi

      nu_coll =  .01 * omgrf
      zieps0 = zi * eps0
      alpha = sqrt(2. * xkt / xm)
      rhol = alpha / omgc
      xkphi = nphi / capr
      omgrfc = omgrf * (1. + zi * dlg_xnuomg)


      xkalp = uxx * xkxsav + uxy * xkysav + uxz * xkphi
      xkbet = uyx * xkxsav + uyy * xkysav + uyz * xkphi
      xkprl = uzx * xkxsav + uzy * xkysav + uzz * xkphi
      xkperp = sqrt(xkalp**2 + xkbet**2)


!     ------------------------------------
!     Optional: leave out upshift in xkprl
!     --------------------------------- --
      if (upshift .eq. 0)         xkprl = uzz * xkphi

      if (xkprl  .eq. 0.0) xkprl  = 1.0e-08
      if (xkperp .eq. 0.0) xkperp = 1.0e-08

      sgn_kprl = sign(1.0, xkprl)
      akprl = abs(xkprl)

!     ---------------------------------
!     Calculate zetal(l) and gammab(l)
!     ---------------------------------

      do l = lmin, lmax
         labs = abs(l)

         zetal(l) = (omgrfc - l * omgc) / (xkprl * alpha)

         gammab(l) = abs(l * omgc / (2.0 * alpha * xkprl**2) &
                                                 * gradprlb / bmod)
         gamma_coll(l) = nu_coll / (akprl * alpha)

         if(xm .eq. xme)gammab(l) = 0.0
         if(abs(gammab(l)) .lt. .01)gammab(l) = .01

      enddo



!     ------------------------------------------------
!     Calculate Brambillas xkrpl_eff using l = 1 only
!     ------------------------------------------------
      y0 = 1.5
      y = y0


      if(sgn_kprl .ge. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            y = y0
            fgam = (sqrt(1. +  4. * gammab(1) * y) - 1.) &
               / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam

      end if


      if(sgn_kprl .lt. 0.0)then
         fgam = 1.0

         if(gammab(1) .gt. 1.0e-05)then
            descrim = 1. - 4. * gammab(1) * y0
            if (descrim .ge. 0.0) y =   y0
            if (descrim .lt. 0.0) y = - y0
            fgam = (1. - sqrt(1. -  4. * gammab(1) * y) ) &
               / (2. * gammab(1) * y)
         endif

         xkprl_eff = xkprl / fgam

      end if


!     Maxwellian distribution
!     -----------------------

    maxwellian: &
    if ( ndist .eq. 0 ) then

        gamma = 0.5 * xkperp**2 * rhol**2
        rgamma = real(gamma)

        if(rgamma .ge. 1.0e-08) &
        call besiexp(gamma, lmax, exil, exilp, lmaxdim, exilovergam)

        if(rgamma .lt. 1.0e-08) &
        call bes_expand(gamma, lmax, exil, exilp, lmaxdim, exilovergam)

        sig0 = 0.0
        sig1 = 0.0
        sig2 = 0.0
        sig3 = 0.0
        sig4 = 0.0
        sig5 = 0.0

        do l = lmin, lmax

            labs = abs(l)

            if(nzfun .eq. 0) call z_approx(sgn_kprl, zetal(l), 0.0, z0, z1, z2)
            if(nzfun .eq. 1) call z_approx(sgn_kprl,zetal(l),gammab(l), z0, z1, z2)
            if(nzfun .eq. 2) call z_smithe(sgn_kprl,zetal(l),gammab(l), z0, z1, z2)
            if(nzfun .eq. 3) call z_table(sgn_kprl,zetal(l),gammab(l), gamma_coll(l), z0, z1, z2)

            al = 1.0 / (xkprl * alpha) * z0
            bl = 1.0 / (xkprl * alpha) * z1
            cl = 1.0 / (xkprl * alpha) * z2

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

    endif maxwellian


    sig1 = sig1 + delta0 * eps0 * omgrfc * xkperp**2 / xk0**2
    sig3 = sig3 + delta0 * eps0 * omgrfc * xkperp**2 / xk0**2


    if (xm .eq. xme) then

        kr = xkperp / xk_cutoff
        step = damping * kr**16 / (1. + kr**16)

        sig3 = sig3 * (1.0 + step)
    end if


!   Swanson's rotation (original):
!   -----------------------------

    sigxx = sig1 + sig0 * xkbet**2
    sigxy = sig2 - sig0 * xkbet * xkalp
    sigxz = sig4 * xkalp + sig5 * xkbet

    sigyx = - sig2 - sig0 * xkbet * xkalp
    sigyy =   sig1 + sig0 * xkalp**2
    sigyz =   sig4 * xkbet - sig5 * xkalp

    sigzx = sig4 * xkalp - sig5 * xkbet
    sigzy = sig4 * xkbet + sig5 * xkalp
    sigzz = sig3

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
  102 format(2i10, 1p8e12.4)
  103 format(4i10, 1p8e12.4)

end subroutine sigma_maxwellian


      subroutine sigmac_stix(i, j, n, m, &
         xm, xn, dlg_xnuomg, &
         xkt, omgc, omgp2, &
         lmin, lmax, nzfun, ibessel, &
         xkxsav, xkysav, nphi, capr, &
         bx,    by,    bz, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         sigxx, sigxy, sigxz, &
         sigyx, sigyy, sigyz, &
         sigzx, sigzy, sigzz, &
         delta0, omgrf, xk0 )

        use constants 

!     ----------------------------------------------------
!     This routine calculates sigma_cold in the Stix frame
!     ----------------------------------------------------

      implicit none

      integer lmin, lmax, nzfun, lmaxdim, l, labs, ibessel, &
          i, j, n, m, nphi

      real xkperp, xkprl, xm, xn, xkt, omgc, omgp2, xkb, akb
      real dakbdkb
      real akprl, gammab, alpha, omgrf
      real cosalp, sinalp, a, b
      real bx, by, bz
      real xkxsav, xkysav, capr
      real xkphi
      real xkalp, xkbet, dlg_xnuomg, xk0, delta0

      complex omgrfc


      complex sig0, sig1, sig2, sig3, sig4, sig5

      complex sigxx, sigxy, sigxz, &
              sigyx, sigyy, sigyz, &
              sigzx, sigzy, sigzz

      real uxx, uxy, uxz, &
           uyx, uyy, uyz, &
           uzx, uzy, uzz


      parameter (lmaxdim = 99)

      complex tpl, tml, tplp, tplpp, zetalp, zetalm, zieps0, &
         bl, gamma





      zieps0 = zi * eps0
      xkphi = nphi / capr
      omgrfc = omgrf * (1. + zi * dlg_xnuomg)


      xkalp = uxx * xkxsav + uxy * xkysav + uxz * xkphi
      xkbet = uyx * xkxsav + uyy * xkysav + uyz * xkphi
      xkprl = uzx * xkxsav + uzy * xkysav + uzz * xkphi



      sig0 = 0.0
      sig1 = zieps0 * omgrfc * omgp2 / (omgrfc**2 - omgc**2)
      sig2 = - eps0 * omgc   * omgp2 / (omgrfc**2 - omgc**2)
      sig3 = zieps0 * omgp2 / omgrfc
      sig4 = 0.0
      sig5 = 0.0


!     -------------------
!     Swanson's rotation:
!     -------------------
      sigxx = sig1 + sig0 * xkbet**2
      sigxy = sig2 - sig0 * xkbet * xkalp
      sigxz = sig4 * xkalp + sig5 * xkbet

      sigyx = - sig2 - sig0 * xkbet * xkalp
      sigyy =   sig1 + sig0 * xkalp**2
      sigyz =   sig4 * xkbet - sig5 * xkalp

      sigzx = sig4 * xkalp - sig5 * xkbet
      sigzy = sig4 * xkbet + sig5 * xkalp
      sigzz = sig3


      return

  101 format(i10, 1p8e12.4)
 1314 format(4i10, 1p9e12.4)
 1312 format(1p9e12.4)
  100 format('ier = ', i5, 'besic failed')
      end subroutine sigmac_stix

!
!***************************************************************************
!


      subroutine besiexp(gamma, lmax, expbes, expbesp, lmaxdim, &
         expbesovergam)

        use bessel_mod

!     ----------------------------------------------------------
!     Calculates exp(-gamma) times the modified bessel functions
!     (and derivatives) of order up to lmax
!     ----------------------------------------------------------

      implicit none

      integer lmax, nmax, ier, l, lmaxdim
      real gammod
      complex gamma, expbes(0: lmaxdim), expbesp(0: lmaxdim), &
         expbesovergam(0: lmaxdim), &
         xil(0: lmaxdim), xilp(0: lmaxdim), exgam

      complex b(100)

      exgam = exp(-gamma)
      gammod = cabs(gamma)

      if(gammod .le. 700.)then
         nmax = lmax + 1
         call besic(gamma, nmax, b, ier)
         if(ier .ne. 0)write(6,100) ier

         do l = 0, lmax
            xil(l) = b(l+1)
         end do

         do l = 0, lmax
           if(l .eq. 0) xilp(0) = xil(1)
           if(l .ne. 0) xilp(l) = xil(l-1) - l / gamma * xil(l)
           expbes(l) = exgam * xil(l)
           expbesp(l) = exgam * xilp(l)
         end do
      end if

      if(gammod .gt. 700.)then
         do l = 0, lmax
            call bes_asym(gamma, l, expbes(l), expbesp(l))
         end do
      end if

      do l = 0, lmax
         expbesovergam(l) = expbes(l) / gamma
      end do

  100 format('ier = ', i5, 'besic failed')
      return
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

      integer lmax, nmax, ier, l, lmaxdim
      real :: factl
      complex gamma, expbes(0:lmaxdim), expbesp(0:lmaxdim), &
         expbesovergam(0:lmaxdim), xilovergam, &
         xil(0:lmaxdim), xilp(0:lmaxdim), exgam

      complex b(100)

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

  100 format(i10, 1p8e12.4)

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
