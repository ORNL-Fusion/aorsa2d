#define ERROR(string) write(*,'(/,a,a,/,a,a,a,i5)') \
    "Error: ", string, "In ",__FILE__,":", __LINE__ ;stop

module zfunction_mod

contains


!
!*******************************************************************************
!

      subroutine z_exact(sgn_kprl, zeta, gamma_, xnu, z0, z1, z2, &
         dz0, dz1, dz2)

!     ----------------------------------------------------
!     z_exact is a wrapper routine for David Smithe's own
!     version of "z_smithe"
!     ----------------------------------------------------
        use ztable_mod

      implicit none

      complex, intent(inout)   :: z0, z1, z2
      complex, intent(out), optional :: dz0, dz1, dz2
      complex, intent(in)    :: zeta
      real, intent(in)       :: sgn_kprl, gamma_, xnu
      complex, dimension(0:4)  :: zfunc

      real alpha_s
      complex zi

      zi = cmplx(0.0, 1.0)

      alpha_s = -2.0 * gamma_

      if ( sgn_kprl .ge. 0.0 ) then

         call ztable_zfunc( real(zeta), alpha_s, xnu, zfunc )
         z0 = zfunc(0)
         z1 = - 0.5 * zfunc(1) - gamma_ / (2.0 * zi) * zfunc(2)
         z2 = 0.5 * zfunc(0) + 0.25 * zfunc(2) &
            + gamma_ / (2.0 * zi) * zfunc(3) &
            - gamma_**2 / 4. * zfunc(4)

        if(present(dz0)) dz0 = zfunc(1)
        if(present(dz1)) dz1 = - 0.5 * zfunc(2) - gamma_ / (2.0 * zi) * zfunc(3)
!         dz2 = 0.5 * zfunc(1) + 0.25 * zfunc(3)
!     .      + gamma_ / (2.0 * zi) * zfunc(4)



      else if ( sgn_kprl .lt. 0.0 ) then

         call ztable_zfunc( - real(zeta), - alpha_s, - xnu, zfunc )
         z0 = - zfunc(0)
         z1 = - 0.5 * zfunc(1) + gamma_ / (2.0 * zi) * zfunc(2)
         z2 = - 0.5 * zfunc(0) - 0.25 * zfunc(2) &
            + gamma_ / (2.0 * zi) * zfunc(3) &
            + gamma_**2 / 4. * zfunc(4)

            if(present(dz0)) dz0 =  zfunc(1)
            if(present(dz1)) dz1 =  0.5 * zfunc(2) - gamma_ / (2.0 * zi) * zfunc(3)
!         dz2 =  0.5 * zfunc(1) + 0.25 * zfunc(3)
!     .      - gamma_ / (2.0 * zi) * zfunc(4)

      endif

      return

      end subroutine z_exact


!
!***************************************************************************
!

subroutine z_approx_dlg ( sgn_kprl, zeta, z, zPrime, zeta_zPrime )

    use hammett, zfun_hammett => zfun

    implicit none

    real, intent(in) :: sgn_kPrl
    complex, intent(in) :: zeta(:)
    complex, intent(out) :: z(:), zPrime(:), zeta_zPrime(:)

    real, parameter :: y0 = 1.5
    real, allocatable :: y(:), descrim(:)
    complex, allocatable :: zetat(:), zFunct(:)
    complex :: fzeta
    complex :: zFunctCheck
    integer :: nL, i

    nL = size ( zeta )

    allocate ( zFunct(nL), y(nL), descrim(nL) )
    y = y0

    ! See pg. 202 of Stix.

    if(sgn_kprl>=0.0)then

      do i=1,nL 
          zFunct(i) = zfun_hammett ( zeta(i) )
      enddo

       z = zfunct 
       zPrime = -2d0 * (1.0 + zeta * z) ! pg. 201 of Stix.
       zeta_zPrime = zeta * zPrime

    else

      do i=1,nL
          zfunct(i) = zfun_hammett (-zeta(i))
      enddo

       z = -zfunct 
       !Not sure about which is correct yet ...
       zPrime = -2d0 * (1.0 - zeta * zfunct) 
       ! or is it
       !zPrime = 2d0 * (1.0 - zeta * z)

       zeta_zPrime = zeta * zPrime 

    endif

    return

end subroutine z_approx_dlg


subroutine z_approx(sgn_kprl, zeta, gamma_, z0, z1, z2)

    ! This routine is for the Brambilla toroidal broadening
    ! method ... NOT the Smithe approach. i.e., gamma_ is 
    ! gamma_brambilla

    use hammett, zfun_hammett => zfun

    implicit none

    real, intent(in) :: sgn_kPrl, gamma_(:)
    complex, intent(in) :: zeta(:)
    complex, intent(out) :: z0(:), z1(:), z2(:)

    real, parameter :: y0 = 1.5
    real, allocatable :: fGam(:), y(:), descrim(:)
    complex, allocatable :: zetat(:), zFunct(:)
    complex :: fzeta
    complex :: zFunctCheck
    integer :: nL, i

    if(any(gamma_<0)) then
        ERROR("gamma_ should be >= 0")
    endif

    nL = size ( zeta )

    allocate ( zFunct(nL), fGam(nL), y(nL), descrim(nL) )
    y = y0
    fGam = 1.0

    if(sgn_kprl>=0.0)then

      where(gamma_>1.0e-5)
          fGam = (sqrt(1.0 +  4.0 * gamma_ * y0) - 1.0) / (2.0 * gamma_ * y0)
      endwhere

      do i=1,nL 
          zFunct(i) = zfun_hammett ( fGam(i) * zeta(i) )
      enddo

       z0 = zfunct * fGam
       z1 = -2d0 * (1.0 + zeta * z0) * fGam
       z2 = zeta * z1 * fGam

    else

      descrim = 1.0 - 4.0 * gamma_ * y0
      where(descrim>=0.0)
              y = y0
      elsewhere
              y = -y0
      endwhere

      where(gamma_>1.0e-5)
          fGam = (1.0 - sqrt(1.0 -  4.0 * gamma_ * y) ) / (2.0 * gamma_ * y)
      endwhere

      do i=1,nL
          zfunct(i) = zfun_hammett (-fGam(i)*zeta(i))
          !zfunct(i) = zfun_hammett (-zeta(i))
      enddo

       z0 = -zfunct * fGam
       z1 = 2d0 * y / abs(y) * (1.0 + zeta * z0) * fGam
       z2 = zeta * z1 * fGam

    endif

    return

end subroutine z_approx

!
!******************************************************************************
!


      subroutine z_approx_im(sgn_kprl, zeta, gamma, z0, z1, z2)

      use ztable_mod
      !use aorsasubs_mod

      implicit none

!     ------------------
!     Finite k_parallel:
!     ------------------
      real gamma, fGam, y0, y, sgn_kprl, descrim
      complex zeta, z0, z1, z2, zfunct, zetat, fzeta


      y0 = 1.5
      y = y0


      if(sgn_kprl .ge. 0.0)then
         fGam = 1.0

         if(gamma .gt. 1.0e-05)then
            y = y0
            fGam = (sqrt(1. +  4. * gamma * y) - 1.) &
               / (2. * gamma * y)
         endif

         zetat = fGam * zeta
!        zfunct = fzeta(zetat)
         call zfun_im (zetat, zfunct)
         z0 = fGam * zfunct
         z1 = fGam * (1.0 + fGam * zeta * zfunct)
         z2 =  fGam**2 * zeta * (1.0 + fGam * zeta * zfunct)
      end if


      if(sgn_kprl .lt. 0.0)then
         fGam = 1.0

         if(gamma .gt. 1.0e-05)then
            descrim = 1. - 4. * gamma * y0
            if (descrim .ge. 0.0) y =   y0
            if (descrim .lt. 0.0) y = - y0
            fGam = (1. - sqrt(1. -  4. * gamma * y) ) &
               / (2. * gamma * y)
         endif


         zetat = - fGam * zeta
!         zfunct = fzeta( zetat)
         call zfun_im (zetat, zfunct)
         z0 = - fGam * zfunct
         z1 = y / abs(y) * fGam * (1.0 - fGam * zeta * zfunct)
         z2 =  fGam**2 * zeta * (1.0 - fGam * zeta * zfunct)
      end if


      return
      end subroutine z_approx_im

!
!******************************************************************************
!

      subroutine z_approx0(sgn_kprl, alpha, xkprl_eff, zeta_eff, &
         al, bl, cl)

      use ztable_mod
      use hammett, zfun_hammett => zfun


      implicit none

!     -----------------------
!     Small k_parallel limit:
!     -----------------------
      real sgn_kprl, alpha, xkprl_eff
      complex zeta_eff, al, bl, cl, zfunct, zetat, fzeta

      if(sgn_kprl .ge. 0.0)then
!        zfunct = fzeta(zeta_eff)
         !call zfun (zeta_eff, zfunct)
         zfunct = zfun_hammett ( zeta_eff )

         al = 1. /(xkprl_eff * alpha) * zfunct
         bl = 1. /(xkprl_eff * alpha) * (1. + zeta_eff * zfunct)
         cl = zeta_eff /(xkprl_eff * alpha) *(1. + zeta_eff * zfunct)
      end if


      if(sgn_kprl .lt. 0.0)then
!        zfunct = fzeta(-zeta_eff)
         !call zfun (-zeta_eff, zfunct)
         zfunct = zfun_hammett ( -zeta_eff )

         al = - 1. /(xkprl_eff * alpha) * zfunct
         bl = - 1. /(xkprl_eff * alpha) * (1. - zeta_eff * zfunct)
         cl = zeta_eff /(xkprl_eff * alpha) * (1. - zeta_eff * zfunct)
            end if


      return
      end subroutine z_approx0


!
!***************************************************************************
!
      subroutine z2pole(zeta, zfunct)

      complex zeta, zfunct, zi, znum
      data pi/3.141592654/
      data zi/(0.0,1.0)/
      data g/1.141592654/
      data sqpi/1.772453851/
      znum = zi*sqpi + g*zeta
      zfunct = znum/(1.0 - zeta*znum)
      return
      end subroutine z2pole
!
!******************************************************************************
!

      subroutine z_smithe(sgn_kprl, zeta, gamma, z0, z1, z2)
        use ztable_mod

      implicit none

      real dx, xmax, xnueff, gamma, sgn_kprl, xmax0
      integer j, npts, nmaxdim
      real x(10000), zeta_mod
      complex z0, z1, z2, zi, f(10000), zeta
      complex cexpx(10000)

      parameter (nmaxdim = 10000)

      zi = cmplx(0.,1.)
      xnueff = 0.0
      zeta_mod = sqrt(conjg(zeta) * zeta)


      if(zeta_mod .gt. 1.0e+02)then

         z0 = - 1.0 / zeta - 0.5 / zeta**3 &
                                + 3.0 * gamma / (zi * zeta**4)
         z1 = - 0.5 / zeta**2 + gamma / (zi * zeta**3) &
                                             - 0.75 / zeta**4
         z2 = - 0.5 / zeta - 0.75 / zeta**3 &
                                + 4.5 * gamma / (zi * zeta**4)
         return

      else

         xmax0 = 7.0
         npts = 10000
         xmax = sgn_kprl * xmax0
         dx = xmax / (npts - 1)


!        ------------
!        calculate Z0
!        ------------
         do  j = 1, npts

            x(j) = (j-1) * dx
            cexpx(j) = exp(zi * zeta * x(j) &
               - x(j)**2 / 4. * (1. + gamma * x(j))**2 &
               - xnueff / 8. * x(j)**3   )

            f(j) = zi * cexpx(j)
         end do

         call dgrat1(x, f, 1, npts, z0, nmaxdim)


!        ------------
!        calculate Z1
!        ------------
         do  j = 1, npts
            f(j) = 0.5 * x(j)* (1.0 + gamma * x(j)) * cexpx(j)
         end do

         call dgrat1(x, f, 1, npts, z1, nmaxdim)


!        ------------
!        calculate Z2
!        ------------
         do  j = 1, npts
            f(j) = zi / 2.0 * (1.0 -  x(j)**2 / 2. &
                 * (1. + gamma * x(j))**2) * cexpx(j)
         end do

         call dgrat1(x, f, 1, npts, z2, nmaxdim)



         return

      end if

      end subroutine z_smithe

!
!*******************************************************************************
!

      subroutine dgrat1(x, f, nx1, nx2, ans, nmaxdim)

      implicit none

      integer n, nx1, nx2, nmaxdim
      real x(nmaxdim), dx
      complex f(nmaxdim), ans

      ans = 0.0

      do n = nx1, nx2 - 1
         dx = x(n+1) - x(n)
         ans = ans + dx * (f(n) + f(n+1)) / 2.0
      end do

      return
      end subroutine dgrat1



      SUBROUTINE ZFUN_IM(Z,FU)

!     ---------------------------------------------------
!     Evaluates only the imaginary part of the Z function
!     ---------------------------------------------------

      implicit none

      complex z, fu, zi
      real pi

      pi = 3.14159
      zi = cmplx (0.0, 1.0)

      fu = zi * sqrt(pi) * cexp (-z*z)

      return
      end SUBROUTINE ZFUN_IM


end module zfunction_mod
