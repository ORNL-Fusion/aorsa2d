module qlsum_myra_gpu_mod
contains

  !
  !*************************************************************************
  !
  subroutine QLSUM(k_uper, b_sum, c_sum, e_sum, f_sum, &
       & sum_wdot, sum_fx0, sum_fy0, W, ZSPEC, ASPEC, BMAG, &
       & lmax, ENORM, UPARMIN, UPARMAX, &
       & NUPAR, NUPER, UPER, UPAR, DFDUPER, DFDUPAR,  &
       & ealphak, ebetak, ebk, nkdim1, nkdim2, mkdim1, mkdim2,   &
       & nkx1, nkx2, nky1, nky2, &
       & uxx, uxy, uxz, &
       & uyx, uyy, uyz, &
       & uzx, uzy, uzz, &
       & nxdim, nydim, xkxsav, xkysav, xkphi, xx, yy, i_global, j_global, &
       & lmaxdim, ndist, nzeta, &
       & gradprlb, bmod, omgc, alpha, xm, upshift, xk_cutoff, rt, nphi)

    implicit none

    integer i_global, j_global, ier, nxdim, nydim, k, lmaxdim, ndist
    integer, intent(IN):: NUPAR, NUPER, lmax
    integer nkx1, nkx2, nky1, nky2, j_upar, k_uper, l
    integer nkdim1, nkdim2, mkdim1, mkdim2
    integer:: NHARM, IHARM, M, N, i, nzeta
    integer i_uprl, upshift
    integer ires, iresmax, nphi

    complex, dimension(:,:), allocatable :: zbeta
    complex, dimension(:,:), allocatable :: zbeta_iharm

    real  y, y0, alpha, xm, akprl, sgn_kprl, omgc, bmod, xkprl_eff,  &
         &   descrim, xme, fgam, gradprlb, gammab, xk_cutoff, rt

    real  uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz
    real  xkphi, sinth(nuper), factc, facte(nuper), factf, sinth_inv(nuper)
    real  xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)
    real  xkperpn, xkperpni, xkrhon, xketan, xkprln, beta
    real, intent(IN):: W, ZSPEC, ASPEC, BMAG
    real, intent(IN):: ENORM, UPARMIN, UPARMAX
    real, dimension(NUPER), intent(IN):: UPER
    real, dimension(NUPAR), intent(IN):: UPAR
    real, dimension(NUPER,NUPAR), intent(IN):: DFDUPER, DFDUPAR
    real:: W2, WCW, RRP, RRM, WC, WCI
    real:: MUT0, SQMUT0, PISQMUT0, SQMUT0I, KPARA1, NPARA1
    real:: ISQ2, SQ2, NWCW, DFACTPAR, DFACTPER, U0(nuper)
    real:: UPAR0, dfdupar0, dfduper0(nuper), du, dui, p 
    real:: time, t1, tmsec, second1, dummy, WI, uperpk(nuper), uperpk2(nuper)
    real:: dzeta, dzetai, zetamax, zetamin, zeta0, zetamax1, zetamax2
    real:: A1, A2, A3, u(nuper)
    real:: temp1, temp2, temp3, factor
    real:: temp1w, temp2w, temp3w

    real :: watchvar

    real, dimension(:),     allocatable :: zetai
    real, dimension(:,:,:),   allocatable :: Jni
    real, dimension(:,:,:,:), allocatable :: Jn
    real, dimension(:,:),   allocatable :: NPARA_sav

    integer, dimension(:),  allocatable :: nres
    integer, dimension(:),  allocatable :: mres

    complex epsx, epsy, epsz
    complex ealphak(nkdim1 : nkdim2, mkdim1 : mkdim2), &
         &        ebetak(nkdim1 : nkdim2, mkdim1 : mkdim2), &
         &           ebk(nkdim1 : nkdim2, mkdim1 : mkdim2)
    complex cexpkx, cexpky, zi, zeta

    real xkperpn_tmp(nkx1 : nkx2, nky1 : nky2), &
         zetamin_tmp(nuper), &
         dzetai_tmp(nuper)

    complex xx(nkdim1 : nkdim2, 1 : nxdim),   &
         &      yy(mkdim1 : mkdim2, 1 : nydim)

    complex cexpn, cexpnp1, cexpnm1, cexp11

    complex cexp1, cexp2, cexp0

    complex b_sum(nuper, nupar), c_sum(nuper, nupar), e_sum(nuper, nupar), f_sum(nuper, nupar)
    complex sum_wdot(nuper), sum_fx0(nuper), sum_fy0(nuper)
    complex b(100)

    real, parameter::  EOVERAMU = 9.58084e7
    real, parameter:: MPC2 = 938271998.38
    real, parameter:: C = 2.99792458e8
    real, parameter:: PI = 3.141592653597932384
    real :: cosbeta_n_m, sinbeta_n_m

    allocate( zbeta(nkx1:nkx2,nky1:nky2) )
    allocate( zbeta_iharm(nkx1:nkx2,nky1:nky2) )

    allocate(zetai(nzeta + 1) )
    allocate(Jni(nuper, -lmaxdim : lmaxdim, nzeta + 1) )
    allocate(NPARA_sav(nkx1 : nkx2, nky1 : nky2) ) 

    allocate(nres(nxdim * nydim) )
    allocate(mres(nxdim * nydim) )

    !   -------------------------------------
    !   initialize allocatable arrays to zero
    !   -------------------------------------
    zbeta = 0.0
    zbeta_iharm = 0.0

    zetai = 0.0
    Jni = 0.0
    NPARA_sav = 0.0

    xme = 9.11e-31
    zi = cmplx(0., 1.)

    do k_uper = 1, nuper

       uperpk = uper(k_uper)
       uperpk2 = uperpk**2

       W2 = W * W
       WI = 1.0 / W

       WCW = BMAG * ZSPEC * EOVERAMU / ASPEC / W
       WC = WCW * W
       WCI = 1.0 /WC

       MUT0 = 0.5 * MPC2 * ASPEC / ENORM
       SQMUT0 = SQRT(MUT0)
       SQMUT0I = 1.0 / SQMUT0
       PISQMUT0 = SQMUT0 * pi

       ISQ2 = SQRT(0.5)
       SQ2 = SQRT(2.0)
       NHARM = lmax

       du = (upar(nupar) - upar(1)) / (nupar - 1)
       dui = 1.0 / du

       ! -------------------------------------- !
       ! ---Interpolate; calculate zeta mesh -- !
       ! -------------------------------------- !

       zetamax = 0.0
       zetamin = 0.0

       do n = nkx1, nkx2
          do m = nky1, nky2
             xkrhon = uxx * xkxsav(n) + uxy * xkysav(m) + uxz * xkphi
             xketan = uyx * xkxsav(n) + uyy * xkysav(m) + uyz * xkphi
             !            xkprln   = uzx * xkxsav(n) + uzy * xkysav(m) + uzz * xkphi

             xkperpn = sqrt(xkrhon**2 + xketan**2) + 1.0e-08

             zeta0 = xkperpn * uperpk(k_uper) * c * sqmut0i * wci

             if (zeta0 .gt. zetamax) zetamax = zeta0
             if (zeta0 .lt. zetamin) zetamin = zeta0
          end do
       end do

       if(zetamax .eq. zetamin)then
          zetamax =  1.0e-06
          zetamin = -1.0e-06
       end if

       dzeta = (zetamax - zetamin) / (nzeta - 1)
       dzetai = 1.0 / dzeta

       zetamin_tmp(k_uper) = zetamin
       dzetai_tmp(k_uper) = dzetai

       ! ------------------------------------------------- !
       ! ---Pre-calculate Bessel functions on zeta mesh -- !
       ! ------------------------------------------------- !

       do i = 1, nzeta + 1
          zetai(i) = zetamin + (i - 1) * dzeta
          zeta = cmplx(zetai(i), 0.0)

          call besjc(zeta, nharm + 2, b, ier)

          !          if(ier .ne. 0) write(6, *) "ier = ", ier

          do iharm = 0, NHARM + 1
             Jni(k_uper, iharm,  i) = b(iharm + 1)
             Jni(k_uper, -iharm, i) = (-1.0)**iharm * b(iharm + 1)
          end do
       end do
    end do


    ! -------------------------------------------- !
    ! ---Prepare to Interpolate Bessel functions-- !
    ! -------------------------------------------- !

    do n = nkx1, nkx2
       do m = nky1, nky2

          xkrhon = uxx * xkxsav(n) + uxy * xkysav(m) + uxz * xkphi
          xketan = uyx * xkxsav(n) + uyy * xkysav(m) + uyz * xkphi
          xkprln = uzx * xkxsav(n) + uzy * xkysav(m) + uzz * xkphi

          ! ------------------------------------
          ! Optional: leave out upshift in xkprl
          ! --------------------------------- --          
          !             if(upshift .eq. 0)xkprln = uzz * xkphi
          if (upshift .eq. 0) xkprln = nphi / rt

          if (upshift .eq. -1) then      
             if (xkperpn  .gt. xk_cutoff) xkprln = uzz * xkphi
          end if

          sgn_kprl = sign(1.0, xkprln)
          akprl = abs(xkprln)

          y0 = 1.5
          y = y0

          l = 1
          if(xkprln.eq. 0)xkprln = 1.0e-06
          gammab = abs(l * omgc / (2.0 * alpha * xkprln**2)  &
               &                                    * gradprlb / bmod)

          if(xm .eq. xme)gammab = 0.0
          !             if(abs(gammab) .gt. 1000.0) gammab = 1000.0
          if(abs(gammab) .lt. .01)gammab = .01


          if(sgn_kprl .ge. 0.0)then
             fgam = 1.0

             if(gammab .gt. 1.0e-05)then
                y = y0
                fgam = (sqrt(1. +  4. * gammab * y) - 1.)   &
                     &               / (2. * gammab * y)
             endif

             xkprl_eff = xkprln / fgam 

          end if


          if(sgn_kprl .lt. 0.0)then
             fgam = 1.0

             if(gammab .gt. 1.0e-05)then
                descrim = 1. - 4. * gammab * y0
                if (descrim .ge. 0.0) y =   y0
                if (descrim .lt. 0.0) y = - y0
                fgam = (1. - sqrt(1. -  4. * gammab * y) )  &
                     &                / (2. * gammab * y)
             endif

             xkprl_eff = xkprln / fgam 

          end if

          if (upshift .ne. 0) xkprln = xkprl_eff

          NPARA_sav(n, m) = xkprln * C * WI

          xkperpn = sqrt(xkrhon**2 + xketan**2) + 1.0e-08
          xkperpni = 1.0 / xkperpn

          cosbeta_n_m = xkrhon * xkperpni
          sinbeta_n_m  = xketan * xkperpni
          zbeta(n,m) = cmplx( cosbeta_n_m , sinbeta_n_m )

          xkperpn_tmp(n, m) = xkperpn * c * sqmut0i * wci
       end do
    end do

    ! ------------------------ !
    ! ---Sum over harmonics--- !
    ! ------------------------ !

    sum_wdot = 0.0
    sum_fx0  = 0.0
    sum_fy0  = 0.0    

    b_sum = 0.0
    c_sum = 0.0
    e_sum = 0.0
    f_sum = 0.0

    call zpow((nkx2-nkx1+1)*(nky2-nky1+1), zbeta, -nharm - 1, zbeta_iharm) 

    call qlsum_gpu(nkx1, nkx2, &
         & nky1, nky2, &
         & nharm, &
         & nuper, nupar, &
         & nres, mres, &
         & wcw, &
         & uparmax, uparmin, &
         & uper, upar, &
         & dfduper, dfdupar, &
         & dui, sqmut0, &
         & xkxsav(nkdim1), xkysav(nkdim1), &
         & npara_sav(nkx1, nky1), Jni(1, -lmaxdim, 1), &
         & xkperpn_tmp(nkx1, nky1), zetamin_tmp(1), &
         & zetai(1), dzetai_tmp(1), &
         & zbeta(nkx1, nky1), zbeta_iharm(nkx1, nky1), &
         & xx(nkx1, i_global), yy(nky1, j_global), &
         & ealphak(nkdim1, mkdim1), ebetak(nkdim1, mkdim1), ebk(nkdim1, mkdim1), &
         & sum_wdot, sum_fx0, sum_fy0, &
         & b_sum, c_sum, e_sum, f_sum)

    deallocate( zbeta )
    deallocate( zbeta_iharm )

    deallocate(nres )
    deallocate(mres )

    deallocate(zetai)
    deallocate(Jni)
    deallocate(NPARA_sav)

    return

  end subroutine QLSUM

  !
  !*************************************************************************
  !

  subroutine zpow(n, z, iharm, zout )
    implicit none
    integer n, iharm
    complex z(n), zout(n)

    integer i,ipow
    complex one, zero
    complex zin

    integer nharm
    logical isodd
    intrinsic mod

    logical use_zdiv
    parameter(use_zdiv=.true.)

    integer nb
    parameter(nb=1024*4*4)
    complex zk(nb)
    integer istart,iend,isize

    one = 1.0d0
    zero = 0.0d0

    if (iharm.eq.0) then
       do i=1,n
          zout(i) = one
       enddo
       return
    endif


    do istart=1,n,nb

       iend = min(n,istart+nb-1)
       isize = iend-istart+1

       do i=1,isize
          zout(istart-1+i) = one
       enddo


       do i=1,isize
          zk(i) = z(istart-1+i)
       enddo

       nharm = abs(iharm)
       do while (nharm .gt. 0)
          isodd = (mod(nharm,2).eq.1)
          if (isodd) then
             do i=1,isize
                zout(istart-1+i) = zout(istart-1+i) * zk(i)
             enddo
          endif
          do i=1,isize
             zk(i) = zk(i) * zk(i)
          enddo
          nharm = int( nharm/2 )
       enddo



       if (iharm.lt.0) then
          if (use_zdiv) then
             do i=1,isize
                zin = zout(istart-1+i)
                call zdiv( zin, zout(istart-1+i) )
             enddo
          else
             do i=1,isize
                zin = zout(istart-1+i)
                zout(istart-1+i) = one/zin
             enddo
          endif
       endif

    enddo

    return
  end subroutine zpow

  !
  !*************************************************************************
  !

  subroutine zdiv( zin, zout )
    implicit none
    complex zin, zout
    real a, b
    real d

    real one
    parameter(one=1.0d0)
    real rd, a_over_b, b_over_a


    a = real(zin)
    b = aimag(zin)

    !       z = (a + i * b)
    !       1/z =  a/(a^2 + b^2) - i * b/(a^2 + b^2)
    !
    !       or    1/(a + (b/a)*b) - i * (b/a) / (a + (b/a)*b)
    !       or    (a/b)/( (a/b)*a + b ) - i * 1/( (a/b)*a + b )
    !        
    if (abs(a).gt.abs(b)) then
       b_over_a = b/a
       d = a + (b_over_a)*b
       rd = one/d
       zout = cmplx( rd, -(b_over_a)*rd )
    else
       a_over_b = a/b
       d = (a_over_b)*a + b
       rd = one/d
       zout = cmplx( (a_over_b)*rd, -rd )
    endif

    return
  end subroutine zdiv



  !
  !*************************************************************************
  !


end module qlsum_myra_gpu_mod

