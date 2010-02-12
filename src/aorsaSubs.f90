module aorsasubs_mod

contains

       subroutine maxwell_dist(u0, NUPAR, NUPER, &
                             UminPara, UmaxPara, &
                             UPERP, UPARA, DFDUPER, DFDUPAR)
       implicit none

       integer, intent(in) ::NUPAR,NUPER
       integer n,m

       real :: UminPara,UmaxPara
       real :: UPERP(NUPER)
       real :: UPARA(NUPAR)
       real, intent(inout) :: DFDUPER(NUPER,NUPAR), DFDUPAR(NUPER,NUPAR)
       real  pi, fnorm, f, u2, u0, u02

       parameter (PI = 3.141592653597932384)

       fnorm = u0**3 / pi**1.5
       u02 = u0**2

       do n = 1, NUPER
          do m = 1, NUPAR
             u2 = UPERP(n)**2 + UPARA(m)**2
             f = fnorm * exp(-u2 * u02)
             DFDUPER(n, m) = f *(-2.0 * UPERP(n) * u02)
             DFDUPAR(n, m) = f *(-2.0 * UPARA(m) * u02)
          end do
       end do


       end subroutine maxwell_dist


!
!***************************************************************************
!


      subroutine besjc (z,n,b,ier)

      use bessel_mod
!
! subroutine besjc (z,n,b,ier)
!
! dimension of           b(n+1)
! arguments
!
! latest revision        october 1978
!
! purpose                to calculate j bessel functions for complex
!                        argument and integer order
!
! usage                  call besjc (z,n,b,ier)
!
! arguments
!
! on input               z
!                          complex argument for which j bessel functions
!                          are to be calculated.  abs(aimag(z)) must be
!                          less than the largest real argument that the
!                          fortran function exp can handle.
!
!                        n
!                          integer, the highest order to be calculated.
!                          n must be greater than or equal to zero.
!
! on output              b
!                          complex vector of length n+1 containing the
!                          bessel function values j-sub-0(z),j-sub-1(z),
!                          ...,j-sub-n(z) in b(1),b(2),...,b(n+1).
!                        ier
!                          an integer error flag.
!                          =0 if all desired orders have been calculated
!                             satisfactorily,
!                          =1 if abs(aimag(z)) is too large,
!                          =2 if n is less than zero,
!                          =2+k if only the first k results are correct.
!                             in the returned values b(m) for m greater
!                             than k, approximately the last
!                             alog10(abs(b(m)/b(k))) significant digits
!                             are in error.
!
! entry points           besjc
!
! special conditions     none
!
! common blocks          none
!
! i/o                    none, except for error messages produced by
!                        calling the error handling routine uliber.
!
! precision              single
!
! specialist             russ rew, ncar, boulder, colorado
!
! accuracy               in tests run on the 7600 with orders from 0
!                        through 10 using random values of the argument
!                        with absolute value less than 50 but
!                        concentrated around the origin, the maximum
!                        relative error (or absolute error when it was
!                        larger) observed was about 8.2e-14.
!
! timing                 on ncar"s control data 7600, besic takes about
!                        .32+.008*n milliseconds when z=(1.0,1.0).
!
! portability            ansi 1966 standard
!
!
!
!
!
!
      complex         z         ! ,b(1)
      real, allocatable :: b(:,:)
      data iorj/0/,xlarge/1000000./
!
      nb = n+1
      allocate ( b(2,nb) )
      call b2slci (real(z),aimag(z),nb,iorj,b,ncalc)
      ier = 0
      if (ncalc .eq. nb) go to 103
      if (ncalc .ge. 0) go to 102
      if (n .ge. 0) go to 101
      ier = 2
      call uliber (ier,'in besjc, n out of range')
      go to 103
  101 ier = 1
      call uliber (ier,'in besjc, x out of range')
      go to 103
  102 ier = 2+ncalc
      call uliber (ier,'in besjc, accuracy lost for some orders')
  103 return
      end subroutine besjc 


!
!********************************************************************
!



!***************************************************************************
!

      subroutine polavg(f, favg, rho, nxdim, nydim, nrhodim, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, r0, vol, fvol)

      implicit none

      integer nxdim, nydim, nrhodim, nnodex, nnodey, nnoderho
      integer n, i, j

      real f(nxdim, nydim), favg(nrhodim), rho(nxdim, nydim), r0, &
         drho, dx, dy, fvol(nrhodim), vol(nrhodim), capr(nxdim)

      logical use_vol_ij
      parameter(use_vol_ij=.false.)

      integer n_ij(nnodex,nnodey)
      real fvol_ij(nnodex,nnodey)
      real vol_ij(nnodex,nnodey)
      logical has_n(nnodex,nnodey)

      do n = 1, nnoderho
          fvol(n) = 0.0
          vol(n) = 0.0
          favg(n) = 0.0
      end do

      if (use_vol_ij) then
         do j=1,nnodey
         do i=1,nnodex
            n = int(rho(i,j) / drho) + 1
            has_n(i,j) = (1.le.n).and.(n.le.nnoderho)
            if(has_n(i,j)) then
               n_ij(i,j) = n
               fvol_ij(i,j) =  dx * dy * f(i,j)
               vol_ij(i,j) =   dx * dy
            endif
         enddo
         enddo

!dir$ novector
         do j=1,nnodey
!dir$ novector
         do i=1,nnodex
           if (has_n(i,j)) then
              n = n_ij(i,j)
              fvol(n) = fvol(n) + fvol_ij(i,j)
              vol(n) = vol(n) + vol_ij(i,j)
           endif
         enddo
         enddo
!dir$ vector

      else

      do i = 1, nnodex
         do j = 1, nnodey
            n = int(rho(i,j) / drho) + 1
            if(n .le. nnoderho)then
               fvol(n) = fvol(n) + dx * dy * f(i,j)
               vol(n) =  vol(n) + dx * dy
            end if
         end do
      end do

      endif

      do n = 1, nnoderho
         favg(n) = 0.0
         if(vol(n) .ne. 0.0)favg(n) = fvol(n) / vol(n)
!         if(favg(n) .eq. 0.0 .and. n .gt. 1)favg(n) = favg(n-1)
      end do

      do n = 2, nnoderho -1
         if(favg(n) .eq. 0.0)then
            favg(n) = (favg(n-1) + favg(n+1)) / 2.0
         end if
      end do

      if (favg(1) .eq. 0.0) favg(1) = favg(2)
      if (favg(nnoderho) .lt. 1.e-03) favg(nnoderho) = favg(nnoderho -1)


  100 format (1i10, 1p8e12.4)
  102 format (2i10)
      return
      end subroutine polavg

!
!***************************************************************************
!

!
!***************************************************************************
!


      SUBROUTINE TRIDI(N,A,B,C,Y,X)

      implicit none

      real A(1),B(1),C(1),Y(1),X(1),D(1001),W(1001)
      real dn
      integer n, nm, j, k
      NM=N-1
      D(1)=C(1)/B(1)
      W(1)=Y(1)/B(1)
      DO 1 J=2,NM
      DN=B(J)-A(J)*D(J-1)
      D(J)=C(J)/DN
    1 W(J)=(Y(J)-A(J)*W(J-1))/DN
      X(N)=(Y(N)-A(N)*W(N-1))/(B(N)-A(N)*D(N-1))
      DO 2 J=1,NM
      K=N-J
    2 X(K)=W(K)-D(K)*X(K+1)
      RETURN
      END SUBROUTINE TRIDI


!
!***************************************************************************
!
      subroutine tridag(a, b, c, r, u, n)

      implicit none

      integer n, nmax
      real a(n), b(n), c(n), r(n), u(n)
      parameter (nmax = 500)
      integer j
      real bet, gam(nmax)

      if(b(1) .eq. 0) then
             write(*,*) 'tridag: rewrite equations'
             stop
     endif
      bet = b(1)
      u(1) = r(1) / bet

      do j = 2, n
         gam(j) = c(j-1) / bet
         bet = b(j) - a(j) * gam(j)
!         if(bet .eq. 0) pause 'tridag failed'

         if(bet .eq. 0) then
          write (6, *) "j = ", j
          write (6, *) "a(j) = ", a(j)
          write (6, *) "b(j) = ", b(j)
          write (6, *) "c(j) = ", c(j)
          write (6, *) "gam(j) = ", gam(j)
          write (6, *) "bet = ", bet
          write (6, *) "tridag failed"
          call exit
       end if


         u(j) = (r(j) - a(j) * u(j-1)) / bet
      end do

      do j = n-1, 1, -1
         u(j) = u(j) - gam(j+1) * u(j+1)
      enddo

      return
      end subroutine tridag

!
!***************************************************************************
!
    subroutine deriv_x ( f, i, j, dx, dfdx, d2fdx2 )
    
        implicit none
        
        integer, intent(in) :: i, j
        real, intent(in) :: f(:,:)
        real, intent(out), optional :: dfdx, d2fdx2
        real, intent(in) :: dx
    
        integer :: nx, ny
    
        nx    = size ( f, dim = 1 )
        ny    = size ( f, dim = 2 )
    
        if(i .ne. 1 .and. i .ne. nx)then

            dfdx = (f(i+1, j) - f(i-1, j)) / (2.0 * dx)
            d2fdx2 = (f(i+1, j) - 2.0 * f(i,j) + f(i-1, j)) / dx**2

        endif
        
        if(i .eq. 1)then

            dfdx = (f(i+1, j) - f(i, j)) / dx
            d2fdx2 = (f(i+2, j) - 2.0 * f(i+1,j) + f(i, j)) / dx**2

        endif
        
        if(i .eq. nx)then

            dfdx = (f(i, j) - f(i-1, j)) / dx
            d2fdx2 = (f(i, j) - 2.0 * f(i-1,j) + f(i-2, j)) / dx**2

        endif
    
    end subroutine deriv_x

!
!***************************************************************************
!
    subroutine deriv_y ( f, i, j, dy, dfdy, d2fdy2 )

        implicit none

        integer, intent(in) :: i, j
        real, intent(in) :: f(:,:)
        real, intent(in) :: dy
        real, intent(out), optional :: dfdy, d2fdy2

        integer :: nx, ny

        nx    = size ( f, dim = 1 )
        ny    = size ( f, dim = 2 )

        if(j .ne. 1 .and. j .ne. ny)then
            dfdy = (f(i, j+1) - f(i, j-1)) / (2.0 * dy)
            d2fdy2 = (f(i, j+1) - 2.0 * f(i,j) + f(i, j-1)) / dy**2
        end if

        if(j .eq. 1)then
            dfdy = (f(i, j+1) - f(i, j)) / dy
            d2fdy2 = (f(i, j+2) - 2.0 * f(i,j+1) + f(i, j)) / dy**2
        end if


        if(j .eq. ny)then
            dfdy = (f(i, j) - f(i, j-1)) / dy
            d2fdy2 = (f(i, j) - 2.0 * f(i,j-1) + f(i, j-2)) / dy**2
        end if

      end subroutine deriv_y

!
!***************************************************************************
!
      subroutine deriv_xy ( f, i, j, dx, dy, d2fdxy)

      implicit none

      integer, intent(in) :: i, j
      real, intent(in) :: f(:,:)
      real, intent(in) :: dx, dy
      real, intent(out), optional :: d2fdxy

      integer :: nx, ny

      nx    = size ( f, dim = 1 )
      ny    = size ( f, dim = 2 )

      d2fdxy = 0.0

      if (i .ne. 1 .and. i .ne. nx .and. &
          j .ne. 1 .and. j .ne. ny) then

         d2fdxy = (f(i+1,j+1) - f(i-1,j+1) - f(i+1,j-1) + f(i-1,j-1) ) &
            / (4.0 * dx * dy)

      endif

      end subroutine deriv_xy

!
!***************************************************************************
!


      subroutine sgrator(x, y, f, nx1, nx2, ny1, ny2, ans, &
         capr, r0, nxmax, nymax)

      implicit none

      integer i, j, nx1, nx2, ny1, ny2, nxmax, nymax
      real f(nxmax, nxmax), x(nxmax), y(nymax), ans, dx, dy
      real capr(nxmax), r0


      ans = 0.0

      do i = nx1, nx2 - 1
         dx = x(i+1) - x(i)

         do j = ny1, ny2 - 1
            dy = y(j+1) - y(j)
            ans = ans + dx * dy* (capr(i)   * f(i, j) &
                                + capr(i+1) * f(i+1, j) &
                                + capr(i)   * f(i, j+1) &
                                + capr(i+1) * f(i+1, j+1) ) / 4.0
         end do
      end do

      ans = ans / r0

      return
      end subroutine sgrator
!
!***************************************************************************
!

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

!
!***************************************************************************
!


      SUBROUTINE ZFUN(Z,FU)
!
!     ROUTINE WHICH EVALUATES THE PLASMA DISPERSION FUNCTION.  USES
!     NUMERICAL INTEGRATION (ABSOLUTE VALUE OF THE COMPLEX ARGUMENT, Z,
!     LESS THAN 5) OR ASYMPTOTIC EXPANSION.
!
      COMPLEX Z,FU,TEMP1,TEMP2,Z2,TPIIOD
      DIMENSION C(21),D(21),E(21),F(21),W(21)
      DATA DELTA/5.0E-1/, YI/-1.0E+0/, N/21/, ITEST/0/
      DATA ZERO/0.0E+0/, OSQPI/5.6418958355E-1/, TPI/6.28318530718E+0/
!
      IF(ITEST .EQ. 1) GO TO 1
      ITEST=1
!***     DEFINE WEIGHTS AND STORE CONSTANTS USED FOR INTEGRATION.
!***     WEIGHTS ARE DERIVED FROM A 3 POINT INTEGRATION SCHEME.
      NM3=N-3
      CONOI=DELTA*OSQPI
      W(1)=CONOI*3.75E-1
      W(2)=CONOI*1.1666666667E+0
      W(3)=CONOI*9.5833333333E-1
      DO 20 I=1,3
        NMIP1=N-I+1
   20 W(NMIP1)=W(I)
      DO 30 I=4,NM3
   30 W(I)=CONOI
      NO2=N/2
      NO2P1=NO2+1
      X=ZERO
      Y=YI
      E(NO2P1)=X
      F(NO2P1)=Y
      TEMP1=CMPLX(X,Y)
      TEMP2=CEXP(-TEMP1*TEMP1)
      C(NO2P1)=REAL(TEMP2)*W(NO2P1)
      D(NO2P1)=AIMAG(TEMP2)*W(NO2P1)
      DO 200 I=1,NO2
      X=DELTA*FLOAT(I)
      NPI=NO2P1+I
      NMI=NO2P1-I
      TEMP1=CMPLX(X,Y)
      TEMP2=CEXP(-TEMP1*TEMP1)
      C(NPI)=REAL(TEMP2)*W(NPI)
      C(NMI)=REAL(TEMP2)*W(NMI)
      D(NPI)=AIMAG(TEMP2)*W(NPI)
      D(NMI)=AIMAG(TEMP2)*W(NMI)*(-1.0E+0)
      E(NMI)=-X
      E(NPI)=X
      F(NPI)=Y
      F(NMI)=Y
  200 CONTINUE
      TPIOD=TPI/DELTA
      TPIIOD=CMPLX(ZERO,TPIOD)
!***     BEGIN CALCULATIONS.
    1 G=REAL(Z)
      YY=AIMAG(Z)
      H=ABS(YY)
      ZA=Z*CONJG(Z)
      IF(ZA .GE. 2.5E+1) GO TO 5
      Z2=Z*Z
!***     NUMERICAL INTEGRATION.
!***     F=1/SQRT(PI)*SUM OF...W(I)*EXP(-X(I)**2)/(X(I)-Z)...I=1,N.
!***     INTEGRATION IS ALONG A LINE X(I) IN THE COMPLEX PLANE, WHERE
!***     THE IMAGINARY PART OF X(I)=YI AND THE DIFFERENCE BETWEEN
!***     SUCCESSIVE REAL PARTS OF X(I)=DELTA.  LIMITS OF INTEGRATION
!***     ARE FROM -DELTA*N/2 TO DELTA*N/2.
!***     COMPUTE THE INTEGRAL BY TAKING THE SUM FROM 1 TO N OF THE
!***     CONSTANTS DIVIDED BY X(I)-Z.  USES REAL ARITHMETIC.
      ZR=0.0E+0
      ZI=0.0E+0
      DO 7 I=1,N
      A=E(I)-G
      B=F(I)-H
      DEN=A*A+B*B
      ODEN=1.0E+0/DEN
      ZR=ZR+(A*C(I)+B*D(I))*ODEN
      ZI=ZI+(A*D(I)-B*C(I))*ODEN
    7 CONTINUE
!***     ADD THE CORRECTION TERM.
      FU=CMPLX(ZR,ZI)+(0.0E+0,-3.5449077018E+0)*CEXP(-Z2-TPIOD* &
         (H-YI)+TPIIOD*G)
      IF(YY .GE. ZERO) GO TO 6
!***     IMAGINARY PART OF ARGUMENT IS NEGATIVE.
      FU=CONJG(FU)+(0.0E+0,3.5449077018E+0)*CEXP(-Z2)
      GO TO 6
!***     MAGNITUDE OF ARGUMENT IS GREATER THAN 5, USE
!***     ASYMPTOTIC EXPANSION.
    5 CALL AEXPAN(Z,G,H,YY,FU)
    6 RETURN
      END SUBROUTINE ZFUN





      SUBROUTINE AEXPAN(Z,G,H,YY,FU)
!
!     ROUTINE WHICH COMPUTES THE PLASMA DISPERSION FUNCTION USING
!     ASYMPTOTIC EXPANSION.  IF THE IMAGINARY PART OF THE ARGUMENT, YY,
!     IS EQUAL TO ZERO REAL ARITHMETIC IS USED.
!
      COMPLEX Z,FU,A,Z2,OZ2
      DATA ZERO/0.0E+0/, N/8/
!
      IF(YY .EQ. ZERO) GO TO 1
!***     COMPLEX ARITHMETIC.
      Z2=Z*Z
      OZ2=1.0E+0/Z2
      FU=-1.0E+0/Z
      A=FU
      EN=5.0E-1
      DO 10 I=1,N
      A=EN*A*OZ2
      FU=FU+A
   10 EN=EN+1.0E+0
      IF(YY .GT. ZERO) GO TO 30
      IF(H .GT. SQRT(G*G+1.72E+2)) GO TO 20
      FU=FU+(0.0E+0,3.5449077018E+0)*CEXP(-Z2)
      GO TO 30
!***     ERROR STOP TO AVOID OVERFLOW.
   20 WRITE (51,100) Z
  100 FORMAT(//1X,'*** ERROR STOP IN FRDCNT ROUTINE, ARGUMENT IS', &
       ' TOO SMALL,  ARG =',1P2E14.6)
      STOP
!***     REAL ARITHMETIC.
    1 X2=G*G
      OX2=1.0E+0/X2
      F=-1.0E+0/G
      B=F
      EN=5.0E-1
      DO 11 I=1,N
      B=EN*B*OX2
      F=F+B
   11 EN=EN+1.0E+0
      C=1.7724538509E+0*EXP(-X2)
      FU=CMPLX(F,C)
   30 RETURN
      END SUBROUTINE AEXPAN
!
!
! TWO POL APPROXIMATION
!
      SUBROUTINE Z2P(Z,F)
      COMPLEX Z,F
      F=CMPLX(.5,.81)/(CMPLX(.51,-.81)-Z)-CMPLX(.5,-.81) &
      /(CMPLX(.51,.81)+Z)
      RETURN
      END SUBROUTINE Z2P

!
! THIS CODE IS SUPERIOR TO PDF AND PDFLL.  NO ERROR AT ZETA=(4.,0.1)
!  disp(z) is the plasma z function
!
      complex function disp(z)
      complex z
      real im
      x=real(z)
      y=aimag(z)
      call zzdisp(x,y,re,im)
      disp=cmplx(re,im)
      return
      end function disp
!
!***************************************************************************
!
!
!  disp1(z) calculates first derivative of disp(z)
!
      complex function disp1(z)
      complex z,dp
      complex zz
      zz=z
      dp=disp(z)
      disp1=-2.0-2.0*zz*dp
      return
      end function disp1
!
!***************************************************************************
!
!
!  calculates plasma z function
!
      subroutine zzdisp(x,y,zzr,zzi)
      x1=abs(x)
      y1=abs(y)
      call wzdisp(x1,y1,wzr1,wzi1)
      if(x.ge.0..and.y.ge.0.) go to 1
      if(x.le.0..and.y.ge.0.) go to 2
      a=2*x1*y1
      b=-(x1*x1-y1*y1)
      abr=2*exp(b)*cos(a)
      abi=-2*exp(b)*sin(a)
      wzr1=abr-wzr1
      wzi1=abi-wzi1
      if(x.le.0..and.y.le.0.) go to 1
    2 wzi1=-wzi1
    1 zzr=-1.7724538509055*wzi1
      zzi=1.7724538509055*wzr1
      return
      end subroutine zzdisp
!
!***************************************************************************
!
!c   needed for z function
!c
      subroutine wzdisp(x,y,re,im)
      equivalence(epsh,epsl,epsy)
      real im,lambda
      integer capn
      logical b
      epsh=1.e-12
      if(y.lt.4.29 .and. x.lt.5.33) go to 10
!
!  (x,y) belongs to q1-r
!
      h=0.
      capn=0
      nu=8
      go to 20
!
!  (x,y) belongs to r
!
   10 s=(1.-y/4.29)*sqrt(1.-x*x/28.41)
      h=1.6*s
      h2=2.*h
      capn=6.+23*s+.5
      nu=9.+21*s+.5
      lambda=h2**capn
   20 b=(h.eq.0. .or. lambda.lt.epsl)
!
!  statement (lambda.lt.epsl) covers the underflow case
!  when h(.gt.0) is very small.
!
      rr=0.
      ri=0.
      sr=0.
      si=0.
      nup=nu+1
      do 100 i=1,nup
      n=nup-i
      np1=n+1
      tr=y+h+np1*rr
      ti=x-np1*ri
      c=.5/(tr*tr+ti*ti)
      rr=c*tr
      ri=c*ti
      if(.not.(h.gt.0..and.n.le.capn)) go to 100
      tr=lambda+sr
      sr=rr*tr-ri*si
      si=ri*tr+rr*si
      lambda=lambda/h2
  100 continue
      cc=1.12837916709551
      if(y.lt.epsy) go to 120
      if(b) go to 110
      re=sr*cc
      go to 130
  110 re=rr*cc
      go to 130
  120 re=exp(-x*x)
  130 if(b) go to 140
      im=si*cc
      go to 150
  140 im=ri*cc
  150 return
      end subroutine wzdisp


!
!***************************************************************************
!
!      real function second1_old(dummy)
!
!      implicit none
!
!      integer :: v(8)
!!      integer mtime, mclock
!      real dummy
!
!!*****Fortran 90 standard for wall clock time
!!      call date_and_time(values=v)
!!      second1_old=(v(5)*3600)+(v(6)*60)+v(7)+0.001d0*v(8)
!
!!*****FALCON:
!      double precision mpi_wtime
!      external mpi_wtime
!      second1_old = MPI_WTIME()
!
!
!!*****EAGLE:
!!      mtime = mclock()
!!      second1_old  = 0.01 * mtime
!
!      return
!      end function second1_old

!
!***************************************************************************
!

      real function second1(dummy)
!     returns a monotonically increasing wall-clock time-stamp in seconds
!     with an attempt to cross multiple days between calls
!     (assumes one day if end-of-the-month is spanned between calls)

      implicit none

      real dummy

      integer :: v(8)

      integer,save :: ld
      logical,save :: elapsed=.false.

!*****Fortran 90 standard for wall clock time
      call date_and_time(values=v)

      if(elapsed)then  ! second has been called before
        if(v(3).ge.ld)then ! days have advanced monotonically
         second1=(v(3)-ld)*86400+(v(5)*3600)+(v(6)*60)+v(7)+0.001d0*v(8)
        else  ! day counter reset past end of the month
          print *,'WARNING: SECOND1 TRIPPED OVER THE END OF THE MONTH'
          print *,'         SECOND1 IS ASSUMING ONLY ONE ELAPSED DAY'
!         proper fix requires full blown calendar program
          ld=0
         second1=(v(3)-ld)*86400+(v(5)*3600)+(v(6)*60)+v(7)+0.001d0*v(8)
        endif
      else  ! first call
        second1=(v(5)*3600)+(v(6)*60)+v(7)+0.001d0*v(8)
        ld=v(3)
        elapsed=.true.
      endif

      return
      end function second1

subroutine midplane(igiven, jgiven, f, fmid, rho,&
    nxdim, nydim, nnodex, nnodey, capr, r0, f0, jmid)

    implicit none

    integer nxdim, nydim, nnodex, nnodey, i, j, jmid, i0
    integer igiven, jgiven, istart

    real f(nxdim, nydim), fmid, rho(nxdim, nydim), r0, capr(nxdim)
    real rhoij, f0

    rhoij = rho(igiven, jgiven)
    fmid =  0.0

    if (rhoij .gt. 1.0) return

    do i = 1, nnodex

        if(capr(i) .ge. r0)then
            istart = i
            go to 200
        endif

    enddo

200 continue

    if(rhoij .ge. rho(istart, jmid)) then

         do i = istart, nnodex - 1
            if(rhoij .ge. rho(i,   jmid) .and. &
              rhoij .lt. rho(i+1, jmid)) i0 = i
         enddo

         fmid = f(i0, jmid) + (f(i0 + 1, jmid) -   f(i0, jmid)) &
                         / (rho(i0 + 1, jmid) - rho(i0, jmid)) &
                         * (rhoij             - rho(i0, jmid))
    endif


    if(rhoij .lt. rho(istart, jmid)) then

         fmid = f0 &
            + (f(istart, jmid) - f0) / (rho(istart, jmid) - 0.0) &
                  * (rhoij - 0.0)
    endif

100 format (1i10, 1p8e12.4)
102 format (2i10)

end subroutine midplane

subroutine flux_to_rz(nnodex, nnodey, profile_in, &
         profile_out, rho, nrho, rho_ij)

      implicit none

      real :: profile_in(nrho), profile_out(nnodex,nnodey)
      real :: rho(nrho), rho_ij(nnodex,nnodey)
      integer :: nnodex, nnodey, nrho, i, j, k

      do i = 1, nnodex
         do j = 1, nnodey
          do k = 1, nrho

             if(rho_ij(i,j) .le. rho(1)) then
             !  need one sided to the left

!              call flush(6)
              profile_out(i,j) = profile_in(1) + &
                  (profile_in(2) - profile_in(1)) / &
                  (rho(2) - rho(1)) * (rho_ij(i,j) - rho(1))

             elseif(rho_ij(i,j) .ge.rho(nrho)) then
             ! need one sidded to the right.

                profile_out(i,j) = profile_in(nrho) + &
                  (profile_in(nrho) - profile_in(nrho - 1)) / &
                  (rho(nrho) - rho(nrho-1)) * (rho_ij(i,j) - rho(nrho))

               elseif(rho_ij(i,j) - rho(k) .ge. 0.0 .and. &
                   rho_ij(i,j) -rho(k+1) .lt. 0.0 ) then
                ! standard  caseprint*, 'last if'
              profile_out(i,j) = profile_in(k)  + &
                  (profile_in(k+1) - profile_in(k))/ &
                  (rho(k+1) - rho(k)) * (rho_ij(i,j) - rho(k))
               endif

             if (profile_out(i,j) .lt. 0.0) profile_out(i,j) = 1.0e-10

          enddo
       end do
      end do

      end subroutine flux_to_rz



end module aorsasubs_mod

