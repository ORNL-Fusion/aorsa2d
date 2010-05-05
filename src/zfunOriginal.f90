module zFunOriginal

contains

      subroutine zfun(z,fu)
!
!     routine which evaluates the plasma dispersion function.  uses
!     numerical integration (absolute value of the complex argument, z,
!     less than 5) or asymptotic expansion.
!
      complex z,fu,temp1,temp2,z2,tpiiod
      dimension c(21),d(21),e(21),f(21),w(21)
      data delta/5.0e-1/, yi/-1.0e+0/, n/21/, itest/0/
      data zero/0.0e+0/, osqpi/5.6418958355e-1/, tpi/6.28318530718e+0/
!
      if(itest .eq. 1) go to 1
      itest=1
!***     define weights and store constants used for integration.
!***     weights are derived from a 3 point integration scheme.
      nm3=n-3
      conoi=delta*osqpi
      w(1)=conoi*3.75e-1
      w(2)=conoi*1.1666666667e+0
      w(3)=conoi*9.5833333333e-1
      do 20 i=1,3
        nmip1=n-i+1
   20 w(nmip1)=w(i)
      do 30 i=4,nm3
   30 w(i)=conoi
      no2=n/2
      no2p1=no2+1
      x=zero
      y=yi
      e(no2p1)=x
      f(no2p1)=y
      temp1=cmplx(x,y)
      temp2=cexp(-temp1*temp1)
      c(no2p1)=real(temp2)*w(no2p1)
      d(no2p1)=aimag(temp2)*w(no2p1)
      do 200 i=1,no2
      x=delta*float(i)
      npi=no2p1+i
      nmi=no2p1-i
      temp1=cmplx(x,y)
      temp2=cexp(-temp1*temp1)
      c(npi)=real(temp2)*w(npi)
      c(nmi)=real(temp2)*w(nmi)
      d(npi)=aimag(temp2)*w(npi)
      d(nmi)=aimag(temp2)*w(nmi)*(-1.0e+0)
      e(nmi)=-x
      e(npi)=x
      f(npi)=y
      f(nmi)=y
  200 continue
      tpiod=tpi/delta
      tpiiod=cmplx(zero,tpiod)
!***     begin calculations.
    1 g=real(z)
      yy=aimag(z)
      h=abs(yy)
      za=z*conjg(z)
      if(za .ge. 2.5e+1) go to 5
      z2=z*z
!***     numerical integration.
!***     f=1/sqrt(pi)*sum of...w(i)*exp(-x(i)**2)/(x(i)-z)...i=1,n.
!***     integration is along a line x(i) in the complex plane, where
!***     the imaginary part of x(i)=yi and the difference between
!***     successive real parts of x(i)=delta.  limits of integration
!***     are from -delta*n/2 to delta*n/2.
!***     compute the integral by taking the sum from 1 to n of the
!***     constants divided by x(i)-z.  uses real arithmetic.
      zr=0.0e+0
      zi=0.0e+0
      do 7 i=1,n
      a=e(i)-g
      b=f(i)-h
      den=a*a+b*b
      oden=1.0e+0/den
      zr=zr+(a*c(i)+b*d(i))*oden
      zi=zi+(a*d(i)-b*c(i))*oden
    7 continue
!***     add the correction term.
      fu=cmplx(zr,zi)+(0.0e+0,-3.5449077018e+0)*cexp(-z2-tpiod* &
         (h-yi)+tpiiod*g)
      if(yy .ge. zero) go to 6
!***     imaginary part of argument is negative.
      fu=conjg(fu)+(0.0e+0,3.5449077018e+0)*cexp(-z2)
      go to 6
!***     magnitude of argument is greater than 5, use
!***     asymptotic expansion.
    5 call aexpan(z,g,h,yy,fu)
    6 return
      end subroutine zfun


      subroutine aexpan(z,g,h,yy,fu)
!
!     routine which computes the plasma dispersion function using
!     asymptotic expansion.  if the imaginary part of the argument, yy,
!     is equal to zero real arithmetic is used.
!
      complex z,fu,a,z2,oz2
      data zero/0.0e+0/, n/8/
!
      if(yy .eq. zero) go to 1
!***     complex arithmetic.
      z2=z*z
      oz2=1.0e+0/z2
      fu=-1.0e+0/z
      a=fu
      en=5.0e-1
      do 10 i=1,n
      a=en*a*oz2
      fu=fu+a
   10 en=en+1.0e+0
      if(yy .gt. zero) go to 30
      if(h .gt. sqrt(g*g+1.72e+2)) go to 20
      fu=fu+(0.0e+0,3.5449077018e+0)*cexp(-z2)
      go to 30
!***     error stop to avoid overflow.
   20 write (51,100) z
  100 format(//1x,'*** error stop in frdcnt routine, argument is', &
       ' too small,  arg =',1p2e14.6)
      stop
!***     real arithmetic.
    1 x2=g*g
      ox2=1.0e+0/x2
      f=-1.0e+0/g
      b=f
      en=5.0e-1
      do 11 i=1,n
      b=en*b*ox2
      f=f+b
   11 en=en+1.0e+0
      c=1.7724538509e+0*exp(-x2)
      fu=cmplx(f,c)
   30 return
      end subroutine aexpan
!
!

end module zfunOriginal
