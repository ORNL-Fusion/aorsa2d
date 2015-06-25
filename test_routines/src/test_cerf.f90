program test_cerf

use constants
use z_erf, only: wofz_f90

implicit none

complex(kind=dbl), allocatable :: x(:), y1(:),y2(:)
real :: xMin, xMax
integer :: nPts, i
real(kind=dbl) :: u, v, u2, v2, u3, v3
logical :: flag
integer :: cerf_wrap, stat
complex :: zIn, zOut

nPts = 21
allocate ( x(nPts), y1(nPts),y2(nPts) )

xMin = -900
xMax =  900

x = (/ (i*(xMax-xMin)/(nPts-1)+xMin,i=0,nPts-1) /) 

do i=1,nPts
      call WOFZ (real(x(i)), aimag(x(i)), u, v, flag)
      if(flag)then
              write(*,*) 'Overflow WILL occur in wofz'
              stop
      endif
      call wofz_f90 (x(i), y1(i), flag)
      if(flag)then
              write(*,*) 'Overflow WILL occur in wofz_f90'
              stop
      endif
      y1(i) = zi*sqrt(pi)*y1(i)
      y2(i) = zi*sqrt(pi)*complex(u,v)
      u3=realpart(y1(i))
      v3=imagpart(y1(i))
      u=realpart(y2(i))
      v=imagpart(y2(i))
      !stat =  cerf_wrap ( real(x(i)), aimag(x(i)), u2, v2 ) 
      write(*,'(6(1x,f8.4))'), &
        realpart(x(i)), imagpart(x(i)), u, v, u3, v3
enddo

stop

end program test_cerf
