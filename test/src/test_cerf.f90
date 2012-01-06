program test_cerf

use constants
use hammett
use z_erfc

implicit none

complex(kind=dbl), allocatable :: x(:), y(:)
real :: xMin, xMax
integer :: nPts, i
real(kind=dbl) :: u, v, u2, v2, u3, v3
logical :: flag
integer :: cerf_wrap, stat

nPts = 11
allocate ( x(nPts), y(nPts) )

xMin = -2
xMax =  2

x = (/ (i*(xMax-xMin)/(nPts-1)+xMin,i=0,nPts-1) /) 

do i=1,nPts
      call WOFZ (real(x(i)), aimag(x(i)), u, v, flag)
      if(flag)then
              write(*,*) 'Overflow WILL occur in wofz'
              stop
      endif
      call wofz_f90 (x(i), y(i), flag)
      if(flag)then
              write(*,*) 'Overflow WILL occur in wofz_f90'
              stop
      endif
      u3=realpart(y(i))
      v3=imagpart(y(i))
      stat =  cerf_wrap ( real(x(i)), aimag(x(i)), u2, v2 ) 
      write(*,'(8(1x,f6.2))'), &
        realpart(x(i)), imagpart(x(i)), u, v, u2, v2, u3, v3
enddo

stop

end program test_cerf
