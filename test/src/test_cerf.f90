program test_cerf

use constants
use hammett

implicit none

complex(kind=dbl), allocatable :: x(:)
real :: xMin, xMax
integer :: nPts, i
real(kind=dbl) :: u, v, u2, v2
logical :: flag
integer :: cerf_wrap, stat

nPts = 11
allocate ( x(nPts) )

xMin = -2
xMax =  2

x = (/ (i*(xMax-xMin)/(nPts-1)+xMin,i=0,nPts-1) /) 

do i=1,nPts
      call WOFZ (real(x(i)), aimag(x(i)), u, v, flag)
      if(flag)then
              write(*,*) 'Overflow WILL occur'
              stop
      endif
      stat =  cerf_wrap ( real(x(i)), aimag(x(i)), u2, v2 ) 
      write(*,'(6(1x,f6.2))'), &
        real(x(i)), aimag(x(i)), u, v, u2, v2
enddo

stop

end program test_cerf
