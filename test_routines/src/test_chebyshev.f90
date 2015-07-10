program test_chebyshev

use chebyshev_mod

implicit none

real, allocatable :: x(:), T(:), U(:)
integer :: nPts, n, i
real :: dx

nPts = 11

allocate ( x(nPts), T(nPts), U(nPts) )

n=25

dx  = 2.0 / (nPts-1)
x   = (/ ((i-1)*dx-1,i=1,nPts) /)

write(*,*) 'Chebyshev T and U results, compare with mathematica worksheet'
do i=1,nPts

    T(i) = chebT(n,x(i))
    U(i) = chebU(n,x(i))

    write(*,*) n, x(i), T(i), U(i)

enddo


end program test_chebyshev
