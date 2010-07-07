module spline_dlg

contains


! See Burden and Faires pg 133.
!
! NOTE: allocate the "spline" variable prior to 
!   function call. Does not work otherwise. Not
!   sure as to the cause, perhaps compiler bug?


function spline ( x, a, fp0, fpn )

    implicit none

    real, intent(in) :: x(:), a(:), fp0, fpn
    real, allocatable :: spline(:,:)

    integer :: n, i
    real, allocatable :: b(:), c(:), d(:), h(:), alpha(:)
    real, allocatable :: u(:), l(:), z(:)

    n = size ( x ) - 1

    allocate ( b(n+1),c(n+1),d(n+1),h(n), alpha(n+1) )
    if(.not.allocated(spline)) allocate ( spline(n+1,4) )
    allocate ( u(n+1), l(n+1), z(n+1) )

    do i=1,n
       h(i) = x(i+1)-x(i)
    enddo 

    alpha(1) = 3 * ( a(2)-a(1) ) / h(1) - 3 * fp0
    alpha(n+1) = 3 * fpn - 3 * ( a(n+1)-a(n) ) / h(n)

    do i=2,n
        alpha(i) = &
            3 * (a(i+1)*h(i-1) - a(i) * ( x(i+1)-x(i-1) ) + a(i-1)*h(i) ) &
                / (h(i-1)*h(i))
    enddo

    l(1) = 2 * h(1)
    u(1) = 0.5
    z(1) = alpha(1) / l(1)

    do i=2,n
        l(i) = 2 * ( x(i+1)-x(i-1) ) - h(i-1)*u(i-1)
        u(i) = h(i) / l(i)
        z(i) = ( alpha(i) - h(i-1)*z(i-1) ) / l(i)
    enddo

    l(n+1) = h(n) * ( 2-u(n) )
    z(n+1) = ( alpha(n+1)-h(n)*z(n) ) / l(n+1)

    c(n+1) = z(n+1)

    do i=n,1,-1

        c(i) = z(i) - u(i)*c(i+1)
        b(i) = (a(i+1)-a(i))/h(i) - h(i)*(c(i+1)+2*c(i))/3
        d(i) = (c(i+1)-c(i))/(3*h(i))

    enddo

    spline(:,1) = a
    spline(:,2) = b
    spline(:,3) = c
    spline(:,4) = d

end function spline

end module spline_dlg
