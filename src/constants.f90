module constants

    integer, parameter :: dbl = selected_real_kind ( p=13 )
    real(kind=dbl), parameter :: eOverAmu = 9.58084e7
    real(kind=dbl), parameter :: mpc2 = 938271998.38
    real(kind=dbl), parameter :: c = 2.99792458e8
    real(kind=dbl), parameter :: pi = 3.141592653597932384
    real(kind=dbl), parameter :: q = 1.60217646e-19 
    real(kind=dbl), parameter :: xme = 9.10938188e-31
    real(kind=dbl), parameter :: xmh = 1.67262158e-27
    real(kind=dbl), parameter :: eps0 = 8.85e-12
    complex, parameter :: zi = cmplx ( 0.0, 1.0 )
    real(kind=dbl), parameter :: xmu0 = 1.2566370614e-06
    real(kind=dbl), parameter ::  clight = 299792458.0 


end module constants
