module constants

    integer, parameter :: dbl = selected_real_kind ( p=13 )
    integer, parameter :: long = selected_real_kind ( p=9 )
    real(kind=dbl), parameter :: eOverAmu = 9.58084d7
    real(kind=dbl), parameter :: mpc2 = 938271998.38
    real(kind=dbl), parameter :: c = 2.99792458d8
    real(kind=dbl), parameter :: pi = 3.141592653597932384
    real(kind=dbl), parameter :: q = 1.60217646d-19 
    real(kind=dbl), parameter :: xme = 9.10938188d-31
    real(kind=dbl), parameter :: xmh = 1.67262158d-27
    real(kind=dbl), parameter :: eps0 = 8.85d-12
    complex, parameter :: zi = cmplx ( 0.0, 1.0 )
    real(kind=dbl), parameter :: mu0 = 1.2566370614d-06
    real(kind=dbl), parameter ::  clight = 299792458.0 


end module constants
