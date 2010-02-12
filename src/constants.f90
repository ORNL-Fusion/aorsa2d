module constants

    real, parameter :: eOverAmu = 9.58084e7
    real, parameter :: mpc2 = 938271998.38
    real, parameter :: c = 2.99792458e8
    real, parameter :: pi = 3.141592653597932384
    real, parameter :: q = 1.60217646e-19 
    real, parameter :: xme = 9.10938188e-31
    real, parameter :: xmh = 1.67262158e-27
    real, parameter :: eps0 = 8.85e-12
    complex, parameter :: zi = cmplx ( 0.0, 1.0 )
    real, parameter :: xmu0 = 1.2566370614e-06
    real, parameter ::  clight = 1.0 / sqrt ( eps0 * xmu0 )
    integer, parameter :: dbl = selected_real_kind ( p=13 )


end module constants
