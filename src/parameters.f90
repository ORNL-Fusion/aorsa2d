module parameters

    implicit none

    integer, parameter :: n_theta_max = 200
    integer, parameter :: n_psi_max = 200
    integer, parameter :: nmodesmax = 288
    integer, parameter :: mmodesmax = 288
    integer, parameter :: nxmx = nmodesmax
    integer, parameter :: nymx = mmodesmax
    integer, parameter :: nrhomax = nmodesmax * 2
    integer, parameter :: nthetamax = nmodesmax * 2
    integer, parameter :: nkdim1 = - nmodesmax / 2
    integer, parameter :: nkdim2 =   nmodesmax / 2
    integer, parameter :: mkdim1 = - mmodesmax / 2
    integer, parameter :: mkdim2 =   mmodesmax / 2
    integer, parameter :: nldim  = nxmx * nymx
    integer, parameter :: nldim3 = 3 * nldim
    integer, parameter :: ndfmax = nldim3
    integer, parameter :: lmaxdim = 25
    integer, parameter :: ndim = 256

end module parameters
