module sigmaInputGeneration


! Interpolation splines

real, allocatable :: zp_omgC(:), zp_omgP2(:)

use sigma, only: spatialSigmaInput_cold

contains

subroutine setupSigmaParameterInterps ( g )

    use grid, only: gridBlock

    implicit none

    type(gridBlock), intent(inout) :: g

    integer :: islpsw, iErr
    real, allocatable :: interpTemp(:)
    real :: zxy11, zxym1, zxy1n, zxymn
    real, allocatable :: zx1(:), zxm(:), zy1(:), zyn(:)
    real :: interpSigma = 0.0
    real, allocatable :: tmpReal(:,:)
 
    ! Setup the splines for all 2D parameters reqd for sigma input 
    ! ------------------------------------------------------------

    allocate( zx1(g%nR), zxm(g%nR), zy1(g%nZ), zyn(g%nZ) ) 
    islpsw  = 255 
    allocate( &
        interpTemp(g%nZ+g%nZ+g%nR), &
        zp_omgc(3*g%nR*g%nZ))

    allocate( tmpReal(g%nR,g%nZ) )
    tmpReal = g%omgc
    call surf1 ( g%nR, g%nZ, g%r, g%z, tmpReal, g%nR, zx1, zxm, &
        zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
        zp_omgc, interpTemp, interpSigma, iErr)
 
end subroutine setupSigmaParameterInterps


function getSigmaInputHere ()

    implicit none

    type(spatialSigmaInput_cold) :: getSigmaInputHere

    psi_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
        psizr, nw, zp_psi, sigma )
 

end function getSigmaInputHere


end module sigmaInputGeneration
