module sigmaInputGeneration

use sigma, only: spatialSigmaInput_cold

contains

subroutine setupSigmaParameterSplines ( g )

    use grid, only: gridBlock
    use fitpack
    use aorsaNamelist, only: nSpec

    implicit none

    type(gridBlock), intent(inout) :: g

    integer :: islpsw, iErr, s
    real, allocatable :: interpTemp(:)
    real :: zxy11, zxym1, zxy1n, zxymn
    real, allocatable :: zx1(:), zxm(:), zy1(:), zyn(:)
    real, allocatable :: tmpReal(:,:)
 
    ! Setup the splines for all 2D parameters reqd for sigma input 
    ! ------------------------------------------------------------

    allocate( zx1(g%nR), zxm(g%nR), zy1(g%nZ), zyn(g%nZ) ) 
    islpsw  = 255 
    allocate( interpTemp(g%nZ+g%nZ+g%nR), tmpReal(g%nR,g%nZ) )

    allocate( g%spline_omgC(3*g%nR*g%nZ,nSpec) )
    do s=1,nSpec

        tmpReal = g%omgc(:,:,s)
        call surf1 ( g%nR, g%nZ, g%r, g%z, tmpReal, g%nR, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            g%spline_omgC(:,s), interpTemp, g%interpSigma, iErr)

    enddo

    allocate( g%spline_omgP2(3*g%nR*g%nZ,nSpec) )
    do s=1,nSpec

        tmpReal = g%omgP2(:,:,s)
        call surf1 ( g%nR, g%nZ, g%r, g%z, tmpReal, g%nR, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            g%spline_omgP2(:,s), interpTemp, g%interpSigma, iErr)

    enddo

    allocate( g%spline_nuOmg(3*g%nR*g%nZ) )
    call surf1 ( g%nR, g%nZ, g%r, g%z, g%nuOmg, g%nR, zx1, zxm, &
        zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
        g%spline_nuOmg, interpTemp, g%interpSigma, iErr)


end subroutine setupSigmaParameterSplines


function getSigmaInputHere ( r, z, g, speciesNo )

    use grid, only: gridBlock
    use fitpack
    use aorsaNamelist, only: nSpec
    use profiles, only: omgRF

    implicit none

    type(gridBlock), intent(in) :: g
    integer, intent(in) :: speciesNo
    real, intent(in) :: r, z
    type(spatialSigmaInput_cold) :: getSigmaInputHere
    real, allocatable :: tmpReal(:,:)
    real :: omgC_here, omgP2_here, nuOmg_here

    allocate( tmpReal(g%nR,g%nZ) )

    tmpReal = g%omgC(:,:,speciesNo)
    omgC_here = surf2 ( r, z, g%nR, g%nZ, g%r, g%z, &
        tmpReal, g%nR, g%spline_omgC(:,speciesNo), g%interpSigma )

    tmpReal = g%omgP2(:,:,speciesNo)
    omgP2_here = surf2 ( r, z, g%nR, g%nZ, g%r, g%z, &
        tmpReal, g%nR, g%spline_omgP2(:,speciesNo), g%interpSigma )

    nuOmg_here = surf2 ( r, z, g%nR, g%nZ, g%r, g%z, &
        g%nuOmg, g%nR, g%spline_nuOmg, g%interpSigma )

    getSigmaInputHere = spatialSigmaInput_cold( &
        omgC_here,&
        omgP2_here,&
        omgRF,&
        nuOmg_here)

end function getSigmaInputHere

end module sigmaInputGeneration
