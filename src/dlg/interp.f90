module interp
    implicit none
    save

    !   interpolation initialisation arrays

    real, allocatable :: zp_bCurv_R(:), zp_bCurv_phi(:), zp_bCurv_z(:), &
        zp_bGrad_R(:), zp_bGrad_phi(:), zp_bGrad_z(:), &
        zp_bR(:), zp_bPhi(:), zp_bz(:), zp_bDotGradB(:), &
        zp_psi(:)
    real :: sigma = 0.0

contains
    subroutine init_interp ()
        !use gc_terms
        use eqdsk_dlg
        use fitpack 

        implicit none
 
        integer :: islpsw, iErr
        real, allocatable :: temp(:)
        real :: zx1(nh), zxm(nh), zy1(nw), zyn(nw), &
            zxy11, zxym1, zxy1n, zxymn
        islpsw  = 255 
        allocate ( zp_bCurv_R(3*nw*nh), &
            zp_bCurv_phi(3*nw*nh), &
            zp_bCurv_z(3*nw*nh),&
            temp(nh+nh+nw), &
            zp_bGrad_R(3*nw*nh), &
            zp_bGrad_phi(3*nw*nh), &
            zp_bGrad_z(3*nw*nh), &
            zp_bDotGradB(3*nw*nh), &
            zp_bR(3*nw*nh), zp_bPhi(3*nw*nh), zp_bz(3*nw*nh), &
            zp_psi(3*nw*nh) )

        !   psi 

        call surf1 ( nw, nh, r, z, psizr, nw, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            zp_psi, temp, sigma, iErr)
     
        !!   b field

        !call surf1 ( nw, nh, r, z, bR, nw, zx1, zxm, &
        !    zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
        !    zp_bR, temp, sigma, iErr)
        !call surf1 ( nw, nh, r, z, bPhi, nw, zx1, zxm, &
        !    zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
        !    zp_bPhi, temp, sigma, iErr)
        !call surf1 ( nw, nh, r, z, bz__, nw, zx1, zxm, &
        !    zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
        !    zp_bz, temp, sigma, iErr)
 
        !!   curvature 
    
        !call surf1 ( nw, nh, r, z, bCurvature_R, nw, zx1, zxm, &
        !    zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
        !    zp_bCurv_R, temp, sigma, iErr)
        !call surf1 ( nw, nh, r, z, bCurvature_phi, nw, zx1, zxm, &
        !    zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
        !    zp_bCurv_phi, temp, sigma, iErr)
        !call surf1 ( nw, nh, r, z, bCurvature_z, nw, zx1, zxm, &
        !    zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
        !    zp_bCurv_z, temp, sigma, iErr)

        !!   gradient

        !call surf1 ( nw, nh, r, z, bGradient_R, nw, zx1, zxm, &
        !    zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
        !    zp_bGrad_R, temp, sigma, iErr)
        !call surf1 ( nw, nh, r, z, bGradient_phi, nw, zx1, zxm, &
        !    zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
        !    zp_bGrad_phi, temp, sigma, iErr)
        !call surf1 ( nw, nh, r, z, bGradient_z, nw, zx1, zxm, &
        !    zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
        !    zp_bGrad_z, temp, sigma, iErr)

        !!   bDotGraB 

        !call surf1 ( nw, nh, r, z, bDotGradB, nw, zx1, zxm, &
        !    zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
        !    zp_bDotGradB, temp, sigma, iErr)
 

    end subroutine init_interp

end module interp
