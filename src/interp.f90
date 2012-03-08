module interp

    use eqdsk_dlg
    use fitpack 

    implicit none

    !   interpolation initialisation arrays

    real, allocatable :: &
        zp_bR(:), zp_bPhi(:), zp_bz(:), zp_psi(:), zp_rho(:)
    real :: sigma = 0.0

    
contains

    subroutine init_interp ()

        implicit none
 
        integer :: islpsw, iErr
        real, allocatable :: temp(:)
        real :: zx1(nh), zxm(nh), zy1(nw), zyn(nw), &
            zxy11, zxym1, zxy1n, zxymn
        islpsw  = 255 
        allocate ( &
            temp(nh+nh+nw), &
            zp_bR(3*nw*nh), zp_bPhi(3*nw*nh), zp_bz(3*nw*nh), &
            zp_psi(3*nw*nh),zp_rho(3*nw*nh) )

        !   psi 

        call surf1 ( nw, nh, r, z, psizr, nw, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            zp_psi, temp, sigma, iErr)
 
        !   rho 

        call surf1 ( nw, nh, r, z, rhoNorm, nw, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            zp_rho, temp, sigma, iErr)
     
        !   b field

        call surf1 ( nw, nh, r, z, bR, nw, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            zp_bR, temp, sigma, iErr)
        call surf1 ( nw, nh, r, z, bPhi, nw, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            zp_bPhi, temp, sigma, iErr)
        call surf1 ( nw, nh, r, z, bz__, nw, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            zp_bz, temp, sigma, iErr)
 
    end subroutine init_interp

    function dlg_interpB ( pos, bMagHere, psiHere, rhoHere )

        implicit none
        
        real :: bR_here, bPhi_here, bz_here, psi_here, rho_here
        real, intent(IN) :: pos(3)
        real :: dlg_interpB(3)
        real, optional, intent(OUT) :: bMagHere, psiHere, rhoHere

        bR_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bR, nw, zp_bR, sigma )
        bPhi_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bPhi, nw, zp_bPhi, sigma )
        bz_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bz__, nw, zp_bz, sigma )
        psi_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            psizr, nw, zp_psi, sigma )
        rho_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            rhoNorm, nw, zp_rho, sigma )


        if ( present (bMagHere) ) &
            bMagHere    = sqrt ( bR_here**2 + bPhi_here**2 + bz_here**2 )

        if ( present (psiHere) ) &
            psiHere    = psi_here 

        if ( present (rhoHere) ) &
            rhoHere    = rho_here 

        dlg_interpB(1)  = bR_here 
        dlg_interpB(2)  = bPhi_here
        dlg_interpB(3)  = bz_here 
    
    end function dlg_interpB

end module interp
