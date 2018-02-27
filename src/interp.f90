module interp

    use eqdsk_dlg
    use fitpack 

    implicit none

    !   interpolation initialisation arrays

    real, allocatable :: &
        zp_bR(:), zp_bPhi(:), zp_bz(:), zp_psi(:), zp_rho(:)
    real :: SplineSigma = 0.0
    real, allocatable :: yp_br(:), yp_bt(:), yp_bz(:)
    
contains

    subroutine init_interp ()

        use aorsaNameList, only: useEqdsk, useAr2Input
        use AR2Input, &
            only:nR_ar2=>nR,nZ_ar2=>nZ, &
            br_ar2=>bR,bt_ar2=>bt,bz_ar2=>bz, &
            r_ar2=>r,z_ar2=>z

        implicit none
 
        integer :: islpsw, iErr
        real, allocatable :: temp(:)
        real, allocatable :: zx1(:), zxm(:), zy1(:), zyn(:)
        real, allocatable :: slp1, slpn
        real :: zxy11, zxym1, zxy1n, zxymn

        islpsw  = 255 

        ar2Input_: &
        if(useAr2Input)then

            if(nZ_ar2.gt.1)then

                allocate ( &
                    temp(nZ_ar2+nZ_ar2+nR_ar2), &
                    zp_bR(3*nR_ar2*nZ_ar2), zp_bPhi(3*nR_ar2*nZ_ar2), zp_bz(3*nR_ar2*nZ_ar2), &
                    zp_psi(3*nR_ar2*nZ_ar2),zp_rho(3*nR_ar2*nZ_ar2) )

                allocate (zx1(nZ_ar2), zxm(nZ_ar2), zy1(nR_ar2), zyn(nR_ar2))

                call surf1 ( nR_ar2, nZ_ar2, r_ar2, z_ar2, br_ar2, nR_ar2, zx1, zxm, &
                    zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
                    zp_bR, temp, SplineSigma, iErr)
                call surf1 ( nR_ar2, nZ_ar2, r_ar2, z_ar2, bt_ar2, nR_ar2, zx1, zxm, &
                    zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
                    zp_bPhi, temp, SplineSigma, iErr)
                call surf1 ( nR_ar2, nZ_ar2, r_ar2, z_ar2, bz_ar2, nR_ar2, zx1, zxm, &
                    zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
                    zp_bz, temp, SplineSigma, iErr)

            else

                islpsw = 3

                allocate(yp_br(nR_ar2),yp_bt(nR_ar2),yp_bz(nR_ar2))
                allocate(temp(nR_ar2))

                call curv1 (nR_ar2,r_ar2,br_ar2,slp1,slpn,islpsw,yp_br,temp,SplineSigma,iErr)
                call curv1 (nR_ar2,r_ar2,bt_ar2,slp1,slpn,islpsw,yp_bt,temp,SplineSigma,iErr)
                call curv1 (nR_ar2,r_ar2,bz_ar2,slp1,slpn,islpsw,yp_bz,temp,SplineSigma,iErr)

            endif

        endif ar2Input_
 
    end subroutine init_interp

    function dlg_interpB ( pos, bMagHere, psiHere, rhoHere )

        use aorsaNameList, only: useEqdsk, useAr2Input
        use AR2Input, &
            only:nR_ar2=>nR,nZ_ar2=>nZ, &
            br_ar2=>bR,bt_ar2=>bt,bz_ar2=>bz,r_ar2=>r,z_ar2=>z

        implicit none
        
        real :: bR_here, bPhi_here, bz_here, psi_here, rho_here
        real, intent(IN) :: pos(3)
        real :: dlg_interpB(3)
        real, optional, intent(OUT) :: bMagHere, psiHere, rhoHere

        if(useAR2Input)then
            if(nZ_ar2.gt.1)then

                bR_here = surf2 ( pos(1), pos(3), nR_ar2, nZ_ar2, r_ar2, z_ar2, &
                    br_ar2, nR_ar2, zp_bR, SplineSigma )
                bPhi_here = surf2 ( pos(1), pos(3), nR_ar2, nZ_ar2, r_ar2, z_ar2, &
                    bt_ar2, nR_ar2, zp_bPhi, SplineSigma )
                bz_here = surf2 ( pos(1), pos(3), nR_ar2, nZ_ar2, r_ar2, z_ar2, &
                    bz_ar2, nR_ar2, zp_bz, SplineSigma )

            else

                bR_here = curv2 (pos(1),nR_ar2,r_ar2,br_ar2,yp_br,SplineSigma)
                bPhi_here = curv2 (pos(1),nR_ar2,r_ar2,bt_ar2,yp_bt,SplineSigma)
                bZ_here = curv2 (pos(1),nR_ar2,r_ar2,bz_ar2,yp_bz,SplineSigma)

            endif
        endif

        if ( present (bMagHere) ) &
            bMagHere    = sqrt ( bR_here**2 + bPhi_here**2 + bz_here**2 )

        if ( present(psiHere).and.useEqdsk) &
            psiHere    = psi_here 

        if ( present (rhoHere).and.useEqdsk ) &
            rhoHere    = rho_here 

    end function dlg_interpB

end module interp
