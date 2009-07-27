module gc_terms
    implicit none
    save
    real, allocatable :: bCurvature_R(:,:), &
        bCurvature_phi(:,:), bCurvature_z(:,:), &
        bGradient_R(:,:), bGradient_phi(:,:), &
        bGradient_z(:,:), bDotGradB(:,:)

contains
    subroutine bCurvature ()
        use eqdsk_dlg
        use dlg
        use constants
        implicit none
       
        integer :: j
        real, dimension(nw,nh) :: bR_B, bPhi_B, bz_B, &
            bDotGradB_R, bDotGradB_phi, bDotGradB_z, &
            bR_B_dR, bR_B_dz, bPhi_B_dR, bPhi_B_dz, &
            bz_B_dR, bz_B_dz, r2D, omega, &
            gradB_R, gradB_z

        do j=1,nh
            r2D(:,j)    = r
        end do

        omega   = q * bMag__ / mi

        bR_B    = bR / bMag__
        bPhi_B  = bPhi / bMag__
        bz_B    = bz__ / bMag__ 

        bR_B_dR = dlg_pDeriv ( bR_B, 1, rStep )
        bR_B_dz = dlg_pDeriv ( bR_B, 2, zStep )
        
        bPhi_B_dR   = dlg_pDeriv ( bPhi_B, 1, rStep )
        bPhi_B_dz   = dlg_pDeriv ( bPhi_B, 2, zStep )
 
        bz_B_dR = dlg_pDeriv ( bz_B, 1, rStep )
        bz_B_dz = dlg_pDeriv ( bz_B, 2, zStep )
    
        bDotGradB_R = bR_B * bR_B_dR + bz_B * bR_B_dz - bPhi_B**2 / r2D
        bDotGradB_phi   = bPhi_B * bR_B / r2D + bR_B * bPhi_B_dR &
            + bz_B * bPhi_B_dz
        bDotGradB_z = bR_B * bz_B_dR + bz_B * bz_B_dz

        allocate ( bCurvature_R(nw,nh), bCurvature_phi(nw,nh), bCurvature_z(nw,nh) )

        bCurvature_R = ( bPhi_B * bDotGradB_z - bz_B * bDotGradB_phi ) / omega
        bCurvature_phi   = -1.0 * ( bR_B * bDotGradB_z - bz_B * bDotGradB_R ) &
            / omega
        bCurvature_z = ( bR_B * bDotGradB_phi - bPhi_B * bDotGradB_R ) / omega

        ! Also in here is the b.Grad b for calculating vPar

        gradB_R = dlg_pDeriv ( bMag__, 1, rStep )
        gradB_z = dlg_pDeriv ( bMag__, 2, zStep )

        allocate ( bDotGradB(nw,nh) )

        bDotGradB   = bR_B * gradB_R + bz_B * gradB_z
       
    end subroutine bCurvature

    subroutine bGradient ()
        use eqdsk_dlg
        use dlg
        use constants
        !use dislin
        implicit none
        real, dimension(nw,nh) :: bR_B, bPhi_B, bz_B, &
            lnB_dR, lnB_dz, lnB, omega

        !!   Plot variables
        !
        !integer :: nLevs, i
        !real, allocatable :: levels(:)
        !real :: levStep

        omega   = q * bMag__ / mi
        lnB = log ( bMag__ )
        lnB_dR  = dlg_pDeriv ( lnB, 1, rStep )
        lnB_dz  = dlg_pDeriv ( lnB, 2, zStep )

        bR_B    = bR / bMag__
        bPhi_B  = bPhi / bMag__
        bz_B    = bz__ / bMag__ 

        allocate ( bGradient_R(nw,nh), bGradient_phi(nw,nh), &
            bGradient_z(nw,nh) )

        bGradient_R = bPhi_B * lnB_dz / ( 2.0 * omega )
        bGradient_phi   = -1.0 * ( bR_B * lnB_dz - bz_B * lnB_dR ) &
            / ( 2.0 * omega )
        bGradient_z = -1.0 * bPhi_B * lnB_dR / ( 2.0 * omega ) 

        !!call setFil ( 'eqdsk.eps' ) 
        !!call setPag ( 'DA4P' )
        !call scrMod ( 'REVERS' )
        !call metaFl ( 'XWIN' )
        !call disIni ()
        !call graf ( 0.0, 3.0, 0.0, 0.5, -2.0, 2.0, -2.0, 0.5 ) 
        !call noClip ()
        !call color ('BLUE')         
        !
        !nLevs   = 101 
        !levStep    = maxVal ( abs ( bz__ ) ) / ( nLevs / 2 )
        !allocate ( levels(nLevs) )
        !levels = (/ (i*levStep,i=-nLevs/2,nLevs/2) /)
        !write (*,*) levels
        !do i=0,nLevs-1
        !    call contur ( r, nw, z, nh, bz__, levels(i) ) 
        !end do
        !call endGrf ()
        !read (*,*) 

    end subroutine bGradient

end module gc_terms
