module bField

implicit none

!real, allocatable, dimension(:,:) :: &
!    g%bMag, bxn, byn, bzn, rho
!real, allocatable, dimension(:,:) :: &
!    g%bR_unit, g%bT_unit, g%bZ_unit

!logical, allocatable :: g%mask(:,:)

contains

    subroutine bFieldEqdsk ( g )

        use aorsa2din_mod, &
        only: r0
        use grid
        use interp
        use eqdsk_dlg

        implicit none

        type(gridBlock), intent(inout) :: g
        integer :: i, j 
        real :: bHere(3)

        allocate ( &
            g%bMag(g%nR,g%nZ), &
            g%bR_unit(g%nR,g%nZ), &
            g%bT_unit(g%nR,g%nZ), &
            g%bZ_unit(g%nR,g%nZ), &
            g%rho(g%nR,g%nZ), &
            g%mask(g%nR,g%nZ) )

        do i=1,g%nR
            do j=1,g%nZ

               bHere = dlg_interpB ( (/g%R(i),0.0,g%z(j)/), &
                            bMagHere = g%bMag(i,j), &
                            rhoHere = g%rho(i,j) )  

               g%bR_unit(i,j) = bHere(1) / g%bMag(i,j)
               g%bT_unit(i,j) = bHere(3) / g%bMag(i,j)
               g%bZ_unit(i,j) = bHere(2) / g%bMag(i,j)

            enddo
        enddo

        r0  = rmaxis__

        !   For regions outside the LCFS apply a correction
        !   if nessecary 

        where(g%rho>1)
            g%rho = 1
        endwhere

        where(g%rho<0)
            g%rho = 0
        endwhere

        g%mask = is_inside_bbbs ( g )
        g%mask = .true.


    end subroutine bFieldEqdsk


    subroutine bFieldAnalytical ( g )

        use aorsa2din_mod, &
        only: r0, b0, &
            bx_frac, by_frac, noPoloidalField
        use grid

        implicit none

        type(gridBlock), intent(inout) :: g

        integer :: i, j
        real :: br_frac, bz_frac

        allocate ( &
            g%bMag(g%nR,g%nZ), &
            g%bR_unit(g%nR,g%nZ), &
            g%bT_unit(g%nR,g%nZ), &
            g%bZ_unit(g%nR,g%nZ), &
            g%rho(g%nR,g%nZ), &
            g%mask(g%nR,g%nZ) )

        do i=1,g%nR
            do j=1,g%nZ

                br_frac = bx_frac
                bz_frac = by_frac 
                 
                g%bT_unit(i,j) = r0 * b0 / g%R(i)  
                g%bR_unit(i,j) = g%bT_unit(i,j) * br_frac
                g%bZ_unit(i,j) = g%bT_unit(i,j) * bz_frac 
 

            enddo
        enddo

        if ( noPoloidalField ) then 

            g%bR_unit   = 0
            g%bZ_unit   = 0

        endif

        g%bMag = sqrt ( g%bR_unit**2 + g%bT_unit**2 + g%bZ_unit**2 )
        
        g%bR_unit = g%bR_unit / g%bMag 
        g%bT_unit = g%bT_unit / g%bMag 
        g%bZ_unit = g%bZ_unit / g%bMag

        allocate ( g%mask(g%nR,g%nZ) )
        g%mask = .true.

    end subroutine bFieldAnalytical


    subroutine bFieldCircular ( g )

        use aorsa2din_mod, &
        only: r0, b0, a, &
            rhoScale, rhoWidth, rhoPower, bPol_frac, noPoloidalField
        use grid
        use dlg

        implicit none

        integer :: i, j
        real :: x, reScaleFac

        type(gridBlock), intent(inout) :: g

        allocate ( &
            g%bMag(g%nR,g%nZ), &
            g%bR_unit(g%nR,g%nZ), &
            g%bT_unit(g%nR,g%nZ), &
            g%bZ_unit(g%nR,g%nZ), &
            g%rho(g%nR,g%nZ), &
            g%mask(g%nR,g%nZ) )

        do i=1,g%nR
            do j=1,g%nZ

                x   = sqrt( (g%R(i)-r0)**2 + (g%Z(j)-0.0)**2  ) / a
                g%rho(i,j) = exp ( - ( x )**rhoPower / ( 2 * rhoWidth**2 ) )

                g%bT_unit(i,j) = r0 * b0 / g%R(i)  

            enddo
        enddo

        g%rho = maxVal(g%rho)-g%rho

        g%bR_unit = -dlg_pDeriv ( g%R, g%z, g%rho, 2 )
        g%bZ_unit = dlg_pDeriv ( g%R, g%z, g%rho, 1 )

        reScaleFac   = maxVal ( sqrt ( g%bR_unit**2 + g%bZ_unit**2 ) )

        g%bR_unit = g%bR_unit / reScaleFac * bPol_frac * b0
        g%bZ_unit = g%bZ_unit / reScaleFac * bPol_frac * b0


        if ( noPoloidalField ) then 

            g%bR_unit   = 0
            g%bZ_unit   = 0

        endif

        g%bMag = sqrt ( g%bR_unit**2 + g%bT_Unit**2 + g%bZ_Unit**2 )

        g%bR_unit = g%bR_unit / g%bMag
        g%bT_unit = g%bT_unit / g%bMag
        g%bZ_unit = g%bZ_unit / g%bMag

        allocate ( g%mask(g%nR,g%nZ) )
        g%mask = .true.

    end subroutine bFieldCircular


    subroutine soloviev ( g )

        use aorsa2din_mod, &
        only: eKappa, r0, b0, a, q0, psiExp, &
            noPoloidalField
        use grid
        use constants

        implicit none

        type(gridBlock), intent(inout) :: g

        real :: q070qa, xiota0, rLim, psi_lim, psi1
        integer :: i, j
        real, allocatable :: psi(:,:), &
            bx(:,:), by(:,:), bz(:,:)
        real :: fRho, gaussian, q07qa

        allocate ( &
            g%bMag(g%nR,g%nZ), &
            bx(g%nR,g%nZ), &
            by(g%nR,g%nZ), &
            bz(g%nR,g%nZ) )


        allocate ( &
            g%bR_unit(g%nR,g%nZ), &
            g%bT_unit(g%nR,g%nZ), &
            g%bZ_unit(g%nR,g%nZ) )

        allocate ( psi(g%nR,g%nZ), &
                    g%rho(g%nR,g%nZ) )

        psi     = 0
        g%rho     = 0

        q07qa   = 0.0
        xiota0  = 1 / q0
        rLim    = r0 + a

        psi_lim = xiota0 * b0 / 2.0 * ( (rLim**2 - r0**2)**2 / 4. / r0**2 )

        do i = 1, g%nR
            do j = 1, g%nZ

                ! Toroidal field
                ! --------------
                
                bz(i,j) = b0 * r0 / g%R(i)
                g%bT_unit(i,j) = b0 * r0 / g%R(i)


                ! Poloidal field
                ! --------------

                psi1 = b0 * xiota0 / 2. * &
                  ((g%R(i) * g%Z(j) / r0 / eKappa)**2 &
                  + (g%R(i)**2 - r0**2)**2 / 4. / r0**2 )

                psi(i,j) = psi1 / psi_lim
                g%rho(i,j) = sqrt(psi(i,j))

                if ( g%rho(i,j) <= 0.0 ) g%rho(i,j) = 1.0e-08
                if ( g%rho(i,j) > 1 ) g%rho(i,j) = 1

                bx(i,j) = -b0 * xiota0 * g%R(i) * g%Z(j) / r0**2 / eKappa**2
                by(i,j) = xiota0 * b0 / 2.0 * &
                      ( 2.0 * g%Z(j)**2 / r0**2 / eKappa**2 &
                       + g%R(i)**2 / r0**2 - 1.0)

                gaussian =  exp(-psi(i,j) / psiExp)

                ! use for regular runs (default)
                ! ------------------------------

                fRho = q07qa + (1. - q07qa) * gaussian
                !fRho = q07qa + (1. - q07qa) * gaussian**(0.5)

                bx(i,j) = bx(i,j) * fRho
                by(i,j) = by(i,j) * fRho

                g%bR_unit(i,j)    = bx(i,j)
                g%bZ_unit(i,j)    = by(i,j)

            enddo
        enddo

        if ( noPoloidalField ) then 

            bx = 0
            by = 0

            g%bR_unit   = 0
            g%bZ_unit   = 0

        endif

        g%bMag = sqrt ( bx**2 + by**2 + bz**2 )

        g%bR_unit = g%bR_unit / g%bMag
        g%bT_unit = g%bT_unit / g%bMag
        g%bZ_unit = g%bZ_unit / g%bMag

        deallocate ( bx, by, bz, psi )

        allocate ( g%mask(g%nR,g%nZ) )
        g%mask = .true.

    end subroutine soloviev

end module bField
