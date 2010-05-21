module bField

implicit none

real, allocatable, dimension(:,:) :: &
    bMod, bxn, byn, bzn, rho
real, allocatable, dimension(:,:) :: &
    brn_, bthn_, bzn_

logical, allocatable :: mask(:,:)

contains

    subroutine bFieldEqdsk ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, r0
        use grid
        use interp
        use eqdsk_dlg

        implicit none

        integer :: i, j 
        real :: bHere(3)

        allocate ( &
            bMod(nPtsX,nPtsY), &
            bxn(nPtsX,nPtsY), &
            byn(nPtsX,nPtsY), &
            bzn(nPtsX,nPtsY), &
            rho(nPtsX,nPtsY) )

        allocate ( &
            brn_(nPtsX,nPtsY), &
            bthn_(nPtsX,nPtsY), &
            bzn_(nPtsX,nPtsY) )

        do i=1,nPtsX
            do j=1,nPtsY

               bHere = dlg_interpB ( (/capR(i),0.0,y(j)/), &
                            bMagHere = bMod(i,j), &
                            rhoHere = rho(i,j) )  
               bxn(i,j) = bHere(1) / bMod(i,j)
               byn(i,j) = bHere(3) / bMod(i,j)
               bzn(i,j) = bHere(2) / bMod(i,j)

               brn_(i,j) = bHere(1) / bMod(i,j)
               bthn_(i,j) = bHere(2) / bMod(i,j)
               bzn_(i,j) = bHere(3) / bMod(i,j)

            enddo
        enddo

        r0  = rmaxis__

        !   For regions outside the LCFS apply a correction
        !   if nessecary 

        where(rho>1)
            rho = 1
        endwhere

        where(rho<0)
            rho = 0
        endwhere

        allocate ( mask(nPtsX,nPtsY) )
        mask = is_inside_bbbs ()


    end subroutine bFieldEqdsk


    subroutine bFieldAnalytical ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, r0, b0, &
            bx_frac, by_frac, noPoloidalField
        use grid

        implicit none

        integer :: i, j
        real :: br_frac, bz_frac

        allocate ( &
            bMod(nPtsX,nPtsY), &
            bxn(nPtsX,nPtsY), &
            byn(nPtsX,nPtsY), &
            bzn(nPtsX,nPtsY) )

        allocate ( &
            brn_(nPtsX,nPtsY), &
            bthn_(nPtsX,nPtsY), &
            bzn_(nPtsX,nPtsY) )

        do i=1,nPtsX
            do j=1,nPtsY

                bzn(i,j) = r0 * b0 / capR(i) 
                bxn(i,j) = bzn(i,j) * bx_frac
                byn(i,j) = bzn(i,j) * by_frac 

                br_frac = bx_frac
                bz_frac = by_frac 
                 
                bthn_(i,j) = r0 * b0 / capR(i)  
                brn_(i,j) = bthn_(i,j) * br_frac
                bzn_(i,j) = bthn_(i,j) * bz_frac 
 

            enddo
        enddo

        if ( noPoloidalField ) then 

            bxn = 0
            byn = 0

            brn_   = 0
            bzn_   = 0

        endif

        bMod = sqrt ( bxn**2 + byn**2 + bzn**2 )
        
        bxn = bxn / bMod
        byn = byn / bMod
        bzn = bzn / bMod

        brn_ = brn_ / bMod
        bthn_ = bthn_ / bMod
        bzn_ = bzn_ / bMod

        allocate ( mask(nPtsX,nPtsY) )
        mask = .true.

    end subroutine bFieldAnalytical


    subroutine bFieldCircular ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, r0, b0, a, &
            rhoScale, rhoWidth, rhoPower, bPol_frac, noPoloidalField
        use grid
        use dlg

        implicit none

        integer :: i, j
        real :: x, reScaleFac

        allocate ( &
            bMod(nPtsX,nPtsY), &
            bxn(nPtsX,nPtsY), &
            byn(nPtsX,nPtsY), &
            bzn(nPtsX,nPtsY), &
            rho(nPtsX,nPtsY) )

        allocate ( &
            brn_(nPtsX,nPtsY), &
            bthn_(nPtsX,nPtsY), &
            bzn_(nPtsX,nPtsY) )

        do i=1,nPtsX
            do j=1,nPtsY

                x   = sqrt( (capR(i)-r0)**2 + (y(j)-0.0)**2  ) / a
                rho(i,j) = exp ( - ( x )**rhoPower / ( 2 * rhoWidth**2 ) )

                bzn(i,j) = r0 * b0 / capR(i) 
                bthn_(i,j) = r0 * b0 / capR(i)  

            enddo
        enddo

        rho = maxVal(rho)-rho

        bxn = -dlg_pDeriv ( capR, y, rho, 2 )
        byn = dlg_pDeriv ( capR, y, rho, 1 )

        reScaleFac   = maxVal ( sqrt ( bxn**2 + byn**2 ) )

        bxn = bxn / reScaleFac * bPol_frac * b0
        byn = byn / reScaleFac * bPol_frac * b0

        brn_ = bxn
        bzn_ = byn

        if ( noPoloidalField ) then 

            bxn = 0
            byn = 0

            brn_   = 0
            bzn_   = 0

        endif

        bMod = sqrt ( bxn**2 + byn**2 + bzn**2 )

        bxn = bxn / bMod
        byn = byn / bMod
        bzn = bzn / bMod

        brn_ = brn_ / bMod
        bthn_ = bthn_ / bMod
        bzn_ = bzn_ / bMod

        allocate ( mask(nPtsX,nPtsY) )
        mask = .true.

    end subroutine bFieldCircular


    subroutine soloviev ()

        use aorsa2din_mod, &
        only: eKappa, r0, b0, a, q0, nPtsX, nPtsY, psiExp, &
            noPoloidalField
        use grid, &
        only: capR, y
        use constants

        implicit none

        real :: q070qa, xiota0, rLim, psi_lim, psi1
        integer :: i, j
        real, allocatable :: psi(:,:), &
            bx(:,:), by(:,:), bz(:,:)
        real :: fRho, gaussian, q07qa

        allocate ( &
            bMod(nPtsX,nPtsY), &
            bxn(nPtsX,nPtsY), &
            byn(nPtsX,nPtsY), &
            bzn(nPtsX,nPtsY), &
            bx(nPtsX,nPtsY), &
            by(nPtsX,nPtsY), &
            bz(nPtsX,nPtsY) )


        allocate ( &
            brn_(nPtsX,nPtsY), &
            bthn_(nPtsX,nPtsY), &
            bzn_(nPtsX,nPtsY) )

        allocate ( psi(nPtsX,nPtsY), &
                    rho(nPtsX,nPtsY) )

        psi     = 0
        rho     = 0

        q07qa   = 0.0
        xiota0  = 1 / q0
        rLim    = r0 + a

        psi_lim = xiota0 * b0 / 2.0 * ( (rLim**2 - r0**2)**2 / 4. / r0**2 )

        do i = 1, nPtsX
            do j = 1, nPtsY

                ! Toroidal field
                ! --------------
                
                bz(i,j) = b0 * r0 / capr(i)
                bthn_(i,j) = b0 * r0 / capr(i)


                ! Poloidal field
                ! --------------

                psi1 = b0 * xiota0 / 2. * &
                  ((capR(i) * y(j) / r0 / eKappa)**2 &
                  + (capR(i)**2 - r0**2)**2 / 4. / r0**2 )

                psi(i,j) = psi1 / psi_lim
                rho(i,j) = sqrt(psi(i,j))

                if ( rho(i,j) <= 0.0 ) rho(i,j) = 1.0e-08
                if ( rho(i,j) > 1 ) rho(i,j) = 1

                bx(i,j) = -b0 * xiota0 * capR(i) * y(j) / r0**2 / eKappa**2
                by(i,j) = xiota0 * b0 / 2.0 * &
                      ( 2.0 * y(j)**2 / r0**2 / eKappa**2 &
                       + capR(i)**2 / r0**2 - 1.0)

                gaussian =  exp(-psi(i,j) / psiExp)

                ! use for regular runs (default)
                ! ------------------------------

                fRho = q07qa + (1. - q07qa) * gaussian
                !fRho = q07qa + (1. - q07qa) * gaussian**(0.5)

                bx(i,j) = bx(i,j) * fRho
                by(i,j) = by(i,j) * fRho

                brn_(i,j)    = bx(i,j)
                bzn_(i,j)    = by(i,j)

            enddo
        enddo

        if ( noPoloidalField ) then 

            bx = 0
            by = 0

            brn_   = 0
            bzn_   = 0

        endif

        bmod = sqrt ( bx**2 + by**2 + bz**2 )

        bxn = bx / bmod
        byn = by / bmod
        bzn = bz / bmod

        brn_ = brn_ / bmod
        bthn_ = bthn_ / bmod
        bzn_ = bzn_ / bmod

        deallocate ( bx, by, bz, psi )

        allocate ( mask(nPtsX,nPtsY) )
        mask = .true.

    end subroutine soloviev

end module bField
