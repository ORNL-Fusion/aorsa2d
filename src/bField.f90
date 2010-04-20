module bField

implicit none

real, allocatable, dimension(:,:) :: &
    bMod, bxn, byn, bzn, rho
real, allocatable, dimension(:,:) :: &
    brn_, bthn_, bzn_

real :: r0__, z0__

contains

    subroutine bFieldEqdsk ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY
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

        r0__  = rmaxis__
        z0__  = zmaxis__


        !   For regions outside the LCFS apply a correction
        !   if nessecary 

        where(rho>1)
            rho = 1
        endwhere

        where(rho<0)
            rho = 0
        endwhere


    end subroutine bFieldEqdsk


    subroutine bFieldAnalytical ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, r0, b0, &
            bx_frac, by_frac
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

                bzn(i,j) = r0 / capR(i) * b0 
                bxn(i,j) = bzn(i,j) * bx_frac
                byn(i,j) = bzn(i,j) * by_frac 

                br_frac = bx_frac
                bz_frac = by_frac 
                 
                bthn_(i,j) = r0 / capR(i) * b0 
                brn_(i,j) = bthn_(i,j) * br_frac
                bzn_(i,j) = bthn_(i,j) * bz_frac 
 
                bMod(i,j) = sqrt ( bxn(i,j)**2 + byn(i,j)**2 + bzn(i,j)**2 )

            enddo
        enddo

        bxn = bxn / bMod
        byn = byn / bMod
        bzn = bzn / bMod

        brn_ = brn_ / bMod
        bthn_ = bthn_ / bMod
        bzn_ = bzn_ / bMod

    end subroutine bFieldAnalytical

end module bField
