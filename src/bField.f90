module bField

implicit none

real, allocatable, dimension(:,:) :: &
    bMod, bxn, byn, bzn

real :: r0__, z0__

contains

    subroutine bFieldEqdsk ()

        use aorsa2din_mod, &
        only: nModesX, nModesY
        use grid
        use interp
        use eqdsk_dlg

        implicit none

        integer :: i, j 
        real :: bHere(3)

        allocate ( &
            bMod(nModesX,nModesY), &
            bxn(nModesX,nModesY), &
            byn(nModesX,nModesY), &
            bzn(nModesX,nModesY) )

        do i=1,nModesX
            do j=1,nModesY

               bHere = dlg_interpB ( (/capR(i),0.0,y(j)/), &
                            bMagHere = bMod(i,j) )  
               bxn(i,j) = bHere(1) / bMod(i,j)
               byn(i,j) = bHere(3) / bMod(i,j)
               bzn(i,j) = bHere(2) / bMod(i,j)

            enddo
        enddo

        r0__  = rmaxis__
        z0__  = zmaxis__

    end subroutine bFieldEqdsk


    subroutine bFieldAnalytical ()

        use aorsa2din_mod, &
        only: nModesX, nModesY
        use grid

        implicit none

        real :: b0
        integer :: i, j

        r0__ = 1.0
        z0__ = 0.0
        b0  = 0.55

        allocate ( &
            bMod(nModesX,nModesY), &
            bxn(nModesX,nModesY), &
            byn(nModesX,nModesY), &
            bzn(nModesX,nModesY) )

        do i=1,nModesX
            do j=1,nModesY

               bxn(i,j) = 0.0 
               byn(i,j) = 0.0
               bzn(i,j) = b0!r0__ / capR(i) * b0 
               bMod(i,j) = sqrt ( bxn(i,j)**2 + byn(i,j)**2 + bzn(i,j)**2 )

            enddo
        enddo

        bxn = bxn / bMod
        byn = byn / bMod
        bzn = bzn / bMod

    end subroutine bFieldAnalytical

end module bField
