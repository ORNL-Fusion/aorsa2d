module bField

implicit none

real, allocatable, dimension(:,:) :: &
    bMod, bxn, byn, bzn

contains

    subroutine bFieldEqdsk ()

        use aorsa2din_mod, &
        only: nModesX, nModesY
        use grid
        use interp

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
               bzn(i,j) = -bHere(2) / bMod(i,j)

            enddo
        enddo

    end subroutine bFieldEqdsk

end module bField
