module bField

implicit none

contains

    subroutine bFieldAR2 ( g )

        use aorsaNamelist, only: r0, noPoloidalField
        use grid
        use AR2Input, &
            only:r_ar2=>r,z_ar2=>z,br_ar2=>br,bt_ar2=>bt,bz_ar2=>bz,bMag_ar2=>bMag

        implicit none

        type(gridBlock), intent(inout) :: g
        integer :: w, i, j 
        real :: bHere(3)

        allocate ( &
            g%bMag(size(g%pt)), &
            g%bR_unit(size(g%pt)), &
            g%bT_unit(size(g%pt)), &
            g%bZ_unit(size(g%pt)), &
            g%rho(size(g%pt)), &
            g%mask(size(g%pt)) )

        do w=1,size(g%pt)

               i = g%pt(w)%i
               j = g%pt(w)%j

               g%bR_unit(w) = br_ar2(i,j) / bMag_ar2(i,j)
               g%bT_unit(w) = bt_ar2(i,j) / bMag_ar2(i,j)
               g%bZ_unit(w) = bz_ar2(i,j) / bMag_ar2(i,j)
               g%bMag(w) = bMag_ar2(i,j)

        enddo

    end subroutine bFieldAR2

end module bField
