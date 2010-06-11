module E_to_lab

contains

subroutine rotate_E_to_lab ( g )

    use solve
    use rotation
    use grid

    implicit none

    type(gridBlock), intent(inout) :: g

    complex :: ELab_RTZ(3)
    integer :: i, j

    allocate ( g%eR(g%nR,g%nZ), &
                g%eTh(g%nR,g%nZ), &
                g%eZ(g%nR,g%nZ) )


    ! Calculate E in the Lab frame and eplus, eminus
    ! ----------------------------------------------

    !isq2 = SQRT(0.5)
    do i = 1, g%nR
        do j = 1, g%nZ

            ELab_RTZ = &
                matMul ( transpose ( g%U_RTZ_to_ABb(i,j,:,:) ), &
                    (/ g%eAlpha(i,j), g%eBeta(i,j), g%eb(i,j) /) )

            g%eR(i,j) = ELab_RTZ(1)
            g%eTh(i,j) = ELab_RTZ(2)
            g%eZ(i,j) = ELab_RTZ(3)

            !eplus(i,j)  = isq2 * (ealpha(i,j) + zi * ebeta(i,j))
            !eminus(i,j) = isq2 * (ealpha(i,j) - zi * ebeta(i,j))

        enddo
     enddo

end subroutine rotate_E_to_lab


end module E_to_lab
