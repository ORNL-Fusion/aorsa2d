module E_to_lab

contains

subroutine rotate_E_to_lab ()

    use solve
    use rotation
    use aorsa2din_mod, &
    only: nPtsX, nPtsY

    implicit none

    complex :: ELab_RTZ(3)
    integer :: i, j

    allocate ( eR(nPtsX,nPtsY), &
                eTh(nPtsX,nPtsY), &
                eZ(nPtsX,nPtsY) )


    ! Calculate E in the Lab frame and eplus, eminus
    ! ----------------------------------------------

    !isq2 = SQRT(0.5)
    do i = 1, nPtsX
        do j = 1, nPtsY

            ELab_RTZ = &
                matMul ( transpose ( U_RTZ_to_ABb(i,j,:,:) ), &
                    (/ eAlpha(i,j), eBeta(i,j), eb(i,j) /) )

            !ELab_RTZ = &
            !    matMul ( U_RTZ_to_ABb(i,j,:,:), &
            !        (/ eAlpha(i,j), eBeta(i,j), eb(i,j) /) )


            eR(i,j) = ELab_RTZ(1)
            eTh(i,j) = ELab_RTZ(2)
            eZ(i,j) = ELab_RTZ(3)

            !eplus(i,j)  = isq2 * (ealpha(i,j) + zi * ebeta(i,j))
            !eminus(i,j) = isq2 * (ealpha(i,j) - zi * ebeta(i,j))

        enddo
     enddo

end subroutine rotate_E_to_lab


end module E_to_lab
