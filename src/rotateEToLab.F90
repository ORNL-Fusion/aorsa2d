module E_to_lab

contains

subroutine rotate_E_to_lab ( g, rhs )

    use solve
    use rotation
    use grid
    use aorsaNamelist, &
        only: nSpec
    use interp, only: dlg_interpB 

    implicit none

    type(gridBlock), intent(inout) :: g
    integer, intent(in) :: rhs

    complex :: ELab_RTZ(3), jPLab_RTZ(3)
    integer :: i, j, s

    real :: mag1, mag2
    real :: R_(3,3)
    real :: bTmp(3), bRu, bTu, bZu, bMagTmp

    if(.not.allocated(g%eR))allocate ( g%eR(g%nR,g%nZ), &
                g%eTh(g%nR,g%nZ), &
                g%eZ(g%nR,g%nZ) )

    if(.not.allocated(g%jP_r))allocate ( g%jP_r(g%nR,g%nZ,nSpec), &
                g%jP_t(g%nR,g%nZ,nSpec), &
                g%jP_z(g%nR,g%nZ,nSpec) )


    ! Rotate fields to the Lab frame and plus/minus
    ! ---------------------------------------------

    do i = 1, g%nR
        do j = 1, g%nZ

            bTmp = dlg_interpB ( (/g%R(i),0.0,g%z(j)/), bMagHere = bMagTmp )  
            bRu = bTmp(1)/bMagTmp
            bTu = bTmp(2)/bMagTmp
            bZu = bTmp(3)/bMagTmp
            R_ = RotMatHere(bRu,bTu,bZu)

#if __noU__==1
            ELab_RTZ = (/ g%eAlpha(i,j), g%eBeta(i,j), g%eb(i,j) /)
#else
            ELab_RTZ = &
                matMul ( transpose ( R_ ), &
                    (/ g%eAlpha(i,j), g%eBeta(i,j), g%eb(i,j) /) )
#endif
            g%eR(i,j) = ELab_RTZ(1)
            g%eTh(i,j) = ELab_RTZ(2)
            g%eZ(i,j) = ELab_RTZ(3)

            do s=1,nSpec

#if __noU__==1
                jPLab_RTZ = (/ g%jAlpha(i,j,s), g%jBeta(i,j,s), g%jb(i,j,s) /)
#else
                jPLab_RTZ = &
                    matMul ( transpose ( R_ ), &
                        (/ g%jAlpha(i,j,s), g%jBeta(i,j,s), g%jb(i,j,s) /) )
#endif 
                g%jP_r(i,j,s) = jPLab_RTZ(1)
                g%jP_t(i,j,s) = jPLab_RTZ(2)
                g%jP_z(i,j,s) = jPLab_RTZ(3)

            enddo    

        enddo
     enddo

end subroutine rotate_E_to_lab


end module E_to_lab
