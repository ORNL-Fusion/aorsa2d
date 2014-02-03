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

    complex :: ELab_RTZ(3), jPLab_RTZ(3), e1,e2,e3
    integer :: i, j, s

    real :: mag1, mag2
    real :: R_(3,3), R_inv(3,3)
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
            R_inv = transpose(R_)
            !write(*,*) 'r: ', g%R(i)
            !write(*,*) 'bu: ', bTmp/bMagTmp
            !write(*,*) 'rot: ', R_

#if __noU__==1
            ELab_RTZ = (/ g%eAlpha(i,j), g%eBeta(i,j), g%eb(i,j) /)
#else
            ELab_RTZ = &
                matMul (  transpose ( R_ ), &
                    (/ g%eAlpha(i,j), g%eBeta(i,j), g%eb(i,j) /) )
            !e1 = g%eAlpha(i,j)
            !e2 = g%eBeta(i,j)
            !e3 = g%eb(i,j)

            !ELab_RTZ(1) = R_inv(1,1)*e1 + R_inv(2,1)*e2 + R_inv(3,1)*e3
            !ELab_RTZ(2) = R_inv(1,2)*e1 + R_inv(2,2)*e2 + R_inv(3,2)*e3
            !ELab_RTZ(3) = R_inv(1,3)*e1 + R_inv(2,3)*e2 + R_inv(3,3)*e3

            !ELab_RTZ(1) = R_inv(1,1)*e1 + R_inv(1,2)*e2 + R_inv(1,3)*e3
            !ELab_RTZ(2) = R_inv(2,1)*e1 + R_inv(2,2)*e2 + R_inv(2,3)*e3
            !ELab_RTZ(3) = R_inv(3,1)*e1 + R_inv(3,2)*e2 + R_inv(3,3)*e3


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
