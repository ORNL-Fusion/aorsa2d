module E_to_lab

contains

subroutine rotate_E_to_lab ( g, rhs )

    use solve
    use rotation
    use grid
    use aorsaNamelist, &
        only: nSpec, UseJpFromFile
    use interp, only: dlg_interpB 
    use read_jp_from_file, only: &
        file_nS=>nS, file_nR=>nR, file_nZ=>nZ, file_r=>r, file_z=>z, &
        file_Jp_r=>Jp_r, file_Jp_t=>Jp_t, file_Jp_z=>Jp_z!, &
        !file_Jp_r_s=>Jp_r_s, file_Jp_t_s=>Jp_t_s, file_Jp_z_s=>Jp_z_s

    implicit none

    type(gridBlock), intent(inout) :: g
    integer, intent(in) :: rhs

    complex :: ELab_RTZ(3), jPLab_RTZ(3), e1,e2,e3, jP1,jP2,jP3
    integer :: i, j, s

    real :: mag1, mag2
    real :: R_rtz_to_abp(3,3), R_abp_to_rtz(3,3)
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
            R_rtz_to_abp = RotMatHere(bRu,bTu,bZu)
            R_abp_to_rtz = transpose(R_rtz_to_abp)

#if __noU__==1
            ELab_RTZ = (/ g%eAlpha(i,j), g%eBeta(i,j), g%eb(i,j) /)
#else
            ELab_RTZ = &
                matMul (  R_abp_to_rtz, (/ g%eAlpha(i,j), g%eBeta(i,j), g%eb(i,j) /) )
#endif
            g%eR(i,j) = ELab_RTZ(1)
            g%eTh(i,j) = ELab_RTZ(2)
            g%eZ(i,j) = ELab_RTZ(3)

            do s=1,nSpec

#if __noU__==1
                jPLab_RTZ = (/ g%jAlpha(i,j,s), g%jBeta(i,j,s), g%jb(i,j,s) /)
#else
                jPLab_RTZ = &
                    matMul ( R_abp_to_rtz, (/ g%jAlpha(i,j,s), g%jBeta(i,j,s), g%jb(i,j,s) /) )
#endif 

                !ReplaceWithJpFromFile: &
                !if(UseJpFromFile)then

                !    g%jP_r(i,j,s) = g%file_JpR_s(i,j,s)
                !    g%jP_t(i,j,s) = g%file_JpT_s(i,j,s)
                !    g%jP_z(i,j,s) = g%file_JpZ_s(i,j,s)

                !else

                g%jP_r(i,j,s) = jPLab_RTZ(1)
                g%jP_t(i,j,s) = jPLab_RTZ(2)
                g%jP_z(i,j,s) = jPLab_RTZ(3)

                !endif ReplaceWithJpFromFile

            enddo    

        enddo
     enddo

end subroutine rotate_E_to_lab


end module E_to_lab
