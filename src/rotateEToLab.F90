module E_to_lab

contains

subroutine rotate_E_to_lab ( g )

    use solve
    use rotation
    use grid
    use aorsaNamelist, &
        only: nSpec
    use interp, only: dlg_interpB 
    use antenna, only: NRHS

    implicit none

    type(gridBlock), intent(inout) :: g

    complex :: ELab_RTZ(3), jPLab_RTZ(3)
    integer :: i, j, s, rhs

    real :: mag1, mag2
    real :: R_(3,3)
    real :: bTmp(3), bRu, bTu, bZu, bMagTmp

    allocate ( g%eR(g%nR,g%nZ,NRHS), &
                g%eTh(g%nR,g%nZ,NRHS), &
                g%eZ(g%nR,g%nZ,NRHS) )

    allocate ( g%jP_r(g%nR,g%nZ,nSpec,NRHS), &
                g%jP_t(g%nR,g%nZ,nSpec,NRHS), &
                g%jP_z(g%nR,g%nZ,nSpec,NRHS) )


    ! Rotate fields to the Lab frame and plus/minus
    ! ---------------------------------------------

    !isq2 = SQRT(0.5)
   
    do rhs=1,NRHS

    do i = 1, g%nR
        do j = 1, g%nZ

            bTmp = dlg_interpB ( (/g%R(i),0.0,g%z(j)/), bMagHere = bMagTmp )  
            bRu = bTmp(1)/bMagTmp
            bTu = bTmp(2)/bMagTmp
            bZu = bTmp(3)/bMagTmp
            R_ = RotMatHere(bRu,bTu,bZu)
            !R_ = g%U_RTZ_to_ABb(i,j,:,:)

#if __noU__==1
            ELab_RTZ = (/ g%eAlpha(i,j,rhs), g%eBeta(i,j,rhs), g%eb(i,j,rhs) /)
#else
            ELab_RTZ = &
                matMul ( transpose ( R_ ), &
                    (/ g%eAlpha(i,j,rhs), g%eBeta(i,j,rhs), g%eb(i,j,rhs) /) )
#endif
            g%eR(i,j,rhs) = ELab_RTZ(1)
            g%eTh(i,j,rhs) = ELab_RTZ(2)
            g%eZ(i,j,rhs) = ELab_RTZ(3)

            !! Check rotated field is equal in magnitude

            !mag1 = abs ( sqrt ( g%eAlpha(i,j)**2 + g%eBeta(i,j)**2 + g%eB(i,j)**2 ))
            !mag2 = abs ( sqrt ( g%eR(i,j)**2 + g%eTh(i,j)**2 + g%eZ(i,j)**2 ) )

            !if(abs(1-mag1/mag2)>1e-4)then
            !    write(*,*) 'ERROR: src/rotateEtoLab.f90 - magntiude was not invariant (E)'
            !    write(*,*) abs(1-mag1/mag2) 
            !    stop
            !endif


            do s=1,nSpec

#if __noU__==1
                jPLab_RTZ = (/ g%jAlpha(i,j,s,rhs), g%jBeta(i,j,s,rhs), g%jb(i,j,s,rhs) /)
#else
                jPLab_RTZ = &
                    matMul ( transpose ( R_ ), &
                        (/ g%jAlpha(i,j,s,rhs), g%jBeta(i,j,s,rhs), g%jb(i,j,s,rhs) /) )
#endif 
                g%jP_r(i,j,s,rhs) = jPLab_RTZ(1)
                g%jP_t(i,j,s,rhs) = jPLab_RTZ(2)
                g%jP_z(i,j,s,rhs) = jPLab_RTZ(3)

                !! Check rotation

                !mag1 = abs ( sqrt ( g%jAlpha(i,j,s)**2 + g%jBeta(i,j,s)**2 + g%jB(i,j,s)**2 ))
                !mag2 = abs ( sqrt ( g%jP_r(i,j,s)**2 + g%jP_t(i,j,s)**2 + g%jP_z(i,j,s)**2 ) )

                !if(abs(1-mag1/mag2)>1e-4)then
                !    write(*,*) 'ERROR: src/rotateEtoLab.f90 - magntiude was not invariant (Jp)'
                !    write(*,*) abs(1-mag1/mag2)
                !    stop 
                !endif

            enddo    

            !eplus(i,j)  = isq2 * (ealpha(i,j) + zi * ebeta(i,j))
            !eminus(i,j) = isq2 * (ealpha(i,j) - zi * ebeta(i,j))

        enddo
     enddo

     enddo ! rhs loop

end subroutine rotate_E_to_lab


end module E_to_lab
