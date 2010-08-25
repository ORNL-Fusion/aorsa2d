module E_to_lab

contains

subroutine rotate_E_to_lab ( g )

    use solve
    use rotation
    use grid
    use aorsa2din_mod, &
    only: nSpec

    implicit none

    type(gridBlock), intent(inout) :: g

    complex :: ELab_RTZ(3), jPLab_RTZ(3)
    integer :: i, j, s

    real :: mag1, mag2

    allocate ( g%eR(g%nR,g%nZ), &
                g%eTh(g%nR,g%nZ), &
                g%eZ(g%nR,g%nZ) )

    allocate ( g%jP_r(g%nR,g%nZ,nSpec), &
                g%jP_t(g%nR,g%nZ,nSpec), &
                g%jP_z(g%nR,g%nZ,nSpec) )


    ! Rotate fields to the Lab frame and plus/minus
    ! ---------------------------------------------

    !isq2 = SQRT(0.5)
    do i = 1, g%nR
        do j = 1, g%nZ

            ELab_RTZ = &
                matMul ( transpose ( g%U_RTZ_to_ABb(i,j,:,:) ), &
                    (/ g%eAlpha(i,j), g%eBeta(i,j), g%eb(i,j) /) )

            g%eR(i,j) = ELab_RTZ(1)
            g%eTh(i,j) = ELab_RTZ(2)
            g%eZ(i,j) = ELab_RTZ(3)

            !! Check rotated field is equal in magnitude

            !mag1 = abs ( sqrt ( g%eAlpha(i,j)**2 + g%eBeta(i,j)**2 + g%eB(i,j)**2 ))
            !mag2 = abs ( sqrt ( g%eR(i,j)**2 + g%eTh(i,j)**2 + g%eZ(i,j)**2 ) )

            !if(abs(1-mag1/mag2)>1e-4)then
            !    write(*,*) 'ERROR: src/rotateEtoLab.f90 - magntiude was not invariant (E)'
            !    write(*,*) abs(1-mag1/mag2) 
            !    stop
            !endif


            do s=1,nSpec

                jPLab_RTZ = &
                    matMul ( transpose ( g%U_RTZ_to_ABb(i,j,:,:) ), &
                        (/ g%jAlpha(i,j,s), g%jBeta(i,j,s), g%jb(i,j,s) /) )
    
                g%jP_r(i,j,s) = jPLab_RTZ(1)
                g%jP_t(i,j,s) = jPLab_RTZ(2)
                g%jP_z(i,j,s) = jPLab_RTZ(3)

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

end subroutine rotate_E_to_lab


end module E_to_lab
