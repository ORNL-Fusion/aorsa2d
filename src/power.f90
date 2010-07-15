module power

contains

subroutine current ( g )

    use grid
    use read_data

    implicit none
   
    type(gridBlock), intent(inout) :: g 

    complex, allocatable :: sigma(:,:,:,:,:,:)
    complex :: ek_nm(3), jVec(3), thisSigma(3,3)

    integer :: i, j, n, m, iStat

    allocate ( &
        g%jAlpha(g%nR,g%nZ), &
        g%jBeta(g%nR,g%nZ), &
        g%jB(g%nR,g%nZ) )

    g%jAlpha = 0
    g%jBeta = 0
    g%jB = 0

    allocate ( &
        sigma(g%nR,g%nZ,g%nMin:g%nMax,g%mMin:g%mMax,3,3), stat = iStat )

    if(iStat/=0)then
            write(*,*) 'ERROR src/power.f90 - allocation failed :('
            stop
    endif 

    call read_sigma ( 'sigma'//g%fNumber//'.nc', sigma = sigma ) 

    do i=1,g%nR
        do j=1,g%nZ
            do n=g%nMin,g%nMax
                do m=g%mMin,g%mMax

                    ek_nm(1) = g%eAlphak(n,m)
                    ek_nm(2) = g%eBetak(n,m)
                    ek_nm(3) = g%eBk(n,m) 

                    thisSigma = sigma(i,j,n,m,:,:)

                    jVec = matMul ( thisSigma, ek_nm )

                    g%jAlpha(i,j) = g%jAlpha(i,j) + jVec(1)
                    g%jBeta(i,j) = g%jBeta(i,j) + jVec(2)
                    g%jB(i,j) = g%jB(i,j) + jVec(3)

                enddo
            enddo
        enddo
    enddo

    deallocate ( sigma )

end subroutine current 


subroutine jDotE ( g )

    use grid

    implicit none
   
    type(gridBlock), intent(inout) :: g 

    integer :: i, j
    complex :: eHere(3), jHere(3)

    allocate ( &
        g%jouleHeating(g%nR,g%nZ) )

    do i=1,g%nR
        do j=1,g%nZ

            eHere = (/ g%eAlpha(i,j), g%eBeta(i,j), g%eB(i,j) /)
            jHere = (/ g%jAlpha(i,j), g%jBeta(i,j), g%jB(i,j) /)

            g%jouleHeating(i,j) = &
                0.5 * realpart ( dot_product ( conjg ( eHere ), jHere ) )

        enddo
    enddo  

end subroutine jDotE

end module power
