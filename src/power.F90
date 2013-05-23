module power

contains

subroutine current ( g )

    use grid
    use read_data
    use aorsaNamelist, &
        only: nSpec, iSigma, fracOfModesInSolution
    use sigma
    use parallel
    use profiles, &
        only: k0, omgrf, mSpec
    use constants
    use antenna, only: NRHS

    implicit none
   
    type(gridBlock), intent(inout) :: g 

    complex :: ek_nm(3), jVec(3), thisSigma(3,3)
    complex :: bFn

    integer :: s, w, rhs
    real :: kr, kz, kVec_stix(3)
    complex, allocatable :: jAlphaTmp(:,:), jBetaTmp(:,:), jBTmp(:,:)

    type(spatialSigmaInput_cold) :: sigmaIn_cold
#if __noU__==1
    real :: R_(3,3)
#endif

    allocate ( &
        g%jAlpha(g%nR,g%nZ,nSpec,NRHS), &
        g%jBeta(g%nR,g%nZ,nSpec,NRHS), &
        g%jB(g%nR,g%nZ,nSpec,NRHS) )

    allocate ( jAlphaTmp(g%nR,g%nZ), jBetaTmp(g%nR,g%nZ), jBTmp(g%nR,g%nZ) )

    g%jAlpha = 0
    g%jBeta = 0
    g%jB = 0

#if __sigma__ != 2
    allocate ( &
        sigma(g%nR,g%nZ,g%nMin:g%nMax,g%mMin:g%mMax,3,3,nSpec), stat = iStat )

    if(iStat/=0)then
            write(*,*) 'ERROR src/power.f90 - allocation failed :('
            stop
    endif 

    call read_sigma ( 'sigma'//g%fNumber//'.nc', sigma = sigmaAll ) 
#endif

    species: &
    do s=1,nSpec

        !do i=1,g%nR
        !    do j=1,g%nZ
        !        do n=g%nS,g%nF
        !            do m=g%mS,g%mF

        workList: &
        do w=1,size(g%wl)

            twoThirdsRule: &
            if(g%wl(w)%m >= g%mMin*fracOfModesInSolution .and. g%wl(w)%m <= g%mMax*fracOfModesInSolution &
                .and. g%wl(w)%n >= g%nMin*fracOfModesInSolution .and. g%wl(w)%n <= g%nMax*fracOfModesInSolution ) then

            
                        bFn = g%xx(g%wl(w)%n,g%wl(w)%i) * g%yy(g%wl(w)%m,g%wl(w)%j)
                        !bFn = xBasis(n,g%rNorm(i)) * yBasis(m,g%zNorm(j))

#if __sigma__ != 2
                        thisSigma = sigmaAll(g%wl(w)%i,g%wl(w)%j,g%wl(w)%n,g%wl(w)%m,:,:,s)
#else
                        if(chebyshevX) then
                            if(g%wl(w)%n>1) then
                                kr = g%wl(w)%n / sqrt ( sin ( pi * (g%rNorm(g%wl(w)%i)+1)/2  ) ) * g%normFacR 
                            else
                                kr = g%wl(w)%n * g%normFacR
                            endif
                        else
                            kr = g%wl(w)%n * g%normFacR
                        endif

                        if(chebyshevY) then
                            if(g%wl(w)%m>1) then
                                kz = g%wl(w)%m / sqrt ( sin ( pi * (g%zNorm(g%wl(w)%j)+1)/2 ) ) * g%normFacZ 
                            else
                                kz = g%wl(w)%m * g%normFacZ
                            endif
                        else
                            kz = g%wl(w)%m * g%normFacZ
                        endif


                        hotPlasma:& 
                        if (iSigma==1 .and. (.not. g%isMetal(g%wl(w)%iPt)) ) then        

                            kVec_stix = matMul( g%U_RTZ_to_ABb(g%wl(w)%iPt,:,:), &
                                (/ kr, g%kPhi(g%wl(w)%i), kz /) ) 

                            thisSigma = sigmaHot_maxwellian &
                                ( mSpec(s), &
                                g%ktSpec(g%wl(w)%iPt,s), &
                                g%omgc(g%wl(w)%iPt,s), &
                                g%omgp2(g%wl(w)%iPt,s), &
                                kVec_stix, g%R(g%wl(w)%i), &
                                omgrf, k0, &
                                g%k_cutoff, s, &
                                g%sinTh(g%wl(w)%iPt), &
                                g%bPol(g%wl(w)%iPt), g%bMag(g%wl(w)%iPt), &
                                g%gradPrlB(g%wl(w)%iPt), &
                                g%nuOmg(g%wl(w)%iPt,s) )

                        endif hotPlasma

                        coldPlasma: &
                        if (iSigma==0 .and. (.not. g%isMetal(g%wl(w)%iPt)) ) then 

                            sigmaIn_cold = spatialSigmaInput_cold( &
                                g%omgc(g%wl(w)%iPt,s), &
                                g%omgp2(g%wl(w)%iPt,s), &
                                omgrf, &
                                g%nuOmg(g%wl(w)%iPt,s) )

                            !thisSigma = sigmaCold_stix &
                            !    ( g%omgc(g%wl(w)%i,g%wl(w)%j,s), &
                            !    g%omgp2(g%wl(w)%i,g%wl(w)%j,s), omgrf, &
                            !    g%nuOmg(g%wl(w)%i,g%wl(w)%j) )
                            thisSigma = sigmaCold_stix ( sigmaIn_cold )


                        endif coldPlasma
#if __noU__==1
                        ! Rotate sigma from alp,bet,prl to r,t,z
                        R_ = g%U_RTZ_to_ABb(g%wl(w)%iPt,:,:)
                        thisSigma = matmul(transpose(R_),matmul(thisSigma,R_))
#endif

                        ! Metal
                        ! -----

                        if (g%isMetal(g%wl(w)%iPt)) then 

                            thisSigma = 0
                            thisSigma(1,1) = metal 
                            thisSigma(2,2) = metal
                            thisSigma(3,3) = metal

                        endif
#endif
                        do rhs=1,NRHS
    
                            ek_nm(1) = g%eAlphak(g%wl(w)%n,g%wl(w)%m,rhs)
                            ek_nm(2) = g%eBetak(g%wl(w)%n,g%wl(w)%m,rhs)
                            ek_nm(3) = g%eBk(g%wl(w)%n,g%wl(w)%m,rhs) 

                            jVec = matMul ( thisSigma, ek_nm ) 

                            g%jAlpha(g%wl(w)%i,g%wl(w)%j,s,rhs) = g%jAlpha(g%wl(w)%i,g%wl(w)%j,s,rhs) + jVec(1) * bFn
                            g%jBeta(g%wl(w)%i,g%wl(w)%j,s,rhs) = g%jBeta(g%wl(w)%i,g%wl(w)%j,s,rhs) + jVec(2) * bFn
                            g%jB(g%wl(w)%i,g%wl(w)%j,s,rhs) = g%jB(g%wl(w)%i,g%wl(w)%j,s,rhs) + jVec(3) * bFn

                        enddo ! rhs loop

                        !g%jAlpha(i,j,s) = g%jAlpha(i,j,s) &
                        !    + ( sigma(i,j,n,m,1,1,s) * g%eAlphak(n,m) &
                        !    + sigma(i,j,n,m,1,2,s) * g%eBetak(n,m) &
                        !    + sigma(i,j,n,m,1,3,s) * g%eBk(n,m) ) * bFn 

                        !g%jBeta(i,j,s) = g%jBeta(i,j,s) &
                        !    + ( sigma(i,j,n,m,2,1,s) * g%eAlphak(n,m) &
                        !    + sigma(i,j,n,m,2,2,s) * g%eBetak(n,m) &
                        !    + sigma(i,j,n,m,2,3,s) * g%eBk(n,m) ) * bFn 

                        !g%jB(i,j,s) = g%jB(i,j,s) &
                        !    + ( sigma(i,j,n,m,3,1,s) * g%eAlphak(n,m) &
                        !    + sigma(i,j,n,m,3,2,s) * g%eBetak(n,m) &
                        !    + sigma(i,j,n,m,3,3,s) * g%eBk(n,m) ) * bFn 

        !            enddo
        !        enddo
        !    enddo
        !enddo

            endif twoThirdsRule

        enddo workList

    enddo species

#if __sigma__ != 2
    deallocate ( sigmaAll )
#else
#ifdef par

    ! Switch to individual 2D arrays for the sum over processors

    do rhs=1,NRHS
    do s=1,nSpec

        jAlphaTmp = g%jAlpha(:,:,s,rhs)
        jBetaTmp = g%jBeta(:,:,s,rhs)
        jBTmp = g%jB(:,:,s,rhs)

        call cGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, jAlphaTmp, g%nR, -1, -1 )
        call cGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, jBetaTmp, g%nR, -1, -1 )
        call cGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, jBTmp, g%nR, -1, -1 )

        call blacs_barrier ( iContext, 'All' ) 

        g%jAlpha(:,:,s,rhs) = jAlphaTmp
        g%jBeta(:,:,s,rhs) = jBetaTmp
        g%jB(:,:,s,rhs) = jBTmp

    enddo
    enddo ! rhs loop
#endif
#endif

end subroutine current 


subroutine jDotE ( g )

    use grid
    use aorsaNamelist, &
        only: nSpec
    use antenna, only: NRHS

    implicit none
   
    type(gridBlock), intent(inout) :: g 

    integer :: i, j, s, rhs

    allocate ( &
        g%jouleHeating(g%nR,g%nZ,nSpec,NRHS) )


    do rhs=1,NRHS
    species: &
    do s=1,nSpec

        do i=1,g%nR
            do j=1,g%nZ

                !eHere = (/ g%eAlpha(i,j), g%eBeta(i,j), g%eB(i,j) /)
                !jHere = (/ g%jAlpha(i,j,s), g%jBeta(i,j,s), g%jB(i,j,s) /)

                !g%jouleHeating(i,j,s) = &
                !    0.5 * realpart ( dot_product ( eHere, jHere ) )

                g%jouleHeating(i,j,s,rhs) = &
                    0.5 * real(real( conjg(g%eAlpha(i,j,rhs)) * g%jAlpha(i,j,s,rhs) &
                                    + conjg(g%eBeta(i,j,rhs)) * g%jBeta(i,j,s,rhs) &
                                    + conjg(g%eB(i,j,rhs)) * g%jB(i,j,s,rhs)  ))

            enddo
        enddo  

    enddo species
    enddo ! rhs loop

end subroutine jDotE

end module power
