module power

contains

subroutine current ( g )

    use grid
    use read_data
    use aorsa2din_mod, &
        only: nSpec, iSigma, fracOfModesInSolution
    use sigma_mod
    use parallel
    use profiles, &
        only: k0, omgrf, mSpec
    use constants

    implicit none
   
    type(gridBlock), intent(inout) :: g 

    complex, allocatable :: sigma(:,:,:,:,:,:,:)
    complex :: ek_nm(3), jVec(3), thisSigma(3,3)
    complex :: bFn

    integer :: i, j, n, m, iStat, s, w
    real :: kr, kz, kVec_stix(3)
    complex, allocatable :: jAlphaTmp(:,:), jBetaTmp(:,:), jBTmp(:,:)

    allocate ( &
        g%jAlpha(g%nR,g%nZ,nSpec), &
        g%jBeta(g%nR,g%nZ,nSpec), &
        g%jB(g%nR,g%nZ,nSpec) )

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

    call read_sigma ( 'sigma'//g%fNumber//'.nc', sigma = sigma ) 
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

                        ek_nm(1) = g%eAlphak(g%wl(w)%n,g%wl(w)%m)
                        ek_nm(2) = g%eBetak(g%wl(w)%n,g%wl(w)%m)
                        ek_nm(3) = g%eBk(g%wl(w)%n,g%wl(w)%m) 

#if __sigma__ != 2
                        thisSigma = sigma(g%wl(w)%i,g%wl(w)%j,g%wl(w)%n,g%wl(w)%m,:,:,s)
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
                        if (iSigma==1 .and. (.not. g%isMetal(g%wl(w)%i,g%wl(w)%j)) ) then        

                            kVec_stix = matMul( g%U_RTZ_to_ABb(g%wl(w)%i,g%wl(w)%j,:,:), &
                                (/ kr, g%kPhi(g%wl(w)%i), kz /) ) 

                            thisSigma = sigmaHot_maxwellian &
                                ( mSpec(s), &
                                g%ktSpec(g%wl(w)%i,g%wl(w)%j,s), &
                                g%omgc(g%wl(w)%i,g%wl(w)%j,s), &
                                g%omgp2(g%wl(w)%i,g%wl(w)%j,s), &
                                kVec_stix, g%R(g%wl(w)%i), &
                                omgrf, k0, &
                                g%k_cutoff, s, &
                                g%sinTh(g%wl(w)%i,g%wl(w)%j), &
                                g%bPol(g%wl(w)%i,g%wl(w)%j), g%bMag(g%wl(w)%i,g%wl(w)%j), &
                                g%gradPrlB(g%wl(w)%i,g%wl(w)%j), &
                                g%nuOmg(g%wl(w)%i,g%wl(w)%j) )

                        endif hotPlasma

                        coldPlasma: &
                        if (iSigma==0 .and. (.not. g%isMetal(g%wl(w)%i,g%wl(w)%j)) ) then 

                            thisSigma = sigmaCold_stix &
                                ( g%omgc(g%wl(w)%i,g%wl(w)%j,s), &
                                g%omgp2(g%wl(w)%i,g%wl(w)%j,s), omgrf, &
                                g%nuOmg(g%wl(w)%i,g%wl(w)%j) )

                        endif coldPlasma

                        ! Metal
                        ! -----

                        if (g%isMetal(g%wl(w)%i,g%wl(w)%j)) then 

                            thisSigma = 0
                            thisSigma(1,1) = metal 
                            thisSigma(2,2) = metal
                            thisSigma(3,3) = metal

                        endif

#endif
                        jVec = matMul ( thisSigma, ek_nm ) 

                        g%jAlpha(g%wl(w)%i,g%wl(w)%j,s) = g%jAlpha(g%wl(w)%i,g%wl(w)%j,s) + jVec(1) * bFn
                        g%jBeta(g%wl(w)%i,g%wl(w)%j,s) = g%jBeta(g%wl(w)%i,g%wl(w)%j,s) + jVec(2) * bFn
                        g%jB(g%wl(w)%i,g%wl(w)%j,s) = g%jB(g%wl(w)%i,g%wl(w)%j,s) + jVec(3) * bFn

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
    deallocate ( sigma )
#else
#ifdef par

    ! Switch to individual 2D arrays for the sum over processors

    do s=1,nSpec

        jAlphaTmp = g%jAlpha(:,:,s)
        jBetaTmp = g%jBeta(:,:,s)
        jBTmp = g%jB(:,:,s)

        call cGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, jAlphaTmp, g%nR, -1, -1 )
        call cGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, jBetaTmp, g%nR, -1, -1 )
        call cGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, jBTmp, g%nR, -1, -1 )

        call blacs_barrier ( iContext, 'All' ) 

        g%jAlpha(:,:,s) = jAlphaTmp
        g%jBeta(:,:,s) = jBetaTmp
        g%jB(:,:,s) = jBTmp

    enddo
#endif
#endif

end subroutine current 


subroutine jDotE ( g )

    use grid
    use aorsa2din_mod, &
    only: nSpec

    implicit none
   
    type(gridBlock), intent(inout) :: g 

    integer :: i, j, s
    complex :: eHere(3), jHere(3)

    allocate ( &
        g%jouleHeating(g%nR,g%nZ,nSpec) )

    species: &
    do s=1,nSpec

        do i=1,g%nR
            do j=1,g%nZ

                !eHere = (/ g%eAlpha(i,j), g%eBeta(i,j), g%eB(i,j) /)
                !jHere = (/ g%jAlpha(i,j,s), g%jBeta(i,j,s), g%jB(i,j,s) /)

                !g%jouleHeating(i,j,s) = &
                !    0.5 * realpart ( dot_product ( eHere, jHere ) )

                g%jouleHeating(i,j,s) = &
                    0.5 * realpart ( conjg(g%eAlpha(i,j)) * g%jAlpha(i,j,s) &
                                    + conjg(g%eBeta(i,j)) * g%jBeta(i,j,s) &
                                    + conjg(g%eB(i,j)) * g%jB(i,j,s)  )

            enddo
        enddo  

    enddo species

end subroutine jDotE

end module power
