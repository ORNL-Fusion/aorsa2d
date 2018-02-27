module power

contains

subroutine current ( g, rhs )

    use grid
    use read_data
    use aorsaNamelist, &
        only: nSpec, iSigma, fracOfModesInSolution, ZeroJp, &
        ZeroJp_rMin, ZeroJp_rMax, ZeroJp_zMin, ZeroJp_zMax, &
        useJpFromFile, ColdIons, nRAll, nZAll, kz_1D
    use sigma
    use parallel
    use profiles, &
        only: omgrf, mSpec
    use constants
    use read_jp_from_file, only: &
        file_nS=>nS, file_nR=>nR, file_nZ=>nZ, file_r=>r, file_z=>z, &
        file_Jp_r=>Jp_r, file_Jp_t=>Jp_t, file_Jp_z=>Jp_z!, &
        !file_Jp_r_s=>Jp_r_s, file_Jp_t_s=>Jp_t_s, file_Jp_z_s=>Jp_z_s

    implicit none
   
    type(gridBlock), intent(inout) :: g 
    integer, intent(in) :: rhs

    complex :: ek_nm(3), jVec(3), thisSigma(3,3)
    complex :: bFn

    integer :: s, w
    real :: kr, kz, kVec_stix(3)
    complex, allocatable :: jAlphaTmp(:,:), jBetaTmp(:,:), jBTmp(:,:)

    type(spatialSigmaInput_cold) :: sigmaIn_cold
#if __noU__==1
    real :: R_(3,3)
#endif
    real :: R, z
    complex(kind=dbl) :: k0

    if (.not.allocated(g%jAlpha)) allocate ( &
        g%jAlpha(g%nR,g%nZ,nSpec), &
        g%jBeta(g%nR,g%nZ,nSpec), &
        g%jB(g%nR,g%nZ,nSpec) )

#if WRITE_SIGMA_TO_OUTPUT > 0 

     if (.not.allocated(g%sig11)) allocate ( &
        g%sig11(g%nR,g%nMin:g%nMax,nSpec), &
        g%sig12(g%nR,g%nMin:g%nMax,nSpec), &
        g%sig13(g%nR,g%nMin:g%nMax,nSpec), &
        g%sig21(g%nR,g%nMin:g%nMax,nSpec), &
        g%sig22(g%nR,g%nMin:g%nMax,nSpec), &
        g%sig23(g%nR,g%nMin:g%nMax,nSpec), &
        g%sig31(g%nR,g%nMin:g%nMax,nSpec), &
        g%sig32(g%nR,g%nMin:g%nMax,nSpec), &
        g%sig33(g%nR,g%nMin:g%nMax,nSpec) )

    if (.not.allocated(g%kr)) allocate( &
        g%kr(g%nMin:g%nMax,g%mMin:g%mMax), g%kz(g%nMin:g%nMax,g%mMin:g%mMax) )
#endif

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

        workList: &
        do w=1,size(g%wl)

            k0  = g%k0(g%wl(w)%iPt)

            !ReplaceWithJpFromFile: &
            !if(useJpFromFile)then

            !else

                twoThirdsRule: &
                if(g%wl(w)%m >= g%mMin*fracOfModesInSolution .and. g%wl(w)%m <= g%mMax*fracOfModesInSolution &
                    .and. g%wl(w)%n >= g%nMin*fracOfModesInSolution .and. g%wl(w)%n <= g%nMax*fracOfModesInSolution ) then

                
                        bFn = g%xx(g%wl(w)%n,g%wl(w)%i) * g%yy(g%wl(w)%m,g%wl(w)%j)

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

                            if(g%nZ.gt.1)then
                                kz = g%wl(w)%m * g%normFacZ
                            else
                                kz = kz_1D
                            endif

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
                            thisSigma(1,1) = 0!metal 
                            thisSigma(2,2) = 0!metal
                            thisSigma(3,3) = 0!metal

                        endif

                        ! This will zero the plasma current in regions
                        ! where we want to specify it manually in the
                        ! RHS.

                        if(ZeroJp)then 
                            R = g%R(g%wl(w)%i)
                            z = g%z(g%wl(w)%j)
                            if(R>=ZeroJp_rMin &
                                    .and.R<=ZeroJp_rMax &
                                    .and.z>=ZeroJp_zMin &
                                    .and.z<=ZeroJp_zMax) then

                                thisSigma = 0

                            endif
                        endif

                        ek_nm(1) = g%eAlphak(g%wl(w)%n,g%wl(w)%m,rhs)
                        ek_nm(2) = g%eBetak(g%wl(w)%n,g%wl(w)%m,rhs)
                        ek_nm(3) = g%eBk(g%wl(w)%n,g%wl(w)%m,rhs) 

#if PRINT_SIGMA_ABP>=1
                        if(g%wl(w)%n.eq.g%nMin.and.s.eq.2)then

                            write(*,*) 
                            write(*,*) g%wl(w)%n, g%wl(w)%n, g%r(g%wl(w)%i)
                            write(*,*) thisSigma(1,1), thisSigma(1,2), thisSigma(1,3)
                            write(*,*) thisSigma(2,1), thisSigma(2,2), thisSigma(2,3)
                            write(*,*) thisSigma(3,1), thisSigma(3,2), thisSigma(3,3)
                        endif
#endif

#if WRITE_SIGMA_TO_OUTPUT > 0

                        g%kr(g%wl(w)%n,g%wl(w)%m) = kr
                        g%kz(g%wl(w)%n,g%wl(w)%m) = kz

                        g%sig11(g%wl(w)%i,g%wl(w)%n,s) = thisSigma(1,1)
                        g%sig12(g%wl(w)%i,g%wl(w)%n,s) = thisSigma(1,2)
                        g%sig13(g%wl(w)%i,g%wl(w)%n,s) = thisSigma(1,3)
                        g%sig21(g%wl(w)%i,g%wl(w)%n,s) = thisSigma(2,1)
                        g%sig22(g%wl(w)%i,g%wl(w)%n,s) = thisSigma(2,2)
                        g%sig23(g%wl(w)%i,g%wl(w)%n,s) = thisSigma(2,3)
                        g%sig31(g%wl(w)%i,g%wl(w)%n,s) = thisSigma(3,1)
                        g%sig32(g%wl(w)%i,g%wl(w)%n,s) = thisSigma(3,2)
                        g%sig33(g%wl(w)%i,g%wl(w)%n,s) = thisSigma(3,3)
#endif

                        thisSigma = transpose(thisSigma) ! this gives consistency with KJ and RS. 
                        jVec = matMul ( thisSigma, ek_nm ) 

                        g%jAlpha(g%wl(w)%i,g%wl(w)%j,s) = g%jAlpha(g%wl(w)%i,g%wl(w)%j,s) + jVec(1) * bFn
                        g%jBeta(g%wl(w)%i,g%wl(w)%j,s) = g%jBeta(g%wl(w)%i,g%wl(w)%j,s) + jVec(2) * bFn
                        g%jB(g%wl(w)%i,g%wl(w)%j,s) = g%jB(g%wl(w)%i,g%wl(w)%j,s) + jVec(3) * bFn

                    endif twoThirdsRule

            !endif ReplaceWithJpFromFile

            ! Where E=0 at the boundary, also set J=0 to avoid
            ! tiny E values combining with sigma to give anomolous
            ! J values in the output.

            if (g%label(g%wl(w)%iPt)==888) then

                g%jAlpha(g%wl(w)%i,g%wl(w)%j,s) = 0
                g%jBeta(g%wl(w)%i,g%wl(w)%j,s) = 0
                g%jB(g%wl(w)%i,g%wl(w)%j,s) = 0

            endif

        enddo workList

    enddo species

#ifdef par

    ! Switch to individual 2D arrays for the sum over processors

    do s=1,nSpec

        jAlphaTmp = g%jAlpha(:,:,s)
        jBetaTmp = g%jBeta(:,:,s)
        jBTmp = g%jB(:,:,s)

        call cGSUM2D ( ICTXT, 'All', ' ', g%nR, g%nZ, jAlphaTmp, g%nR, -1, -1 )
        call cGSUM2D ( ICTXT, 'All', ' ', g%nR, g%nZ, jBetaTmp, g%nR, -1, -1 )
        call cGSUM2D ( ICTXT, 'All', ' ', g%nR, g%nZ, jBTmp, g%nR, -1, -1 )

        call blacs_barrier ( ICTXT, 'All' ) 

        g%jAlpha(:,:,s) = jAlphaTmp
        g%jBeta(:,:,s) = jBetaTmp
        g%jB(:,:,s) = jBTmp

    enddo

        deallocate(jAlphaTmp,jBetaTmp,jBTmp)

#endif
#endif

end subroutine current 


subroutine jDotE ( g, rhs)

    use grid
    use aorsaNamelist, &
        only: nSpec

    implicit none
   
    type(gridBlock), intent(inout) :: g 
    integer, intent(in) :: rhs

    integer :: i, j, s

    if(.not.allocated(g%jouleHeating))allocate ( g%jouleHeating(g%nR,g%nZ,nSpec) )

    species: &
    do s=1,nSpec

        do i=1,g%nR
            do j=1,g%nZ

                g%jouleHeating(i,j,s) = &
                    0.5 * real(real( conjg(g%eAlpha(i,j)) * g%jAlpha(i,j,s) &
                                    + conjg(g%eBeta(i,j)) * g%jBeta(i,j,s) &
                                    + conjg(g%eB(i,j)) * g%jB(i,j,s)  ))

            enddo
        enddo  

    enddo species

end subroutine jDotE

end module power
