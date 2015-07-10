module profiles

use constants

implicit none

real(kind=dbl) :: omgrf, k0
real(kind=dbl), allocatable, dimension(:) :: mSpec, qSpec, tSpec, dSpec
real(kind=dbl), allocatable, dimension(:) :: zSpec, amuSpec
real(kind=dbl), allocatable, dimension(:,:,:) :: &
    omgc, omgp2
real, allocatable, dimension(:) :: &
    tLim, dLim, dAlpha, dBeta, tAlpha, tBeta
!real, allocatable :: nuOmg2D(:,:)

contains

    subroutine init_AR2_Profiles( g )

        use AR2Input, only: &
            ar2_nS=>nS,ar2_amu=>amu,ar2_AtomicZ=>AtomicZ, &
            ar2_nR=>nR, ar2_nZ=>nZ, ar2_r=>r, ar2_z=>z, &
            Density_m3, Temp_eV, nuOmg, jant_r, jant_t, jant_z
        use grid
        use aorsaNameList, only: &
            nSpec,freqcy
        use parallel, only: iAm

        use fitpack

        implicit none

        type(gridBlock), intent(inout) :: g

        real, allocatable :: Tmp2DArr(:,:)

        integer :: islpsw, iErr, i, j, s, w
        real, allocatable :: temp(:)
        real, allocatable :: zx1(:), zxm(:), zy1(:), zyn(:)
        real :: zxy11, zxym1, zxy1n, zxymn

        real, allocatable :: zp_tmp(:)
        real :: sigma_dp
        real, allocatable :: ar2_r_dp(:), ar2_z_dp(:)
        real :: ThisR,ThisZ

        sigma_dp = 1.0
        islpsw  = 255 

        allocate(temp(ar2_nZ+ar2_nZ+ar2_nR))
        allocate(zp_tmp(3*ar2_nZ*ar2_nR) )
        allocate(zx1(ar2_nZ), zxm(ar2_nZ), zy1(ar2_nR), zyn(ar2_nR))
        allocate(Tmp2DArr(ar2_nR,ar2_nZ))
        allocate(ar2_r_dp(ar2_nR),ar2_z_dp(ar2_nZ))

        ar2_r_dp = ar2_r
        ar2_z_dp = ar2_z

        ! At present these two lines are duplicates of the
        ! init_profiles() routine below. This is not nice
        ! and I should change it to get done somewhere else.

        omgrf = 2.0 * pi * freqcy
        k0 = omgrf / clight
 
        ! Overwrite some nameList varibles with those from 
        ! the ar2Input file.

        nSpec = ar2_nS

        allocate ( &
            mSpec(nSpec), zSpec(nSpec), &
            qSpec(nSpec), amuSpec(nSpec) )

        amuSpec = ar2_amu
        zSpec = ar2_AtomicZ

        if(iAm==0)then
            write(*,*) 'Overwriting namelist inputs for nSpec and profiles'
            write(*,*) '    nSpec: ', nSpec
            write(*,*) '    zSpec: ', zSpec
            write(*,*) '    amuSpec: ', amuSpec
        endif

        mSpec       = amuSpec * xmh
        mSpec(1)    = xme  
        qSpec       = zSpec * q 

        ! Interpolate the ar2Input data to the grid.

        !allocate ( &
        !    g%densitySpec ( g%nR, g%nZ, nSpec ), & 
        !    g%ktSpec ( g%nR, g%nZ, nSpec ), &
        !    g%nuOmg(g%nR,g%nZ) )
        allocate ( &
            g%densitySpec (size(g%pt), nSpec ), & 
            g%ktSpec (size(g%pt), nSpec ), &
            g%nuOmg(size(g%pt), nSpec) )
 
        do s=1,nSpec

            ! Density interpolation

            Tmp2DArr = Density_m3(:,:,s) 

            call surf1 ( ar2_nR, ar2_nZ, ar2_r_dp, ar2_z_dp, Tmp2DArr, ar2_nR, zx1, zxm, &
                zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
                zp_tmp, temp, sigma_dp, iErr)

            do w=1,size(g%pt)
                i=g%pt(w)%i
                j=g%pt(w)%j

                ThisR = g%r(i)
                ThisZ = g%z(j)

                g%DensitySpec(w,s) = surf2 ( ThisR, ThisZ, ar2_nR, ar2_nZ, ar2_r_dp, ar2_z_dp, &
                    Tmp2DArr, ar2_nR, zp_tmp, sigma_dp )

            enddo

            ! Temp interpolation

            Tmp2DArr = Temp_eV(:,:,s) 

            call surf1 ( ar2_nR, ar2_nZ, ar2_r_dp, ar2_z_dp, Tmp2DArr, ar2_nR, zx1, zxm, &
                zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
                zp_tmp, temp, sigma_dp, iErr)

            do w=1,size(g%pt)
                i=g%pt(w)%i
                j=g%pt(w)%j

                ThisR = g%r(i)
                ThisZ = g%z(j)

                g%kTSpec(w,s) = surf2 ( ThisR, ThisZ, ar2_nR, ar2_nZ, ar2_r_dp, ar2_z_dp, &
                    Tmp2DArr, ar2_nR, zp_tmp, sigma_dp ) * q
            
            enddo

            ! nuOmg (artificial absorption) interpolation

            Tmp2DArr = nuOmg(:,:,s) 

            call surf1 ( ar2_nR, ar2_nZ, ar2_r_dp, ar2_z_dp, Tmp2DArr, ar2_nR, zx1, zxm, &
                zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
                zp_tmp, temp, sigma_dp, iErr)

            do w=1,size(g%pt)
                i=g%pt(w)%i
                j=g%pt(w)%j

                ThisR = g%r(i)
                ThisZ = g%z(j)

                g%nuOmg(w,s) = surf2 ( ThisR, ThisZ, ar2_nR, ar2_nZ, ar2_r_dp, ar2_z_dp, &
                    Tmp2DArr, ar2_nR, zp_tmp, sigma_dp )
            enddo
        enddo

        ! JAnt interpolation

        if(.not.allocated(g%jR))allocate ( &
            g%jR(g%nR,g%nZ), &
            g%jT(g%nR,g%nZ), &
            g%jZ(g%nR,g%nZ) )

        Tmp2DArr = jant_r(:,:) 
        call surf1 ( ar2_nR, ar2_nZ, ar2_r_dp, ar2_z_dp, Tmp2DArr, ar2_nR, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            zp_tmp, temp, sigma_dp, iErr)

        do i=1,g%nR
            do j=1,g%nZ
                ThisR = g%r(i)
                ThisZ = g%z(j)
                g%jR(i,j) = surf2 ( ThisR, ThisZ, ar2_nR, ar2_nZ, ar2_r_dp, ar2_z_dp, &
                    Tmp2DArr, ar2_nR, zp_tmp, sigma_dp )
            enddo
        enddo

        Tmp2DArr = jant_t(:,:) 
        call surf1 ( ar2_nR, ar2_nZ, ar2_r_dp, ar2_z_dp, Tmp2DArr, ar2_nR, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            zp_tmp, temp, sigma_dp, iErr)

        do i=1,g%nR
            do j=1,g%nZ
                ThisR = g%r(i)
                ThisZ = g%z(j)
                g%jT(i,j) = surf2 ( ThisR, ThisZ, ar2_nR, ar2_nZ, ar2_r_dp, ar2_z_dp, &
                    Tmp2DArr, ar2_nR, zp_tmp, sigma_dp )
            enddo
        enddo

        Tmp2DArr = jant_z(:,:) 
        call surf1 ( ar2_nR, ar2_nZ, ar2_r_dp, ar2_z_dp, Tmp2DArr, ar2_nR, zx1, zxm, &
            zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
            zp_tmp, temp, sigma_dp, iErr)

         do i=1,g%nR
            do j=1,g%nZ
                ThisR = g%r(i)
                ThisZ = g%z(j)
                g%jZ(i,j) = surf2 ( ThisR, ThisZ, ar2_nR, ar2_nZ, ar2_r_dp, ar2_z_dp, &
                    Tmp2DArr, ar2_nR, zp_tmp, sigma_dp )
            enddo
        enddo

        ! Santiy checking

        if(count(g%kTSpec<=0)>0)then
                write(*,*) 'ERROR: -ve temp'
                stop
        endif

        if(count(g%DensitySpec<=0)>0)then
                write(*,*) 'ERROR: -ve density'
                write(*,*) '    This is possibly due to a bad interpolation.'
                write(*,*) '    Try creating the ar2Input file with higher res.'
                stop
        endif

        if(count(g%nuOmg<0)>0)then
                where(g%nuOmg<0.0) g%nuOmg = 0.0
        endif

        call omega_freqs ( g )

    end subroutine init_AR2_Profiles

    !subroutine init_profiles ( nSpec )

    !    use fitpack_dp 

    !    use aorsaNamelist, &
    !    only: freqcy, zSpecIn, amuSpecIn, &
    !        tSpecIn, dSpecIn, tLimIn, dLimIn, &
    !        dAlphaIn, dBetaIn, tAlphaIn, tBetaIn
    !    use grid

    !    implicit none

    !    integer, intent(in) :: nSpec

    !    omgrf = 2.0 * pi * freqcy
    !    k0 = omgrf / clight
    !    
    !    allocate ( &
    !        mSpec(nSpec), zSpec(nSpec), &
    !        qSpec(nSpec), amuSpec(nSpec), &
    !        dSpec(nSpec), tSpec(nSpec), &
    !        dLim(nSpec), tLim(nSpec), &
    !        dAlpha(nSpec), dBeta(nSpec), &
    !        tAlpha(nSpec), tBeta(nSpec) )
    !    
    !    zSpec       = zSpecIn(1:nSpec)
    !    amuSpec     = amuSpecIn(1:nSpec) 
    !    tSpec       = tSpecIn(1:nSpec) ! [eV]
    !    dSpec       = dSpecIn(1:nSpec)
    !    tLim        = tLimIn(1:nSpec)
    !    dLim        = dLimIn(1:nSpec)
    !    dAlpha      = dAlphaIn(1:nSpec)
    !    dBeta       = dBetaIn(1:nSpec)
    !    tAlpha      = tAlphaIn(1:nSpec)
    !    tBeta       = tBetaIn(1:nSpec)
    !    
    !    mSpec       = amuSpec * xmh
    !    mSpec(1)    = xme  
    !    qSpec       = zSpec * q 
 

    !end subroutine init_profiles


    !subroutine flat_profiles ( g, parabolic )

    !    use aorsaNamelist, &
    !    only: nSpec, xNuOmg, r0, a
    !    use bField
    !    use grid

    !    implicit none

    !    logical, intent(in), optional :: parabolic
    !    type(gridBlock), intent(inout) :: g

    !    integer :: i, j, s
    !    real :: scaleFacT, scaleFacD

    !    ! create profile
    !    ! --------------
    !  
    !    allocate ( &
    !        g%densitySpec ( g%nR, g%nZ, nSpec ), & 
    !        g%ktSpec ( g%nR, g%nZ, nSpec ), &
    !        g%nuOmg(g%nR,g%nZ) )
    !   
    !    ! ions
    !    ! ---- 

    !    do s=2,nSpec
    !    

    !        ! Flat profiles
    !        ! -------------

    !        g%ktSpec(:,:,s)       = tSpec(s) * q 
    !        g%densitySpec(:,:,s)  = dSpec(s) 


    !        ! Parabolic profiles
    !        ! ------------------

    !        if(present(parabolic))then
    !            if(parabolic)then

    !                !if(r0+a>g%

    !                do i=1,g%nR
    !                    scaleFacT = (tSpec(s)-tLim(s))/a**2
    !                    scaleFacD = (dSpec(s)-dLim(s))/a**2
    !                    g%ktSpec(i,:,s)       = ( tSpec(s)-scaleFacT*(g%R(i)-r0)**2 ) * q  
    !                    g%densitySpec(i,:,s)  = ( dSpec(s)-scaleFacD*(g%R(i)-r0)**2 ) 
    !                enddo

    !                if(any(g%ktSpec(:,:,s)<0))then
    !                    write(*,*) 'ERROR - src/profiles.f90: ktSpec < 0'
    !                    stop
    !                endif

    !                if(any(g%densitySpec(:,:,s)<=0))then
    !                    write(*,*) 'ERROR - src/profiles.f90: densitySpec <= 0'
    !                    stop
    !                endif

    !            endif
    !        endif
    !    
    !    enddo

    !    ! electrons
    !    ! ---------

    !    g%ktSpec(:,:,1) = tSpec(1) * q

    !    if(present(parabolic))then
    !        if(parabolic)then
    !            scaleFacT = (tSpec(1)-tLim(1))/a**2
    !            do i=1,g%nR
    !                g%ktSpec(i,:,1) = ( tSpec(1)-scaleFacT*(g%R(i)-r0)**2 ) * q  
    !            enddo
    !        endif
    !    endif


    !    do i=1,g%nR
    !        do j=1,g%nZ

    !            g%densitySpec(i,j,1)  = sum ( g%densitySpec(i,j,2:nSpec)*zSpec(2:nSpec) )

    !        enddo
    !    enddo

    !    g%nuOmg = xNuOmg

    !    call omega_freqs ( g )


    !end subroutine


    !subroutine circular_profiles ( g )

    !    use aorsaNamelist, &
    !    only: nSpec, xNuOmgOutside, a, r0, xNuOmg
    !    use bField
    !    use parallel
    !    use grid

    !    implicit none

    !    type(gridBlock), intent(inout) :: g

    !    integer :: s, i, j
    !    real, allocatable :: distance(:,:)

    !    allocate ( &
    !        g%densitySpec ( g%nR, g%nZ, nSpec ), & 
    !        g%ktSpec ( g%nR, g%nZ, nSpec ), &
    !        distance(g%nR,g%nZ), &
    !        g%nuOmg(g%nR,g%nZ) )
    ! 
    !    do i=1,g%nR
    !        do j=1,g%nZ

    !            distance(i,j)   = sqrt( (g%R(i)-r0)**2 + (g%Z(j)-0.0)**2  ) / a

    !        enddo
    !    enddo

    !    do s=1,nSpec 
    !            
    !        g%ktSpec(:,:,s) = ( exp ( - ( distance )**tBeta(s) / ( 2 * tAlpha(s)**2 ) ) &
    !            * (tSpec(s)-tLim(s)) + tLim(s) ) * q
    !    enddo

    !    do s=2,nSpec 
    !        g%densitySpec(:,:,s)  = &
    !            dLim(s) + ( dSpec(s)-dLim(s) ) &
    !            * exp ( - ( distance )**dBeta(s) / ( 2 * dAlpha(s)**2 ) )
    !    enddo

    !    do i=1,g%nR
    !        do j=1,g%nZ

    !            g%densitySpec(i,j,1)  = sum ( g%densitySpec(i,j,2:nSpec)*zSpec(2:nSpec) )

    !        enddo
    !    enddo

    !    call omega_freqs ( g )

    !    g%nuOmg = xNuOmg
    !    where ( distance/a>=0.99)
    !        g%nuOmg = xNuOmgOutside
    !    endwhere

    !end subroutine circular_profiles

    !
    !subroutine flux_profiles ( g )

    !    use aorsaNamelist, &
    !    only: nSpec, xNuOmgOutside, xNuOmg
    !    use bField
    !    use parallel
    !    use grid

    !    implicit none

    !    type(gridBlock), intent(inout) :: g

    !    integer :: s, i, j, sWidthX, sWidthY, nS
    !    real, allocatable :: smoothingTmp(:,:), smoothingDen(:,:)

    !    allocate ( &
    !        g%densitySpec ( g%nR, g%nZ, nSpec ), & 
    !        g%ktSpec ( g%nR, g%nZ, nSpec ), &
    !        g%nuOmg(g%nR,g%nZ) )

    !    ! Catch bad rho values

    !    if ( any ( g%rho > 1 ) .or. any ( g%rho < 0 ) ) then 

    !        write(*,*) 'profiles.f90: ERROR - bad rho values'
    !        stop

    !    endif


    !    ! ions and e for temp since they can 
    !    ! be specified independantly
    !    ! ----------------------------------

    !    do s=1,nSpec 
    !        g%ktSpec(:,:,s) = ( &
    !            tLim(s) + ( tSpec(s)-tLim(s) ) * (1d0 - g%rho**tBeta(s))**tAlpha(s) ) * q 
    !    enddo


    !    ! ions only for density, e density calculated
    !    ! for quasi neutrality
    !    ! -------------------------------------------

    !    do s=2,nSpec 
    !        g%densitySpec(:,:,s)  = &
    !            dLim(s) + ( dSpec(s)-dLim(s) ) * (1d0 - g%rho**dBeta(s))**dAlpha(s)
    !    enddo

    !    !   Enforce limits just in case. I would hope
    !    !   this is not nessecary but I have not double checked

    !    do i=1,g%nR
    !        do j=1,g%nZ
    !            do s=2,nSpec
    !        
    !                if (g%ktSpec(i,j,s)<tLim(s)*q) &
    !                    g%ktSpec(i,j,s) = tLim(s)*q

    !                if (g%ktSpec(i,j,s)>tSpec(s)*q) &
    !                    g%ktSpec(i,j,s) = tSpec(s)*q

    !                if (g%densitySpec(i,j,s)<dLim(s)) &
    !                    g%densitySpec(i,j,s) = dLim(s)

    !                if (g%densitySpec(i,j,s)>dSpec(s)) &
    !                    g%densitySpec(i,j,s) = dSpec(s)

    !            enddo
    !        enddo
    !    enddo

    !    ! electrons density for quasi neutrality
    !    ! --------------------------------------

    !    do i=1,g%nR
    !        do j=1,g%nZ

    !            g%densitySpec(i,j,1)  = sum ( g%densitySpec(i,j,2:nSpec)*zSpec(2:nSpec) )

    !        enddo
    !    enddo


    !    ! smooth densities
    !    ! ----------------

    !    sWidthX = g%nR/20 
    !    sWidthY = g%nZ/20
    !    allocate ( smoothingTmp(1-sWidthX:g%nR+sWidthX,1-sWidthY:g%nZ+sWidthY) )
    !    allocate ( smoothingDen(1-sWidthX:g%nR+sWidthX,1-sWidthY:g%nZ+sWidthY) )

    !    do s=1,nSpec
    !        smoothingTmp = minVal(g%ktSpec(:,:,s))
    !        smoothingTmp(1:g%nR,1:g%nZ) = g%ktSpec(:,:,s)
    !        smoothingDen = minVal(g%densitySpec(:,:,s))
    !        smoothingDen(1:g%nR,1:g%nZ) = g%densitySpec(:,:,s)

    !        do nS = 1, 5
    !            do i=1,g%nR
    !                do j=1,g%nZ

    !                    smoothingTmp(i,j)  = &
    !                        sum ( smoothingTmp(i-sWidthX:i+sWidthX,j-sWidthY:j+sWidthY) ) &
    !                        / ((sWidthX*2+1) * (sWidthY*2+1))
    !                    smoothingDen(i,j)  = &
    !                        sum ( smoothingDen(i-sWidthX:i+sWidthX,j-sWidthY:j+sWidthY) ) &
    !                        / ((sWidthX*2+1) * (sWidthY*2+1))

    !                enddo
    !            enddo
    !        enddo
    !        g%densitySpec(:,:,s) = smoothingDen(1:g%nR,1:g%nZ)
    !        g%ktSpec(:,:,s) = smoothingTmp(1:g%nR,1:g%nZ)
    !    enddo

    !    deallocate ( smoothingTmp, smoothingDen )

    !    call omega_freqs ( g )


    !    ! create a 2D map of collisional damping parameter
    !    ! for application outside the last closed flux surface
    !    ! ----------------------------------------------------
    !    
    !    g%nuOmg = xNuOmg
    !    where ( g%rho>=0.99)
    !        g%nuOmg = xNuOmgOutside
    !    endwhere

    !end subroutine flux_profiles

    
    subroutine omega_freqs ( g )

        use aorsaNamelist, &
        only: nSpec
        use bField
        use grid

        implicit none

        type(gridBlock), intent(inout) :: g

        integer :: i, j, w

        !   calculate the cyclotron and plasma freqs
        !   ----------------------------------------
       
        !allocate ( &
        !    g%omgc ( g%nR, g%nZ, nSpec ), &
        !    g%omgp2 ( g%nR, g%nZ, nSpec ) )
        allocate ( &
            g%omgc ( size(g%pt), nSpec ), &
            g%omgp2 ( size(g%pt), nSpec ) )
        
        do w=1,size(g%pt)
            i=g%pt(w)%i
            j=g%pt(w)%j
        !do i=1,g%nR
        !    do j=1,g%nZ
        
                g%omgc(w,:) = qSpec * g%bMag(w) / mSpec
                g%omgp2(w,:)    = g%densitySpec(w,:) * qSpec**2 / ( eps0 * mSpec )
        
        !    enddo
        !enddo
        enddo

    end subroutine omega_freqs

end module profiles
