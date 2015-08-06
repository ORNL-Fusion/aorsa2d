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

contains

    function interp1d(xArr,yArr,x)

        implicit none

        real, intent(in) :: xArr(:), x 
        complex, intent(in) :: yArr(:) 

        integer :: n, iL, iR
        real :: i_
        complex :: interp1d

        n = size(xArr)

        i_ = (x-xArr(1))/(xArr(n)-xArr(1))*real(n-1.0)+1.0
        iL = floor(i_)
        iR = ceiling(i_)
        if(iL.lt.1)then 
                write(*,*) 'ERROR 1', i_, iL, iR
        endif
        if(iR.gt.n)then 
                write(*,*) 'ERROR 2', i_, iL, iR
        endif

        if(iL.eq.iR)then
            interp1d = yArr(iR) 
        else
            interp1d = yArr(iL) + (yArr(iR)-yArr(iL)) * (x - xArr(iL)) / (xArr(iR)-xArr(iL)) 
        endif

    end function interp1d


    subroutine init_JpFromFile( g )

        use read_jp_from_file, only: &
            jp_nS=>nS, jp_nR=>nR, jp_nZ=>nZ, jp_r=>r, jp_z=>z, &
            file_Jp_r=>Jp_r, file_Jp_t=>Jp_t, file_Jp_z=>Jp_z
        use grid
        use parallel, only: iAm
        use fitpack

        implicit none

        type(gridBlock), intent(inout) :: g

        real, allocatable :: Tmp2DArr_re(:,:),Tmp2DArr_im(:,:)

        integer :: islpsw, iErr, i, j, s, w
        real, allocatable :: temp(:), temp_im(:)
        real, allocatable :: zx1(:), zxm(:), zy1(:), zyn(:)
        real :: zxy11, zxym1, zxy1n, zxymn

        real, allocatable :: zx1_im(:), zxm_im(:), zy1_im(:), zyn_im(:)
        real :: zxy11_im, zxym1_im, zxy1n_im, zxymn_im

        real, allocatable :: zp_tmp(:), zp_tmp_im(:)
        real :: sigma_dp
        real, allocatable :: Jp_r_dp(:), Jp_z_dp(:)
        real :: ThisR,ThisZ
        real :: re_, im_, i_
        integer :: iL, iR
        
        if(.not.allocated(g%file_JpR))allocate ( &
            g%file_JpR(g%nR,g%nZ), &
            g%file_JpT(g%nR,g%nZ), &
            g%file_JpZ(g%nR,g%nZ) )
    
        if(Jp_nZ.gt.1)then

            sigma_dp = 1.0
            islpsw  = 255 

            allocate(temp(Jp_nZ+Jp_nZ+Jp_nR))
            allocate(zp_tmp(3*Jp_nZ*Jp_nR) )
            allocate(zx1(Jp_nZ), zxm(Jp_nZ), zy1(Jp_nR), zyn(Jp_nR))

            allocate(temp_im(Jp_nZ+Jp_nZ+Jp_nR))
            allocate(zp_tmp_im(3*Jp_nZ*Jp_nR) )
            allocate(zx1_im(Jp_nZ), zxm_im(Jp_nZ), zy1_im(Jp_nR), zyn_im(Jp_nR))

            allocate(Tmp2DArr_re(Jp_nR,Jp_nZ), Tmp2DArr_im(Jp_nR,Jp_nZ))
            allocate(Jp_r_dp(Jp_nR),Jp_z_dp(Jp_nZ))

            Jp_r_dp = Jp_r
            Jp_z_dp = Jp_z

            Tmp2DArr_re = real(real(file_Jp_r)) 
            call surf1 ( Jp_nR, Jp_nZ, Jp_r, Jp_z, Tmp2DArr_re, Jp_nR, zx1, zxm, &
                zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
                zp_tmp, temp, sigma_dp, iErr)

            Tmp2DArr_im = real(aimag(file_Jp_r)) 
            call surf1 ( Jp_nR, Jp_nZ, Jp_r_dp, Jp_z_dp, Tmp2DArr_im, Jp_nR, zx1_im, zxm_im, &
                zy1_im, zyn_im, zxy11_im, zxym1_im, zxy1n_im, zxymn_im, islpsw, &
                zp_tmp_im, temp_im, sigma_dp, iErr)

            do i=1,g%nR
                do j=1,g%nZ

                    ThisR = g%r(i)
                    ThisZ = g%z(j)

                    re_ = surf2 ( ThisR, ThisZ, Jp_nR, Jp_nZ, Jp_r, Jp_z, &
                        Tmp2DArr_re, Jp_nR, zp_tmp, sigma_dp )
                    im_ = surf2 ( ThisR, ThisZ, Jp_nR, Jp_nZ, Jp_r_dp, Jp_z_dp, &
                        Tmp2DArr_im, Jp_nR, zp_tmp_im, sigma_dp )

                    g%file_JpR(i,j) = cmplx(re_,im_)

                    write(*,*) ThisR, ThisZ, g%file_JpR(i,j), re_, im_

                enddo
            enddo

        else

            do i=1,g%nR

                    ThisR = g%r(i)

                    g%file_JpR(i,1) = interp1d(Jp_r,file_Jp_r(:,1),ThisR)
                    g%file_JpT(i,1) = interp1d(Jp_r,file_Jp_t(:,1),ThisR)
                    g%file_JpZ(i,1) = interp1d(Jp_r,file_Jp_z(:,1),ThisR)
      
            enddo

        endif

    end subroutine init_JpFromFile


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
                write(*,*) '    This is possibly due to a bad interpolation.'
                write(*,*) '    Try creating the ar2Input file with higher res.'
                write(*,*) '    Really need to fix this crap.'
                stop
        endif

        if(count(g%DensitySpec<=0)>0)then
                write(*,*) 'ERROR: -ve density'
                write(*,*) '    This is possibly due to a bad interpolation.'
                write(*,*) '    Try creating the ar2Input file with higher res.'
                write(*,*) '    Really need to fix this crap.'
                stop
        endif

        if(count(g%nuOmg<0)>0)then
                where(g%nuOmg<0.0) g%nuOmg = 0.0
        endif

        call omega_freqs ( g )

    end subroutine init_AR2_Profiles

    
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
       
        allocate ( &
            g%omgc ( size(g%pt), nSpec ), &
            g%omgp2 ( size(g%pt), nSpec ) )
        
        do w=1,size(g%pt)
            i=g%pt(w)%i
            j=g%pt(w)%j
        
                g%omgc(w,:) = qSpec * g%bMag(w) / mSpec
                g%omgp2(w,:)    = g%densitySpec(w,:) * qSpec**2 / ( eps0 * mSpec )
        
        enddo

    end subroutine omega_freqs

end module profiles
