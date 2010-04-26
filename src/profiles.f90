module profiles

use constants

implicit none

real :: omgrf, k0
real, allocatable, dimension(:) :: mSpec, qSpec, tSpec, dSpec
integer, allocatable, dimension(:) :: zSpec, amuSpec
real, allocatable, dimension(:,:,:) :: &
    omgc, omgp2, densitySpec, ktSpec
real, allocatable, dimension(:) :: &
    tLim, dLim, dAlpha, dBeta, tAlpha, tBeta

contains

    subroutine init_profiles ()

        use aorsa2din_mod, &
        only: freqcy, nSpec, zSpecIn, amuSpecIn, &
            tSpecIn, dSpecIn, tLimIn, dLimIn, &
            dAlphaIn, dBetaIn, tAlphaIn, tBetaIn

        implicit none

        omgrf = 2.0 * pi * freqcy
        k0 = omgrf / clight
        
        allocate ( &
            mSpec(nSpec), zSpec(nSpec), &
            qSpec(nSpec), amuSpec(nSpec), &
            dSpec(nSpec), tSpec(nSpec), &
            dLim(nSpec), tLim(nSpec), &
            dAlpha(nSpec), dBeta(nSpec), &
            tAlpha(nSpec), tBeta(nSpec) )
        
        zSpec       = zSpecIn(1:nSpec)
        amuSpec     = amuSpecIn(1:nSpec) 
        tSpec       = tSpecIn(1:nSpec) ! [eV]
        dSpec       = dSpecIn(1:nSpec)
        tLim        = tLimIn(1:nSpec)
        dLim        = dLimIn(1:nSpec)
        dAlpha      = dAlphaIn(1:nSpec)
        dBeta       = dBetaIn(1:nSpec)
        tAlpha      = tAlphaIn(1:nSpec)
        tBeta       = tBetaIn(1:nSpec)
        
        mSpec       = amuSpec * xmh
        mSpec(1)    = xme  
        qSpec       = zSpec * q 
 
    end subroutine init_profiles


    subroutine flat_profiles ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nSpec
        use bField

        implicit none

        integer :: i, j, s

        ! create profile
        ! --------------
       
        allocate ( &
            densitySpec ( nPtsX, nPtsY, nSpec ), & 
            ktSpec ( nPtsX, nPtsY, nSpec ) )
       
        ! ions
        ! ---- 

        do s=2,nSpec
        
            ktSpec(:,:,s)       = tSpec(s) * q 
            densitySpec(:,:,s)  = dSpec(s) 
        
        enddo

        ! electrons
        ! ---------

        ktSpec(:,:,1) = tSpec(1) * q

        do i=1,nPtsX
            do j=1,nPtsY

                densitySpec(i,j,1)  = sum ( densitySpec(i,j,2:nSpec)*zSpec(2:nSpec) )

            enddo
        enddo

        call omega_freqs ()

    end subroutine

    
    subroutine flux_profiles ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nSpec
        use bField
        use parallel

        implicit none

        integer :: s, i, j

        allocate ( &
            densitySpec ( nPtsX, nPtsY, nSpec ), & 
            ktSpec ( nPtsX, nPtsY, nSpec ) )
     
        ! Catch bad rho values

        if ( any ( rho > 1 ) .or. any ( rho < 0 ) ) then 

            write(*,*) 'profiles.f90: ERROR - bad rho values'
            stop

        endif

        ! ions
        ! ----

        do s=2,nSpec
        
            ktSpec(:,:,s)       = ( &
                tLim(s) + ( tSpec(s)-tLim(s) ) * (1d0 - rho**tBeta(s))**tAlpha(s) ) * q 
            densitySpec(:,:,s)  = &
                dLim(s) + ( dSpec(s)-dLim(s) ) * (1d0 - rho**dBeta(s))**dAlpha(s)
      
        enddo

        !   Enforce limits just in case. I would hope
        !   this is not nessecary but I have not double checked

        do i=1,nPtsX
            do j=1,nPtsY
                do s=2,nSpec
            
                    if (ktSpec(i,j,s)<tLim(s)*q) &
                        ktSpec(i,j,s) = tLim(s)*q

                    if (ktSpec(i,j,s)>tSpec(s)*q) &
                        ktSpec(i,j,s) = tSpec(s)*q

                    if (densitySpec(i,j,s)<dLim(s)) &
                        densitySpec(i,j,s) = dLim(s)

                    if (densitySpec(i,j,s)>dSpec(s)) &
                        densitySpec(i,j,s) = dSpec(s)

                enddo
            enddo
        enddo

        ! electrons
        ! ---------

        ktSpec(:,:,1) = tSpec(1) * q

        do i=1,nPtsX
            do j=1,nPtsY

                densitySpec(i,j,1)  = sum ( densitySpec(i,j,2:nSpec)*zSpec(2:nSpec) )

            enddo
        enddo


        call omega_freqs ()

    end subroutine flux_profiles

    
    subroutine omega_freqs ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nSpec
        use bField

        implicit none

        integer :: i, j

        !   calculate the cyclotron and plasma freqs
        !   ----------------------------------------
       
        allocate ( &
            omgc ( nPtsX, nPtsY, nSpec ), &
            omgp2 ( nPtsX, nPtsY, nSpec ) )
        
        do i=1,nPtsX
            do j=1,nPtsY
        
                omgc(i,j,:) = qSpec * bMod(i,j) / mSpec
                omgp2(i,j,:)    = densitySpec(i,j,:) * qSpec**2 / ( eps0 * mSpec )
        
            enddo
        enddo

    end subroutine omega_freqs

end module profiles
