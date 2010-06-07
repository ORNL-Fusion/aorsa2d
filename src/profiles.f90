module profiles

use constants

implicit none

real(kind=dbl) :: omgrf, k0
real, allocatable, dimension(:) :: mSpec, qSpec, tSpec, dSpec
integer, allocatable, dimension(:) :: zSpec, amuSpec
real(kind=dbl), allocatable, dimension(:,:,:) :: &
    omgc, omgp2, densitySpec, ktSpec
real, allocatable, dimension(:) :: &
    tLim, dLim, dAlpha, dBeta, tAlpha, tBeta
real, allocatable :: nuOmg2D(:,:)

contains

    subroutine init_profiles ()

        use aorsa2din_mod, &
        only: freqcy, nSpec, zSpecIn, amuSpecIn, &
            tSpecIn, dSpecIn, tLimIn, dLimIn, &
            dAlphaIn, dBetaIn, tAlphaIn, tBetaIn, &
            nPtsX, nPtsY, xNuOmg

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
 
        allocate ( nuOmg2D(nPtsX,nPtsY) )

        nuOmg2D = xNuOmg


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


    subroutine circular_profiles ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nSpec, xNuOmgOutside, a, r0
        use bField
        use parallel
        use grid, &
        only: capR, y

        implicit none

        integer :: s, i, j
        real, allocatable :: distance(:,:)

        allocate ( &
            densitySpec ( nPtsX, nPtsY, nSpec ), & 
            ktSpec ( nPtsX, nPtsY, nSpec ), &
            distance(nPtsX,nPtsY) )
     
        do i=1,nPtsX
            do j=1,nPtsY

                distance(i,j)   = sqrt( (capR(i)-r0)**2 + (y(j)-0.0)**2  ) / a

            enddo
        enddo

        do s=1,nSpec 
                
            ktSpec(:,:,s) = ( exp ( - ( distance )**tBeta(s) / ( 2 * tAlpha(s)**2 ) ) &
                * (tSpec(s)-tLim(s)) + tLim(s) ) * q
        enddo

        do s=2,nSpec 
            densitySpec(:,:,s)  = &
                dLim(s) + ( dSpec(s)-dLim(s) ) &
                * exp ( - ( distance )**dBeta(s) / ( 2 * dAlpha(s)**2 ) )
        enddo

        do i=1,nPtsX
            do j=1,nPtsY

                densitySpec(i,j,1)  = sum ( densitySpec(i,j,2:nSpec)*zSpec(2:nSpec) )

            enddo
        enddo

        call omega_freqs ()

        where ( distance/a>=0.99)
            nuOmg2D = xNuOmgOutside
        endwhere

    end subroutine circular_profiles

    
    subroutine flux_profiles ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nSpec, xNuOmgOutside
        use bField
        use parallel

        implicit none

        integer :: s, i, j, sWidthX, sWidthY, nS
        real, allocatable :: smoothingTmp(:,:), smoothingDen(:,:)

        allocate ( &
            densitySpec ( nPtsX, nPtsY, nSpec ), & 
            ktSpec ( nPtsX, nPtsY, nSpec ) )
     
        ! Catch bad rho values

        if ( any ( rho > 1 ) .or. any ( rho < 0 ) ) then 

            write(*,*) 'profiles.f90: ERROR - bad rho values'
            stop

        endif


        ! ions and e for temp since they can 
        ! be specified independantly
        ! ----------------------------------

        do s=1,nSpec 
            ktSpec(:,:,s) = ( &
                tLim(s) + ( tSpec(s)-tLim(s) ) * (1d0 - rho**tBeta(s))**tAlpha(s) ) * q 
        enddo


        ! ions only for density, e density calculated
        ! for quasi neutrality
        ! -------------------------------------------

        do s=2,nSpec 
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

        ! electrons density for quasi neutrality
        ! --------------------------------------

        do i=1,nPtsX
            do j=1,nPtsY

                densitySpec(i,j,1)  = sum ( densitySpec(i,j,2:nSpec)*zSpec(2:nSpec) )

            enddo
        enddo


        ! smooth densities
        ! ----------------

        sWidthX = nPtsX/20 
        sWidthY = nPtsY/20
        allocate ( smoothingTmp(1-sWidthX:nPtsX+sWidthX,1-sWidthY:nPtsY+sWidthY) )
        allocate ( smoothingDen(1-sWidthX:nPtsX+sWidthX,1-sWidthY:nPtsY+sWidthY) )

        do s=1,nSpec
            smoothingTmp = minVal(ktSpec(:,:,s))
            smoothingTmp(1:nPtsX,1:nPtsY) = ktSpec(:,:,s)
            smoothingDen = minVal(densitySpec(:,:,s))
            smoothingDen(1:nPtsX,1:nPtsY) = densitySpec(:,:,s)

            do nS = 1, 5
                do i=1,nPtsX
                    do j=1,nPtsY

                        smoothingTmp(i,j)  = &
                            sum ( smoothingTmp(i-sWidthX:i+sWidthX,j-sWidthY:j+sWidthY) ) &
                            / ((sWidthX*2+1) * (sWidthY*2+1))
                        smoothingDen(i,j)  = &
                            sum ( smoothingDen(i-sWidthX:i+sWidthX,j-sWidthY:j+sWidthY) ) &
                            / ((sWidthX*2+1) * (sWidthY*2+1))

                    enddo
                enddo
            enddo
            densitySpec(:,:,s) = smoothingDen(1:nPtsX,1:nPtsY)
            ktSpec(:,:,s) = smoothingTmp(1:nPtsX,1:nPtsY)
        enddo

        deallocate ( smoothingTmp, smoothingDen )

        call omega_freqs ()


        ! create a 2D map of collisional damping parameter
        ! for application outside the last closed flux surface
        ! ----------------------------------------------------

        where ( rho>=0.99)
            nuOmg2D = xNuOmgOutside
        endwhere

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
