module profiles

implicit none

real :: omgrf, xk0
real, allocatable, dimension(:) :: mSpec, qSpec, tSpec, dSpec
integer, allocatable, dimension(:) :: zSpec, amuSpec
real, allocatable, dimension(:,:,:) :: omgc, omgp2, densitySpec, ktSpec

contains

    subroutine init_profiles ()

        use aorsa2din_mod, &
        only: freqcy, nSpec, zSpecIn, amuSpecIn, tSpecIn, dSpecIn
        use constants

        implicit none

        omgrf = 2.0 * pi * freqcy
        xk0 = omgrf / clight
        
        allocate ( &
            mSpec(nSpec), zSpec(nSpec), &
            qSpec(nSpec), amuSpec(nSpec), &
            dSpec(nSpec), tSpec(nSpec) )
        
        zSpec       = zSpecIn(1:nSpec)
        amuSpec     = amuSpecIn(1:nSpec) 
        tSpec       = tSpecIn(1:nSpec) ! [eV]
        dSpec       = dSpecIn(1:nSpec)
        
        mSpec       = amuSpec * xmh
        mSpec(1)    = xme  
        qSpec       = zSpec * q 
 
    end subroutine init_profiles


    subroutine flat_profiles ()

        use aorsa2din_mod, &
        only: nModesX, nModesY, nSpec
        use constants
        use bField

        implicit none

        integer :: i, j, s

        !   create profile
        !   --------------
       
        allocate ( &
            densitySpec ( nModesX, nModesY, nSpec ), & 
            ktSpec ( nModesX, nModesY, nSpec ) )
        
        do s=1,nSpec
        
            ktSpec(:,:,s)       = tSpec(s) * q 
            densitySpec(:,:,s)  = dSpec(s) 
        
        enddo
 
        !   calculate the cyclotron and plasma freqs
        !   ----------------------------------------
        
        allocate ( &
            omgc ( nModesX, nModesY, nSpec ), &
            omgp2 ( nModesX, nModesY, nSpec ) )
        
        do i=1,nModesX
            do j=1,nModesY
        
                omgc(i,j,:) = qSpec * bMod(i,j) / mSpec
                omgp2(i,j,:)    = densitySpec(i,j,:) * qSpec**2 / ( eps0 * mSpec )
        
            enddo
        enddo

    end subroutine

end module profiles
