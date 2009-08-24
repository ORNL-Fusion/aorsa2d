module read_particle_f
implicit none

real, allocatable :: f_rzvv ( :,:,:,: ), &
            vPerp_binCenters ( : ), &
            vPar_binCenters ( : ), &
            R_binCenters ( : ), &
            z_binCenters ( : ), &
            R_binEdges ( : ), &
            z_binEdges ( : ), &
            vPerp_binEdges ( : ), &
            vPar_binEdges ( : )
integer :: vperp_nbins, vpar_nbins, &
           r_nbins, z_nbins
real :: vPerp_binSize, vPar_binSize

contains
    subroutine dlg_particle_f ( capR, y, uperp, upara, nuper, nupar, &
            mi, eNorm_keV, p_f_rzvv, p_dfduPerp, p_dfduPar, &
            aorsa_rho, aorsa_dfduPer_rho, aorsa_dfduPar_rho )

        use netcdf
        use dlg
        use eqdsk_dlg
        use interp

        implicit none
       
        !   capR, y are the R, z values from AORSA, i.e., capR(i) and y(j)
 
        integer :: nuper, nupar, status_,  &
            vper_index, vpar_index
        integer ::  nr, nz, k, l
        real :: enorm, e_, c_, u_, mi, vc_, pi_
        character(len=100) :: ncfilename
        real :: uperp(:), upara(:), &
            uPerp_binSize, uPar_binSize
        real, intent(IN) :: capR, y!, vc_mks
        real, intent(in) :: eNorm_keV 
        real, allocatable :: f_vv(:,:), &
            uPerp_binCenters(:), &
            uPar_binCenters(:), uPerp_2D(:,:), &
            uPerp_binEdges(:), &
            uPar_binEdges(:), &
            dfduPer(:,:), &
            dfduPar(:,:), &
            uPar_2D(:,:)
 
        real, intent(OUT), dimension(nuper,nupar), optional :: &
            p_f_rzvv, p_dfduPerp, p_dfduPar
        integer :: R_index, z_index

        logical :: interpS = .false.
        logical :: interpV = .false.

        !   Interpolation variables for fitpack.f

        integer :: i, m, n, islpsw, iErr, j
        real, allocatable :: zx1(:), zxm(:), zy1(:), zyn(:), zp(:), &
            temp(:), zp_duPer(:), zp_duPar(:)
        real :: zxy11, zxym1, zxy1n, zxymn, surf2
        real, allocatable :: yp_dfduPer(:), yp_dfduPar(:)
        real :: spl1 = 0.0, spln = 0.0, curv2

        !   Flux surface average variables
        
        integer :: q1, q2, q3, q4
        real, allocatable :: R_all(:), z_all(:), &
            psi_all(:), rho_all(:), dfduPer_rho(:,:,:), psiNorm_all(:), &
            dfduPar_rho(:,:,:), f_all(:,:,:), f_rho_vv(:,:,:)
        integer :: nFluxAvg
        real, allocatable :: rho_binEdges(:), rho_binCenters(:)
        real :: dRho, psiRange
        integer :: rhoCnt, rho_nBins
        real, intent(IN), optional :: aorsa_rho(:)
        real, intent(OUT), optional :: &
            aorsa_dfduPer_rho(:,:,:), aorsa_dfduPar_rho(:,:,:)
        logical, allocatable :: maskRho(:), maskVPer(:), maskVPar(:)
        integer, allocatable :: rhoII(:), vPerII(:), vParII(:)
        integer, allocatable :: rhoII_all(:), vPerII_all(:), vParII_all(:)
        integer :: minGoodRhoII

        allocate ( &
            uPerp_binCenters ( vPerp_nBins ), &
            uPar_binCenters ( vPar_nBins ), &
            uPerp_2D ( vPerp_nBins, vPar_nBins ), &
            uPar_2D ( vPerp_nBins, vPar_nBins ), &
            uPerp_binEdges ( vPerp_nBins + 1 ), &
            uPar_binEdges ( vPar_nBins + 1 ), &
            dfduPer ( vPerp_nBins, vPar_nBins ), &
            dfduPar ( vPerp_nBins, vPar_nBins ) )

        uPerp_binCenters    = 0.0
        uPar_binCenters = 0.0
        uPerp_binEdges = 0.0
        uPar_binEdges = 0.0

        !   Normalise the velocity grid, v -> u
        !   u = v / vc_

        e_  = 1.609e-19
        c_  = 3.0e+08
        pi_ = 3.14159265

        enorm   = enorm_keV * 1000.0

        u_  =  mi * c_**2 / ( 2.0 * e_ * enorm )
        vc_ = c_ / sqrt ( u_ )

        uPerp_binCenters    = vPerp_binCenters / vc_
        uPar_binCenters = vPar_binCenters / vc_
        uPerp_binSize   = vPerp_binSize / vc_
        uPar_binSize    = vPar_binSize / vc_
        uPerp_binEdges    = vPerp_binEdges / vc_
        uPar_binEdges = vPar_binEdges / vc_

        do i=1,vPar_nBins
            uPerp_2D(:,i) = uPerp_binCenters
        end do
        do i=1,vPerp_nBins
            uPar_2D(i,:) = uPar_binCenters
        end do
 
        cartesian: &
        if ( .not. present ( aorsa_rho ) ) then 
 
            p_f_rzvv    = 0.0
            p_dfduPerp  = 0.0
            p_dfduPar   = 0.0
    
            !   Either interpolate (slow) or just use the nearest R,z location
            !   from p2f.f90
    
            allocate ( f_vv ( vPerp_nBins, vPar_nBins ) )
            f_vv = 0.0
    
            if ( interpS ) then
    
                !   Interpolate f_rzvv to the AORSA R,z grid point, i.e., capR(i), y(j)
        
                
                m   = R_nBins
                n   = z_nBins
                allocate ( zx1(m), zxm(m), zy1(n), zyn(n), zp(3*m*n), temp(n+n+m) )
                
                islpsw  = 255
              
                !   Yes, I know this is a poor 4D interpolation method. Feel free
                !   to improve it.
         
                do i=1, vPerp_nBins
                    do j=1, vPar_nBins
                        
                        call surf1 ( m, n, R_binCenters, z_binCenters, f_rzvv(:,:,i,j), &
                           m, zx1, zxm, zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
                           zp, temp, sigma, iErr)
                        
                        f_vv(i,j) = surf2 ( capR, y, m, n, &
                           R_binCenters, z_binCenters, f_rzvv(:,:,i,j), m, zp )
                
                    end do
                end do
                
                deallocate ( zx1, zxm, zy1, zyn, zp, temp )
        
            else
               
                !   Find R,z index of nearest p2f grid point
     
                R_index = int ( (capR-minVal(R_binEdges)) / &
                        (maxVal(R_binEdges)-minVal(R_binEdges))*R_nBins+1 )
                z_index =int ( (y-minVal(z_binEdges)) / &
                        (maxVal(z_binEdges)-minVal(z_binEdges))*z_nBins+1 )
                
                if ( R_index > 0 .and. R_index < R_nBins+1 &
                    .and. z_index > 0 .and. z_index < z_nBins+1 ) then
    
                    !write(*,*) 'DLG: ', capR, y, int(R_index), int(z_index), R_nBins, z_nBins
                    f_vv  = f_rzvv(R_index,z_index,:,:)
                
                else
    
                    write(*,*) 'DLG: R,z point outside p2f grid'
                    write(*,'(6(1x,f5.2))') capR, y, minVal ( R_binCenters ), maxVal ( R_binCenters ), &
                        minVal ( z_binCenters ), maxVal ( z_binCenters )
                    f_vv = 0.0
    
                end if
    
            endif
    
            !   Normalise f to uperp/upar integral = 1.0
            !   This will have to change when we alter the way
            !   the density profiles are dealt with in AORSA.
    
            !   While this integration method is crude, it works OK.
    
            if ( sum ( f_vv ) .ne. 0.0 ) then
        
                f_vv    = f_vv / &
                    sum ( f_vv * uPerp_binSize * &
                    uPar_binSize * 2.0 * pi_ * uPerp_2D ) 
        
                !   If f are all zeros then the above integral will give
                !   NaN since there is a divide by zero. This is OK, since we can
                !   just set f and its derivatives to zero here. This usually happens
                !   out near, or eve outside the last closed flux sturface.
            
            else 
    
                f_vv = 0.0
    
            end if
    
    
            !   Calculate uPerp and uPar partial derivatives on the p2f grid
    
            dfduPer    = dlg_pDeriv ( f_vv, 1, uPerp_binSize )
            dfduPar    = dlg_pDeriv ( f_vv, 2, uPar_binSize )
   
            !!   Transform to spherical coords for use in constructing
            !!   the bounce average gradient terms in spherical coords
            !!   for the bounce averaged QL Wdot calculation in AORSA.
            !!   i.e., the same thing they get from CQL3D
    
            !dfdu    = ( uPerp_2D * dfduPer + uPar_2D * dfduPar ) &
            !             / sqrt ( uPerp_2D**2 + uPar_2D**2 )
            !dfdTh    = ( -uPar_2D * dfduPer + uPerp_2D * dfduPar ) &
            !             / ( uPerp_2D**2 + uPar_2D**2 )
 
            !   Now put f_vv onto the AORSA uPerp/uPar grid
    
            m   = vPerp_nBins
            n   = vPar_nBins
    
            if ( interpV ) then ! This interp gives problems, do not use.
    
                allocate ( zx1(m), zxm(m), zy1(n), zyn(n), zp(3*m*n), temp(n+n+m) )
                zx1 = 0.0
                zxm = 0.0
                zy1 = 0.0
                zp  = 0.0
                temp    = 0.0
                islpsw  = 255
    
                call surf1 ( m, n, uPerp_binCenters, uPar_binCenters, &
                    f_vv, &
                    m, zx1, zxm, &
                    zy1, zyn, zxy11, zxym1, zxy1n, zxymn, islpsw, &
                    zp, temp, sigma, iErr )
    
                do i=1,nuper
                    do j=1,nupar
           
                        if ( uperp(i) <= maxVal(uPerp_binCenters) .and. &
                                upara(j) >= minval(uPar_binCenters) .and. &
                                upara(j) <= maxVal(uPar_binCenters) ) then
     
                            p_f_rzvv(i,j) = surf2 ( uperp(i), upara(j), m, n, &
                                uPerp_binCenters, uPar_binCenters, f_vv, m, zp )
                
                            !   Force positive, should only be important down in the 
                            !   noise. Not sure why there are NaNs appearing sometimes 
                            !   though.
    
                            if ( p_f_rzvv(i,j) < 0 ) then 
                                !write (*,*) 'DLG: correcting -ve f'
                                p_f_rzvv(i,j) = 0.0
                            end if
    
                            if ( p_f_rzvv(i,j) .ne. p_f_rzvv(i,j) ) then
                                write(*,*) 'DLG: fvv interp gives NaN!? - fixing'
                                p_f_rzvv(i,j) = 0.0
                            end if
    
                        else
                           
                            p_f_rzvv(i,j)   = 0.0
    
                        end if
        
                    end do
                end do 
        
                deallocate ( zx1, zxm, zy1, zyn, zp, temp )
    
            else ! Use nearest neighbour as value
    
                do i=1,nuper
                    do j=1,nupar
         
                    vPer_index = int ( (uperp(i)-minVal(uPerp_binEdges)) / &
                            (maxVal(uPerp_binEdges)-minVal(uPerp_binEdges))*vPerp_nBins+1 )
                    vPar_index =int ( (upara(j)-minVal(uPar_binEdges)) / &
                            (maxVal(uPar_binEdges)-minVal(uPar_binEdges))*vPar_nBins+1 )
                   
                     if ( vPer_index > 0 .and. vPer_index < vPerp_nBins+1 &
                        .and. vPar_index > 0 .and. vPar_index < vPar_nBins+1 ) then
                        
                        p_f_rzvv(i,j)  = f_vv(vPer_index, vPar_index)
                        p_dfduPerp(i,j)  = dfduPer(vPer_index, vPar_index)
                        p_dfduPar(i,j)  = dfduPar(vPer_index, vPar_index)
                       
                        !if ( 0.5*mi*((uperp(i)*vc_)**2+(upara(j)*vc_)**2)/e_*1d-3 < 3 ) then

                        !    write(*,*) 'dlg: energy below 3keV, setting derivs to 0', &
                        !        0.5*mi*((uperp(i)*vc_)**2+(upara(j)*vc_)**2)/e_*1d-3

                        !    p_f_rzvv(i,j)  = 0 
                        !    p_dfduPerp(i,j)  = 0 
                        !    p_dfduPar(i,j)  = 0 
 
                        !endif

                        sanityCheck: &
                        if ( p_f_rzvv(i,j) .ne. p_f_rzvv(i,j) .or. &
                            p_f_rzvv(i,j) < 0 .or. &
                            p_dfduPerp(i,j) .ne. p_dfduPerp(i,j) .or. &
                            p_dfduPar(i,j) .ne. p_dfduPar(i,j) ) then 

                            write(*,*) 'dlg: ERROR - read particle f failed sanity check.'
                            write(*,*) p_f_rzvv(i,j)
                            write(*,*) p_dfduPerp(i,j)
                            write(*,*) p_dfduPar(i,j)

                        endif sanityCheck
   
                    else
    
                        !write(*,*) 'DLG: vPer,vPar point outside p2f grid'
                        !write(*,*) vPer_index, vPar_index, vPerp_nBins, vPar_nBins
                        !write(*,*) uperp(i), minVal (uPerp_binEdges), maxVal ( uPerp_binEdges )
                        !write(*,*) upara(j), minVal (uPar_binEdges), maxVal ( uPar_binEdges )
    
                        p_f_rzvv(i,j) = 0.0
                        p_dfduPerp(i,j)  = 0.0 
                        p_dfduPar(i,j)  = 0.0
    
                    end if
    
    
                    end do 
                end do
    
            end if
    
            deallocate ( f_vv )
   
        endif cartesian


        !   Create flux surface averaged dfdu & dfdTh

        flux_surf_avg: &
        if ( present ( aorsa_rho ) ) then

            !write(*,*) 'DLG: Creating flux surface average dfdu*'

            allocate ( f_all(R_nBins * z_nBins, vPerp_nBins, vPar_nBins), &
                   R_all(R_nBins * z_nBins), &
                   z_all(R_nBins * z_nBins), &
                   psi_all(R_nBins * z_nBins), & 
                   psiNorm_all(R_nBins * z_nBins), & 
                   rho_all(R_nBins * z_nBins) )

            f_all   = 0
            R_all   = 0
            z_all   = 0
            nFluxAvg    = 0
            psi_all = 0
            psiNorm_all = 0
            rho_all = 0

            do i = 1, R_nBins
                do j = 1, z_nBins
    
                q1  = count ( (R_binCenters(i) - rbbbs) > 0 .and. (z_binCenters(j) - zbbbs) > 0 ) 
                q2  = count ( R_binCenters(i) - rbbbs > 0 .and. z_binCenters(j) - zbbbs <= 0 ) 
                q3  = count ( R_binCenters(i) - rbbbs <= 0 .and. z_binCenters(j) - zbbbs > 0 ) 
                q4  = count ( R_binCenters(i) - rbbbs <= 0 .and. z_binCenters(j) - zbbbs <= 0 ) 
    
                if ( q1 > 0 .and. q2 > 0 .and. q3 > 0 .and. q4 > 0 ) then 
    
                    nFluxAvg    = nFluxAvg + 1         
                    f_all(nFluxAvg,:,:) = f_rzvv(i,j,:,:) * R_binCenters(i)
                    R_all(nFluxAvg) = R_binCenters(i)
                    z_all(nFluxAvg) = z_binCenters(j)
    
                endif
    
                enddo
            enddo 
    
            do i = 1, nFluxAvg
    
                psi_all(i)  = surf2 ( r_all(i), z_all(i), nw, nh, r, z, &
                    psizr, nw, zp_psi, sigma )
    
            enddo

    
            psiRange    = abs ( siMag - siBry )
            psiNorm_all = ( psi_all - siMag ) / psiRange
            rho_all = sqrt ( psiNorm_all )

            deallocate ( psi_all )
            deallocate ( psiNorm_all )
    
            !   Bin by rho coord.
    
            rho_nBins   = R_nBins/2
   
            allocate ( dfduPer_rho(rho_nBins, vPerp_nBins, vPar_nBins), &
                   dfduPar_rho(rho_nBins, vPerp_nBins, vPar_nBins) )

            dfduPer_rho = 0
            dfduPar_rho = 0

            allocate ( rho_binEdges(rho_nBins+1), rho_binCenters(rho_nBins) )
    
            rho_binEdges  = (/ (i*1.0,i=0,rho_nBins) /) / ( rho_nBins )
            dRho    = abs(rho_binEdges(1)-rho_binEdges(2))
            rho_binCenters  = rho_binEdges(2:) - dRho/2.0

            allocate ( f_rho_vv ( rho_nBins, vPerp_nBins, vPar_nBins ) )
            allocate ( maskRho ( R_nBins * z_nBins ) )
            maskRho = .false.
   
            minGoodRhoII    = 0 
            do i = 1, rho_nBins 
   
                maskRho(1:nFluxAvg) =  rho_all(1:nFluxAvg) >= rho_binCenters(i) &
                            .and.  rho_all(1:nFluxAvg) < rho_binCenters(i)+dRho
                rhoCnt  = count ( maskRho )
             
                !   There is a problem here, need to fill in the case where
                !   rhoCnt<0 or the AORSA results will have zeros at rho less
                !   than some value ;-)
                if ( rhoCnt > 0 ) then
                    if ( minGoodRhoII .eq. 0 ) minGoodRhoII = i 
                    f_rho_vv(i,:,:) = sum ( f_all(pack((/(i,i=1,nFluxAvg)/),maskRho),:,:), 1 )
                    f_rho_vv(i,:,:) = f_rho_vv(i,:,:) / rhoCnt
                    !write(*,*) pack((/(i,i=1,nFluxAvg)/),maskRho)

                endif

                !   Normalise to 1

                if ( sum ( f_rho_vv(i,:,:) ) .ne. 0.0 ) then
    
                    f_rho_vv(i,:,:)    = f_rho_vv(i,:,:) / &
                        sum ( f_rho_vv(i,:,:) * uPerp_binSize * &
                        uPar_binSize * 2.0 * pi_ * uPerp_2D ) 
    
                else 

                    f_rho_vv(i,:,:) = 0.0

                end if

                !   Take derivatives
   
                dfduPer_rho(i,:,:)    = &
                    dlg_pDeriv ( f_rho_vv(i,:,:), 1, uPerp_binSize )
                dfduPar_rho(i,:,:)    = &
                    dlg_pDeriv ( f_rho_vv(i,:,:), 2, uPar_binSize )
 
            enddo

            !   Fill in those small rho points that did not have any R,z values

            do i = 1, rho_nBins 
   
                maskRho(1:nFluxAvg) =  rho_all(1:nFluxAvg) >= rho_binCenters(i) &
                            .and.  rho_all(1:nFluxAvg) < rho_binCenters(i)+dRho
                rhoCnt  = count ( maskRho )

                if ( rhoCnt .eq. 0 ) then

                    f_rho_vv(i,:,:) = f_rho_vv(minGoodRhoII,:,:)
                    dfduPer_rho(i,:,:) = dfduPer_rho(minGoodRhoII,:,:)
                    dfduPar_rho(i,:,:) = dfduPar_rho(minGoodRhoII,:,:)

                endif

            enddo

 
            call dlg_write_rho_f ( f_rho_vv, dfduPer_rho, dfduPar_rho, f_all )

            deallocate ( r_all, z_all )
            deallocate ( f_rho_vv )
            deallocate ( f_all )
            deallocate ( maskRho )
            deallocate ( rho_all )
   
            !   Put dfduPer_rho and dfduPar_rho on the aorsa_rho grid
            !   Nearest neighbour in rho for now ... :-(

            allocate ( rhoII_all(size(aorsa_rho)) )
            allocate ( vPerII_all ( nuper ), vParII_all ( nupar ) )

            rhoII_all = int( (aorsa_rho-minVal(rho_binEdges)) / &
                    (maxVal(rho_binEdges)-minVal(rho_binEdges)) * rho_nBins + 1 )

            if ( aorsa_rho(size(aorsa_rho)) .eq. 1.0 ) rhoII_all(size(aorsa_rho)) = rho_nBins
            where ( rhoII_all .eq. 0 ) rhoII_all = 1

            vPerII_all = int ( (uperp-minVal(uPerp_binEdges)) / &
                    (maxVal(uPerp_binEdges)-minVal(uPerp_binEdges))*vPerp_nBins+1 )
            vParII_all =int ( (upara-minVal(uPar_binEdges)) / &
                    (maxVal(uPar_binEdges)-minVal(uPar_binEdges))*vPar_nBins+1 )

            allocate ( maskRho(size(aorsa_rho)),maskVPer(nuper),maskVPar(nupar) )

            maskRho     = rhoII_all > 0 .and. rhoII_all < rho_nBins + 1
            maskVPer    = vPerII_all > 0 .and. vPerII_all < vPerp_nBins+1 
            maskVPar    = vParII_all > 0 .and. vParII_all < vPar_nBins+1

            allocate ( rhoII(count(maskRho)), vPerII(count(maskVPer)), vParII(count(maskVPar)) )

            rhoII   = pack ( rhoII_all, maskRho )
            vPerII   = pack ( vPerII_all, maskVPer )
            vParII   = pack ( vParII_all, maskVPar )

            !write(*,*) rhoII_all
            !write(*,*) rhoII
            !write(*,*) maskRho
            !write(*,*) aorsa_rho
            !allocate ( aorsa_dfduPer_rho ( size(aorsa_rho), nuper, nupar ) )
            !allocate ( aorsa_dfduPar_rho ( size(aorsa_rho), nuper, nupar ) )

            aorsa_dfduPer_rho   = 0
            aorsa_dfduPar_rho   = 0

            !aorsa_dfduPer_rho( rhoII,vPerII,vParII )   = dfduPer_rho( rhoII,vPerII,vParII )
            !aorsa_dfduPar_rho( rhoII,vPerII,vParII )   = dfduPar_rho( rhoII,vPerII,vParII )

            do i=1,size(aorsa_rho)
                do j=1,nuper
                    do k=1,nupar

                        if ( maskRho(i) .and. maskVPer(j) .and. maskVPar(k) ) then
                            aorsa_dfduPer_rho(i,j,k) = &
                                dfduPer_rho ( rhoII_all(i), vPerII_all(j), vParII_all(k) )
                            aorsa_dfduPar_rho(i,j,k) = &
                                dfduPar_rho ( rhoII_all(i), vPerII_all(j), vParII_all(k) )
                        endif

                        !if ( 0.5*mi*((uperp(j)*vc_)**2+(upara(k)*vc_)**2)/e_*1d-3 < 3 ) then

                        !    write(*,*) 'dlg: energy below 3keV, setting derivs to 0[rho]', &
                        !        0.5*mi*((uperp(j)*vc_)**2+(upara(k)*vc_)**2)/e_*1d-3

                        !    aorsa_dfduPer_rho(i,j,k)  = 0 
                        !    aorsa_dfduPar_rho(i,j,k)  = 0 
 
                        !endif

                        sanityCheck2: &
                        if ( aorsa_dfduPer_rho(i,j,k) .ne.  aorsa_dfduPer_rho(i,j,k) .or. &
                                aorsa_dfduPar_rho(i,j,k) .ne. aorsa_dfduPar_rho(i,j,k) ) then
                            write(*,*) 'dlg: ERROR - rho averaged f failed sanity check.'
                            write(*,*) aorsa_dfduPer_rho(i,j,k)
                            write(*,*) aorsa_dfduPar_rho(i,j,k)
                            read(*,*)
                        endif sanityCheck2

                    enddo
                enddo
            enddo 

            deallocate ( rhoII, vPerII, vParII, maskRho, maskVPer, maskVPar )
            deallocate ( rhoII_all, vPerII_all, vParII_all )
            deallocate ( dfduPer_rho, dfduPar_rho )
  
            !allocate ( temp(nFluxAvg), yp_dfduPer(nFluxAvg), yp_dfduPar(nFluxAvg) ) 
    
            !!   curv1 initialises the spline (fitpack.f)
    
            !call curv1 ( nFluxAvg, rho_binCenters, dfduPer_rho, &
            !                spl1, spln, 3, yp_dfduPer, temp, sigma, iErr )
            !call curv1 ( nFluxAvg, rho_binCenters, dfduPar_rho, &
            !                spl1, spln, 3, yp_dfduPar, temp, sigma, iErr )
    
            !allocate ( aorsa_dfduPer_rho(size(aorsa_rho)), &
            !            aorsa_dfduPar_rho(size(aorsa_rho)) )
    
            !do i = 1, size ( aorsa_rho )
            !    
            !    !   curv2 evaluates the spline (fitpack.f)
            !    
            !    aorsa_dfduPer_rho(i) = curv2 ( aorsa_rho(i), nFluxAvg, &
            !        rho_binCenters, dfduPer_rho, &
            !        yp_dfduPer, sigma )
            !    aorsa_dfduPar_rho(i) = curv2 ( aorsa_rho(i), nFluxAvg, &
            !        rho_binCenters, dfduPar_rho, &
            !        yp_dfduPar, sigma )
    
            !enddo
     
            !deallocate ( temp, yp_dfduPer, yp_dfduPar ) 

       
        endif flux_surf_avg

        deallocate ( uPerp_binCenters, &
            uPar_binCenters, &
            uPerp_2D, &
            uPar_2D, &
            uPerp_binEdges, &
            uPar_binEdges, &
            dfduPer, &
            dfduPar )


        !   Write a diagnostic f file to check magnitude etc

        if ( capR > 1.6 .and. capR < 1.7 .and. abs ( y ) < 0.2 ) then

            !write(*,*) 'DLG: writing p_f_aorsa_grid.nc' 
            !call dlg_write_p_f_aorsa_grid ( enorm, nuper, nupar, capR, y, uPerp, uPara, &
            !    p_f_rzvv, p_dfduPerp, p_dfduPar )

        end if

    end subroutine dlg_particle_f

    function dlg_getDensity ( capR, y, ncFileNameIn )
        use netcdf
        use dlg
        implicit none
       
        !   capR, y are the R, z values from AORSA, i.e., capR(i) and y(j)
 
        integer :: nc_id, density_dim_ids(2), &
            R_nBins, z_nBins, &
            R_binCenters_id, z_binCenters_id, &
            R_binEdges_id, z_binEdges_id, &
            density_id
        integer ::  nR, nz, k, l
        character(len=*), intent(in), optional :: ncFileNameIn
        character(len=100) ::  ncFileName
        real, intent(IN) :: capR, y
        real, allocatable :: density(:,:), &
            R_bincenters(:), z_binCenters(:), &
            R_binEdges(:), z_binEdges(:)
        integer :: R_index, z_index, insaneCnt
        real :: dlg_getDensity
        real :: minR, minz, RStep, zStep, RRange, zRange

        !   Read particle f netCdf file
       
        if ( present ( ncFileNameIn ) ) then 
            ncFileName = ncFileNameIn
        else 
            ncFileName  = 'fdis.dav.nc'
        endif

        call dlg_check ( nf90_open ( path = ncFileName, mode = nf90_nowrite, ncid = nc_id ) )
        call dlg_check ( nf90_inq_varId ( nc_id, 'density', density_id ) )
        call dlg_check ( nf90_inq_varId ( nc_id, 'R_binCenters', &
            R_binCenters_id ) )
        call dlg_check ( nf90_inq_varId ( nc_id, 'z_binCenters', &
            z_binCenters_id ) )
        !call dlg_check ( nf90_inq_varId ( nc_id, 'R_binEdges', &
        !    R_binEdges_id ) )
        !call dlg_check ( nf90_inq_varId ( nc_id, 'z_binEdges', &
        !    z_binEdges_id ) )
        
        call dlg_check ( nf90_inquire_variable ( &
            nc_id, density_id, dimIds = density_dim_ids ) )
        
        call dlg_check ( nf90_inquire_dimension ( nc_id, density_dim_ids(1), &
            len = R_nBins ) ) 
        call dlg_check ( nf90_inquire_dimension ( nc_id, density_dim_ids(2), &
            len = z_nBins ) ) 

        allocate ( density ( R_nBins, z_nBins ), &
            R_binCenters ( R_nBins ), &
            z_binCenters ( z_nBins ) )!, &
        !    R_binEdges ( R_nBins+1 ), &
        !    z_binEdges ( z_nBins+1 ) )

        density = 0.0

        call dlg_check ( nf90_get_var ( nc_id, density_id, density ) )
        call dlg_check ( nf90_get_var ( nc_id, R_binCenters_id, &
            R_binCenters ) )
        call dlg_check ( nf90_get_var ( nc_id, z_binCenters_id, &
            z_binCenters ) )
        !call dlg_check ( nf90_get_var ( nc_id, R_binEdges_id, &
        !    R_binEdges ) )
        !call dlg_check ( nf90_get_var ( nc_id, z_binEdges_id, &
        !    z_binEdges ) )
        
        call dlg_check ( nf90_close ( nc_id ) )

        !   Either interpolate (slow) or just use the nearest R,z location
        !   from p2f.f90
          
        !   Find R,z index of nearest p2f grid point

        RStep   = R_binCenters(2) - R_binCenters(1)
        zStep   = z_binCenters(2) - z_binCenters(1)

        RRange  = ( maxVal ( R_binCenters ) - minVal ( R_binCenters ) ) + RStep
        zRange  = ( maxVal ( z_binCenters ) - minVal ( z_binCenters ) ) + zStep

        minR    = minVal ( R_binCenters ) - RStep / 2.0
        minz    = minVal ( z_binCenters ) - zStep / 2.0

        !R_index = int ( (capR-minVal(R_binEdges)) / &
        !        (maxVal(R_binEdges)-minVal(R_binEdges))*R_nBins+1 )
        !z_index =int ( (y-minVal(z_binEdges)) / &
        !        (maxVal(z_binEdges)-minVal(z_binEdges))*z_nBins+1 )
        R_index = nint ( ( capR - minR ) / RRange * R_nBins+1 )
        z_index = nint ( ( y - minz ) / zRange * z_nBins+1 )
        
        if ( R_index > 0 .and. R_index < R_nBins+1 &
            .and. z_index > 0 .and. z_index < z_nBins+1 ) then

            !write(*,*) 'DLG: ', capR, y, density(R_index,z_index)
            dlg_getDensity  = density(R_index,z_index)
        
        else

            !write(*,*) 'DLG: BAD ', capR, y, R_index,z_index, R_nBins, z_nBins
            dlg_getDensity = 0.0

        end if

        sanity_check: &
        if ( dlg_getDensity /= dlg_getDensity &
            .or. dlg_getDensity < 0 &
            .or. dlg_getDensity * 0 /= 0 ) then

            write(*,*) 'ERROR: sanity failure on reading density from netcdf file, NaN, Inf or -ve value detected'
            stop
        endif sanity_check

        deallocate ( density, &
            R_binCenters, &
            z_binCenters )!, &
            !R_binEdges, &
            !z_binEdges)

        return

    end function dlg_getDensity

    subroutine dlg_write_rho_f ( &
                    f_rho_vv, &
                    dfduPer, &
                    dfduPar, &
                    f_all )

        use netcdf 
        implicit none

        real, intent(IN) :: &
            f_rho_vv(:,:,:), &
            dfduPer(:,:,:), &
            dfduPar(:,:,:), &
            f_all(:,:,:)

        character(len=100) :: ncFileName

        integer :: nc_id, &
            f_rho_vv_id, &
            dfduPer_id, &
            dfduPar_id, &
            scalar_id, &
            vPer_nBins_id, &
            vPar_nBins_id, &
            rho_nBins_id, &
            nRz_id, &
            f_all_id


        ncFileName = 'output/p_f_rho.nc'
        
        call dlg_check ( nf90_create ( ncFileName, nf90_clobber, nc_id ) )
        call dlg_check ( nf90_def_dim ( nc_id, "vPer_nBins", size(f_rho_vv,2) , vPer_nBins_id ) )
        call dlg_check ( nf90_def_dim ( nc_id, "vPar_nBins", size(f_rho_vv,3) , vPar_nBins_id ) )
        call dlg_check ( nf90_def_dim ( nc_id, "rho_nBins", size(f_rho_vv,1) , rho_nBins_id ) )
        call dlg_check ( nf90_def_dim ( nc_id, "scalar", 1 , scalar_id ) )
        call dlg_check ( nf90_def_dim ( nc_id, "nRz", size(f_all,1) , nRz_id ) )

        call dlg_check ( nf90_def_var ( nc_id, "f_rho_vv", NF90_REAL, &
          (/ rho_nBins_id, vPer_nBins_id, vPar_nBins_id /), f_rho_vv_id ) )
        call dlg_check ( nf90_def_var ( nc_id, "dfduPer", NF90_REAL, &
          (/ rho_nBins_id, vPer_nBins_id, vPar_nBins_id /), dfduPer_id ) )
        call dlg_check ( nf90_def_var ( nc_id, "dfduPar", NF90_REAL, &
          (/ rho_nBins_id, vPer_nBins_id, vPar_nBins_id /), dfduPar_id ) )
        call dlg_check ( nf90_def_var ( nc_id, "f_all", NF90_REAL, &
          (/ nRz_id, vPer_nBins_id, vPar_nBins_id /), f_all_id ) )
       
        call dlg_check ( nf90_enddef ( nc_id ) )
        
        call dlg_check ( nf90_put_var ( nc_id, f_rho_vv_id, f_rho_vv ) )
        call dlg_check ( nf90_put_var ( nc_id, dfduPer_id, dfduPer ) )
        call dlg_check ( nf90_put_var ( nc_id, dfduPar_id, dfduPar ) )
        call dlg_check ( nf90_put_var ( nc_id, f_all_id, f_all ) )
   
        call dlg_check ( nf90_close ( nc_id ) )

    end subroutine dlg_write_rho_f



    subroutine dlg_write_p_f_aorsa_grid ( enorm, nuper, nupar, capR, y, uPerp, uPara, &
        p_f_vv, p_dfduPerp, p_dfduPar )
        use netcdf 
        implicit none

        integer, intent(IN) :: nuper, nupar 
        real, intent(IN) :: uPerp(:), uPara(:), capR, y, &
            p_f_vv(:,:), p_dfduPerp(:,:), p_dfduPar(:,:), enorm
        character(len=100) :: ncFileName
        integer :: nc_id, n_uper_id, n_upar_id, &
            p_f_vv_id, p_dfduPer_id, p_dfduPar_id, &
            uPerp_id, uPara_id, capR_id, y_id, nR_id, nz_id, &
            scalar_id, p_enorm_id


        ncFileName = 'output/p_f_aorsa_grid.nc'
        
        call dlg_check ( nf90_create ( ncFileName, nf90_clobber, nc_id ) )
        call dlg_check ( nf90_def_dim ( nc_id, "n_uper", nuper , n_uper_id ) )
        call dlg_check ( nf90_def_dim ( nc_id, "n_upar", nupar , n_upar_id ) )
        !call dlg_check ( nf90_def_dim ( nc_id, "n_R", nR , nR_id ) )
        !call dlg_check ( nf90_def_dim ( nc_id, "n_z", nz , nz_id ) )
        call dlg_check ( nf90_def_dim ( nc_id, "scalar", 1 , scalar_id ) )

        
        call dlg_check ( nf90_def_var ( nc_id, "p_f_vv", NF90_REAL, &
          (/ n_uper_id, n_upar_id /), p_f_vv_id ) )
        call dlg_check ( nf90_def_var ( nc_id, "p_dfduPer", NF90_REAL, &
          (/ n_uper_id, n_upar_id /), p_dfduPer_id ) )
        call dlg_check ( nf90_def_var ( nc_id, "p_dfduPar", NF90_REAL, &
          (/ n_uper_id, n_upar_id /), p_dfduPar_id ) )
        call dlg_check ( nf90_def_var ( nc_id, "p_enorm", NF90_REAL, &
          scalar_id, p_enorm_id ) )
        
        call dlg_check ( nf90_def_var ( nc_id, "uPerp", NF90_REAL, &
          (/ n_uper_id /), uPerp_id ) )
        call dlg_check ( nf90_def_var ( nc_id, "uPara", NF90_REAL, &
          (/ n_upar_id /), uPara_id ) )
         call dlg_check ( nf90_def_var ( nc_id, "capR", NF90_REAL, &
          (/ scalar_id /), capR_id ) )
        call dlg_check ( nf90_def_var ( nc_id, "y", NF90_REAL, &
          (/ scalar_id /), y_id ) )
        
        call dlg_check ( nf90_enddef ( nc_id ) )
        
        call dlg_check ( nf90_put_var ( nc_id, uPerp_id, uperp ) )
        call dlg_check ( nf90_put_var ( nc_id, uPara_id, upara ) )
        call dlg_check ( nf90_put_var ( nc_id, p_f_vv_id, p_f_vv ) )
        call dlg_check ( nf90_put_var ( nc_id, p_dfduPer_id, p_dfduPerp ) )
        call dlg_check ( nf90_put_var ( nc_id, p_dfduPar_id, p_dfduPar ) )
        call dlg_check ( nf90_put_var ( nc_id, capR_id, capR ) )
        call dlg_check ( nf90_put_var ( nc_id, y_id, y ) )
        call dlg_check ( nf90_put_var ( nc_id, p_enorm_id, enorm ) )
       
        call dlg_check ( nf90_close ( nc_id ) )

    end subroutine dlg_write_p_f_aorsa_grid

!    subroutine dlg_read_p_f_aorsa_grid ( nR, nz, nuPer, nuPar, p_f_rzvv, p_dfduPerp, p_dfduPar )
!        use netcdf
!        implicit none
!        
!        character(len=100) :: ncFileName
!        integer :: nc_id, p_f_vvid, p_dfduPerp_id, p_dfduPar_id, &
!            p_f_vvdim_ids(4), nuPerp_, nuPar_, nuPer, nuPar, nR, nz, &
!            nR_, nz_
!        real, dimension (:,:,:,:) :: p_f_rzvv, p_dfduPerp, p_dfduPar
!
!
!        ncFileName  = 'output/p_f_aorsa_grid.nc'
!
!        call dlg_check ( nf90_open ( path = ncFileName, mode = nf90_nowrite, ncid = nc_id ) )
!        call dlg_check ( nf90_inq_varId ( nc_id, 'p_f_rzvv', p_f_vvid ) )
!        call dlg_check ( nf90_inq_varId ( nc_id, 'p_dfduPer', &
!            p_dfduPerp_id ) )
!        call dlg_check ( nf90_inq_varId ( nc_id, 'p_dfduPar', &
!            p_dfduPar_id ) )
!
!        call dlg_check ( nf90_inquire_variable ( &
!            nc_id, p_f_vvid, dimIds = p_f_vvdim_ids ) )
!        call dlg_check ( nf90_inquire_dimension ( nc_id, p_f_vvdim_ids(3), &
!            len = nuPerp_ ) ) 
!        call dlg_check ( nf90_inquire_dimension ( nc_id, p_f_vvdim_ids(4), &
!            len = nuPar_ ) ) 
!        call dlg_check ( nf90_inquire_dimension ( nc_id, p_f_vvdim_ids(1), &
!            len = nR_ ) ) 
!        call dlg_check ( nf90_inquire_dimension ( nc_id, p_f_vvdim_ids(2), &
!            len = nz_ ) ) 
!
!        if ( nuPer /= nuPerp_ .or. nuPar /= nuPar_ .or. nR_ /= nR .or. nz_ /= nz ) then
!            write(*,*) 'DLG: ERROR ( dlg_read_p_f_aorsa_grid )'
!            read (*,*)
!        end if
! 
!        !allocate ( p_f_vv ( nuPerp_, nuPar_ ), &
!        !    p_dfduPerp ( nuPerp_, nuPar_ ), &
!        !    p_dfduPar ( nuPerp_, nuPar_ ) )
!
!        call dlg_check ( nf90_get_var ( nc_id, p_f_vvid, p_f_rzvv ) )
!        call dlg_check ( nf90_get_var ( nc_id, p_dfduPerp_id, p_dfduPerp ) )
!        call dlg_check ( nf90_get_var ( nc_id, p_dfduPar_id, p_dfduPar ) )
!
!        call dlg_check ( nf90_close ( nc_id ) )
!
!    end subroutine dlg_read_p_f_aorsa_grid
!    
    subroutine dlg_check ( status_ )
        use netcdf
        integer, intent ( in) :: status_
      
        if(status_ /= nf90_noerr) then 
            write(*,*) trim(nf90_strerror(status_))
            stop "Stopped"
        end if

    end subroutine dlg_check
 

    subroutine init_particleFile ( mpi_pId )
        use netcdf
        implicit none
        integer :: nc_id, f_vv_dim_ids(2), &
            f_vv_id, vperp_bincenters_id, &
            vpar_bincenters_id, vperp_binsize_id, vpar_binsize_id, &
            f_vvid, f_vvdim_ids(4), &
            r_bincenters_id, z_bincenters_id, &
            r_binedges_id, z_binedges_id, &
            vperp_binedges_id, vpar_binedges_id
        integer :: k, l, mpi_pId
        character(len=100) :: ncfilename

        !   Read particle f netCdf file

        if ( mpi_pId .eq. 0 ) &
            write(*,*) 'DLG: Reading fdis.dav.nc ...' 

        ncFileName  = 'fdis.dav.nc'

        call dlg_check ( nf90_open ( path = ncFileName, mode = nf90_nowrite, ncid = nc_id ) )
        call dlg_check ( nf90_inq_varId ( nc_id, 'f_rzvv', f_vvid ) )
        call dlg_check ( nf90_inq_varId ( nc_id, 'vPerp_binCenters', &
            vPerp_binCenters_id ) )
        call dlg_check ( nf90_inq_varId ( nc_id, 'vPar_binCenters', &
            vPar_binCenters_id ) )
        call dlg_check ( nf90_inq_varId ( nc_id, 'R_binCenters', &
            R_binCenters_id ) )
        call dlg_check ( nf90_inq_varId ( nc_id, 'z_binCenters', &
            z_binCenters_id ) )
        call dlg_check ( nf90_inq_varId ( nc_id, 'vPerp_binSize', &
            vPerp_binSize_id ) )
        call dlg_check ( nf90_inq_varId ( nc_id, 'vPar_binSize', &
            vPar_binSize_id ) )
         call dlg_check ( nf90_inq_varId ( nc_id, 'R_binEdges', &
            R_binEdges_id ) )
        call dlg_check ( nf90_inq_varId ( nc_id, 'z_binEdges', &
            z_binEdges_id ) )
        call dlg_check ( nf90_inq_varId ( nc_id, 'vPerp_binEdges', &
            vPerp_binEdges_id ) )
        call dlg_check ( nf90_inq_varId ( nc_id, 'vPar_binEdges', &
            vPar_binEdges_id ) )
        
        call dlg_check ( nf90_inquire_variable ( &
            nc_id, f_vvid, dimIds = f_vvdim_ids ) )
        
        call dlg_check ( nf90_inquire_dimension ( nc_id, f_vvdim_ids(1), &
            len = R_nBins ) ) 
        call dlg_check ( nf90_inquire_dimension ( nc_id, f_vvdim_ids(2), &
            len = z_nBins ) ) 
        call dlg_check ( nf90_inquire_dimension ( nc_id, f_vvdim_ids(3), &
            len = vPerp_nBins ) ) 
        call dlg_check ( nf90_inquire_dimension ( nc_id, f_vvdim_ids(4), &
            len = vPar_nBins ) ) 

        allocate ( f_rzvv ( R_nBins, z_nBins, vPerp_nBins, vPar_nBins ), &
            vPerp_binCenters ( vPerp_nBins ), &
            vPar_binCenters ( vPar_nBins ), &
            R_binCenters ( R_nBins ), &
            z_binCenters ( z_nBins ), &
            R_binEdges ( R_nBins+1 ), &
            z_binEdges ( z_nBins+1 ), &
            vPerp_binEdges ( vPerp_nBins + 1 ), &
            vPar_binEdges ( vPar_nBins + 1 ) )

        vPerp_binCenters = 0.0
        vPar_binCenters = 0.0
        vPerp_binEdges = 0.0
        vPar_binEdges = 0.0

        call dlg_check ( nf90_get_var ( nc_id, f_vvid, f_rzvv ) )
        call dlg_check ( nf90_get_var ( nc_id, vPerp_binCenters_id, &
            vPerp_binCenters ) )
        call dlg_check ( nf90_get_var ( nc_id, vPar_binCenters_id, &
            vPar_binCenters ) )
        call dlg_check ( nf90_get_var ( nc_id, vPerp_binSize_id, &
            vPerp_binSize ) )
        call dlg_check ( nf90_get_var ( nc_id, vPar_binSize_id, &
            vPar_binSize ) )
        call dlg_check ( nf90_get_var ( nc_id, R_binCenters_id, &
            R_binCenters ) )
        call dlg_check ( nf90_get_var ( nc_id, z_binCenters_id, &
            z_binCenters ) )
        call dlg_check ( nf90_get_var ( nc_id, R_binEdges_id, &
            R_binEdges ) )
        call dlg_check ( nf90_get_var ( nc_id, z_binEdges_id, &
            z_binEdges ) )
        call dlg_check ( nf90_get_var ( nc_id, vPerp_binEdges_id, &
            vPerp_binEdges ) )
        call dlg_check ( nf90_get_var ( nc_id, vPar_binEdges_id, &
            vPar_binEdges ) )
        
        call dlg_check ( nf90_close ( nc_id ) )

        if ( mpi_pId .eq. 0 ) &
            write(*,*) 'DLG: DONE' 

    end subroutine init_particleFile 

end module read_particle_f
