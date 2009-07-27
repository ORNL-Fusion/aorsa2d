module gc_integrate

contains
    function dlg_gc_velocity ( pos, vPerp, vPar )
        use gc_terms
        use interp
        use eqdsk_dlg

        implicit none

        real :: dlg_gc_velocity(3)
        real, intent(IN) :: pos(3), vPerp, vPar
        real :: grad_R, grad_phi, grad_z, &
            curv_R, curv_phi, curv_z, &
            unitb_R, unitb_phi, unitb_z
        real :: surf2, bMagHere, bHere(3), &
            vgc_R, vgc_phi, vgc_z

        grad_R  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bGradient_R, nw, zp_bGrad_R, sigma )
        grad_phi  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bGradient_phi, nw, zp_bGrad_phi, sigma )
        grad_z  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bGradient_z, nw, zp_bGrad_z, sigma )

        curv_R  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bCurvature_R, nw, zp_bCurv_R, sigma )
        curv_phi  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bCurvature_phi, nw, zp_bCurv_phi, sigma )
        curv_z  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bCurvature_z, nw, zp_bCurv_z, sigma )

        bHere   = dlg_interpB ( pos, bMagHere = bMagHere )
        unitb_R = bHere(1) / bMagHere
        unitb_phi   = bHere(2) / bMagHere
        unitb_z = bHere(3) / bMagHere

        vgc_R   = vPar * unitb_R + vPerp**2 * grad_R + vPar**2 * curv_R 
        vgc_phi   = vPar * unitb_phi + vPerp**2 * grad_phi + vPar**2 * curv_phi
        vgc_z   = vPar * unitb_z + vPerp**2 * grad_z + vPar**2 * curv_z

        dlg_gc_velocity = (/ vgc_R, vgc_phi, vgc_z /)

    end function dlg_gc_velocity

    function dlg_vPerp ( pos, u )
        use constants
        implicit none

        real :: bHere(3),bMagHere
        real, intent(IN) :: pos(3), u
        real :: dlg_vPerp

        bHere   = dlg_interpB ( pos, bMagHere = bMagHere )
        dlg_vPerp   = sqrt ( 2.0 * u * bMagHere / mi )

    end function dlg_vPerp

    function dlg_vPar ( pos, u )
        use constants
        use gc_terms
        use interp
        use eqdsk_dlg
        implicit none
        
        real :: dlg_vPar, bDotGradB_here, surf2
        real, intent(IN) :: pos(3), u

        bDotGradB_here  = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bDotGradB, nw, zp_bDotGradB, sigma )

        dlg_vPar    = -u / mi * bDotGradB_here 

    end function dlg_vPar

    function dlg_interpB ( pos, bMagHere )
        use eqdsk_dlg
        use interp
        implicit none
        
        real :: bR_here, bPhi_here, bz_here
        real, intent(IN) :: pos(3)
        real :: dlg_interpB(3)
        real, optional, intent(OUT) :: bMagHere
        real :: surf2

        bR_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bR, nw, zp_bR, sigma )
        bPhi_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bPhi, nw, zp_bPhi, sigma )
        bz_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
            bz__, nw, zp_bz, sigma )

        if ( present (bMagHere) ) &
            bMagHere    = sqrt ( bR_here**2 + bPhi_here**2 + bz_here**2 )

        dlg_interpB(1)  = bR_here 
        dlg_interpB(2)  = bPhi_here
        dlg_interpB(3)  = bz_here 
    
    end function dlg_interpB

    subroutine gc_orbit ( start_R, start_z, start_vPerp, start_vPar, weight, &
        R_nBins, z_nBins, R_range, z_range, rBottom, zBottom, orbit_2D, stillIn, plot )
        use eqdsk_dlg
        use gc_terms
        use constants

#if USE_DISLIN     
        use dislin 
#endif

        use interp
        !use rzvv_grid
        !use init_mpi
        implicit none
       
        integer, parameter :: r4 = selected_real_kind(p=6,r=37)
        real, intent(IN) :: start_R, start_z, start_vPerp, start_vPar, &
            weight
        integer, parameter :: maxSteps = 5000
        integer :: stepCnt, var_dt, i, j, k, l, ii
        real, dimension (maxSteps) :: rTrack, zTrack, &
            vPerpTrack, vParTrack, distance, dtArray, &
            R_index, z_index, &
            vPerp_index, vPar_index
        real :: tau, dTau, pos(3), &
            u, bMagHere, bHere(3), vPerp, vPar, &
            vgc(3), dt, dtMin, dtMax
        logical :: stillIn
        logical :: firstOrbit
        logical, optional, intent(IN) :: plot
        real :: k1_vPar, k1_vgc(3), k2_vPar, k2_vgc(3), &
            k3_vPar, k3_vgc(3), k4_vPar, k4_vgc(3), lastStep, &
            psi_here, surf2
        
        real :: R_sigma, z_sigma, vPerp_sigma, vPar_sigma

        integer, dimension (maxSteps) :: vPerp_start, vPerp_stop, &
            vPar_start, vPar_stop, R_start, R_stop, &
            z_start, z_stop
        integer :: sR

        !   Additional variables for the AORSA interface

        integer :: R_nBins, z_nBins
        real :: R_range, z_range, rBottom, &
            zBottom, orbit_2D(R_nBins,z_nBins)
        
        !real :: f_rzvv_update(R_nBins,z_nBins,vPerp_nBins,vPar_nBins)

#if USE_DISLIN

        !   Plotting variables
    
        integer :: nLevs, nxPag, nyPag
        real :: levStep
        real, allocatable :: levels(:)

#endif
        
        !   Initialize variables
        
        stepCnt = 1
        firstOrbit  = .true.
        stillIn = .true.
        dTau    = 0.0
        tau = 0.0
        var_dt  = 1

        !! Interpolation testing variables

        !real, allocatable :: bR_interp(:,:)
        !real :: b_interp(3)
        !integer :: j

10      continue

        pos(1)  = start_R
        pos(2)  = 0.0
        pos(3)  = start_z

        vPerp   = start_vPerp
        vPar    = start_vPar

        bHere   = dlg_interpB ( pos, bMagHere = bMagHere )

        u   = mi * vPerp**2 / ( 2.0 * bMagHere ) 

        dtMin   = 0.005e-7
        if ( var_dt == 1 ) then
            dtMax   = 5.0e-7
        else 
            dtMin   = 0.0005e-7
            dtMax   = 0.05e-7
        end if

        dt  = dtMin

        !!   Test the surf1/2 interpolations

        !allocate ( bR_interp(nw,nh) )
        !do i=1,nw
        !    do j=1,nh
        !        b_interp    = dlg_interpB ( (/ r(i), 0.0, z(j) /) )
        !        bR_interp(i,j)  = b_interp(1)
        !    end do
        !end do

#if USE_DISLIN
        if ( present ( plot ) ) then 
            if ( plot ) then  
           
                !write(*,*) 'DLG GC: ', start_R, start_z
 
                call scrMod ( 'REVERS' )
                call setPag ( 'DA4P' )! nxPag = 2100, nyPag = 2970
                call metaFl ( 'XWIN' )
                call disIni ()
                call winMod ( 'NONE' )
                call erase ()
                call noChek ()
                call unit ( 0 )
                call getPag ( nxPag, nyPag )
                call axsPos ( 300, 1700 )
                call axsLen ( 1500, 1400 )
                call graf ( real(0.0,r4), real(3.0,r4), &
                    real(0.0,r4), real(0.5,r4), real(-2.0,r4), &
                    real(2.0,r4), real(-2.0,r4), real(0.5,r4) ) 
                call noClip ()
                call color ('BLUE')         
                
                nLevs   = 101 
                levStep    = 0.05!maxVal ( abs ( bPhi ) ) / ( nLevs / 2 )
                if ( .not. allocated ( levels ) ) &
                    allocate ( levels(nLevs) )
                levels = (/ (i*levStep,i=-nLevs/2,nLevs/2) /)
                do i=1,nLevs
                    call contur ( real(r,r4), nw, real(z,r4), nh, &
                        real(psizr,r4), real(levels(i),r4) ) 
                end do
                call color ( 'MAGENTA' )
                call curve ( real(rbbbs,r4), real(zbbbs,r4), nbbbs )

                call endGrf ()
                call color ('RED' )
                !call setGrf ( 'NONE', 'NONE', 'NONE', 'NONE' ) 

            end if 
        end if 
#endif

        do 

            vPerp   = dlg_vPerp ( pos, u ) 
            vgc = dlg_gc_velocity ( pos, vPerp, vPar )
            k1_vPar   = dt * dlg_vPar ( pos, u ) 
            k1_vgc  = dt * vgc

            vPerp   = dlg_vPerp ( pos + k1_vgc / 2.0, u ) 
            vgc = dlg_gc_velocity ( pos + k1_vgc / 2.0, vPerp, vPar + k1_vPar / 2.0 )
            k2_vPar   = dt * dlg_vPar ( pos + k1_vgc / 2.0, u ) 
            k2_vgc  = dt * vgc

            vPerp   = dlg_vPerp ( pos + k2_vgc / 2.0, u ) 
            vgc = dlg_gc_velocity ( pos + k2_vgc / 2.0, vPerp, vPar + k2_vPar / 2.0 )
            k3_vPar   = dt * dlg_vPar ( pos + k2_vgc / 2.0, u ) 
            k3_vgc  = dt * vgc

            vPerp   = dlg_vPerp ( pos + k3_vgc, u ) 
            vgc = dlg_gc_velocity ( pos + k3_vgc, vPerp, vPar + k3_vPar )
            k4_vPar   = dt * dlg_vPar ( pos + k3_vgc, u ) 
            k4_vgc  = dt * vgc

            vPar    = vPar + ( k1_vPar + 2.0 * k2_vPar + 2.0 * k3_vPar &
                + k4_vPar ) / 6.0
            pos   = pos + ( k1_vgc + 2.0 * k2_vgc + 2.0 * k3_vgc + k4_vgc ) / 6.0

            !   Interpolate to get psi at new pos
            !   and if psi(pos) is outside the last 
            !   closed flux surface then dump the particle
            !   The 0.98 is a factor suggested by EFJ to keep
            !   away boundary of the flux grid.

            psi_here = surf2 ( pos(1), pos(3), nw, nh, r, z, &
                psizr, nw, zp_psi, sigma )

            if ( psi_here > sibry * 0.98 ) stillIn = .false.

            if ( (.not. stillIn) .or. (.not. firstOrbit) .or. stepCnt+1 >= maxSteps ) then

                if ( stepCnt+1 >= maxSteps ) then
                    if ( var_dt == 1 ) then 
                        
                        !   The algorithm seems to trace the particles well,
                        !   including finding the orbit end points. Howevere,
                        !   there seem to be two exceptions that I have found so
                        !   far, these being particles whose orbits track very
                        !   near the plasma boundary and those that are on the
                        !   trapped passing boundary. Reducing the step size
                        !   does not seem to fix the trajectories near the
                        !   plasma boundary while the trapped/passing boundary
                        !   orbits are rare.
                           
                        !write(*,*) 'Suspected bad trace [long], adjusting dt range'
                        stepCnt = 1
                        var_dt = 0
                        !if ( present ( plot ) ) then
                        !    if ( plot ) call disFin ()
                        !endif
                        go to 10 
                    else
                        write(*,*) 'DLG: bad'
                        !nP_bad  = nP_bad + 1
                    end if
                end if
                
                if ( .not. stillIn ) write(*,*) 'DLG: wall' 
                !write(*,*) 'DLG: Orbit complete'
                exit

            end if

            rTrack(stepCnt) = pos(1)
            zTrack(stepCnt) = pos(3)
            vPerpTrack(stepCnt) = vPerp
            vParTrack(stepCnt)  = vPar
            tau = tau + dt
            dtArray(stepCnt)    = dt
           
            distance(stepCnt)   = sqrt ( ( start_R - pos(1) )**2 &
                + ( start_z - pos(3) )**2 )
      
            !   Detect end of orbit
             
            if ( stepCnt > 50 ) then

                if ( ( distance(stepCnt) - distance(stepCnt-1) ) > 0 .AND. &
                     ( distance(stepCnt-1) - distance(stepCnt-2) ) < 0 .AND. &
                     distance(stepCnt-1) < 0.0006  ) firstOrbit = .false.

            end if

            if ( firstOrbit .and. stepCnt > 1 ) then 

                !   Adjust time step 

                lastStep    = sqrt ( ( rTrack(stepCnt) - rTrack(stepCnt-1) )**2 &
                    + ( zTrack(stepCnt) - zTrack(stepCnt-1) )**2 )

                if ( distance(stepCnt) > 0.1 ) then
        
                    dt  = 0.01 / lastStep * dt
           
                else
   
                    if ( distance(stepCnt) > 0.01 ) then
                        dt = 0.001 / lastStep * dt
                    else
                        dt = 0.0001 / lastStep * dt
                    end if
 
                end if
                
                if ( dt < dtMin ) dt = dtMin
                if ( dt > dtMax ) dt = dtMax
 
            end if 

            stepCnt = stepCnt + 1

#if USE_DISLIN           
            !   Plot track

            if ( present ( plot ) ) then
                if ( plot ) then  

                    if ( stepCnt > 2 ) then

                        call axsPos ( 300, 1700 )
                        call axsLen ( 1500, 1400 )
                        call color ( 'WHITE' )
                        call graf ( real(0.0,r4), real(3.0,r4), &
                            real(0.0,r4), real(0.5,r4), real(-2.0,r4), &
                            real(2.0,r4), real(-2.0,r4), real(0.5,r4) ) 
                        call curve ( real(rTrack(1:stepCnt-1),r4), &
                            real(zTrack(1:stepCnt-1),r4), stepCnt-1 )
                        call endGrf ()
                        !call axsPos ( 300, 2600 )
                        !call axsLen ( 1500, 700 )
                        !call graf ( start_R-0.01, start_R+0.01, start_R-0.01,&
                        !    0.005, start_z-0.01, start_z+0.01,start_z-0.01, 0.005 ) 
                        !call color ( 'RED' )
                        !call curve ( rTrack(1:stepCnt-1), zTrack(1:stepCnt-1), stepCnt-1 )
                        !call endGrf ()

                    end if

                endif
            end if
#endif
                
       end do

        !call disFin ()
        
        !   Calculate f_rzvv indices and add this orbit
        !   to the f_rzvv particle count. At this point the
        !   f_rzvv function will only be the number of particles
        !   but since its a regular grid we can divide by the
        !   volume element after counting all the particles ;-)

        if ( stillIn ) then 
        
            R_index = ( rTrack - rBottom ) / R_range * R_nBins + 1
            z_index = ( zTrack - zBottom ) / z_range * z_nBins + 1
            !vPerp_index = vPerpTrack / vPerp_range * vPerp_nBins + 1
            !vPar_index  = ( vParTrack + vPar_range ) / ( 2.0 * vPar_range ) &
            !    * vPar_nBins + 1

            do i=1,stepCnt
                
                if ( int(R_index(i)) <= R_nBins .and. int(R_index(i)) >= 1 .and. &
                    int(z_index(i)) <= z_nBins .and. int(z_index(i)) >= 1 ) then
                    
                   !write(*,*) 'Adding orbit data', &
                   !     int(R_index(i)),int(z_index(i))
                   orbit_2D(int(R_index(i)),int(z_index(i))) &
                       = orbit_2D(int(R_index(i)),int(z_index(i))) &
                            + dtArray(i) / tau * weight
                else
                    !write(*,*) 'Not adding orbit data'
                end if

            end do

            !!   Try adding a particle of some finite size 
            !!   described by a gaussian in 4D space

            !R_sigma = 0.01
            !z_sigma = 0.01
            !vPerp_sigma = 1e-4*c
            !vPar_sigma = 1e-4*c
            !
            !sR  = 2

            !!   Loop only over a few cells around the particle location, 
            !!   otherwise this is prohibitavely slow. Its already a factor
            !!   of 10 slower than using the delta functions above :-(

            !vPerp_start = int(vPerp_index)-sR
            !where (vPerp_start < 1) vPerp_start = 1
            !vPerp_stop  = int(vPerp_index)+sR
            !where (vPerp_stop > vPerp_nBins) vPerp_stop = vPerp_nBins
 
            !vPar_start = int(vPar_index)-sR
            !where (vPar_start < 1) vPar_start = 1
            !vPar_stop  = int(vPar_index)+sR
            !where (vPar_stop > vPar_nBins) vPar_stop = vPar_nBins
 
            !R_start = int(R_index)-sR
            !where (R_start < 1) R_start = 1
            !R_stop  = int(R_index)+sR
            !where (R_stop > R_nBins) R_stop = R_nBins
 
            !z_start = int(z_index)-sR
            !where (z_start < 1) z_start = 1
            !z_stop  = int(z_index)+sR
            !where (z_stop > z_nBins) z_stop = z_nBins

            !do ii=1,stepCnt          

            !do i=vPerp_start(ii),vPerp_stop(ii)
            !    do j=vPar_start(ii),vPar_stop(ii)
            !        do l=z_start(ii),z_stop(ii)
            !            do k=R_start(ii),R_stop(ii)
            !               
            !                !   The * vPerp_binCenters at the end of this
            !                !   command is to ensure the particle has a gaussian
            !                !   in real phase space ;-)
 
            !                f_rzvv_update(k,l,i,j) =  dtArray(ii) / tau * exp ( -( &
            !                    ( R_binCenters(k) - rTrack(ii) )**2 / ( 2.0 * R_sigma**2 ) &
            !                    + ( z_binCenters(l) - zTrack(ii) )**2 / ( 2.0 * z_sigma**2 ) &
            !                    + ( vPerp_binCenters(i) - vPerpTrack(ii) )**2 / ( 2.0 * vPerp_sigma**2 ) &
            !                    + ( vPar_binCenters(j) - vParTrack(ii) )**2 / ( 2.0 * vPar_sigma**2 ) ) ) &
            !                    * vPerp_binCenters(i)

            !            end do
            !        end do
            !    end do
            !end do
           
            !f_rzvv(R_start(ii):R_stop(ii),z_start(ii):z_stop(ii),&
            !        vPerp_start(ii):vPerp_stop(ii),vPar_start(ii):vPar_stop(ii))  = &
            !    f_rzvv(R_start(ii):R_stop(ii),z_start(ii):z_stop(ii),&
            !        vPerp_start(ii):vPerp_stop(ii),vPar_start(ii):vPar_stop(ii)) &
            !    + f_rzvv_update(R_start(ii):R_stop(ii),z_start(ii):z_stop(ii),&
            !                vPerp_start(ii):vPerp_stop(ii),vPar_start(ii):vPar_stop(ii))
   
            !end do 

        
        end if

    end subroutine gc_orbit 

end module gc_integrate
