module eqdsk_dlg
    implicit none
    save    
    character (len=10) :: case_ ( 6 )
    integer :: idum, nw, nh, nbbbs, limitr
    real :: rdim, zdim, rcentr, rleft, zmid, &
        rmaxis__, zmaxis__, simag, sibry, bcentr, & 
        current__, xdum, rStep, zStep, fStep
    
    real, allocatable :: fpol (:), psizr (:,:), &
        pres (:), ffprim (:), pprime (:), qpsi (:), &
        rbbbs (:), zbbbs (:), rlim__ (:), zlim (:), &
        r (:), z (:), bR(:,:), bPhi(:,:), bz__(:,:), &
        fluxGrid (:), fpolRZ(:,:), bMag__(:,:)

contains
    subroutine read_geqdsk ( eqdsk_fileName, plot )
        !use dislin
        use dlg
 
        implicit none
        
        integer :: i, j, iErr
        character(len=*), intent(IN) :: eqdsk_fileName
        logical, intent(IN), optional :: plot
        real, allocatable :: yp(:), xp(:), temp(:), ss(:), yp_c(:)
        real :: spl1 = 0.0, spln = 0.0, t, curv2, sigma = 0.0

        !!   Plotting variables
    
        !integer :: nLevs
        !real :: levStep, lev

        !   Read in variables from geqdsk file

        open ( unit = 8, file = eqdsk_fileName, status = 'OLD' )

        read ( 8, 2000 ) ( case_ (i), i=1, 6 ), idum, nw, nh 
        read ( 8, 2020 ) rdim,zdim,rcentr,rleft,zmid 
        read ( 8, 2020 ) rmaxis__,zmaxis__,simag,sibry,bcentr 
        read ( 8, 2020 ) current__,simag,xdum,rmaxis__,xdum 
        read ( 8, 2020 ) zmaxis__,xdum,sibry,xdum,xdum 
       
        write(*,*) idum, nw, nh 
        allocate ( fpol ( nw ), pres ( nw ), ffprim ( nw ), &
            pprime ( nw ), psizr ( nw, nh ), qpsi ( nw ), &
            r ( nw ), z ( nh ), fluxGrid ( nw ) )
        
        read ( 8, 2020 ) ( fpol (i), i=1, nw ) 
        read ( 8, 2020 ) ( pres (i), i=1, nw ) 
        read ( 8, 2020 ) ( ffprim (i), i=1, nw ) 
        read ( 8, 2020 ) ( pprime (i), i=1, nw ) 
        read ( 8, 2020 ) ( ( psizr (i,j), i=1, nw ), j=1, nh ) 
        read ( 8, 2020 ) ( qpsi (i), i=1, nw ) 
        
        read ( 8, 2022 ) nbbbs,limitr 
        
        allocate ( rbbbs ( nbbbs ), zbbbs ( nbbbs ), &
            rlim__ ( limitr ), zlim ( limitr ) )
        
        read ( 8, 2020 ) ( rbbbs (i), zbbbs (i), i=1,nbbbs ) 
        read ( 8, 2020 ) ( rlim__ (i), zlim (i), i=1,limitr ) 
         
        2000 format (6a8,3i4) 
        2020 format (5e16.9)
        2022 format (2i5) 
        
        close ( unit = 8 )
        
        !   Calculate other required variables
       
        rStep   = rdim / ( nw - 1 )
        zStep   = zdim / ( nh - 1 )
        fStep   = ( sibry - simag ) / ( nw - 1 )

        r   = (/ (i,i=0,nw-1) /) * rStep + rleft
        z   = (/ (i,i=0,nh-1) /) * zStep + zmid - zdim / 2.0

        fluxGrid    = (/ (i,i=0,nW-1) /) * fStep + simag

        allocate ( bR ( nw, nh ), bz__ ( nw, nh ), &
            bPhi ( nw, nh ), bMag__(nw,nh) )

        bR  = dlg_pDeriv ( psizr, 2, zStep )
        bz__  = -1.0 * dlg_pDeriv ( psizr, 1, rStep )

        !   Remember psi = - R * A

        do i=1,nw
            do j=1,nh
                bR(i,j)  = -bR(i,j) / r(i)
                bz__(i,j)  = -bz__(i,j) / r(i)
            enddo
        enddo
       
        allocate ( temp(4*nw), fpolRZ(nw,nh), yp_c(nw) ) 

        flux_grid_direction: &
        if ( fluxGrid(1) > fluxGrid(2) ) then 

            call curv1 ( nw, -fluxGrid, fpol, spl1, spln, 3, yp_c, temp, sigma, iErr )
            if ( iErr .ne. 0 ) then 
                    write(*,*) 'eqdsk_dlg.f90 [103]: curv1 error', iErr
                    stop
            endif
            do i=1,nw
                do j=1,nh
                    
                    !   curv2 evaluates the spline (fitpack.f)
                    !t   =  ( psizr(i,j) - simag ) / ( sibry - simag )
                    fPolRZ(i,j) = curv2 ( -psizr(i,j), nw, -fluxGrid, fpol, yp_c, sigma )
                    bPhi(i,j)   = fpolRZ(i,j) / r(i)

                end do
            end do

        else

            call curv1 ( nw, fluxGrid, fpol, spl1, spln, 3, yp_c, temp, sigma, iErr )
            if ( iErr .ne. 0 ) then 
                    write(*,*) 'eqdsk_dlg.f90 [103]: curv1 error', iErr
                    stop
            endif
            do i=1,nw
                do j=1,nh

                    fPolRZ(i,j) = curv2 ( psizr(i,j), nw, fluxGrid, fpol, yp_c, sigma )
                    bPhi(i,j)   = fpolRZ(i,j) / r(i)

                end do
            end do

        endif flux_grid_direction
 
        deallocate ( temp, yp_c ) 
     
!        !   Test the fitpack interpolation
!
!        do i=1,nw
!            
!            t   =  ( fluxGrid(i) - simag ) / ( sibry - simag )
!!            call kurv2 ( t, xs, ys, nw, fluxGrid, fpol, xp, yp, ss, 0.0 )
!            write (*,*) i, t, fluxGrid(i), fpol(i), xs, ys, &
!                curv2 ( fluxGrid(i), nw, fluxGrid, fpol, yp_c, 0.0 )
!         
!        end do

        bMag__    = sqrt ( bR**2 + bPhi**2 + bz__**2 )

        !!   Plotting

        !if ( present ( plot ) ) then
        !    if ( plot ) then 
        !        
        !        call setFil ( 'eqdsk.eps' ) 
        !        call setPag ( 'DA4P' )
        !        call metaFl ( 'EPS' )
        !        call disIni ()
        !        call graf ( 0.0, 3.0, 0.0, 0.5, -2.0, 2.0, -2.0, 0.5 ) 
        !        nLevs   = 101 
        !        levStep    = ( maxVal ( psizr ) - minVal ( psizr ) ) / nLevs
        !        do i=0,nLevs-1
        !            lev    = i * levStep - nLevs/2*levStep
        !            call contur ( r, nw, z, nh, psizr, lev ) 
        !        end do
        !        call endGrf ()  
        !        call disFin ()

        !    end if 
        !end if

    end subroutine read_geqdsk

    function is_inside_bbbs ( rIn, zIn )

        implicit none
        logical :: is_inside_bbbs
        real, intent(in) :: rIn, zIn
        integer :: q1, q2, q3, q4

        q1  = count ( rIn - rbbbs > 0 .and. zIn - zbbbs > 0 )
        q2  = count ( rIn - rbbbs > 0 .and. zIn - zbbbs .le. 0 )
        q3  = count ( rIn - rbbbs .le. 0 .and. zIn - zbbbs > 0 )
        q4  = count ( rIn - rbbbs .le. 0 .and. zIn - zbbbs .le. 0 )

        if ( q1 > 0 .and. q2 > 0 .and. q3 > 0 .and. q4 > 0 ) then

           is_inside_bbbs    = .true. 

        else

           is_inside_bbbs   = .false.

        endif
       
        return

    end function is_inside_bbbs

    function is_inside_lim ( rIn, zIn )

        implicit none
        logical :: is_inside_lim
        real, intent(in) :: rIn, zIn
        integer :: q1, q2, q3, q4
        real, allocatable :: newR(:), newZ(:)
        integer :: nInterp, nLim, i, j, cnt
        real :: m, b, d, dStep

        ! create an interpolated limiter boundary

        nInterp = 100
        nLim    = size ( rLim__ )
        cnt = 1

        allocate ( newR(nInterp*nLim+nLim), &
            newZ(nInterp*nLim+nLim) )
    
        newR(1:nLim)    = rLim__
        newZ(1:nLim)    = zLim

        do i=1,nLim-1

            !   get slope

            dStep   = (rLim__(i+1) - rLim__(i)) / nInterp

            if (abs(dStep) .gt. 0 ) then 

                m   = ( zLim(i+1)-zLim(i) ) &
                        / ( rLim__(i+1) - rLim__(i) )
                b   = zLim(i) - m * rLim__(i)

                !   distance

                d   = sqrt ( ( rLim__(i+1) - rLim__(i) )**2 &
                        + ( zLim(i+1) - zLim(i) )**2 )

                do j = 1, nInterp

                        newR(nLim+cnt)  = rLim__(i) + dStep*j
                        newZ(nLim+cnt)  = m * (rLim__(i) + dStep*j) + b 
                        cnt = cnt + 1

                enddo

            endif

        enddo

        !   interpolated eqdsk limite boundary

        q1  = count ( rIn - newR(1:nLim+cnt-1) .ge. 0 .and. zIn - newZ(1:nLim+cnt-1) .ge. 0 )
        q2  = count ( rIn - newR(1:nLim+cnt-1) .ge. 0 .and. zIn - newZ(1:nLim+cnt-1) .le. 0 )
        q3  = count ( rIn - newR(1:nLim+cnt-1) .le. 0 .and. zIn - newZ(1:nLim+cnt-1) .ge. 0 )
        q4  = count ( rIn - newR(1:nLim+cnt-1) .le. 0 .and. zIn - newZ(1:nLim+cnt-1) .le. 0 )

        !!   coarse eqdsk limiter boundary

        !q1  = count ( rIn - rLim__ .ge. 0 .and. zIn - zLim .ge. 0 )
        !q2  = count ( rIn - rLim__ .ge. 0 .and. zIn - zLim .le. 0 )
        !q3  = count ( rIn - rLim__ .le. 0 .and. zIn - zLim .ge. 0 )
        !q4  = count ( rIn - rLim__ .le. 0 .and. zIn - zLim .le. 0 )

        if ( q1 > 0 .and. q2 > 0 .and. q3 > 0 .and. q4 > 0 ) then

           is_inside_lim    = .true. 

        else

           is_inside_lim   = .false.

        endif
      
        deallocate ( newR, newZ )

        return

    end function is_inside_lim

    function is_inside_domain ( rIn, zIn )

        use aorsa2din_mod, only: rwLeft, rwRight, yTop, yBottom

        implicit none
        logical :: is_inside_domain
        real, intent(in) :: rIn, zIn
        integer :: q1, q2, q3, q4
        real, allocatable :: newR(:), newZ(:)
        integer :: nInterp, nLim, i, j, cnt
        real :: m, b, d, dStep
        real :: rDomain(5), zDomain(5)

        ! create an interpolated limiter boundary

        rDomain = (/ rwLeft, rwRight, rwRight, rwLeft, rwLeft /)
        zDomain = (/ yBottom, yBottom, yTop, yTop, yBottom /)

        nInterp = 100
        nLim    = size ( rDomain )
        cnt = 1

        allocate ( newR(nInterp*nLim+nLim), &
            newZ(nInterp*nLim+nLim) )
    
        newR(1:nLim)    = rDomain 
        newZ(1:nLim)    = zDomain

        do i=1,nLim-1

            !   get slope

            dStep   = (rDomain(i+1) - rDomain(i)) / nInterp

            if (abs(dStep) .gt. 0 ) then 

                m   = ( zDomain(i+1)-zDomain(i) ) &
                        / ( rDomain(i+1) - rDomain(i) )
                b   = zDomain(i) - m * rDomain(i)

                !   distance

                d   = sqrt ( ( rDomain(i+1) - rDomain(i) )**2 &
                        + ( zDomain(i+1) - zDomain(i) )**2 )

                do j = 1, nInterp

                        newR(nLim+cnt)  = rDomain(i) + dStep*j
                        newZ(nLim+cnt)  = m * (rDomain(i) + dStep*j) + b 
                        cnt = cnt + 1

                enddo

            endif

        enddo

        !   interpolated eqdsk limite boundary

        q1  = count ( rIn - newR(1:nLim+cnt-1) .ge. 0 .and. zIn - newZ(1:nLim+cnt-1) .ge. 0 )
        q2  = count ( rIn - newR(1:nLim+cnt-1) .ge. 0 .and. zIn - newZ(1:nLim+cnt-1) .le. 0 )
        q3  = count ( rIn - newR(1:nLim+cnt-1) .le. 0 .and. zIn - newZ(1:nLim+cnt-1) .ge. 0 )
        q4  = count ( rIn - newR(1:nLim+cnt-1) .le. 0 .and. zIn - newZ(1:nLim+cnt-1) .le. 0 )

        if ( q1 > 0 .and. q2 > 0 .and. q3 > 0 .and. q4 > 0 ) then

           is_inside_domain = .true. 

        else

           is_inside_domain = .false.

        endif
      
        deallocate ( newR, newZ )

        return

    end function is_inside_domain



end module eqdsk_dlg
