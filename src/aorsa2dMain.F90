program aorsa2dMain
    
    use constants
    use eqdsk_dlg
    use aorsasubs_mod
    use sigma_mod
    use aorsa2din_mod
    use interp
    use fourier_mod
    use write_data
        
    implicit none

!   Variable list
!   -------------

    real :: omgrf, xk0
    real, allocatable, dimension(:) :: mSpec, qSpec, tSpec, dSpec
    integer, allocatable, dimension(:) :: zSpec, amuSpec
    real, allocatable, dimension(:,:,:) :: omgc, omgp2, densitySpec, ktSpec
    integer :: nSpec
    real, allocatable, dimension(:,:) :: btau
    real :: sqx
    real, allocatable, dimension(:,:) :: uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz
    real :: dx, dy, xRange, yRange
    real, allocatable, dimension(:) :: capR, xkphi
    real, allocatable, dimension(:) :: y
    integer :: i, j, m, n, s, iRow, iCol, p
    integer :: nkx, nky, nRow, nCol
    real, allocatable, dimension(:) :: xkxsav, xkysav
    real :: xk_cutOff
    real :: sinTh
    real, allocatable, dimension(:,:) :: gradPrlB
    real, allocatable, dimension(:,:) :: &
        dxuxx, dxxuxx, dxuxy, dxxuxy, dxuxz, dxxuxz, &
        dxuyx, dxxuyx, dxuyy, dxxuyy, dxuyz, dxxuyz, &
        dxuzx, dxxuzx, dxuzy, dxxuzy, dxuzz, dxxuzz
    real, allocatable, dimension(:,:) :: &
        dyuxx, dyyuxx, dyuxy, dyyuxy, dyuxz, dyyuxz, &
        dyuyx, dyyuyx, dyuyy, dyyuyy, dyuyz, dyyuyz, &
        dyuzx, dyyuzx, dyuzy, dyyuzy, dyuzz, dyyuzz
    real, allocatable, dimension(:,:) :: &
            dxyuxx, dxyuxy, dxyuxz, &
            dxyuyx, dxyuyy, dxyuyz, &
            dxyuzx, dxyuzy, dxyuzz  
    integer :: kxL, kxR, kyL, kyR
    complex, allocatable, dimension(:,:) :: xx, xx_inv, yy, yy_inv
    complex :: cexpkxky
    complex :: &
        sigxx, sigxy, sigxz, &
        sigyx, sigyy, sigyz, &
        sigzx, sigzy, sigzz
    complex :: &
        sigxxTmp, sigxyTmp, sigxzTmp, &
        sigyxTmp, sigyyTmp, sigyzTmp, &
        sigzxTmp, sigzyTmp, sigzzTmp
    complex :: &
        xkxx, xkxy, xkxz, &
        xkyx, xkyy, xkyz, &
        xkzx, xkzy, xkzz
    real :: rnx, rny, rnPhi 
    complex :: &
        dxx, dxy, dxz, &
        dyx, dyy, dyz, &
        dzx, dzy, dzz 
    complex :: &
        fdk, fek, ffk, &
        fgk, fak, fpk, &
        frk, fqk, fsk
    real, allocatable, dimension(:,:) :: &
        bMod, bxn, byn, bzn
    real :: bHere(3)
    complex, allocatable, dimension(:) :: &
        sss, ttt, qqq
    complex, allocatable :: aMat(:,:), brhs(:)
    real :: antSigX, antSigY
    complex, allocatable, dimension(:,:) :: &
        xjx, xjy, xjz
    integer :: info
    integer, allocatable, dimension(:) :: ipiv
    complex, allocatable, dimension(:,:) :: &
        ealphak, ebetak, eBk
    complex, allocatable, dimension(:,:) :: &
       ealpha, ealphax, ealphay, &
       ebeta, ebetax, ebetay, &
       eB, eBx, eBy


!   read namelist input data
!   ------------------------

    write(*,*) 'Reading namelist'
    call read_nameList ()


!   define x mesh: capr(i)
!   --------------------------------------------

    allocate ( &
        capR ( nModesX ), &
        xkphi ( nModesX ) )

    rwLeft   = 0.2
    rwRight  = 1.7
    xRange  = rwRight - rwLeft
    dx = xRange / nModesX
    do i = 1, nModesX

        capr(i) = (i-1) * dx + rwLeft 
        xkphi(i) = nphi / capr(i)

    enddo


!   define y mesh: y(j), yprime(j)
!---------------------------------

    allocate ( y ( nModesY ) )

    yTop    =  1.1
    yBot    = -1.1
    yRange  = yTop - yBot
    dy = yRange / nModesY
    do j = 1, nModesY

        y(j) = (j-1) * dy + yBot

    enddo


!   read g-eqdsk file
!   -----------------

    write(*,*) 'Reading eqdsk'
    call read_geqdsk ( eqdsk, plot = .false. )
    call init_interp ()

    allocate ( &
        bMod(nModesX,nModesY), &
        bxn(nModesX,nModesY), &
        byn(nModesX,nModesY), &
        bzn(nModesX,nModesY) )

    do i=1,nModesX
        do j=1,nModesY

           bHere = dlg_interpB ( (/capR(i),0.0,y(j)/), &
                        bMagHere = bMod(i,j) )  
           bxn(i,j) = bHere(1) / bMod(i,j)
           byn(i,j) = bHere(3) / bMod(i,j)
           bzn(i,j) = -bHere(2) / bMod(i,j)

        enddo
    enddo


!   calculate the cyclotron and plasma freqs
!   ----------------------------------------

    omgrf = 2.0 * pi * freqcy
    xk0 = omgrf / clight

    nSpec   = 2

    allocate ( &
        mSpec(nSpec), zSpec(nSpec), &
        qSpec(nSpec), amuSpec(nSpec), &
        dSpec(nSpec), tSpec(nSpec) )

    zSpec       = (/ -1, 2 /)
    amuSpec     = (/ 0, 4 /) 
    tSpec       = (/ 400.0, 400.0 /) ! [eV]
    dSpec       = (/ 0.6e18, 0.3e18 /)

    mSpec       = amuSpec * xmh
    mSpec(1)    = xme  
    qSpec       = zSpec * q 

    allocate ( &
        densitySpec ( nModesX, nModesY, nSpec ), & 
        ktSpec ( nModesX, nModesY, nSpec ) )

    do s=1,nSpec

        ktSpec(:,:,s)       = tSpec(s) * q 
        densitySpec(:,:,s)  = dSpec(s) 

    enddo

    allocate ( &
        omgc ( nModesX, nModesY, nSpec ), &
        omgp2 ( nModesX, nModesY, nSpec ) )

    do i=1,nModesX
        do j=1,nModesY

            omgc(i,j,:) = qSpec * bMod(i,j) / mSpec
            omgp2(i,j,:)    = densitySpec(i,j,:) * qSpec**2 / ( eps0 * mSpec )

        enddo
    enddo


!   calculate rotation matrix U
!   ---------------------------

    write(*,*) 'Building rotation matrix U'
    allocate ( btau ( nModesX, nModesY ) )

    allocate ( &
        uxx( nModesX, nModesY ), & 
        uxy( nModesX, nModesY ), &
        uxz( nModesX, nModesY ), &
        uyx( nModesX, nModesY ), & 
        uyy( nModesX, nModesY ), &
        uyz( nModesX, nModesY ), &
        uzx( nModesX, nModesY ), & 
        uzy( nModesX, nModesY ), &
        uzz( nModesX, nModesY ) )

    do i = 1, nModesX
        do j = 1, nModesY

            btau(i,j) = sqrt(bxn(i,j)**2 + byn(i,j)**2)

            sqx = sqrt(1.0 - bxn(i,j)**2)

            uxx(i, j) =   sqx
            uxy(i, j) = - bxn(i, j) * byn(i, j) / sqx
            uxz(i, j) = - bxn(i, j) * bzn(i, j) / sqx
            uyx(i, j) =   0.0
            uyy(i, j) =   bzn(i, j) / sqx
            uyz(i, j) = - byn(i, j) / sqx
            uzx(i, j) =   bxn(i, j)
            uzy(i, j) =   byn(i, j)
            uzz(i, j) =   bzn(i, j)

        enddo
    enddo


!   Calculate kx and ky values 
!   --------------------------

    kxL = -nModesX/2+1
    kxR =  nModesX/2
    kyL = -nModesY/2+1
    kyR =  nModesY/2

    allocate ( &
        xkxsav (kxL:kxR), &
        xkysav (kyL:kyR) )

    nkx = size ( xkxsav, 1 )
    nky = size ( xkysav, 1 )

    do n = kxL, kxR 
        xkxsav(n) = 2.0 * pi * n / xRange
    enddo

    do m = kyL, kyR 
        xkysav(m) = 2.0 * pi * m / yRange
    enddo

    xk_cutoff   = sqrt ( xkxsav( kxR )**2 &
                            + xkysav( kyR )**2 ) * xkperp_cutoff


!   Take numerical derivatives of the rotation matrix U
!   ---------------------------------------------------

    write(*,*) 'Calculating U derivatives'

    allocate ( gradPrlB (nModesX,nModesY) )

    allocate ( &  
        dxuxx(nModesX,nModesY), dxxuxx(nModesX,nModesY), &
        dxuxy(nModesX,nModesY), dxxuxy(nModesX,nModesY), &
        dxuxz(nModesX,nModesY), dxxuxz(nModesX,nModesY), &
        dxuyx(nModesX,nModesY), dxxuyx(nModesX,nModesY), &
        dxuyy(nModesX,nModesY), dxxuyy(nModesX,nModesY), &
        dxuyz(nModesX,nModesY), dxxuyz(nModesX,nModesY), &
        dxuzx(nModesX,nModesY), dxxuzx(nModesX,nModesY), &
        dxuzy(nModesX,nModesY), dxxuzy(nModesX,nModesY), &
        dxuzz(nModesX,nModesY), dxxuzz(nModesX,nModesY) )

    allocate ( &  
        dyuxx(nModesX,nModesY), dyyuxx(nModesX,nModesY), &
        dyuxy(nModesX,nModesY), dyyuxy(nModesX,nModesY), &
        dyuxz(nModesX,nModesY), dyyuxz(nModesX,nModesY), &
        dyuyx(nModesX,nModesY), dyyuyx(nModesX,nModesY), &
        dyuyy(nModesX,nModesY), dyyuyy(nModesX,nModesY), &
        dyuyz(nModesX,nModesY), dyyuyz(nModesX,nModesY), &
        dyuzx(nModesX,nModesY), dyyuzx(nModesX,nModesY), &
        dyuzy(nModesX,nModesY), dyyuzy(nModesX,nModesY), &
        dyuzz(nModesX,nModesY), dyyuzz(nModesX,nModesY) )

    allocate ( &
        dxyuxx(nModesX,nModesY), dxyuxy(nModesX,nModesY), dxyuxz(nModesX,nModesY), &
        dxyuyx(nModesX,nModesY), dxyuyy(nModesX,nModesY), dxyuyz(nModesX,nModesY), &
        dxyuzx(nModesX,nModesY), dxyuzy(nModesX,nModesY), dxyuzz(nModesX,nModesY) )
                
    do i = 1, nModesX
        do j = 1, nModesY

            !   Brambilla approximation:
            !   -----------------------

            sinTh = y(j) / sqrt ( (capR(i)-rmaxis__)**2 + (y(j)-zmaxis__)**2 )
            gradprlb(i,j) = bMod(i,j) / capr(i) * abs ( btau(i,j) * sinTh )

            if (nzfun == 0) gradPrlB(i,j) = 1.0e-10

            call deriv_x(uxx, i, j, dx, dfdx = dxuxx(i,j), d2fdx2 = dxxuxx(i,j))
            call deriv_x(uxy, i, j, dx, dfdx = dxuxy(i,j), d2fdx2 = dxxuxy(i,j))
            call deriv_x(uxz, i, j, dx, dfdx = dxuxz(i,j), d2fdx2 = dxxuxz(i,j))

            call deriv_x(uyx, i, j, dx, dfdx = dxuyx(i,j), d2fdx2 = dxxuyx(i,j))
            call deriv_x(uyy, i, j, dx, dfdx = dxuyy(i,j), d2fdx2 = dxxuyy(i,j))
            call deriv_x(uyz, i, j, dx, dfdx = dxuyz(i,j), d2fdx2 = dxxuyz(i,j))

            call deriv_x(uzx, i, j, dx, dfdx = dxuzx(i,j), d2fdx2 = dxxuzx(i,j))
            call deriv_x(uzy, i, j, dx, dfdx = dxuzy(i,j), d2fdx2 = dxxuzy(i,j))
            call deriv_x(uzz, i, j, dx, dfdx = dxuzz(i,j), d2fdx2 = dxxuzz(i,j))

            call deriv_y(uxx, i, j, dy, dfdy = dyuxx(i,j), d2fdy2 = dyyuxx(i,j))
            call deriv_y(uxy, i, j, dy, dfdy = dyuxy(i,j), d2fdy2 = dyyuxy(i,j))
            call deriv_y(uxz, i, j, dy, dfdy = dyuxz(i,j), d2fdy2 = dyyuxz(i,j))

            call deriv_y(uyx, i, j, dy, dfdy = dyuyx(i,j), d2fdy2 = dyyuyx(i,j))
            call deriv_y(uyy, i, j, dy, dfdy = dyuyy(i,j), d2fdy2 = dyyuyy(i,j))
            call deriv_y(uyz, i, j, dy, dfdy = dyuyz(i,j), d2fdy2 = dyyuyz(i,j))

            call deriv_y(uzx, i, j, dy, dfdy = dyuzx(i,j), d2fdy2 = dyyuzx(i,j))
            call deriv_y(uzy, i, j, dy, dfdy = dyuzy(i,j), d2fdy2 = dyyuzy(i,j))
            call deriv_y(uzz, i, j, dy, dfdy = dyuzz(i,j), d2fdy2 = dyyuzz(i,j))

            call deriv_xy(uxx, i, j, dx, dy, d2fdxy = dxyuxx(i,j))
            call deriv_xy(uxy, i, j, dx, dy, d2fdxy = dxyuxy(i,j))
            call deriv_xy(uxz, i, j, dx, dy, d2fdxy = dxyuxz(i,j))

            call deriv_xy(uyx, i, j, dx, dy, d2fdxy = dxyuyx(i,j))
            call deriv_xy(uyy, i, j, dx, dy, d2fdxy = dxyuyy(i,j))
            call deriv_xy(uyz, i, j, dx, dy, d2fdxy = dxyuyz(i,j))

            call deriv_xy(uzx, i, j, dx, dy, d2fdxy = dxyuzx(i,j))
            call deriv_xy(uzy, i, j, dx, dy, d2fdxy = dxyuzy(i,j))
            call deriv_xy(uzz, i, j, dx, dy, d2fdxy = dxyuzz(i,j))

       enddo
    enddo


!   precompute basis functions xx(n,i), yy(m,j)
!   -------------------------------------------

    allocate ( &
        xx(kxL:kxR,nModesX), xx_inv(kxL:kxR,nModesX), &
        yy(kyL:kyR,nModesY), yy_inv(kyL:kyR,nModesY) )

    do i = 1, nModesX
        do n = kxL, kxR 
            xx(n, i) = exp(zi * xkxsav(n) * capR(i))
            xx_inv(n,i) = 1.0/xx(n,i)
        enddo
    enddo

    do j = 1, nModesY
        do m = kyL, kyR 
            yy(m,j) = exp(zi * xkysav(m) * y(j))
            yy_inv(m,j) = 1.0/yy(m,j)
        enddo
    enddo



!   Load x, y and z equations for spatial point (i,j) and mode number (n,m)
!   ------------------------------------------------------------------------

    write(*,*) 'Filling aMat, size: ', &
        nModesX*nModesY*3*nkx*nky*3*2*8.0 / 1024.0**2

    allocate ( &
        sss(nModesX*nModesY*3), &
        ttt(nModesX*nModesY*3), &
        qqq(nModesX*nModesY*3) )

    allocate ( aMat(nModesX*nModesY*3,nkx*nky*3) ) 
    allocate ( brhs(nModesX*nModesY*3) )

    i_loop: &
    do i=1,nModesX

        !   progress indicator
        !   ------------------

        do p=1,7 
            write(*,'(a)',advance='no') char(8)
        enddo
        write(*,'(1x,f5.1,a)',advance='no') real(i)/nModesX*100, '%'

        j_loop: &
        do j=1,nModesY

            iRow = (j-1) * 3 + (i-1) * nModesY * 3 + 1

            n_loop: &
            do n=kxL,kxR
                m_loop: &
                do m=kyL,kyR

                    !write(*,*) i, j, m, n

                    iCol = (n-kxL) * 3 * nModesY + (m-kyL) * 3 + 1

                    cexpkxky = xx(n, i) * yy(m, j)

                    sigxx = 0.0
                    sigxy = 0.0
                    sigxz = 0.0
                    sigyx = 0.0
                    sigyy = 0.0
                    sigyz = 0.0
                    sigzx = 0.0
                    sigzy = 0.0
                    sigzz = 0.0


                    !   interior plasma region:
                    !   ----------------------

                    do s=1,nSpec

                        if (iSigma==1) & ! hot plasma        
                        call sigmaHot_maxwellian(i, j, n, m, &
                            gradprlb(i,j), bmod(i,j), &
                            mSpec(s), densitySpec(i,j,s), xnuomg, &
                            ktSpec(i,j,s), omgc(i,j,s), omgp2(i,j,s), &
                            -lmax, lmax, nzfun, &
                            xkxsav(n), xkysav(m), nphi, capr(i), &
                            uxx(i,j), uxy(i,j), uxz(i,j), &
                            uyx(i,j), uyy(i,j), uyz(i,j), &
                            uzx(i,j), uzy(i,j), uzz(i,j), &
                            sigxxTmp, sigxyTmp, sigxzTmp, &
                            sigyxTmp, sigyyTmp, sigyzTmp, &
                            sigzxTmp, sigzyTmp, sigzzTmp, &
                            delta0, 0, omgrf, xk0, &
                            upshift, damping, xk_cutoff )
                        
                        if (iSigma==0) & ! cold plasma 
                        call sigmaCold_stix( &
                            xnuomg, &
                            omgc(i,j,s), omgp2(i,j,s), &
                            xkxsav(n), xkysav(m), nphi, capr(i), &
                            uxx(i,j), uxy(i,j), uxz(i,j), &
                            uyx(i,j), uyy(i,j), uyz(i,j), &
                            uzx(i,j), uzy(i,j), uzz(i,j), &
                            sigxxTmp, sigxyTmp, sigxzTmp, &
                            sigyxTmp, sigyyTmp, sigyzTmp, &
                            sigzxTmp, sigzyTmp, sigzzTmp, &
                            omgrf )

                        sigxx = sigxx + sigxxTmp 
                        sigxy = sigxy + sigxyTmp 
                        sigxz = sigxz + sigxzTmp 
                                       
                        sigyx = sigyx + sigyxTmp 
                        sigyy = sigyy + sigyyTmp 
                        sigyz = sigyz + sigyzTmp 
                                       
                        sigzx = sigzx + sigzxTmp 
                        sigzy = sigzy + sigzyTmp 
                        sigzz = sigzz + sigzzTmp 

                    enddo

                    xkxx = 1.0 + zi / (eps0 * omgrf) * sigxx
                    xkxy =       zi / (eps0 * omgrf) * sigxy
                    xkxz =       zi / (eps0 * omgrf) * sigxz

                    xkyx =       zi / (eps0 * omgrf) * sigyx
                    xkyy = 1.0 + zi / (eps0 * omgrf) * sigyy
                    xkyz =       zi / (eps0 * omgrf) * sigyz

                    xkzx =       zi / (eps0 * omgrf) * sigzx
                    xkzy =       zi / (eps0 * omgrf) * sigzy
                    xkzz = 1.0 + zi / (eps0 * omgrf) * sigzz

                    rnx = xkxsav(n) / xk0
                    rny = xkysav(m) / xk0
                    rnPhi = xkphi(i) / xk0

                    dxx = (xkxx - rny**2 - rnphi**2) * uxx(i,j) &
                        +  xkyx * uyx(i,j) &
                        +  xkzx * uzx(i,j) &
                        + rnx * (rny * uxy(i,j) + rnphi * uxz(i,j)) &
                        - zi * rnphi / xk0 * &
                                 (uxz(i,j) / capr(i) + dxuxz(i,j)) &
                        - zi * rny / xk0 * (dxuxy(i,j) - 2. * dyuxx(i,j)) &
                        - zi * rnx / xk0 * dyuxy(i,j) &
                        + 1. / xk0**2 * (dyyuxx(i,j) - dxyuxy(i,j))

                    dxy =  xkxy * uxx(i,j) &
                        + (xkyy - rny**2 - rnphi**2) * uyx(i,j) &
                        +  xkzy * uzx(i,j) &
                        + rnx * (rny * uyy(i,j) + rnphi * uyz(i,j)) &
                        - zi * rnphi / xk0 * &
                                 (uyz(i,j) / capr(i) + dxuyz(i,j)) &
                        - zi * rny / xk0 * (dxuyy(i,j) - 2. * dyuyx(i,j)) &
                        - zi * rnx / xk0 * dyuyy(i,j) &
                        + 1. / xk0**2 * (dyyuyx(i,j) - dxyuyy(i,j))

                    dxz =  xkxz * uxx(i,j) &
                        +  xkyz * uyx(i,j) &
                        + (xkzz - rny**2 - rnphi**2) * uzx(i,j) &
                        + rnx * (rny * uzy(i,j) + rnphi * uzz(i,j)) &
                        - zi * rnphi / xk0 * &
                                 (uzz(i,j) / capr(i) + dxuzz(i,j)) &
                        - zi * rny / xk0 * (dxuzy(i,j) - 2. * dyuzx(i,j)) &
                        - zi * rnx / xk0 * dyuzy(i,j) &
                        + 1. / xk0**2 * (dyyuzx(i,j) - dxyuzy(i,j))

                    dyx = (xkxx - rnx**2 - rnphi**2) * uxy(i,j) &
                        +  xkyx * uyy(i,j) &
                        +  xkzx * uzy(i,j) &
                        + rny * (rnx * uxx(i,j) + rnphi * uxz(i,j)) &
                        - zi * rny / xk0 * &
                                 (dxuxx(i,j) + uxx(i,j) / capr(i)) &
                        - zi * rnphi / xk0 * dyuxz(i,j) &
                        - zi * rnx / xk0 * &
                        (dyuxx(i,j) - uxy(i,j) / capr(i) - 2.* dxuxy(i,j)) &
                        + 1. / xk0**2 * (dxxuxy(i,j) - dxyuxx(i,j) &
                        - dyuxx(i,j)/ capr(i) + dxuxy(i,j) / capr(i))

                    dyy =  xkxy * uxy(i,j) &
                        + (xkyy - rnx**2 - rnphi**2) * uyy(i,j) &
                        +  xkzy * uzy(i,j) &
                        + rny * (rnx * uyx(i,j) + rnphi * uyz(i,j)) &
                        - zi * rny / xk0 * &
                                 (dxuyx(i,j) + uyx(i,j) / capr(i)) &
                        - zi * rnphi / xk0 * dyuyz(i,j) &
                        - zi * rnx / xk0 * &
                        (dyuyx(i,j) - uyy(i,j) / capr(i) - 2.* dxuyy(i,j)) &
                        + 1. / xk0**2 * (dxxuyy(i,j) - dxyuyx(i,j) &
                        - dyuyx(i,j)/ capr(i) + dxuyy(i,j) / capr(i))

                    dyz =  xkxz * uxy(i,j) &
                        +  xkyz * uyy(i,j) &
                        + (xkzz - rnx**2 - rnphi**2) * uzy(i,j) &
                        + rny * (rnx * uzx(i,j) + rnphi * uzz(i,j)) &
                        - zi * rny / xk0 * &
                                 (dxuzx(i,j) + uzx(i,j) / capr(i)) &
                        - zi * rnphi / xk0 * dyuzz(i,j) &
                        - zi * rnx / xk0 * &
                        (dyuzx(i,j) - uzy(i,j) / capr(i) - 2.* dxuzy(i,j)) &
                        + 1. / xk0**2 * (dxxuzy(i,j) - dxyuzx(i,j) &
                        - dyuzx(i,j)/ capr(i) + dxuzy(i,j) / capr(i))

                    dzx = (xkxx - rnx**2 - rny**2) * uxz(i,j) &
                        +  xkyx * uyz(i,j) &
                        +  xkzx * uzz(i,j) &
                        + rnphi * (rny * uxy(i,j) + rnx * uxx(i,j)) &
                        + zi * rny / xk0 * 2. * dyuxz(i,j) &
                        - zi * rnphi / xk0 * &
                           (dyuxy(i,j) + dxuxx(i,j) - uxx(i,j) / capr(i)) &
                        + zi * rnx / xk0 * &
                           (uxz(i,j) / capr(i) + 2.* dxuxz(i,j)) &
                        - 1. / (xk0**2 * capr(i)) * &
                           (uxz(i,j)/ capr(i) - dxuxz(i,j)) &
                        + 1. / xk0**2  * (dxxuxz(i,j) + dyyuxz(i,j))

                    dzy =  xkxy * uxz(i,j) &
                        + (xkyy - rnx**2 - rny**2) * uyz(i,j) &
                        +  xkzy * uzz(i,j) &
                        + rnphi * (rny * uyy(i,j) + rnx * uyx(i,j)) &
                        + zi * rny / xk0 * 2. * dyuyz(i,j) &
                        - zi * rnphi / xk0 * &
                           (dyuyy(i,j) + dxuyx(i,j) - uyx(i,j) / capr(i)) &
                        + zi * rnx / xk0 * &
                           (uyz(i,j) / capr(i) + 2.* dxuyz(i,j)) &
                        - 1. / (xk0**2 * capr(i)) * &
                           (uyz(i,j)/ capr(i) - dxuyz(i,j)) &
                        + 1. / xk0**2  * (dxxuyz(i,j) + dyyuyz(i,j))

                    dzz =  xkxz * uxz(i,j) &
                        +  xkyz * uyz(i,j) &
                        + (xkzz - rnx**2 - rny**2) * uzz(i,j) &
                        + rnphi * (rny * uzy(i,j) + rnx * uzx(i,j)) &
                        + zi * rny / xk0 * 2. * dyuzz(i,j) &
                        - zi * rnphi / xk0 * &
                           (dyuzy(i,j) + dxuzx(i,j) - uzx(i,j) / capr(i)) &
                        + zi * rnx / xk0 * &
                           (uzz(i,j) / capr(i) + 2.* dxuzz(i,j)) &
                        - 1. / (xk0**2 * capr(i)) * &
                           (uzz(i,j)/ capr(i) - dxuzz(i,j)) &
                        + 1. / xk0**2  * (dxxuzz(i,j) + dyyuzz(i,j))


                    fdk = dxx * cexpkxky
                    fek = dxy * cexpkxky
                    ffk = dxz * cexpkxky

                    fgk = dyx * cexpkxky
                    fak = dyy * cexpkxky
                    fpk = dyz * cexpkxky

                    frk = dzx * cexpkxky
                    fqk = dzy * cexpkxky
                    fsk = dzz * cexpkxky

                    !   boundary conditions
                    !   -------------------

                    if ( i==1 .or. i==nModesX &
                            .or. j==1 .or. j==nModesY ) then

                        fdk = cExpKxKy
                        fek = 0
                        ffk = 0

                        fgk = 0
                        fak = cExpKxKy
                        fpk = 0

                        frk = 0
                        fqk = 0
                        fsk = cExpKxKy
                    
                    endif

                    sss(iCol)   = fdk
                    sss(iCol+1) = fek
                    sss(iCol+2) = ffk

                    ttt(iCol)   = fgk
                    ttt(iCol+1) = fak
                    ttt(iCol+2) = fpk

                    qqq(iCol)   = frk
                    qqq(iCol+1) = fqk
                    qqq(iCol+2) = fsk

                enddo m_loop
            enddo n_loop 

        
            aMat(iRow,:)    = sss
            aMat(iRow+1,:)  = ttt
            aMat(iRow+2,:)  = qqq


        enddo j_loop
    enddo i_loop 
    write(*,*)

!   Antenna current
!   ---------------

    write(*,*) 'Building antenna current (brhs)'

    allocate ( &
        xjx(nModesX,nModesY), &
        xjy(nModesX,nModesY), &
        xjz(nModesX,nModesY) )

    !   note curden is in Amps per meter of toroidal length (2.*pi*rt).

    antSigX = 0.1
    antSigY = 0.2

    do i = 1, nModesX
        do j = 1, nModesY

            xjx(i,j) = 0.0
            xjy(i,j) = 1.0 / dx &
                * exp ( &
                -( (capR(i)-rant)**2/antSigX**2 + (y(j)-0.0)**2/antSigY**2 ) &
                      )
            xjz(i,j) = 0.0

            !   boundary conditions
            !   -------------------

            if ( i==1 .or. i==nModesX &
                    .or. j==1 .or. j==nModesY ) then
                xjx(i,j)    = 0
                xjy(i,j)    = 0
                xjz(i,j)    = 0
            endif

       enddo
    enddo

    xjx = -zi / omgrf / eps0 * xjx
    xjy = -zi / omgrf / eps0 * xjy
    xjz = -zi / omgrf / eps0 * xjz

    do i = 1, nModesX
        do j = 1, nModesY

            iRow = (j-1) * 3 + (i-1) * nModesY * 3 + 1

            brhs(iRow)      = xjx(i,j)
            brhs(iRow+1)    = xjy(i,j)
            brhs(iRow+2)    = xjz(i,j)

        enddo
    enddo

!   Write the run input data to disk
!   --------------------------------

    write(*,*) 'Writing run input data to file'
    call write_runData ( 'runData.nc', &
        capR, y, bxn, byn, bzn, bmod, xjy )


!   Solve complex, linear system
!   ----------------------------

    write(*,*) 'Solving complex linear system'
    
    nRow    = nModesX*nModesY*3
    nCol    = nkx * nky * 3

    allocate ( ipiv ( nRow ) )

    call cgesv ( nRow, 1, aMat, nRow, ipiv, brhs, nRow, info )
            
    write(*,*) '    LAPACK status: ', info


!   Extract the k coefficents from the solution
!   for each field component
!   -------------------------------------------

    allocate ( &
        ealphak(kxL:kxR,kyL:kyR), &
        ebetak(kxL:kxR,kyL:kyR), &
        eBk(kxL:kxR,kyL:kyR) )

    do n=kxL,kxR
        do m=kyL,kyR

            iRow = (m-kyL) * 3 + (n-kxL) * nModesY * 3 + 1
    
            ealphak(n,m)    = brhs(iRow)
            ebetak(n,m)     = brhs(iRow+1)
            eBk(n,m)        = brhs(iRow+2)

        enddo
    enddo


!   Inverse Fourier transform the k coefficents
!   to real space
!   -------------------------------------------

    write(*,*) 'Inverse Fourier transforming the k coeffs'

    call sftinv2d ( xkxsav, xkysav, xx, yy, &
       ealphak, f = ealpha, fx = ealphax, fy = ealphay )
    call sftinv2d ( xkxsav, xkysav, xx, yy, &
       ebetak, f = ebeta, fx = ebetax, fy = ebetay )
    call sftinv2d ( xkxsav, xkysav, xx, yy, &
       eBk, f = eB, fx = eBx, fy = eBy )


!   Write data to file
!   ------------------

    write(*,*) 'Writing solution to file'
    call write_solution ( 'solution.nc', ealpha, ebeta, eB )



!
!!     ----------------------------------------------
!!     Calculate E in the Lab frame and eplus, eminus
!!     ----------------------------------------------
!      isq2 = SQRT(0.5)
!      do i = 1, nModesX
!         do j = 1, nModesY
!
!            ex(i,j)   = uxx(i,j) * ealpha(i,j) &
!                      + uyx(i,j) * ebeta(i,j) &
!                      + uzx(i,j) * eb(i,j)
!
!            ey(i,j)   = uxy(i,j) * ealpha(i,j) &
!                      + uyy(i,j) * ebeta(i,j) &
!                      + uzy(i,j) * eb(i,j)
!
!            ez(i,j)   = uxz(i,j) * ealpha(i,j) &
!                      + uyz(i,j) * ebeta(i,j) &
!                      + uzz(i,j) * eb(i,j)
!
!            eplus(i,j)  = isq2 * (ealpha(i,j) + zi * ebeta(i,j))
!          eminus(i,j) = isq2 * (ealpha(i,j) - zi * ebeta(i,j))
!
!         end do
!      end do
!
end program aorsa2dMain


