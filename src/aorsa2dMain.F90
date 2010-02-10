program aorsa2dMain
    
    use constants
    use eqdsk_dlg
    use aorsasubs_mod
    use sigma_mod
    use aorsa2din_mod
    use interp
        
    implicit none

!   Variable list
!   -------------

    real :: omgrf, xk0
    real, allocatable, dimension(:) :: mSpec, qSpec
    integer, allocatable, dimension(:) :: zSpec, amuSpec
    real, allocatable, dimension(:,:,:) :: omgc, omgp2, densitySpec, ktSpec
    integer :: nSpec
    real, allocatable, dimension(:,:) :: btau
    real :: sqx
    real, allocatable, dimension(:,:) :: uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz
    real :: dx, dy, xRange, yRange
    real, allocatable, dimension(:) :: capR, xkphi
    real, allocatable, dimension(:) :: y
    integer :: i, j, m, n, s
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


!   read namelist input data
!   ------------------------

    call read_nameList ()


!   define x mesh: capr(i)
!   --------------------------------------------

    allocate ( &
        capR ( nmodesx ), &
        xkphi ( nmodesx ) )

    rwLeft   = 1.0
    rwRight  = 2.0
    xRange  = rwRight - rwLeft
    dx = xRange / nModesX
    do i = 1, nModesX

        capr(i) = (i-1) * dx + xLeft 
        xkphi(i) = nphi / capr(i)

    enddo


!   define y mesh: y(j), yprime(j)
!---------------------------------

    allocate ( y ( nModesY ) )

    yTop    =  1.0
    yBot    = -1.0
    yRange  = yTop - yBot
    dy = yRange / nmodesy
    do j = 1, nModesY

        y(j) = (j-1) * dy + yBot

    enddo


!   read g-eqdsk file
!   -----------------

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
           bxn(i,j) = bHere(1)
           byn(i,j) = bHere(2)
           bzn(i,j) = bHere(3)

        enddo
    enddo

!   calculate the ion cyclotron freqs
!   ---------------------------------

    omgrf = 2.0 * pi * freqcy
    xk0 = omgrf / clight

    nSpec   = 3

    allocate ( &
        mSpec(nSpec), zSpec(nSpec), &
        qSpec(nSpec), amuSpec(nSpec) )

    amuSpec     = (/ 0, amu1, amu2 /) 
    mSpec       = amuSpec * xmh
    mSpec(1)    = xme  
    zSpec       = (/ -1, z1, z2 /)
    qSpec       = zSpec * q 

    allocate ( &
        densitySpec ( nModesX, nModesY, nSpec ), & 
        ktSpec ( nModesX, nModesY, nSpec ) )

    ktSpec(:,:,1)   = te0 * q 
    ktSpec(:,:,2)   = ti01 * q 
    ktSpec(:,:,3)   = ti02 * q 

    densitySpec(:,:,1)   = xn0
    densitySpec(:,:,2)   = xn1 
    densitySpec(:,:,3)   = xn2 

    allocate ( &
        omgc ( nModesX, nModesY, nSpec ), &
        omgp2 ( nModesX, nModesY, nSpec ) )

    do i=1,nModesX
        do j=1,nModesY

            omgc(i,j,:) = qSpec * b0 / mSpec
            omgp2(i,j,:)    = densitySpec(i,j,:) * qSpec**2 / ( eps0 * mSpec )

        enddo
    enddo
    
!   calculate rotation matrix U
!   ---------------------------

    allocate ( btau ( nmodesx, nmodesy ) )

    allocate ( &
        uxx( nmodesx, nmodesy ), & 
        uxy( nmodesx, nmodesy ), &
        uxz( nmodesx, nmodesy ), &
        uyx( nmodesx, nmodesy ), & 
        uyy( nmodesx, nmodesy ), &
        uyz( nmodesx, nmodesy ), &
        uzx( nmodesx, nmodesy ), & 
        uzy( nmodesx, nmodesy ), &
        uzz( nmodesx, nmodesy ) )

    do i = 1, nmodesx
        do j = 1, nmodesy

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



!   Calculate the integrated volume on even mesh:
!   --------------------------------------------

    kxL = -nModesx/2
    kxR =  nModesx/2
    kyL = -nModesY/2+1
    kyR =  nModesY/2

    allocate ( &
        xkxsav (kxL:kxR), &
        xkysav (kyL:kyR) )

    do m = kxL, kxR 
        xkxsav(m) = 2.0 * pi * m / xRange
    enddo

    do n = kyL, kyR 
        xkysav(n) = 2.0 * pi * n / yRange
    enddo

    xk_cutoff   = sqrt ( xkxsav( nmodesx/2 )**2 &
                            + xkysav( nmodesy/2 )**2 ) * xkperp_cutoff


!   Take numerical derivatives
!   --------------------------

    allocate ( gradPrlB (nmodesx,nmodesy) )

    allocate ( &  
        dxuxx(nmodesx,nmodesy), dxxuxx(nmodesx,nmodesy), &
        dxuxy(nmodesx,nmodesy), dxxuxy(nmodesx,nmodesy), &
        dxuxz(nmodesx,nmodesy), dxxuxz(nmodesx,nmodesy), &
        dxuyx(nmodesx,nmodesy), dxxuyx(nmodesx,nmodesy), &
        dxuyy(nmodesx,nmodesy), dxxuyy(nmodesx,nmodesy), &
        dxuyz(nmodesx,nmodesy), dxxuyz(nmodesx,nmodesy), &
        dxuzx(nmodesx,nmodesy), dxxuzx(nmodesx,nmodesy), &
        dxuzy(nmodesx,nmodesy), dxxuzy(nmodesx,nmodesy), &
        dxuzz(nmodesx,nmodesy), dxxuzz(nmodesx,nmodesy) )

    allocate ( &  
        dyuxx(nmodesx,nmodesy), dyyuxx(nmodesx,nmodesy), &
        dyuxy(nmodesx,nmodesy), dyyuxy(nmodesx,nmodesy), &
        dyuxz(nmodesx,nmodesy), dyyuxz(nmodesx,nmodesy), &
        dyuyx(nmodesx,nmodesy), dyyuyx(nmodesx,nmodesy), &
        dyuyy(nmodesx,nmodesy), dyyuyy(nmodesx,nmodesy), &
        dyuyz(nmodesx,nmodesy), dyyuyz(nmodesx,nmodesy), &
        dyuzx(nmodesx,nmodesy), dyyuzx(nmodesx,nmodesy), &
        dyuzy(nmodesx,nmodesy), dyyuzy(nmodesx,nmodesy), &
        dyuzz(nmodesx,nmodesy), dyyuzz(nmodesx,nmodesy) )

    allocate ( &
        dxyuxx(nmodesx,nmodesy), dxyuxy(nmodesx,nmodesy), dxyuxz(nmodesx,nmodesy), &
        dxyuyx(nmodesx,nmodesy), dxyuyy(nmodesx,nmodesy), dxyuyz(nmodesx,nmodesy), &
        dxyuzx(nmodesx,nmodesy), dxyuzy(nmodesx,nmodesy), dxyuzz(nmodesx,nmodesy) )
                
    do i = 1, nmodesx
        do j = 1, nmodesy

            !   Brambilla approximation:
            !   -----------------------

            sinTh = y(j) / sqrt ( (capR(i)-rmaxis__)**2 + y(j)**2 )
            gradprlb(i,j) = bmod(i,j) / capr(i) * abs ( btau(i,j) * sinTh )

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


!   precompute xx(n,i), yy(m,j)
!   ---------------------------

    allocate ( &
        xx(kxL:kxR,nModesX), xx_inv(kxL:kxR,nModesX), &
        yy(kyL:kyR,nModesY), yy_inv(kyL:kyR,nModesY) )

    do i = 1, nmodesx
        do m = kxL, kxR 
            xx(m, i) = exp(zi * xkxsav(m) * capR(i))
            xx_inv(m,i) = 1.0/xx(m,i)
        enddo
    enddo

    do j = 1, nmodesy
        do n = kyL, kyR 
            yy(n,j) = exp(zi * xkysav(n) * y(j))
            yy_inv(n,j) = 1.0/yy(n,j)
        enddo
    enddo



!   Load x, y and z equations for spatial point (i,j) and mode number (n,m)
!   ------------------------------------------------------------------------

    i_loop: &
    do i=1,nmodesx
        j_loop: &
        do j=1,nmodesy
            m_loop: &
            do m=kxL,kxR
                n_loop: &
                do n=kyL,kyR

                write(*,*) i, j, m, n

                cexpkxky = xx(m, i) * yy(n, j)

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

                hot_plasma: &
                if (isigma .eq. 1)then

                    do s=1,nSpec

                        call sigma_maxwellian(i, j, m, n, &
                            gradprlb(i,j), bmod(i,j), &
                            mSpec(s), densitySpec(i,j,s), xnuomg, &
                            ktSpec(i,j,s), omgc(i,j,s), omgp2(i,j,s), &
                            -lmax, lmax, nzfun, &
                            xkxsav(m), xkysav(n), nphi, capr(i), &
                            bxn(i,j), byn(i,j), bzn(i,j), &
                            uxx(i,j), uxy(i,j), uxz(i,j), &
                            uyx(i,j), uyy(i,j), uyz(i,j), &
                            uzx(i,j), uzy(i,j), uzz(i,j), &
                            sigxxTmp, sigxyTmp, sigxzTmp, &
                            sigyxTmp, sigyyTmp, sigyzTmp, &
                            sigzxTmp, sigzyTmp, sigzzTmp, &
                            delta0, 0, omgrf, xk0, &
                            upshift, damping, xk_cutoff )

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

                endif hot_plasma

                cold_plasma: &
                if (isigma .eq. 0)then 

                    do s=1,nSpec

                        call sigmac_stix(i, j, m, n, &
                            mSpec(s), densitySpec(i,j,s), xnuomg, &
                            ktSpec(i,j,s), omgc(i,j,s), omgp2(i,j,s), &
                            -lmax, lmax, nzfun, ibessel, &
                            xkxsav(m), xkysav(n), nphi, capr(i), &
                            bxn(i,j), byn(i,j), bzn(i,j), &
                            uxx(i,j), uxy(i,j), uxz(i,j), &
                            uyx(i,j), uyy(i,j), uyz(i,j), &
                            uzx(i,j), uzy(i,j), uzz(i,j), &
                            sigxxTmp, sigxyTmp, sigxzTmp, &
                            sigyxTmp, sigyyTmp, sigyzTmp, &
                            sigzxTmp, sigzyTmp, sigzzTmp, &
                            delta0, omgrf, xk0 )

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

                endif cold_plasma



                xkxx = 1.0 + zi / (eps0 * omgrf) * sigxx
                xkxy =       zi / (eps0 * omgrf) * sigxy
                xkxz =       zi / (eps0 * omgrf) * sigxz

                xkyx =       zi / (eps0 * omgrf) * sigyx
                xkyy = 1.0 + zi / (eps0 * omgrf) * sigyy
                xkyz =       zi / (eps0 * omgrf) * sigyz

                xkzx =       zi / (eps0 * omgrf) * sigzx
                xkzy =       zi / (eps0 * omgrf) * sigzy
                xkzz = 1.0 + zi / (eps0 * omgrf) * sigzz

                rnx = xkxsav(m) / xk0
                rny = xkysav(n) / xk0
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

                enddo n_loop
            enddo m_loop 
        enddo j_loop
    enddo i_loop 


!!   Antenna current
!!   ---------------
!
!    xant    = rant - rmaxis
!    iant    = int((rant - rwleft) / dx) + 1
!
!    !   note curden is in Amps per meter of toroidal length (2.*pi*rt).
!
!    xjantx = curdnx / dx
!    xjanty = curdny / dx
!    xjantz = curdnz / dx
!    xjant=sqrt(xjantx**2 + xjanty**2 + xjantz**2)
!
!    theta_antr = theta_ant / 180. * pi
!
!    dpsiant = dpsiant0
!    dthetant = dthetant0 / 360. * 2.0 * pi
!    yant_max = antlen / 2.0
!
!    do i = 1, nmodesx
!        do j = 1, nmodesy
!
!            delta_theta = theta0(i,j) - theta_antr
!            gaussantth = exp(-delta_theta**2 / dthetant**2)
!            gausspsi = exp(-(psi(i,j) - psiant)**2 / dpsiant**2)
!
!            shapey = 0.0
!            deltay = y(j) - yant
!
!            if(capr(i) .gt. rmaxis) then
!
!                !   if(i_antenna .eq. 1) antenna current is cos(ky * y)  (DEFAULT)
!                !   ------------------------------------------------ --------------
!
!                if (i_antenna .eq. 1 .and. abs(deltay) .lt. yant_max)then
!                    shapey = cos(xk0 * antlc * deltay)
!                endif
!
!            endif
!
!            xjx(i,j) = 0.0
!            xjy(i,j) = xjanty * shapey  * gausspsi
!            xjz(i,j) = 0.0
!
!       enddo
!    enddo
!
!
!
!!     Fourier transform the Stix electric fields:
!!     ------------------------------------------
!
!      call sft2d_parallel(ealpha, xmax, ymax, nmodesx, nmodesy, &
!         nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2, &
!         -nmodesx /2, nmodesx /2, -nmodesy /2 + 1, nmodesy /2, xx_inv, yy_inv, ealphak, dx, dy, &
!         myid, nproc, icontxt)
!
!      call sft2d_parallel(ebeta, xmax, ymax, nmodesx, nmodesy, &
!         nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2, &
!         -nmodesx /2, nmodesx /2, -nmodesy /2 + 1, nmodesy /2, xx_inv, yy_inv, ebetak, dx, dy, &
!         myid, nproc, icontxt)
!
!      call sft2d_parallel(eb, xmax, ymax, nmodesx, nmodesy, &
!         nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2, &
!         -nmodesx /2, nmodesx /2, -nmodesy /2 + 1, nmodesy /2, xx_inv, yy_inv, ebk, dx, dy, &
!         myid, nproc, icontxt)
!
!
!
!      do n = -nmodesx /2, nmodesx /2
!         do m = -nmodesy /2 + 1, nmodesy /2
!
!            ealphakmod(n, m) = sqrt(conjg(ealphak(n, m))* ealphak(n, m))
!            ebetakmod(n, m) = sqrt(conjg(ebetak(n, m)) * ebetak(n, m))
!            ebkmod(n, m) = sqrt(conjg(ebk(n, m)) * ebk(n, m))
!
!         end do
!      end do
!
!
!
!!     ----------------------------------------------
!!     Calculate E in the Lab frame and eplus, eminus
!!     ----------------------------------------------
!      isq2 = SQRT(0.5)
!      do i = 1, nmodesx
!         do j = 1, nmodesy
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


