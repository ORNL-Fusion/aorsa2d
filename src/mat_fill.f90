module mat_fill

implicit none

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
complex, allocatable, dimension(:) :: &
    sss, ttt, qqq
complex, allocatable :: aMat(:,:)

contains

    subroutine amat_fill

        use aorsa2din_mod, &
        only: nModesX, nModesY, upShift, damping, &
            lmax, delta0, nPhi, nSpec, nZFun, &
            iSigma, xnuomg
        use grid
        use sigma_mod
        use rotation
        use constants
        use profiles
        use bField

        implicit none

        integer :: iRow, iCol, i, j, n, m, p, s


        write(*,*) 'Filling aMat, size: ', &
            nModesX*nModesY*3*nkx*nky*3*2*8.0 / 1024.0**2

        allocate ( &
            sss(nModesX*nModesY*3), &
            ttt(nModesX*nModesY*3), &
            qqq(nModesX*nModesY*3) )

        allocate ( aMat(nModesX*nModesY*3,nkx*nky*3) ) 

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

                        !write(*,*) i, j, m, n, kxL, kxR, kyL, kyR

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
                                mSpec(s), densitySpec(i,j,s), xnuomg, &
                                ktSpec(i,j,s), omgc(i,j,s), omgp2(i,j,s), &
                                -lmax, lmax, nzfun, &
                                xkxsav(n), xkysav(m), nphi, capr(i), &
                                sigxxTmp, sigxyTmp, sigxzTmp, &
                                sigyxTmp, sigyyTmp, sigyzTmp, &
                                sigzxTmp, sigzyTmp, sigzzTmp, &
                                delta0, 0, omgrf, xk0, &
                                upshift, damping, xk_cutoff )
                            
                            if (iSigma==0) & ! cold plasma 
                            call sigmaCold_stix(i, j, &
                                xnuomg, &
                                omgc(i,j,s), omgp2(i,j,s), &
                                xkxsav(n), xkysav(m), nphi, capr(i), &
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

                        !if ( i==1 .or. i==nModesX &
                        !        .or. j==1 .or. j==nModesY ) then

                        !    sigxx = complex ( 0.0, 1e18 ) 
                        !    sigxy = 0 
                        !    sigxz = 0 
                        !                   
                        !    sigyx = 0 
                        !    sigyy = complex ( 0.0, 1e18 )
                        !    sigyz = 0 
                        !                   
                        !    sigzx = 0 
                        !    sigzy = 0 
                        !    sigzz = complex ( 0.0, 1e18 )

                        !endif

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


                        fdk = cexpkxky !* dxx  
                        fek = cexpkxky !* dxy  
                        ffk = cexpkxky !* dxz  
                                       !       
                        fgk = cexpkxky !* dyx  
                        fak = cexpkxky !* dyy  
                        fpk = cexpkxky !* dyz  
                                       !       
                        frk = cexpkxky !* dzx  
                        fqk = cexpkxky !* dzy  
                        fsk = cexpkxky !* dzz  

                        !!   boundary conditions
                        !!   -------------------

                        !if ( i==1 .or. i==nModesX &
                        !        .or. j==1 .or. j==nModesY ) then

                        !    fdk = cExpKxKy
                        !    fek = 0
                        !    ffk = 0

                        !    fgk = 0
                        !    fak = cExpKxKy
                        !    fpk = 0

                        !    frk = 0
                        !    fqk = 0
                        !    fsk = cExpKxKy
                        !
                        !endif

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

    end subroutine amat_fill

end module mat_fill
