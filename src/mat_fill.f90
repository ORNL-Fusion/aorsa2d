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
            iSigma, xnuomg, nPtsX, nPtsY, npRow, npCol
        use grid
        use sigma_mod
        use rotation
        use constants
        use profiles
        use bField
        use parallel

        implicit none

        integer :: iRow, iCol, i, j, n, m, p, s, ii, jj

        !   scalapack indicies
        !   see http://www.netlib.org/scalapack/slug/node76.html

        integer :: l_sp, m_sp, pr_sp, pc_sp, x_sp, y_sp
        integer :: localRow, localCol


        if (iAm == 0) &
        write(*,100), &
            nPtsX*nPtsY*3*nModesX*nModesY*3*2*8.0 / 1024.0**2, &
            nRowLocal*nColLocal*2*8.0 / 1024.0**2
        100 format (' Filling aMat [global size: ',f8.1,' MB, local size: ',f8.1' MB]')

        !allocate ( aMat(nPtsX*nPtsY*3,nModesX*nModesY*3) ) 
        allocate ( aMat(nRowLocal,nColLocal) )
        aMat    = 0

        i_loop: &
        do i=1,nPtsX

            !!   progress indicator
            !!   ------------------

            !do p=1,7 
            !    write(*,'(a)',advance='no') char(8)
            !enddo
            !write(*,'(1x,f5.1,a)',advance='no') real(i)/nPtsX*100, '%'

            j_loop: &
            do j=1,nPtsY

                iRow = (i-1) * 3 * nPtsY + (j-1) * 3 + 1

                n_loop: &
                do n=kxL,kxR
                    m_loop: &
                    do m=kyL,kyR

                        !write(*,*) i, j, m, n, kxL, kxR, kyL, kyR

                        iCol = (n-kxL) * 3 * nModesY + (m-kyL) * 3 + 1

                        ii_loop: &
                        do ii=0,2
                            jj_loop: &
                            do jj=0,2

                                ! scalapack indicies for 2D block cyclic data format
                                ! see http://www.netlib.org/scalapack/slug/node76.html
                                !       and
                                ! http://acts.nersc.gov/scalapack/hands-on/example4.html
 
                                pr_sp   = mod ( rowStartProc + (iRow-1+ii)/rowBlockSize, npRow )
                                pc_sp   = mod ( colStartProc + (iCol-1+jj)/colBlockSize, npCol )

                                myProc: &
                                if ( myRow==pr_sp .and. myCol==pc_sp ) then

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
                                      
                                        !if(iAm==0) &
                                        !write(*,*) xnuomg, omgc(i,j,s), omgp2(i,j,s), &
                                        !    xkxsav(n), xkysav(m), nphi, capr(i), omgrf

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

                                    !if(iAm==0) then
                                    !write(*,*) sigxx, sigxy, sigxz
                                    !write(*,*) sigyx, sigyy, sigyz
                                    !write(*,*) sigzx, sigzy, sigzz
                                    !endif

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

                                    !if (ii==0 .and. jj==0) &
                                    dxx = (xkxx - rny**2 - rnphi**2) * uxx(i,j) &
                                        +  xkyx * uyx(i,j) &
                                        +  xkzx * uzx(i,j) &
                                        + rnx * (rny * uxy(i,j) + rnphi * uxz(i,j)) &
                                        - zi * rnphi / xk0 * &
                                                 (uxz(i,j) / capr(i) + dxuxz(i,j)) &
                                        - zi * rny / xk0 * (dxuxy(i,j) - 2. * dyuxx(i,j)) &
                                        - zi * rnx / xk0 * dyuxy(i,j) &
                                        + 1. / xk0**2 * (dyyuxx(i,j) - dxyuxy(i,j))

                                    !if (ii==1 .and. jj==0) &
                                    dxy =  xkxy * uxx(i,j) &
                                        + (xkyy - rny**2 - rnphi**2) * uyx(i,j) &
                                        +  xkzy * uzx(i,j) &
                                        + rnx * (rny * uyy(i,j) + rnphi * uyz(i,j)) &
                                        - zi * rnphi / xk0 * &
                                                 (uyz(i,j) / capr(i) + dxuyz(i,j)) &
                                        - zi * rny / xk0 * (dxuyy(i,j) - 2. * dyuyx(i,j)) &
                                        - zi * rnx / xk0 * dyuyy(i,j) &
                                        + 1. / xk0**2 * (dyyuyx(i,j) - dxyuyy(i,j))

                                    !if (ii==2 .and. jj==0) &
                                    dxz =  xkxz * uxx(i,j) &
                                        +  xkyz * uyx(i,j) &
                                        + (xkzz - rny**2 - rnphi**2) * uzx(i,j) &
                                        + rnx * (rny * uzy(i,j) + rnphi * uzz(i,j)) &
                                        - zi * rnphi / xk0 * &
                                                 (uzz(i,j) / capr(i) + dxuzz(i,j)) &
                                        - zi * rny / xk0 * (dxuzy(i,j) - 2. * dyuzx(i,j)) &
                                        - zi * rnx / xk0 * dyuzy(i,j) &
                                        + 1. / xk0**2 * (dyyuzx(i,j) - dxyuzy(i,j))

                                    !if (ii==0 .and. jj==1) &
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

                                    !if (ii==1 .and. jj==1) &
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

                                    !if (ii==2 .and. jj==1) &
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

                                    !if (ii==0 .and. jj==2) &
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

                                    !if (ii==1 .and. jj==2) &
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

                                    !if (ii==2 .and. jj==2) &
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

                                    l_sp    = ( iRow-1+ii ) / ( npRow * rowBlockSize )
                                    m_sp    = ( iCol-1+jj ) / ( npCol * colBlockSize )

                                    pr_sp   = mod ( rowStartProc + (iRow-1+ii)/rowBlockSize, npRow )
                                    pc_sp   = mod ( colStartProc + (iCol-1+jj)/colBlockSize, npCol )

                                    x_sp    = mod ( iRow-1+ii, rowBlockSize ) + 1
                                    y_sp    = mod ( iCol-1+jj, colBlockSize ) + 1

                                    localRow    = l_sp*rowBlockSize+x_sp
                                    localCol    = m_sp*colBlockSize+y_sp

                                    !write(*,*) iRow+ii, iCol+jj, l_sp, m_sp, pr_sp, pc_sp, x_sp, y_sp, localRow, localCol

                                    if (ii==0 .and. jj==0) aMat(localRow,localCol) = cexpkxky * dxx  
                                    if (ii==1 .and. jj==0) aMat(localRow,localCol) = cexpkxky * dxy  
                                    if (ii==2 .and. jj==0) aMat(localRow,localCol) = cexpkxky * dxz  
                               
                                    if (ii==0 .and. jj==1) aMat(localRow,localCol) = cexpkxky * dyx  
                                    if (ii==1 .and. jj==1) aMat(localRow,localCol) = cexpkxky * dyy  
                                    if (ii==2 .and. jj==1) aMat(localRow,localCol) = cexpkxky * dyz  
                               
                                    if (ii==0 .and. jj==2) aMat(localRow,localCol) = cexpkxky * dzx  
                                    if (ii==1 .and. jj==2) aMat(localRow,localCol) = cexpkxky * dzy  
                                    if (ii==2 .and. jj==2) aMat(localRow,localCol) = cexpkxky * dzz  

                                    !!   boundary conditions
                                    !!   -------------------

                                    !if ( i==1 .or. i==nPtsX &
                                    !        .or. j==1 .or. j==nPtsY ) then

                                    !    aMat(iRow+0,iCol+0) = cexpkxky  
                                    !    aMat(iRow+1,iCol+0) = 0  
                                    !    aMat(iRow+2,iCol+0) = 0 
                                    !                    
                                    !    aMat(iRow+0,iCol+1) = 0  
                                    !    aMat(iRow+1,iCol+1) = cexpkxky   
                                    !    aMat(iRow+2,iCol+1) = 0  
                                    !                  
                                    !    aMat(iRow+0,iCol+2) = 0  
                                    !    aMat(iRow+1,iCol+2) = 0  
                                    !    aMat(iRow+2,iCol+2) = cexpkxky   
                                    !
                                    !endif

                                endif myProc

                            enddo jj_loop
                        enddo ii_loop

                    enddo m_loop
                enddo n_loop 

            enddo j_loop
        enddo i_loop 

        !if(iAm==0) then
        !write(*,*) 
        !do i=1,nRowLocal
        !    write(*,*) real ( aMat(i,:) )
        !enddo
        !endif
    end subroutine amat_fill

end module mat_fill
