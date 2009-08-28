!
!***************************************************************************
!

      subroutine current_cd(xjpx, xjpy, xjpz, &
         nxdim, nydim, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xm, q, xn, xkt, omgc, omgp2, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         exk, eyk, ezk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndist, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper, dfdupar, &
         UminPara, UmaxPara, UPERP, UPARA, &
         vc_mks, df_cql_uprp, df_cql_uprl, rho, rho_a, nbessj, &
         zeffcd, clight, x, y, rt, b0, ftrap, omgrf, &
         xjpx_ehst, xjpy_ehst, xjpz_ehst, nkperp, &
         zi, eps0, v0i, xk0, kperp_max, i_sav, j_sav, upshift, &
         damping, xk_cutoff,zLoc, eNormIN, mask)

!-----------------------------------------------------------------------
!     This subroutine calculates the plasma current for a single species
!-----------------------------------------------------------------------

        use size_mod

      implicit none

      logical ismine
      real :: eNormIN
      complex zi
      real eps0, v0i, xk0, kperp_max
      integer i_sav, j_sav, upshift
      real damping, xk_cutoff
      !DLG:
      real :: zLoc(nydim)
      integer, intent(in) :: mask(nmodesmax,mmodesmax)

      integer nproc, myid, ngrid, id, ndist, nbessj, nkperp
      integer  i, j, n, m, lmax, nzfun, ibessel, nphi
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx, &
         icontxt, nboundary
      integer dlen_, desc_amat(dlen_)


      integer nxdim, nydim, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, isigma

      integer ftrap
      real zeffcd, clight, rt, b0, damp, r1, xm1, ceta, epsa, xmut2, &
         eta0, xnexp, xjtild, signkz, cfit, c1, afit, yt, bmaxa, vphase, &
         omgrf, c_ehst, akprl, xkprl, xkphi, a, rmaxa, rmina, &
         wphase, xlnlam, vth
      real x(nxdim), y(nydim)


      real xm, omgc(nxdim, nydim), omgp2(nxdim, nydim), &
           psi(nxdim, nydim), psilim
      real capr(nxdim), xnuomg, delta0
      real gradprlb(nxdim, nydim), bmod(nxdim, nydim), &
                               bmod_mid(nxdim, nydim)

      complex xx(nkdim1 : nkdim2, 1 : nxdim), &
              yy(mkdim1 : mkdim2, 1 : nydim)
      complex xjpx(nxdim, nydim), xjpy(nxdim, nydim), xjpz(nxdim, nydim)

      complex xjpx_ehst(nxdim, nydim), &
              xjpy_ehst(nxdim, nydim), &
              xjpz_ehst(nxdim, nydim)


      complex sigxx, sigxy, sigxz, &
              sigyx, sigyy, sigyz, &
              sigzx, sigzy, sigzz


      complex cexpkxky
      complex exk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              eyk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              ezk(nkdim1 : nkdim2, mkdim1 : mkdim2)

      real xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)
      real q, xn(nxdim, nydim), xkt(nxdim, nydim)
      real bxn(nxdim, nydim), byn(nxdim, nydim), bzn(nxdim, nydim)

      real uxx(nxdim, nydim), uxy(nxdim, nydim), uxz(nxdim,nydim), &
           uyx(nxdim, nydim), uyy(nxdim, nydim), uyz(nxdim,nydim), &
           uzx(nxdim, nydim), uzy(nxdim, nydim), uzz(nxdim,nydim)

      real rho(nxdim, nydim)

      integer  :: n_psi_dim, nuper, nupar, n_psi

      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
      real :: UminPara,UmaxPara
      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: vc_mks, rho_a(n_psi_dim)



      if(nphi .gt. 0) signkz =  1.0
      if(nphi .lt. 0) signkz = -1.0
      if(nphi .eq. 0) signkz =  1.0


      ceta = 8.27
      xnexp = 2.48
      cfit = .0987
      afit = 12.3
      damp = 23.83/(0.678 + zeffcd)


      xjpx(:,:) = 0.0
      xjpy(:,:) = 0.0
      xjpz(:,:) = 0.0

      xjpx_ehst(:,:) = 0.0
      xjpy_ehst(:,:) = 0.0
      xjpz_ehst(:,:) = 0.0


!--Loop over mode numbers and mesh:
      do i = 1, nnodex
         xkphi = nphi / capr(i)

         do j = 1, nnodey

            ngrid = (i - 1) * nnodey + j
            id = mod(ngrid, nproc)
            if(id .eq. myid)then

            xjpx(i,j) = 0.0
            xjpy(i,j) = 0.0
            xjpz(i,j) = 0.0

            xjpx_ehst(i,j) = 0.0
            xjpy_ehst(i,j) = 0.0
            xjpz_ehst(i,j) = 0.0

            !if(psi(i,j) .le. psilim .and. nboundary .eq. 1)then
            if(mask(i,j) == 1 .and. nboundary .eq. 1)then


               do n = nkx1, nkx2
                  do m = nky1, nky2

!                    ----------------
!                    Calculate C_ehst
!                    ----------------

!                     xkprl = uzx(i, j) * xkxsav(n)
!     .                    + uzy(i, j) * xkysav(m)
!     .                    + uzz(i, j) * xkphi

                     xkprl = uzz(i, j) * xkphi


                     akprl = abs(xkprl)
                     if(akprl .lt. 1.0e-05)akprl = 1.0e-05

                     vphase = omgrf / akprl
                     if (vphase .ge. clight)vphase = clight

                     c_ehst = 0.0

                     if(xkt(i, j) .ne. 0.0) then

                        xlnlam = 24. - alog(sqrt(xn(i, j) &
                           /1.e+06) / (xkt(i, j)/ abs(q)))
                        vth = sqrt(xkt(i, j) / xm)
                        wphase = vphase / vth

                        if (wphase .lt. 10.0 .and. wphase .gt. 0.1)then

!                       -----------------------------
!                       Assume circular flux surfaces
!                       -----------------------------
                        a = sqrt(x(i)**2 + y(j)**2)
                        rmaxa = rt + a
                        rmina = rt - a
                        bmaxa = b0 * rt / rmina
!                       ------------------------------------
!                       End circular flux surface assumption
!                       ------------------------------------


                        epsa = (rmaxa - rmina)/(rmaxa + rmina)


                        xmut2 = 1. - bmod(i, j) / bmaxa
                        eta0 = 8. * wphase**2 / (5. + zeffcd) + &
                           ceta / zeffcd**.707 + damp / wphase
                        r1 = 1. - epsa**.77 * sqrt(3.5**2 + wphase**2) &
                           /(3.5 * epsa**.77 + wphase)
                        xm1 = 1.0
                        c1 = 1.0

                        if(ftrap .ne. 0)then
                           if (xmut2 .gt. 1.0E-06) then
                              xm1 = 1.0 + afit * &
                                              (xmut2**1.5 / wphase**3)
                              yt = (1. - xmut2) * wphase**2 / xmut2
                              if(yt.ge.0.0)c1 = 1. - &
                                                exp(-(cfit*yt)**xnexp)
                           endif
                        endif

                        xjtild = c1 * xm1 * eta0 * r1


                        c_ehst = -19.19E+15 / xlnlam * &
                           (xkt(i, j) &
                           / abs(q))/ xn(i, j) * &
                           xjtild * signkz
!                        c_ehst = 1.0
                        end if
                     end if

                     if (vphase .ge. clight .or. akprl .lt. 1.0e-05) &
                        c_ehst = 0.0



                     if (isigma .eq. 1) &
                        call sigmad_cql3d(i, j, n, m, rho(i,j),rho_a, &
                            gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
                            xm, q, xn(i,j), xnuomg, &
                            xkt(i,j), omgc(i,j), omgp2(i,j), &
                            -lmax, lmax, nzfun, ibessel, &
                            xkxsav(n), xkysav(m), nphi, capr(i), &
                            bxn(i,j), byn(i,j), bzn(i,j), &
                            uxx(i,j), uxy(i,j), uxz(i,j), &
                            uyx(i,j), uyy(i,j), uyz(i,j), &
                            uzx(i,j), uzy(i,j), uzz(i,j), &
                            sigxx, sigxy, sigxz, &
                            sigyx, sigyy, sigyz, &
                            sigzx, sigzy, sigzz, &
                            delta0, ndist, &
                            nupar, nuper, n_psi, &
                            n_psi_dim, dfduper, dfdupar, &
                            UminPara, UmaxPara, UPERP, UPARA, &
                            vc_mks, df_cql_uprp, df_cql_uprl, nbessj, &
                            nkperp, zi, eps0, v0i, omgrf, xk0, &
                            kperp_max, i_sav, j_sav, upshift, &
                            damping, xk_cutoff,zLoc(j),eNormIN)

                     if (isigma .eq. 0) &
                        call sigmac_stix(i, j, n, m, &
                            xm, q, xn(i,j), xnuomg, &
                            xkt(i,j), omgc(i,j), omgp2(i,j), &
                            -lmax, lmax, nzfun, ibessel, &
                            xkxsav(n), xkysav(m), nphi, capr(i), &
                            bxn(i,j), byn(i,j), bzn(i,j), &
                            uxx(i,j), uxy(i,j), uxz(i,j), &
                            uyx(i,j), uyy(i,j), uyz(i,j), &
                            uzx(i,j), uzy(i,j), uzz(i,j), &
                            sigxx, sigxy, sigxz, &
                            sigyx, sigyy, sigyz, &
                            sigzx, sigzy, sigzz, &
                            delta0, zi, eps0, v0i, omgrf, xk0, &
                            kperp_max, i_sav, j_sav)


                     cexpkxky = xx(n, i) * yy(m, j)

                     xjpx(i,j) = xjpx(i,j) + (sigxx * exk(n,m) &
                                     + sigxy * eyk(n,m) &
                                     + sigxz * ezk(n,m)) * cexpkxky
                     xjpy(i,j) = xjpy(i,j) + (sigyx * exk(n,m) &
                                     + sigyy * eyk(n,m) &
                                     + sigyz * ezk(n,m)) * cexpkxky
                     xjpz(i,j) = xjpz(i,j) + (sigzx * exk(n,m) &
                                     + sigzy * eyk(n,m) &
                                     + sigzz * ezk(n,m)) * cexpkxky


                     xjpx_ehst(i,j) = xjpx_ehst(i,j) + (sigxx * exk(n,m) &
                                 + sigxy * eyk(n,m) &
                                 + sigxz * ezk(n,m)) * cexpkxky * c_ehst
                     xjpy_ehst(i,j) = xjpy_ehst(i,j) + (sigyx * exk(n,m) &
                                 + sigyy * eyk(n,m) &
                                 + sigyz * ezk(n,m)) * cexpkxky * c_ehst
                     xjpz_ehst(i,j) = xjpz_ehst(i,j) + (sigzx * exk(n,m) &
                                 + sigzy * eyk(n,m) &
                                 + sigzz * ezk(n,m)) * cexpkxky * c_ehst

                  end do
               end do

            end if

            end if

         end do
      end do

      call blacs_barrier(icontxt, 'All')

      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpx, &
         nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpy, &
         nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpz, &
         nxdim, -1, -1)


      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpx_ehst, &
         nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpy_ehst, &
         nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpz_ehst, &
         nxdim, -1, -1)


      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (10i10)

      end



!
!***************************************************************************
!

      subroutine current(xjpx, xjpy, xjpz, nxdim, nydim, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xm, q, xn, xkt, omgc, omgp2, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         exk, eyk, ezk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndist, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper, dfdupar, &
         UminPara, UmaxPara, UPERP, UPARA, &
         vc_mks, df_cql_uprp, df_cql_uprl, rho, rho_a, nbessj, nkperp, &
         zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift, &
         damping, xk_cutoff,zLoc,eNormIN)

!-----------------------------------------------------------------------
!     This subroutine calculates the plasma current for a single species
!-----------------------------------------------------------------------

      implicit none

      logical ismine
      real :: eNormIN
      complex zi
      real eps0, v0i, omgrf, xk0, kperp_max
      integer i_sav, j_sav, upshift
      real damping, xk_cutoff
      !DLG:
      real :: zLoc(nydim)

      integer nproc, myid, ngrid, id, ndist, nbessj, nkperp
      integer  i, j, n, m, lmax, nzfun, ibessel, nphi
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx, &
         icontxt, nboundary
      integer dlen_, desc_amat(dlen_)


      integer nxdim, nydim, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, isigma

      real xm, omgc(nxdim, nydim), omgp2(nxdim, nydim), &
           psi(nxdim, nydim), psilim
      real capr(nxdim), xnuomg, delta0
      real gradprlb(nxdim, nydim), bmod(nxdim, nydim), &
                               bmod_mid(nxdim, nydim)

      complex xx(nkdim1 : nkdim2, 1 : nxdim), &
              yy(mkdim1 : mkdim2, 1 : nydim)
      complex xjpx(nxdim, nydim), xjpy(nxdim, nydim), xjpz(nxdim, nydim)

      complex sigxx, sigxy, sigxz, &
              sigyx, sigyy, sigyz, &
              sigzx, sigzy, sigzz


      complex cexpkxky
      complex exk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              eyk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              ezk(nkdim1 : nkdim2, mkdim1 : mkdim2)

      real xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)
      real q, xn(nxdim, nydim), xkt(nxdim, nydim)
      real bxn(nxdim, nydim), byn(nxdim, nydim), bzn(nxdim, nydim)

      real uxx(nxdim, nydim), uxy(nxdim, nydim), uxz(nxdim,nydim), &
           uyx(nxdim, nydim), uyy(nxdim, nydim), uyz(nxdim,nydim), &
           uzx(nxdim, nydim), uzy(nxdim, nydim), uzz(nxdim,nydim)

      real rho(nxdim, nydim)

      integer  :: n_psi_dim, nuper, nupar, n_psi

      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
      real :: UminPara,UmaxPara
      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: vc_mks, rho_a(n_psi_dim)


      xjpx(:,:) = 0.0
      xjpy(:,:) = 0.0
      xjpz(:,:) = 0.0


!--Loop over mode numbers and mesh:
      do i = 1, nnodex
         do j = 1, nnodey

            ngrid = (i - 1) * nnodey + j
            id = mod(ngrid, nproc)
            if(id .eq. myid)then

            xjpx(i,j) = 0.0
            xjpy(i,j) = 0.0
            xjpz(i,j) = 0.0

            if(psi(i,j) .le. psilim .and. nboundary .eq. 1)then

               do n = nkx1, nkx2
                  do m = nky1, nky2

                     if (isigma .eq. 1) &
                        call sigmad_cql3d(i, j, n, m, rho(i,j),rho_a, &
                            gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
                            xm, q, xn(i,j), xnuomg, &
                            xkt(i,j), omgc(i,j), omgp2(i,j), &
                            -lmax, lmax, nzfun, ibessel, &
                            xkxsav(n), xkysav(m), nphi, capr(i), &
                            bxn(i,j), byn(i,j), bzn(i,j), &
                            uxx(i,j), uxy(i,j), uxz(i,j), &
                            uyx(i,j), uyy(i,j), uyz(i,j), &
                            uzx(i,j), uzy(i,j), uzz(i,j), &
                            sigxx, sigxy, sigxz, &
                            sigyx, sigyy, sigyz, &
                            sigzx, sigzy, sigzz, &
                            delta0, ndist, &
                            nupar, nuper, n_psi, &
                            n_psi_dim, dfduper, dfdupar, &
                            UminPara, UmaxPara, UPERP, UPARA, &
                            vc_mks, df_cql_uprp, df_cql_uprl, nbessj, &
                            nkperp, zi, eps0, v0i, omgrf, xk0, &
                            kperp_max, i_sav, j_sav, upshift, &
                            damping, xk_cutoff,zLoc(j),eNormIN)

                     if (isigma .eq. 0) &
                        call sigmac_stix(i, j, n, m, &
                            xm, q, xn(i,j), xnuomg, &
                            xkt(i,j), omgc(i,j), omgp2(i,j), &
                            -lmax, lmax, nzfun, ibessel, &
                            xkxsav(n), xkysav(m), nphi, capr(i), &
                            bxn(i,j), byn(i,j), bzn(i,j), &
                            uxx(i,j), uxy(i,j), uxz(i,j), &
                            uyx(i,j), uyy(i,j), uyz(i,j), &
                            uzx(i,j), uzy(i,j), uzz(i,j), &
                            sigxx, sigxy, sigxz, &
                            sigyx, sigyy, sigyz, &
                            sigzx, sigzy, sigzz, &
                            delta0, zi, eps0, v0i, omgrf, xk0, &
                            kperp_max, i_sav, j_sav)


                     cexpkxky = xx(n, i) * yy(m, j)

                     xjpx(i,j) = xjpx(i,j) + (sigxx * exk(n,m) &
                                     + sigxy * eyk(n,m) &
                                     + sigxz * ezk(n,m)) * cexpkxky
                     xjpy(i,j) = xjpy(i,j) + (sigyx * exk(n,m) &
                                     + sigyy * eyk(n,m) &
                                     + sigyz * ezk(n,m)) * cexpkxky
                     xjpz(i,j) = xjpz(i,j) + (sigzx * exk(n,m) &
                                     + sigzy * eyk(n,m) &
                                     + sigzz * ezk(n,m)) * cexpkxky

                  end do
               end do

            end if

            end if

         end do
      end do

      call blacs_barrier(icontxt, 'All')

      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpx, &
         nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpy, &
         nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpz, &
         nxdim, -1, -1)

      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (10i10)

      end


!
!***************************************************************************
!


      subroutine current_1(xjpx, xjpy, xjpz, nxdim, nydim, &
         nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xm, q, xn, xkt, omgc, omgp2, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         exk, eyk, ezk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndist, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper, dfdupar, &
         UminPara, UmaxPara, UPERP, UPARA, &
         vc_mks, df_cql_uprp, df_cql_uprl, rho, rho_a, nbessj, nkperp, &
         zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift, &
         damping, xk_cutoff, zLoc,eNormIN)

!-----------------------------------------------------------------------
!     This subroutine calculates the plasma current for a single species
!-----------------------------------------------------------------------

      implicit none

      logical ismine
      real :: eNormIN
      complex zi
      real eps0, v0i, omgrf, xk0, kperp_max
      integer i_sav, j_sav, upshift
      real damping, xk_cutoff
      !DLG:
      real :: zLoc(nydim)

      integer nproc, myid, ngrid, id, ndist, nbessj, nkperp
      integer  i, j, n, m, lmax, nzfun, ibessel, nphi
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx, &
         icontxt, nboundary
      integer dlen_, desc_amat(dlen_)


      integer nxdim, nydim, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, isigma

      real xm, omgc(nxdim, nydim), omgp2(nxdim, nydim), &
           psi(nxdim, nydim), psilim
      real capr(nxdim), xnuomg, delta0
      real gradprlb(nxdim, nydim), bmod(nxdim, nydim), &
                               bmod_mid(nxdim, nydim)

      complex xx(nkdim1 : nkdim2, 1 : nxdim), &
              yy(mkdim1 : mkdim2, 1 : nydim)
      complex xjpx(nxdim, nydim), xjpy(nxdim, nydim), xjpz(nxdim, nydim)

      complex sigxx, sigxy, sigxz, &
              sigyx, sigyy, sigyz, &
              sigzx, sigzy, sigzz


      complex cexpkxky
      complex exk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              eyk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              ezk(nkdim1 : nkdim2, mkdim1 : mkdim2)

      real xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)
      real q, xn(nxdim, nydim), xkt(nxdim, nydim)
      real bxn(nxdim, nydim), byn(nxdim, nydim), bzn(nxdim, nydim)

      real uxx(nxdim, nydim), uxy(nxdim, nydim), uxz(nxdim,nydim), &
           uyx(nxdim, nydim), uyy(nxdim, nydim), uyz(nxdim,nydim), &
           uzx(nxdim, nydim), uzy(nxdim, nydim), uzz(nxdim,nydim)

      real rho(nxdim, nydim)

      integer  :: n_psi_dim, nuper, nupar, n_psi

      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
      real :: UminPara,UmaxPara
      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: vc_mks, rho_a(n_psi_dim)


      xjpx(:,:) = 0.0
      xjpy(:,:) = 0.0
      xjpz(:,:) = 0.0


!--Loop over mode numbers and mesh:
      do i = 1, nnodex
         do j = 1, nnodey

            ngrid = (i - 1) * nnodey + j
            id = mod(ngrid, nproc)
            if(id .eq. myid)then

            xjpx(i,j) = 0.0
            xjpy(i,j) = 0.0
            xjpz(i,j) = 0.0

            if(psi(i,j) .le. psilim .and. nboundary .eq. 1)then

               do n = nkx1, nkx2
                  do m = nky1, nky2

                     if (isigma .eq. 1) &
                        call sigmad_cql3d(i, j, n, m, rho(i,j),rho_a, &
                            gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
                            xm, q, xn(i,j), xnuomg, &
                            xkt(i,j), omgc(i,j), omgp2(i,j), &
                            -lmax, lmax, nzfun, ibessel, &
                            xkxsav(n), xkysav(m), nphi, capr(i), &
                            bxn(i,j), byn(i,j), bzn(i,j), &
                            uxx(i,j), uxy(i,j), uxz(i,j), &
                            uyx(i,j), uyy(i,j), uyz(i,j), &
                            uzx(i,j), uzy(i,j), uzz(i,j), &
                            sigxx, sigxy, sigxz, &
                            sigyx, sigyy, sigyz, &
                            sigzx, sigzy, sigzz, &
                            delta0, ndist, &
                            nupar, nuper, n_psi, &
                            n_psi_dim, dfduper, dfdupar, &
                            UminPara, UmaxPara, UPERP, UPARA, &
                            vc_mks, df_cql_uprp, df_cql_uprl, nbessj, &
                            nkperp, zi, eps0, v0i, omgrf, xk0, &
                            kperp_max, i_sav, j_sav, upshift, &
                            damping, xk_cutoff,zLoc(j),eNormIN)

                     if (isigma .eq. 0) &
                        call sigmac_stix(i, j, n, m, &
                            xm, q, xn(i,j), xnuomg, &
                            xkt(i,j), omgc(i,j), omgp2(i,j), &
                            -lmax, lmax, nzfun, ibessel, &
                            xkxsav(n), xkysav(m), nphi, capr(i), &
                            bxn(i,j), byn(i,j), bzn(i,j), &
                            uxx(i,j), uxy(i,j), uxz(i,j), &
                            uyx(i,j), uyy(i,j), uyz(i,j), &
                            uzx(i,j), uzy(i,j), uzz(i,j), &
                            sigxx, sigxy, sigxz, &
                            sigyx, sigyy, sigyz, &
                            sigzx, sigzy, sigzz, &
                            delta0, zi, eps0, v0i, omgrf, xk0, &
                            kperp_max, i_sav, j_sav)


                     cexpkxky = xx(n, i) * yy(m, j)

                     xjpx(i,j) = xjpx(i,j) + (sigxx * exk(n,m) &
                                     + sigxy * eyk(n,m) &
                                     + sigxz * ezk(n,m)) * cexpkxky
                     xjpy(i,j) = xjpy(i,j) + (sigyx * exk(n,m) &
                                     + sigyy * eyk(n,m) &
                                     + sigyz * ezk(n,m)) * cexpkxky
                     xjpz(i,j) = xjpz(i,j) + (sigzx * exk(n,m) &
                                     + sigzy * eyk(n,m) &
                                     + sigzz * ezk(n,m)) * cexpkxky

                  end do
               end do

            end if

            end if

         end do
      end do

      call blacs_barrier(icontxt, 'All')

      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpx, &
         nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpy, &
         nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpz, &
         nxdim, -1, -1)

      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (10i10)

      end


!
!***************************************************************************
!


      subroutine current_2(xjpx, xjpy, xjpz, nxdim, nydim, &
         nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xm, q, xn, xkt, omgc, omgp2, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         exk, eyk, ezk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndist, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper, dfdupar, &
         UminPara, UmaxPara, UPERP, UPARA, &
         vc_mks, df_cql_uprp, df_cql_uprl, rho, rho_a, nbessj, nkperp, &
         zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift, &
         damping, xk_cutoff, zLoc, eNormIN)

!-----------------------------------------------------------------------
!     This subroutine calculates the plasma current for a single species
!-----------------------------------------------------------------------

      implicit none

      logical ismine

      real :: eNormIN
      complex zi
      real eps0, v0i, omgrf, xk0, kperp_max
      integer i_sav, j_sav, upshift
      real damping, xk_cutoff
      !DLG:
      real :: zLoc(nydim)

      integer nproc, myid, ngrid, id, ndist, nbessj, nkperp
      integer  i, j, n, m, lmax, nzfun, ibessel, nphi
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx, &
         icontxt, nboundary
      integer dlen_, desc_amat(dlen_)


      integer nxdim, nydim, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, isigma

      real xm, omgc(nxdim, nydim), omgp2(nxdim, nydim), &
           psi(nxdim, nydim), psilim
      real capr(nxdim), xnuomg, delta0
      real gradprlb(nxdim, nydim), bmod(nxdim, nydim), &
                               bmod_mid(nxdim, nydim)

      complex xx(nkdim1 : nkdim2, 1 : nxdim), &
              yy(mkdim1 : mkdim2, 1 : nydim)
      complex xjpx(nxdim, nydim), xjpy(nxdim, nydim), xjpz(nxdim, nydim)

      complex sigxx, sigxy, sigxz, &
              sigyx, sigyy, sigyz, &
              sigzx, sigzy, sigzz


      complex cexpkxky
      complex exk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              eyk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              ezk(nkdim1 : nkdim2, mkdim1 : mkdim2)

      real xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)
      real q, xn(nxdim, nydim), xkt(nxdim, nydim)
      real bxn(nxdim, nydim), byn(nxdim, nydim), bzn(nxdim, nydim)

      real uxx(nxdim, nydim), uxy(nxdim, nydim), uxz(nxdim,nydim), &
           uyx(nxdim, nydim), uyy(nxdim, nydim), uyz(nxdim,nydim), &
           uzx(nxdim, nydim), uzy(nxdim, nydim), uzz(nxdim,nydim)

      real rho(nxdim, nydim)

      integer  :: n_psi_dim, nuper, nupar, n_psi

      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
      real :: UminPara,UmaxPara
      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: vc_mks, rho_a(n_psi_dim)


      xjpx(:,:) = 0.0
      xjpy(:,:) = 0.0
      xjpz(:,:) = 0.0


!--Loop over mode numbers and mesh:
      do i = 1, nnodex
         do j = 1, nnodey

            ngrid = (i - 1) * nnodey + j
            id = mod(ngrid, nproc)
            if(id .eq. myid)then

            xjpx(i,j) = 0.0
            xjpy(i,j) = 0.0
            xjpz(i,j) = 0.0

            if(psi(i,j) .le. psilim .and. nboundary .eq. 1)then

               do n = nkx1, nkx2
                  do m = nky1, nky2

                     if (isigma .eq. 1) &
                        call sigmad_cql3d(i, j, n, m, rho(i,j),rho_a, &
                            gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
                            xm, q, xn(i,j), xnuomg, &
                            xkt(i,j), omgc(i,j), omgp2(i,j), &
                            -lmax, lmax, nzfun, ibessel, &
                            xkxsav(n), xkysav(m), nphi, capr(i), &
                            bxn(i,j), byn(i,j), bzn(i,j), &
                            uxx(i,j), uxy(i,j), uxz(i,j), &
                            uyx(i,j), uyy(i,j), uyz(i,j), &
                            uzx(i,j), uzy(i,j), uzz(i,j), &
                            sigxx, sigxy, sigxz, &
                            sigyx, sigyy, sigyz, &
                            sigzx, sigzy, sigzz, &
                            delta0, ndist, &
                            nupar, nuper, n_psi, &
                            n_psi_dim, dfduper, dfdupar, &
                            UminPara, UmaxPara, UPERP, UPARA, &
                            vc_mks, df_cql_uprp, df_cql_uprl, nbessj, &
                            nkperp, zi, eps0, v0i, omgrf, xk0, &
                            kperp_max, i_sav, j_sav, upshift, &
                            damping, xk_cutoff,zLoc(j),eNormIN)

                     if (isigma .eq. 0) &
                        call sigmac_stix(i, j, n, m, &
                            xm, q, xn(i,j), xnuomg, &
                            xkt(i,j), omgc(i,j), omgp2(i,j), &
                            -lmax, lmax, nzfun, ibessel, &
                            xkxsav(n), xkysav(m), nphi, capr(i), &
                            bxn(i,j), byn(i,j), bzn(i,j), &
                            uxx(i,j), uxy(i,j), uxz(i,j), &
                            uyx(i,j), uyy(i,j), uyz(i,j), &
                            uzx(i,j), uzy(i,j), uzz(i,j), &
                            sigxx, sigxy, sigxz, &
                            sigyx, sigyy, sigyz, &
                            sigzx, sigzy, sigzz, &
                            delta0, zi, eps0, v0i, omgrf, xk0, &
                            kperp_max, i_sav, j_sav)


                     cexpkxky = xx(n, i) * yy(m, j)

                     xjpx(i,j) = xjpx(i,j) + (sigxx * exk(n,m) &
                                     + sigxy * eyk(n,m) &
                                     + sigxz * ezk(n,m)) * cexpkxky
                     xjpy(i,j) = xjpy(i,j) + (sigyx * exk(n,m) &
                                     + sigyy * eyk(n,m) &
                                     + sigyz * ezk(n,m)) * cexpkxky
                     xjpz(i,j) = xjpz(i,j) + (sigzx * exk(n,m) &
                                     + sigzy * eyk(n,m) &
                                     + sigzz * ezk(n,m)) * cexpkxky

                  end do
               end do

            end if

            end if

         end do
      end do

      call blacs_barrier(icontxt, 'All')

      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpx, &
         nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpy, &
         nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpz, &
         nxdim, -1, -1)

      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (10i10)

      end


!
!***************************************************************************
!


      subroutine cur_slo(xjpx, xjpy, xjpz, nxdim, nydim, &
         nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xm, q, xn, eslow, omgc, omgp2, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         exk, eyk, ezk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         xkte, zeffcd, myid, nproc, &
         zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav)

      implicit none

      complex zi
      real eps0, v0i, omgrf, xk0, kperp_max
      integer i_sav, j_sav

      logical ismine


      integer nproc, myid, ngrid, id
      integer  i, j, n, m, lmax, nzfun, ibessel, nphi
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx, &
         icontxt, nboundary
      integer dlen_, desc_amat(dlen_)


      integer nxdim, nydim, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, isigma

      real xm, omgc(nxdim, nydim), omgp2(nxdim, nydim), &
           psi(nxdim, nydim), psilim, zeffcd
      real capr(nxdim), xnuomg

      complex xx(nkdim1 : nkdim2, 1 : nxdim), &
              yy(mkdim1 : mkdim2, 1 : nydim)
      complex xjpx(nxdim, nydim), xjpy(nxdim, nydim), xjpz(nxdim, nydim)

      complex sigxx, sigxy, sigxz, &
              sigyx, sigyy, sigyz, &
              sigzx, sigzy, sigzz


      complex cexpkxky
      complex exk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              eyk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              ezk(nkdim1 : nkdim2, mkdim1 : mkdim2)

      real xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)
      real q, xn(nxdim, nydim), eslow, xkte(nxdim, nydim)
      real bxn(nxdim, nydim), byn(nxdim, nydim), bzn(nxdim, nydim)

      real uxx(nxdim, nydim), uxy(nxdim, nydim), uxz(nxdim,nydim), &
           uyx(nxdim, nydim), uyy(nxdim, nydim), uyz(nxdim,nydim), &
           uzx(nxdim, nydim), uzy(nxdim, nydim), uzz(nxdim,nydim)


!--Loop over mode numbers and mesh:
      do i = 1, nnodex
         do j = 1, nnodey

            ngrid = (j - 1) * nnodex + i
            id = mod(ngrid, nproc)
            if(id .eq. myid)then

            xjpx(i,j) = 0.0
            xjpy(i,j) = 0.0
            xjpz(i,j) = 0.0

            if(psi(i,j) .le. psilim .and. nboundary .eq. 1)then

               do n = nkx1, nkx2
                  do m = nky1, nky2

                     call sigmah_slow(i, j, n, m, &
                         xm, q, xn(i,j), xnuomg, &
                         eslow, omgc(i,j), omgp2(i,j), &
                         -lmax, lmax, nzfun, ibessel, &
                         xkxsav(n), xkysav(m), nphi, capr(i), &
                         bxn(i,j), byn(i,j), bzn(i,j), &
                         uxx(i,j), uxy(i,j), uxz(i,j), &
                         uyx(i,j), uyy(i,j), uyz(i,j), &
                         uzx(i,j), uzy(i,j), uzz(i,j), &
                         sigxx, sigxy, sigxz, &
                         sigyx, sigyy, sigyz, &
                         sigzx, sigzy, sigzz, &
                         xkte(i,j), zeffcd, zi, eps0, v0i, omgrf, xk0, &
                         kperp_max, i_sav, j_sav)


                     cexpkxky = xx(n, i) * yy(m, j)

                     xjpx(i,j) = xjpx(i,j) + (sigxx * exk(n,m) &
                                     + sigxy * eyk(n,m) &
                                     + sigxz * ezk(n,m)) * cexpkxky
                     xjpy(i,j) = xjpy(i,j) + (sigyx * exk(n,m) &
                                     + sigyy * eyk(n,m) &
                                     + sigyz * ezk(n,m)) * cexpkxky
                     xjpz(i,j) = xjpz(i,j) + (sigzx * exk(n,m) &
                                     + sigzy * eyk(n,m) &
                                     + sigzz * ezk(n,m)) * cexpkxky

                  end do
               end do

            end if
            end if

         end do
      end do

      call blacs_barrier(icontxt, 'All')

      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpx, &
         nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpy, &
         nxdim, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, xjpz, &
         nxdim, -1, -1)

      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (10i10)

      end

!
!***************************************************************************
!

