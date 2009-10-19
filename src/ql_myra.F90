module ql_myra_mod

      use qlsum_myra_mod
      use read_particle_f
      use gc_integrate
      use netcdf
      use dlg_p_check
      use size_mod
      !use mpi
      include 'mpif.h'
#     include <pnetcdf.inc>
      contains

!
!***************************************************************************
!

      subroutine ql_myra_write(bqlavg, cqlavg, eqlavg, fqlavg, &
         wdot_inout, fx0_inout, fy0_inout, fz0_inout, &
         vol, nrhodim, nnoderho, drho, &
         dx, dy, r0, nxdim, nydim, nnodex, nnodey, &
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
         UminPara_cql, UmaxPara_cql, UPERP_cql, UPARA_cql, UPERP, UPARA, &
         vc_mks_cql, df_cql_uprp, df_cql_uprl, rho, rho_a, nbessj, &
         zeffcd, clight, x, y, rt, b0, ftrap, omgrf, &
         nkperp, lmaxdim, nzeta_wdot, theta_, &
         n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg, &
         n_bin, upshift, i_write, xk_cutoff, zLoc, eNormIN, mask)

!----------------------------------------------------------------------
!     This subroutine calculates wdot and the quasi-linear operator for
!     a single species
!----------------------------------------------------------------------

      use dlg
      use aorsa2din_mod, only: nModesX, nModesY, eNorm_factor, &
        ana_maxwellian
      use write_pql
      implicit none

      integer, parameter :: DBL = selected_real_kind ( 13, 300 )
      integer, intent(IN) :: mask(nmodesmax,nmodesmax)
      real :: eNormIN
      integer splitrank, nloops, partition, start, finish
      integer pstart, pfinish, status(MPI_STATUS_SIZE), ierr
      integer recv_size, remainder
      real, allocatable :: bql_store(:,:,:)
      real, allocatable :: cql_store(:,:,:)
      real, allocatable :: eql_store(:,:,:)
      real, allocatable :: fql_store(:,:,:)

      real, allocatable :: bql4D(:,:,:,:)
      real, allocatable :: cql4D(:,:,:,:)
      real, allocatable :: eql4D(:,:,:,:)
      real, allocatable :: fql4D(:,:,:,:)

      real, allocatable :: bql4D_cyl(:,:,:,:)
      real, allocatable :: cql4D_cyl(:,:,:,:)
      real, allocatable :: eql4D_cyl(:,:,:,:)
      real, allocatable :: fql4D_cyl(:,:,:,:)

      real :: ql_th, vTmp

      logical iam_root
      integer left_neighbor, right_neighbor
      real token
      real upara_test
      integer mi_max, mi_min

      logical ismine

      real u2, fnorm, f_cql

      real, dimension(:,:), allocatable :: DFDUPER0, DFDUPAR0

      integer n_theta_max, n_psi_max, i_psi_eq, n_bin, upshift
      integer n_theta_(n_psi_max), n_theta, ntheta_giv, i_write

      integer nproc, myid, ngrid, id, ndist, nbessj, nkperp, lmaxdim
      integer  i, j, n, m, lmax, nzfun, ibessel, nphi
      integer rsrc, csrc, myrow, mycol, nprow, npcol, lrindx, lcindx, &
         icontxt, nboundary, nzeta_wdot
      integer dlen_, desc_amat(dlen_)
      integer nxdim, nydim, nnodex, nnodey, nrhodim, nnoderho, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, isigma
      integer ftrap, k, ni0, mi0

      real dfdth, dfdupar_check, dfduper_check, dfdth_check
      real uperp0_grid, upara0_grid, zeta, eta, ai, bi, ci, di
      real dfduper0_intplt, dfdupar0_intplt, xk_cutoff

      real theta_(n_theta_max, n_psi_max), thegiv, dthetag
      real deriv, dtheta, dtau_ratio_giv, tau_bounce_giv
      real uperp0, upara0, xlamda, derivb, dtheta0, factor, upara_mi
      real duperp, dupara

      real zeffcd, clight, rt, b0, damp, r1, xm1, ceta, epsa, xmut2, &
         eta0, xnexp, xjtild, signkz, cfit, c1, afit, yt, bmaxa, vphase, &
         omgrf, c_ehst, akprl, xkprl, xkphi, a, rmaxa, rmina, &
         wphase, xlnlam, vth
      real x(nxdim), y(nydim), argd, bratio, drho


      real xm, omgc(nxdim, nydim), omgp2(nxdim, nydim), &
           psi(nxdim, nydim), psilim, alpha
!DLG:   Define variables
      real capr(nxdim), xnuomg(:,:), delta0, dx, dy, r0
      real gradprlb(nxdim, nydim), bmod(nxdim, nydim), &
                               bmod_mid(nxdim, nydim)

      complex xx(nkdim1 : nkdim2, 1 : nxdim), &
              yy(mkdim1 : mkdim2, 1 : nydim)

      complex wdoti, fx0i, fy0i


      real bqlavg(nuper, nupar, nnoderho)
      real cqlavg(nuper, nupar, nnoderho)
      real eqlavg(nuper, nupar, nnoderho)
      real fqlavg(nuper, nupar, nnoderho)

      real wdot_inout(nxdim, nydim)
      real  fx0_inout(nxdim, nydim)
      real  fy0_inout(nxdim, nydim)
      real  fz0_inout(nxdim, nydim)

      real wdot(nnodex, nnodey)
      real fx0(nnodex, nnodey)
      real fy0(nnodex, nnodey)
      real fz0(nnodex, nnodey)

!      real count(0 : 10000, 1), sum_count

      real vol(nrhodim), dldbavg(nrhodim)



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

      integer  :: n_psi_dim, nuper, nupar, n_psi, mi, ni

      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: UPERP_cql(NUPER), UPARA_cql(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)

      real bql, cql, eql, fql
      real bqlTmp, cqlTmp, eqlTmp, fqlTmp

      complex, dimension(:),  allocatable :: b_sum
      complex, dimension(:),  allocatable :: c_sum
      complex, dimension(:),  allocatable :: e_sum
      complex, dimension(:),  allocatable :: f_sum

      complex, dimension(:),  allocatable :: wdot_sum
      complex, dimension(:),  allocatable :: sum_fx0
      complex, dimension(:),  allocatable :: sum_fy0

      real, dimension(:,:,:), allocatable :: factvol
      real, dimension(:,:), allocatable :: factvol2d

      real, dimension(:,:,:), allocatable :: bqlvol
      real, dimension(:,:), allocatable :: bqlvol2d

      real, dimension(:,:,:), allocatable :: cqlvol
      real, dimension(:,:), allocatable :: cqlvol2d

      real, dimension(:,:,:), allocatable :: eqlvol
      real, dimension(:,:), allocatable :: eqlvol2d

      real, dimension(:,:,:), allocatable :: fqlvol
      real, dimension(:,:), allocatable :: fqlvol2d

      real :: W, ENORM, ZSPEC, ASPEC, BMAG
      real :: UminPara,UmaxPara
      real :: UminPara_cql,UmaxPara_cql
      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)
      real :: vc_mks, vc_mks_cql, rho_a(n_psi_dim)
      real :: eps0, pi, emax, u0, u_0, u_, costh, costh0, psic

      parameter (eps0 = 8.85e-12)
      parameter (PI = 3.141592653597932384)
      real, parameter :: e = 1.6e-19


      integer :: nwork, ip
      logical :: has_work
      integer, dimension(nnodex*nnodey) :: i_table, j_table

!DLG: Define more variables
      real, dimension(nuPer,nuPar) :: p_f_rzvv_, &
       p_dfduPerp_,p_dfduPar_
      real :: zLoc(nnodey)
      real :: startR, startz,vPer,vPar
      complex :: wdot_sum_dlg(nupar)
      integer :: wdot_sum_res(nupar)
      complex, allocatable :: wdot_sum_UU(:,:)
      integer :: wdot_sum_res_UU(nuper,nupar)
      complex, allocatable :: wdoti_RZ(:,:)
      real :: wdot__(nnodex, nnodey)
      real, allocatable :: wdot_orbit(:,:)
      real, allocatable :: wdot_orbit_(:,:)
      logical :: plotVar, goodOrbit
      integer :: nc_id, nR_id, nz_id, scalar_id, wdot_id
      character(len=100) :: ncFileName, pncFileName
      integer :: orbit_id, capR_id, zLoc_id, nRMax_id
      logical :: ba_wdot
      integer :: dumpQL
      integer :: ql_b_id, ql_c_id, ql_e_id, ql_f_id
      integer :: nuper_id, nupar_id, uper_id, upar_id
      integer :: pnc_id, mpi_iErr, mpi_ret
      integer :: pnR_id, pnz_id, pnuper_id, pnupar_id
      integer :: pql_b_id, pql_c_id, pql_e_id, pql_f_id
      integer :: pql_vPer_id, pql_vPar_id, pql_vc_mks_id
      integer(KIND=MPI_OFFSET_KIND) :: pnc_start(4), pnc_cnt(4)
      integer(KIND=MPI_OFFSET_KIND) :: pnc_nR, pnc_nz, &
       pnc_nuper, pnc_nupar, pnc_scalar
      integer :: what, ivalue
      integer :: dimIds(4)
      real :: tmpArray(nuper,nupar)
      real :: ql_tmp(nuper,nupar), vMin, vMax
      integer :: pql_maxV_id, pql_minV_id, &
       qlMax(2), qlMin(2), pql_R_id, pql_z_id
      integer(KIND=MPI_OFFSET_KIND) :: bM_start(2), bM_cnt(2)
      real :: vPar_binSize, vPar_stripSize
      real :: dvPer
      integer :: npRowOut, npColOut, myRowDLG, myColDLG
      integer :: iStart, iFinish, jStart, jFinish

      allocate( dfduper0(nuper, nupar) )
      allocate( dfdupar0(nuper, nupar) )

      allocate(b_sum(NUPAR) )
      allocate(c_sum(NUPAR) )
      allocate(e_sum(NUPAR) )
      allocate(f_sum(NUPAR) )

!DLG: Allocate
      allocate(wdot_sum_UU(nuper,nupar))
      allocate(wdot_orbit(nnodex,nnodey))
      allocate(wdot_orbit_(nnodex,nnodey))
      allocate(wdoti_RZ(nnodex,nnodey))

      allocate(wdot_sum(NUPER) )
      allocate(sum_fx0(NUPER) )
      allocate(sum_fy0(NUPER) )

      allocate(factvol(nuper, nupar, nnoderho) )
      allocate(factvol2d(nuper, nupar) )

      allocate(bqlvol(nuper, nupar, nnoderho) )
      allocate(bqlvol2d(nuper, nupar) )

      allocate(cqlvol(nuper, nupar, nnoderho) )
      allocate(cqlvol2d(nuper, nupar) )

      allocate(eqlvol(nuper, nupar, nnoderho) )
      allocate(eqlvol2d(nuper, nupar) )

      allocate(fqlvol(nuper, nupar, nnoderho) )
      allocate(fqlvol2d(nuper, nupar) )

!     -------------------------------------
!efd  initialize allocatable arrays to zero
!     -------------------------------------
      wdot = 0.0
      fx0 = 0.0
      fy0 = 0.0
      fz0 = 0.0

      b_sum = 0.0
      c_sum = 0.0
      e_sum = 0.0
      f_sum = 0.0

!DLG: Initialise variables to zero
      wdot__ = 0.0
      wdot_sum_dlg = 0.0
      wdot_sum_UU = 0.0
      wdot_orbit = 0.0
      wdot_orbit_ = 0.0
      wdoti_RZ = 0.0
      wdot_sum = 0.0
      sum_fx0 = 0.0
      sum_fy0 = 0.0

!      count = 0.0

!      factvol = 0.0
!      factvol2d = 0.0

      bqlvol = 0.0
      bqlvol2d = 0.0

      cqlvol = 0.0
      cqlvol2d = 0.0

      eqlvol = 0.0
      eqlvol2d = 0.0

      fqlvol = 0.0
      fqlvol2d = 0.0

      !DLG:   Dump QL coeffs to netCdf file for plotting.

      dumpQL = i_write



      W = omgrf
      ZSPEC = q / 1.6e-19

      if(ndist .eq. 0)then   !--Maxwellian--!

         do ni = 1, NUPER
            UPERP(ni) = (real(ni-1)/real(NUPER-1))
         end do

         UminPara = -1.0
         UmaxPara =  1.0

         do mi = 1, NUPAR
            UPARA(mi) = (-1.0 + 2. * (real(mi-1) / real(NUPAR-1)))
         end do

      else   !--non-Maxwellian--!

            vc_mks = vc_mks_cql
            UminPara = UminPara_cql
            UmaxPara = UmaxPara_cql

            do ni = 1, nuper
               uperp(ni) = uperp_cql(ni)
            end do

            do mi = 1, nupar
               upara(mi) = upara_cql(mi)
            end do

      end if


      do n = 1, nnoderho

         do mi = 1, nupar
            do ni = 1, nuper


               bqlvol(ni, mi, n) = 0.0
               cqlvol(ni, mi, n) = 0.0
               eqlvol(ni, mi, n) = 0.0
               fqlvol(ni, mi, n) = 0.0


               bqlavg(ni, mi, n) = 0.0
               cqlavg(ni, mi, n) = 0.0
               eqlavg(ni, mi, n) = 0.0
               fqlavg(ni, mi, n) = 0.0

            end do
         end do
      end do



      nwork = 0
      do j=1,nnodey
         do i=1,nnodex
           !has_work = (psi(i,j) .le. psilim .and. nboundary .eq. 1)
           has_work = (mask(i,j) == 1 .and. nboundary .eq. 1 .and. bmod_mid(i,j) > 0)
           if (has_work) then
              nwork = nwork + 1
              i_table(nwork) = i
              j_table(nwork) = j
           endif
         enddo
      enddo


!     -----------------------
!     Loop over spatial mesh:
!     -----------------------

      nloops = nwork
      partition=int(nloops/nproc)
      remainder=mod(nloops,nproc)
      splitrank=nproc-remainder
      if(myid.lt.splitrank)then
         start=1+myid*partition
         finish=start+partition-1
      else
         start=1+myid*partition+myid-splitrank
         finish=start+partition
      endif

!      call blacs_gridInfo ( &
!           iContxt, npRowOut, npColOut, myRowDLG, myColDLG )
!
!      call blacs_barrier(icontxt, 'All')
!
!!DLG:   Adjust if statement for ndist >= 1
!
!      if(i_write .eq. 1) then
!!             allocate(bql_store(start:finish,1:nuper,1:nupar),
!!     .         cql_store(start:finish,1:nuper,1:nupar),
!!     .         eql_store(start:finish,1:nuper,1:nupar),
!!     .         fql_store(start:finish,1:nuper,1:nupar))
!
!      allocate(bql4D(nModesX/npRowOut,nModesY/npColOut,1:nuPer,1:nuPar))
!      allocate(cql4D(nModesX/npRowOut,nModesY/npColOut,1:nuPer,1:nuPar))
!      allocate(eql4D(nModesX/npRowOut,nModesY/npColOut,1:nuPer,1:nuPar))
!      allocate(fql4D(nModesX/npRowOut,nModesY/npColOut,1:nuPer,1:nuPar))
!
!      allocate(bql4D_cyl(nModesX/npRowOut, &
!                nModesY/npColOut,1:nuPer,1:nuPar))
!      !allocate(cql4D_cyl(nModesX/npRowOut,nModesY/npColOut,1:nuPer,1:nu
!      !allocate(eql4D_cyl(nModesX/npRowOut,nModesY/npColOut,1:nuPer,1:nu
!      !allocate(fql4D_cyl(nModesX/npRowOut,nModesY/npColOut,1:nuPer,1:nu
!
!
!      endif

!spatial_loop: &
!do ip = start,finish
!iStart    = nModesX / npRowOut * myRowDLG
!jStart    = nModesY / npColOut * myColDLG
!iFinish   = nModesX / npRowOut * ( myRowDLG + 1 )
!jFinish   = nModesY / npColOut * ( myColDLG + 1 )

ip_loop: do ip = start, finish
!i_loop: do i = iStart + 1, iFinish
!  j_loop: do j = jStart + 1, jFinish
!      inMask: &
!      if (mask(i,j) == 1) then
      i = i_table(ip)
      j = j_table(ip)

      xkphi = nphi / capr(i)
      if(xkphi .eq. 0.0)xkphi = 1.0e-05

      alpha = sqrt(2.0 * xkt(i,j) / xm)
      n = int(rho(i,j) / drho) + 1

      !DLG: Set maxwellian v grid range using eNorm_factor
      if (ndist .eq. 0) then

          if ( eNormIN < 0 ) then
              vc_mks =  3.0 * alpha    !--Maxwellian only--!
          else
              vc_mks = vc_mks_cql!sqrt(2d0*1d3*eNormIN*e/xm)
          endif

          vc_mks_cql = vc_mks
          UminPara_cql = -1.0
          UmaxPara_cql = 1.0

      endif

      !has_work = (psi(i,j) .le. psilim .and. nboundary .eq. 1)
      !has_work = (mask(i,j) == 1 .and. nboundary .eq. 1 .and. bmod_mid(i,j) > 0)

      !hasWork: &
      !if (has_work) then

      u0 = vc_mks / alpha

      Emax = 0.5 * xm * vc_mks**2
      Enorm = Emax / 1.6e-19
      ASPEC = xm / 1.67e-27
      BMAG = omgc(i,j) * (xm / q)

    duperp = (uperp(nuper) - uperp(1)) / (nuper - 1)
      dupara = (upara(nupar) - upara(1)) / (nupar - 1)

    bratio = bmod_mid(i,j) / bmod(i,j)
    if (bratio .gt. 1.0) bratio = 1.0
    psic = 1.0 / bratio

!           -----------------------------------------------------
!           get CQL3D distribution function on the midplane:
!           used for Wdot only - not the quasilinear coefficients
!           -----------------------------------------------------

        if(ndist .eq. 0)then   !--Maxwellian--!

         call maxwell_dist(u0, NUPAR, NUPER, &
                UminPara, UmaxPara, &
                UPERP, UPARA, DFDUPER, DFDUPAR)

        else   !--non-Maxwellian--!

           call cql3d_dist(nupar, nuper, n_psi, &
                n_psi_dim, rho_a, rho(i,j), &
                UminPara,UmaxPara, &
                df_cql_uprp, df_cql_uprl, &
                UPERP, UPARA, DFDUPER0, DFDUPAR0)

!              ---------------------------------------------------------
!              map CQL3D distribution function off the midplane for Wdot
!              ---------------------------------------------------------
               ifBratio: &
               if(bratio .gt. 0.0)then

                dfduper = 0.0
                dfdupar = 0.0

                  ni_loop: do ni = 1, nuper
                     mi_loop: do mi = 1, nupar

                        argd =  uperp(ni)**2 * (1. - bratio) &
                                                      + upara(mi)**2
                        if (argd .le. 0.0) argd = 1.0e-06

                    uperp0 = uperp(ni) * sqrt(bratio)
                        upara0  = sign(1.0, upara(mi)) * sqrt(argd)

                    dfduper(ni, mi) = 0.0
                        dfdupar(ni, mi) = 0.0

                    if(upara0 .ge. upara(1) .and. &
                                         upara0 .le. upara(nupar)) then

                            ni0 = int((uperp0 - uperp(1)) / duperp) + 1
                     mi0 = int((upara0 - upara(1)) / dupara) + 1

                     dfduper0_intplt = dfduper0(ni0, mi0)
                     dfdupar0_intplt = dfdupar0(ni0, mi0)

                       if (ni0 .lt. nuper .and. mi0 .lt. nupar) then

                           uperp0_grid = uperp(1) + (ni0 - 1) * duperp
                     upara0_grid = upara(1) + (mi0 - 1) * dupara

                     zeta = (uperp0 - uperp0_grid) / duperp
                     eta  = (upara0 - upara0_grid) / dupara

                           ai = dfduper0(ni0, mi0)
                           bi = dfduper0(ni0+1,mi0) - dfduper0(ni0,mi0)
                           ci = dfduper0(ni0,mi0+1) - dfduper0(ni0,mi0)
                           di = dfduper0(ni0+1,mi0+1)+ dfduper0(ni0,mi0) &
                              - dfduper0(ni0+1,mi0)- dfduper0(ni0,mi0+1)

                     dfduper0_intplt = ai + bi * zeta &
                                           + ci * eta + di * zeta * eta

                           ai = dfdupar0(ni0, mi0)
                           bi = dfdupar0(ni0+1,mi0) - dfdupar0(ni0,mi0)
                           ci = dfdupar0(ni0,mi0+1) - dfdupar0(ni0,mi0)
                           di = dfdupar0(ni0+1,mi0+1)+ dfdupar0(ni0,mi0) &
                              - dfdupar0(ni0+1,mi0)- dfdupar0(ni0,mi0+1)

                     dfdupar0_intplt = ai + bi * zeta &
                                           + ci * eta + di * zeta * eta

                           endif


                     if (upara0 .ne. 0.0)then

                        dfdupar(ni, mi) = dfdupar0_intplt * &
                                 upara(mi) / upara0

                              dfduper(ni, mi) = dfduper0_intplt * &
                                 sqrt(bratio) + dfdupar0_intplt * &
                                 uperp(ni) / upara0 * (1.0 - bratio)

!                             dfdth = upara(mi) * dfduper(ni, mi)
!     .                             - uperp(ni) * dfdupar(ni, mi)

                           endif

                  endif


            if (.not. ana_maxwellian) go to 5000
!                       ----------------------------
!                       optional analytic Maxwellian
!                       ----------------------------

                  alpha = sqrt(2.0 * xkt(i, j) / xm)
!                        vc_mks = 3.5 * alpha
                        u0 = vc_mks / alpha

                      fnorm = u0**3 / pi**1.5

                    u2 = uperp(ni)**2 + upara(mi)**2

                        f_cql = exp(-u2 * u0**2) * fnorm
                      dfduper(ni, mi) = -f_cql * 2.* uperp(ni) * u0**2
                        dfdupar(ni, mi) = -f_cql * 2.* upara(mi) * u0**2
 5000                   continue


                     enddo mi_loop
                  enddo ni_loop
          endif ifBratio

            endif

!DLG:   Overwrite the df terms with those from the particle list for ql_myra

    if ( ndist .eq. 2 .and. (.not. ana_maxwellian) ) then

        call dlg_particle_f ( capR(i), zLoc(j), &
            uPerp, uPara, &
            nuper, nupar, xm, eNormIN, &
            p_f_rzvv_, p_dfduPerp_, p_dfduPar_ )
        
        dfdupar   = 0.0
        dfduper = 0.0
        
        dfdupar = p_dfduPar_
        dfduper = p_dfduPerp_
    
    end if

!           ----------------------------------
!           Loop over perpendicular velocities
!           ----------------------------------

            ni_loop2: do ni = 1, nuper

          if (ndist .eq. 0) &
                call QLSUM_MAXWELLIAN(ni, b_sum, c_sum, e_sum, f_sum, &
                     wdot_sum(ni), sum_fx0(ni), sum_fy0(ni), W, ZSPEC, &
                     ASPEC, BMAG, lmax, ENORM, UminPara, UmaxPara, &
                     NUPAR, NUPER, UPERP, UPARA, DFDUPER, DFDUPAR, &
                     exk, eyk, ezk, nkdim1, nkdim2, mkdim1, mkdim2, &
                     nkx1, nkx2, nky1, nky2, &
                     uxx(i,j), uxy(i,j), uxz(i,j), &
                     uyx(i,j), uyy(i,j), uyz(i,j), &
                     uzx(i,j), uzy(i,j), uzz(i,j), &
                     nxdim, nydim, xkxsav, xkysav, xkphi, xx, yy, i, j, &
                     lmaxdim, ndist, nzeta_wdot, &
                     gradprlb(i,j), bmod(i,j), omgc(i,j), alpha, xm, &
                     upshift, xk_cutoff )

!DLG:   Adjust if statement for ndist >= 1
          if (ndist .ge. 1) &
               call QLSUM_NON_MAXWELLIAN(ni, b_sum, c_sum, e_sum, f_sum, &
                     wdot_sum(ni), sum_fx0(ni), sum_fy0(ni), W, ZSPEC, &
                     ASPEC, BMAG, lmax, ENORM, UminPara, UmaxPara, &
                     NUPAR, NUPER, UPERP, UPARA, DFDUPER, DFDUPAR, &
                     exk, eyk, ezk, nkdim1, nkdim2, mkdim1, mkdim2, &
                     nkx1, nkx2, nky1, nky2, &
                     uxx(i,j), uxy(i,j), uxz(i,j), &
                     uyx(i,j), uyy(i,j), uyz(i,j), &
                     uzx(i,j), uzy(i,j), uzz(i,j), &
                     nxdim, nydim, xkxsav, xkysav, xkphi, xx, yy, i, j, &
                     lmaxdim, ndist, nzeta_wdot, &
                     gradprlb(i,j), bmod(i,j), omgc(i,j), alpha, xm, &
                     upshift, xk_cutoff, wdot_sum_dlg, wdot_sum_res )


!!DLG:   Retain vPar wdot sum in array
!        wdot_sum_UU(ni,:) = wdot_sum_dlg
!        wdot_sum_res_UU(ni,:) = wdot_sum_res
!              -----------------------------
!              Loop over parallel velocities
!              -----------------------------

               upara_test = sqrt(1.0 - (UPERP(ni))**2)
             mi_max = ceiling(upara_test*(nupar-1)/2 + (nupar+1)/2)
             mi_min = floor(-1.0*upara_test*(nupar-1)/2 &
                              + (nupar+1)/2)

               mi_loop2: do mi = 1, nupar

                bql = 0.0
                  cql = 0.0
                  eql = 0.0
                  fql = 0.0

                if(mi.ge.mi_min .and. mi.le.mi_max)then

                  bql = 1.0 / (8. * emax * dupara) &
                        * eps0 * omgp2(i,j) / omgrf * real(b_sum(mi))

                  cql = 1.0 / (8. * emax * dupara) &
                        * eps0 * omgp2(i,j) / omgrf * real(c_sum(mi))

                  eql = 1.0 / (8. * emax * dupara) &
                        * eps0 * omgp2(i,j) / omgrf * real(e_sum(mi))

                  fql = 1.0 / (8. * emax * dupara) &
                        * eps0 * omgp2(i,j) / omgrf * real(f_sum(mi))



!                  if(bql .ne. 0.0)count(myid, 1) = count(myid, 1) + 1.0


                  if(n .le. nnoderho)then

!                    ----------------------
!                    calculate midplane us
!                    ----------------------
                     argd = uperp(ni)**2 * (1.-bratio) + upara(mi)**2
                     if (argd .le. 0.0) argd = 1.0e-06

                     uperp0 = uperp(ni) * sqrt(bratio)
                     upara0 = sign(1.0, upara(mi)) * sqrt(argd)

                     u_  = sqrt(uperp(ni)**2 + upara(mi)**2)
                     if (u_  .eq. 0.0) u_  = 1.0e-08
                     u_0 = u_


                     ni0 = int((uperp0 - uperp(1)) / duperp) + 1
                     mi0 = int((upara0 - upara(1)) / dupara) + 1

!                    --------------------------------------
!                    bounce average and map to midplane us
!                    --------------------------------------
                     if(ni0 .ge. 1 .and. ni0 .le. nuper .and. &
                           mi0 .ge. 1 .and. mi0 .le. nupar)then



                       costh0 = upara0 / u_0
                  costh = upara(mi) / u_

                  if(costh0 .eq. 0.0)costh0 = 1.0e-08

                  upara_mi = upara(mi)
                  if(upara_mi .eq. 0.0) &
                              upara_mi = (upara(mi) + upara(mi+1)) / 2.0


!                        factor = abs(psic * upara0 / upara_mi)
!                       factor = abs(sqrt(psic))
                        factor = 1.0


!                        factvol(ni0, mi0, n) = factvol(ni0, mi0, n)
!     .                      + dx * dy * capr(i) / r0 * factor

                        bqlvol(ni0, mi0, n) = bqlvol(ni0, mi0, n) &
                            + dx * dy * capr(i) / r0 * bql * factor

                        cqlvol(ni0, mi0, n) = cqlvol(ni0, mi0, n) &
                            + dx * dy * capr(i) / r0 * cql * factor &
                            * costh / costh0 / sqrt(psic)

                        eqlvol(ni0, mi0, n) = eqlvol(ni0, mi0, n) &
                            + dx * dy * capr(i) / r0 * eql * factor &
                            * costh / costh0 / psic

                        fqlvol(ni0, mi0, n) = fqlvol(ni0, mi0, n) &
                            + dx * dy * capr(i) / r0 * fql * factor &
                            * (costh / costh0)**2 / psic**1.5

                     endif



              endif

                  else
                  endif
!DLG:   Adjust if statement for ndist >= 1
        if(i_write.eq.1)then

            !bql_store(ip,ni,mi) = bql
            !cql_store(ip,ni,mi) = cql
            !eql_store(ip,ni,mi) = eql
            !fql_store(ip,ni,mi) = fql

            if ( dumpql .eq. 1 ) then
            bql4D(i-iStart,j-jStart,ni,mi) = bql/xn(i,j)
            cql4D(i-iStart,j-jStart,ni,mi) = cql/xn(i,j)
            eql4D(i-iStart,j-jStart,ni,mi) = eql/xn(i,j)
            fql4D(i-iStart,j-jStart,ni,mi) = fql/xn(i,j)

            ql_th   = atan2 ( uPerp(ni), uPara(mi) )
            vTmp    = sqrt ( uPerp(ni)**2 + uPara(mi)**2 ) &
                * vc_mks

            bqlTmp  = bql/xn(i,j)*vc_mks**4/vTmp**2
            cqlTmp  = cql/xn(i,j)*vc_mks**3/vTmp
            eqlTmp  = eql/xn(i,j)*vc_mks**3/(vTmp*sin(ql_th))
            fqlTmp  = fql/xn(i,j)*vc_mks**2/sin(ql_th)

            if ( bqlTmp .ne. bqlTmp .or. bqlTmp*0 .ne. 0 ) bqlTmp = 0
            if ( cqlTmp .ne. cqlTmp .or. cqlTmp*0 .ne. 0 ) cqlTmp = 0
            if ( eqlTmp .ne. eqlTmp .or. eqlTmp*0 .ne. 0 ) eqlTmp = 0
            if ( fqlTmp .ne. fqlTmp .or. fqlTmp*0 .ne. 0 ) fqlTmp = 0

            bql4D_cyl(i-iStart,j-jStart,ni,mi)  = &
                 0.5 * ( bqlTmp + fqlTmp + (-bqlTmp + fqlTmp ) &
                 * cos ( 2d0 * ql_th ) &
                 + ( cqlTmp + eqlTmp ) * sin ( 2d0 * ql_th ) )
            endif

        endif


               enddo mi_loop2

            enddo ni_loop2

!DLG: Overwrite vPerp = 0 row with vPerp+1 row

        if ( dumpQL .eq. 1 ) then
            bql4D_cyl(i-iStart,j-jStart,1,:) = &
                bql4D_cyl(i-iStart,j-jStart,2,:)

        endif

          wdoti = 0.0
          fx0i = 0.0
          fy0i = 0.0

!!DLG:   Calculate R,z orbit trajectory for each uPerp/uPar that absorbed power
!!       Then distribute the power (wdot) at that uPerp/uPar
!!       along that R,z orbit to give a 2D (R,z) wdot for each
!!       uPerp/uPar that absorbs power.
!
!        if ( ndist .ge. 1 ) then
!            ba_wdot = .false.
!        else
!            ba_wdot = .false.
!        endif
!
!        bounceAveWdot: &
!        if ( ba_wdot ) then
!
!        wdoti_RZ = 0.0
!
!        mi_loop3: do mi = 1, nuper
!            ni_loop3: do ni = 1, nupar
!
!                ! Only trace orbit if power is absorbed here
!                !if ( wdot_sum_res_UU(mi,ni) .gt. 0 ) then
!                if ( abs(real(wdot_sum_UU(mi,ni))) .gt. 0.1 ) then
!
!                    startR  = capR(i)
!                    startz  = zLoc(j)
!                    vPer    = uPerp(mi) * vc_mks
!                    vPar    = uPara(ni) * vc_mks
!
!                    if ( myid .eq. 1 ) then
!                        plotVar = .false.
!                    else
!                        plotVar = .false.
!                    endif
!
!                    wdot_orbit = 0.0
!
!                    call gc_orbit ( startR, startz, &
!                            vPer, vPar, &
!                            1.0, nnodex,nnodey, &
!                            capR(nnodex)-capR(1), &
!                            zLoc(nnodey)-zLoc(1), &
!                            capR(1), zLoc(1), wdot_orbit, &
!                            goodOrbit, plot = plotVar )
!
!                    if ( goodOrbit ) then
!
!                        wdoti_RZ = wdoti_RZ + &
!                          wdot_sum_UU(mi,ni) * wdot_orbit &
!                           * duperp
!
!!                        write(*,*) maxVal(abs(real(wdoti_RZ))),
!!     .                    maxVal(wdot_orbit),duperp,
!!     .                    maxVal(capR)-minVal(capR),
!!     .                    abs(maxVal(zLoc)-minVal(zLoc))/2.0,
!!     .                    minVal(capR), size(capR), size(zLoc)
!!
!!                        write(*,*) maxVal(wdot_orbit)
!!
!                       wdot_orbit_ = wdot_orbit
!
!                    endif
!
!                endif
!
!            enddo ni_loop3
!        enddo mi_loop3
!        endif bounceAveWdot

!           --------------------------------------------
!           integrate wdot over perpendicular velocities
!           --------------------------------------------

          ni_loop4: do ni = 2, nuper - 1
             wdoti = wdoti + wdot_sum(ni) * duperp
             fx0i = fx0i + sum_fx0(ni) * duperp
             fy0i = fy0i + sum_fy0(ni) * duperp
          enddo ni_loop4


          wdoti = wdoti + 0.5 * wdot_sum(1)     * duperp &
                          + 0.5 * wdot_sum(nuper) * duperp

          fx0i = fx0i + 0.5 * sum_fx0(1)     * duperp &
                        + 0.5 * sum_fx0(nuper) * duperp

          fy0i = fy0i + 0.5 * sum_fy0(1)     * duperp &
                        + 0.5 * sum_fy0(nuper) * duperp


!!DLG:   Use bounce averaged wdot, i.e., power deposition, not absorption.
!        if ( ba_wdot ) then
!         wdot__    = wdot__ - pi / 2.0 * eps0 * omgp2 &
!                    / omgrf * real ( wdoti_RZ )
!        endif

        wdot(i,j) = - pi / 2.0 * eps0 * omgp2(i,j) &
                    / omgrf * real(wdoti)

            fx0(i,j)  = - pi / 2.0 * eps0 * omgp2(i,j) &
                            / omgrf * real(fx0i) / (2.0 * omgrf)
            fy0(i,j)  = - pi / 2.0 * eps0 * omgp2(i,j) &
                             / omgrf * real(fy0i)/ (2.0 * omgrf)

            fz0(i,j)  = xkphi / omgrf * wdot(i,j)


!     --------------------------------
!     end loop over i,j spatial points
!     --------------------------------
        !enddo spatial_loop
        !endif hasWork
        !endif inMask
        !enddo j_loop
        !enddo i_loop
        enddo ip_loop

!DLG:   Overwrite wdot with the bounce averaged version
!        if ( ba_wdot ) wdot = wdot__

!DLG:   Adjust if statement for ndist >= 1

      if( i_write .eq. 99)then

         iam_root = (myid .eq. 0)
         left_neighbor = mod( myid - 1 + nproc, nproc)
         right_neighbor = mod( myid + 1, nproc )
         token = 1.0

         if(iam_root)then

           open(unit=43,file='out_orbitrf.coef',  status='replace', &
                                                   form='formatted')

         do ip=start,finish
            i=i_table(ip)
            j=j_table(ip)
            do ni=1,nuper
              do mi=1,nupar
              if(abs(bql_store(ip,ni,mi)) .gt. 10**(-10)) then

                 write(43,102) i, j, ni, mi, &
                       bql_store(ip,ni,mi), cql_store(ip,ni,mi), &
                       eql_store(ip,ni,mi), fql_store(ip,ni,mi)
                endif
            enddo
            enddo
           enddo

         close(43)

         call MPI_SEND(token,1, MPI_REAL,right_neighbor, &
                   2, MPI_COMM_WORLD, ierr)
           call MPI_RECV(token,1, MPI_REAL,left_neighbor, &
                   2, MPI_COMM_WORLD, status, ierr)
         else

           call MPI_RECV(token,1, MPI_REAL,left_neighbor, &
                   2, MPI_COMM_WORLD, status, ierr)

           open(unit=43,file='out_orbitrf.coef',  status='old', &
                          form='formatted', position='append')

         do ip=start,finish
            i=i_table(ip)
            j=j_table(ip)
            do ni=1,nuper
              do mi=1,nupar
              if(abs(bql_store(ip,ni,mi)) .gt. 10**(-10)) then

                 write(43,102) i, j, ni, mi, &
                       bql_store(ip,ni,mi), cql_store(ip,ni,mi), &
                       eql_store(ip,ni,mi), fql_store(ip,ni,mi)
                endif
            enddo
            enddo
           enddo

         close(43)

         call MPI_SEND(token,1, MPI_REAL,right_neighbor, &
                   2, MPI_COMM_WORLD, ierr)

         endif
      endif

      call blacs_barrier(icontxt, 'All')

!     -------------------
!     Sum over processors
!     -------------------

!      call dgsum2d(icontxt, 'All', ' ', 5001, 1, count,
!     .      5001, -1, -1)

      call dgsum2d(icontxt, 'All', ' ', nnodex, nnodey, wdot, &
            nnodex, -1, -1)

      call dgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fx0, &
            nnodex, -1, -1)

      call dgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fy0, &
            nnodex, -1, -1)

      call dgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fz0, &
            nnodex, -1, -1)

      do n = 1, nnoderho

         do mi = 1, nupar
            do ni = 1, nuper
!              factvol2d(ni, mi) = factvol(ni, mi, n)
               bqlvol2d(ni, mi) = bqlvol(ni, mi, n)
               cqlvol2d(ni, mi) = cqlvol(ni, mi, n)
               eqlvol2d(ni, mi) = eqlvol(ni, mi, n)
               fqlvol2d(ni, mi) = fqlvol(ni, mi, n)
            enddo
         enddo

!         call dgsum2d(icontxt, 'All', ' ', nuper, nupar, factvol2d,
!     .      nuper, -1, -1)
         call dgsum2d(icontxt, 'All', ' ', nuper, nupar, bqlvol2d, &
            nuper, -1, -1)
         call dgsum2d(icontxt, 'All', ' ', nuper, nupar, cqlvol2d, &
            nuper, -1, -1)
         call dgsum2d(icontxt, 'All', ' ', nuper, nupar, eqlvol2d, &
            nuper, -1, -1)
         call dgsum2d(icontxt, 'All', ' ', nuper, nupar, fqlvol2d, &
            nuper, -1, -1)




         do mi = 1, nupar
            do ni = 1, nuper
!             factvol(ni, mi, n) = factvol2d(ni, mi)
               bqlvol(ni, mi, n) = bqlvol2d(ni, mi)
               cqlvol(ni, mi, n) = cqlvol2d(ni, mi)
               eqlvol(ni, mi, n) = eqlvol2d(ni, mi)
               fqlvol(ni, mi, n) = fqlvol2d(ni, mi)
            enddo
         enddo


      enddo


!     ------------------------------------
!     Divide by volume element (passed in)
!     ------------------------------------
      do n = 1, nnoderho
         do mi = 1, nupar
            do ni = 1, nuper

!           --------------
!           bounce average
!           --------------
            if (vol(n) .ne. 0.0)then
            bqlavg(ni, mi, n) = bqlvol(ni, mi, n) / vol(n) * dldbavg(n)
            cqlavg(ni, mi, n) = cqlvol(ni, mi, n) / vol(n) * dldbavg(n)
            eqlavg(ni, mi, n) = eqlvol(ni, mi, n) / vol(n) * dldbavg(n)
            fqlavg(ni, mi, n) = fqlvol(ni, mi, n) / vol(n) * dldbavg(n)
          end if

          if (bqlavg(ni, mi, n) .lt. 0.0)bqlavg(ni, mi, n) = 0.0
          if (cqlavg(ni, mi, n) .lt. 0.0)cqlavg(ni, mi, n) = 0.0
          if (eqlavg(ni, mi, n) .lt. 0.0)eqlavg(ni, mi, n) = 0.0
          if (fqlavg(ni, mi, n) .lt. 0.0)fqlavg(ni, mi, n) = 0.0

            enddo
         enddo
      enddo

!DLG:   Try pNetCDF to write the whole QL operator ;-)

        if ( i_write .eq. 1 ) then

            call blacs_barrier(icontxt, 'All')

            what    = 10
            call blacs_get ( icontxt, what, ivalue )

            if (myId .eq. 0) write(*,*) 'pre: ', vc_mks
            call write_pql_dlg ( iValue, myId, nNodeX, nNodeY, &
                nuPer, nuPar, &
                vc_mks, capR, zLoc, uPerp, uPara, &
                npRowOut, npColout, myRow, myCol, &
                bql4d_cyl, eNormIN )

            call blacs_barrier(icontxt, 'All')

        end if

    !   Finished p-netCDF ql dump.


      deallocate( dfduper0 )
      deallocate( dfdupar0 )


      deallocate(b_sum)
      deallocate(c_sum)
      deallocate(e_sum)
      deallocate(f_sum)

!DLG:   Deallocate
      deallocate(wdot_sum_UU)
      deallocate(wdot_orbit)
      deallocate(wdot_orbit_)
      deallocate(wdoti_RZ)

      deallocate(wdot_sum)
      deallocate(sum_fx0)
      deallocate(sum_fy0)

      deallocate(factvol)
      deallocate(factvol2d)

      deallocate(bqlvol)
      deallocate(bqlvol2d)

      deallocate(cqlvol)
      deallocate(cqlvol2d)

      deallocate(eqlvol)
      deallocate(eqlvol2d)

      deallocate(fqlvol)
      deallocate(fqlvol2d)
!DLG:   Adjust if statement for ndist >= 1
      if(i_write.ge.1) then
         !deallocate(bql_store, cql_store, eql_store,  fql_store)
         deallocate(bql4D, cql4D, eql4D,  fql4D)
         deallocate(bql4D_cyl)

      endif

      wdot_inout(1:nnodex,1:nnodey) = wdot(1:nnodex,1:nnodey)
      fx0_inout(1:nnodex,1:nnodey) = fx0(1:nnodex,1:nnodey)
      fy0_inout(1:nnodex,1:nnodey) = fy0(1:nnodex,1:nnodey)
      fz0_inout(1:nnodex,1:nnodey) = fz0(1:nnodex,1:nnodey)

!      if (myid .eq. 0)then
!         sum_count = 0.0
!         do n = 0, 5000
!            write(43, 101)n, count(n, 1)
!          sum_count = sum_count + count(n, 1)
!         end do
!       write(43, *)"sum_count = ", sum_count
!      end if

      return

 1311 format(1p9e12.4)
  100 format (1p8e12.4)
  101 format (1i6, 1p8e12.4)
  102 format (4i6, 1p8e12.4)

      end subroutine ql_myra_write


!
!***************************************************************************
!


      subroutine wdot_qlcheck(wdot_check, &
         nnoderho, nrhodim, &
         bqlavg, cqlavg, xm, omgrf, xktavg, ndist, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper, dfdupar, &
         UminPara_cql, UmaxPara_cql, UPERP_cql, UPARA_cql, UPERP, UPARA, &
         vc_mks_cql, df_cql_uprp, df_cql_uprl, rhon, rho_a, myid, &
         dldbavg, eNormIN)

!     ------------------------------------------------------------------
!     This subroutine calculates the flux averaged wdot for checking QL
!     ------------------------------------------------------------------

      use read_particle_f
      use aorsa2din_mod, only: eNorm_factor
      implicit none

      real :: eNormIN
      integer nnoderho, nrhodim, n, m, ndist, myid
      real bqlavg(nuper, nupar, nnoderho)
      real cqlavg(nuper, nupar, nnoderho)
      real xktavg(nnoderho), alpha, e
      real ans, dfdth, dfdu, u, xm, omgrf, jacobian, upara_mi
      real, dimension(:,:), allocatable :: wdot_int
      real wdot_check(nrhodim), rhon(nrhodim)
      real dldbavg(nrhodim)

      integer  :: n_psi_dim, nuper, nupar, n_psi, mi0, ni0

      real :: UPERP(NUPER), UPARA(NUPAR)
      real :: UPERP_cql(NUPER), UPARA_cql(NUPAR)
      real :: DFDUPER(NUPER,NUPAR),DFDUPAR(NUPER,NUPAR)
      real :: DFDUPER0, DFDUPAR0, UPARA0
      real :: W, ENORM
      real :: UminPara,UmaxPara
      real :: UminPara_cql,UmaxPara_cql

      real :: df_cql_uprp(NUPER, NUPAR, n_psi_dim)
      real :: df_cql_uprl(NUPER, NUPAR, n_psi_dim)


      real :: vc_mks, vc_mks_cql, rho_a(n_psi_dim)
      real :: eps0, pi, emax, u0, dfdu0, dfdth0, u_0

      parameter (eps0 = 8.85e-12)
      parameter (PI = 3.141592653597932384)

!DLG:   Define variables
      real, allocatable :: aorsa_dfduPer_rho(:,:,:), &
       aorsa_dfduPar_rho(:,:,:)

!DLG:   Allocate arrays
      allocate ( aorsa_dfduPer_rho (nnoderho,nuper,nupar) )
      allocate ( aorsa_dfduPar_rho (nnoderho,nuper,nupar) )

      allocate(wdot_int(nuper, nupar) )

!     ------------------------------------
!efd  initialize allocatable array to zero
!     ------------------------------------

      wdot_int = 0.0

      e = 1.6e-19
      W = omgrf

      if(ndist .eq. 0)then   !--Maxwellian--!

         do n = 1, NUPER
            UPERP(n) = (real(n-1)/real(NUPER-1))
         end do

         UminPara = -1.0
         UmaxPara =  1.0

         do m = 1, NUPAR
            UPARA(m) = (-1.0 + 2. * (real(m-1) / real(NUPAR-1)))
         end do

      else   !--non-Maxwellian--!

         vc_mks = vc_mks_cql
         UminPara = UminPara_cql
         UmaxPara = UmaxPara_cql

         do n = 1, nuper
            uperp(n) = uperp_cql(n)
         end do

         do m = 1, nupar
            upara(m) = upara_cql(m)
         end do

      end if

!DLG:   Read in the flux surfaced averaged dfduPer and dfduPar for
!       the flux surface average QL calculation

      if ( ndist .eq. 2 ) then
        call dlg_particle_f ( 0.0, 0.0, &
          uPerp, uPara, &
          nuper, nupar, xm, eNormIN, &
          aorsa_rho = rhon(1:nnoderho), &
          aorsa_dfduPer_rho = aorsa_dfduPer_rho, &
          aorsa_dfduPar_rho = aorsa_dfduPar_rho  )

        write(*,*) 'ql_myra.F:1302', sum(aorsa_dfduPer_rho), &
            sum(aorsa_dfduPar_rho)

      end if

!     -------------------
!     Loop over rho mesh:
!     -------------------
      do n = 1, nnoderho

         alpha = sqrt(2.0 * xktavg(n) / xm)
!DLG: Set maxwellian v grid with enorm_factor as well
         if (ndist .eq. 0) then

          if ( eNormIN < 0 ) then
            vc_mks= 3.0*alpha    !--Maxwellian only--!
          else
            vc_mks = vc_mks_cql
          endif

            endif

         u0 = vc_mks / alpha

         Emax = 0.5 * xm * vc_mks**2
         Enorm = Emax / 1.6e-19

!        ------------------------------------------------
!        get CQL3D distribution function on the midplane
!        ------------------------------------------------

         if(ndist .eq. 0)then   !--Maxwellian--!

            call maxwell_dist(u0, NUPAR, NUPER, &
                       UminPara, UmaxPara, &
                       UPERP, UPARA, DFDUPER, DFDUPAR)

         else   !--non-Maxwellian--!

            call cql3d_dist(nupar, nuper, n_psi, &
                       n_psi_dim, rho_a, rhon(n), &
                       UminPara,UmaxPara, &
                       df_cql_uprp, df_cql_uprl, &
                       UPERP, UPARA, DFDUPER, DFDUPAR)

         end if

!DLG:   Overwrite dfduPer and dfduPar in wdot_qlcheck
         if ( ndist .eq. 2 ) then

            dfduPer = aorsa_dfduPer_rho(n,:,:)
            dfduPar = aorsa_dfduPar_rho(n,:,:)

         endif

!        -----------------------------
!        loop over MIDPLANE velocities
!        -----------------------------

         do ni0 = 1, nuper
            do mi0 = 1, nupar

               dfdupar0 = dfdupar(ni0, mi0)
               dfduper0 = dfduper(ni0, mi0)

               u_0 = sqrt(uperp(ni0)**2 + upara(mi0)**2)
               if (u_0 .eq. 0.0) u_0 = 1.0e-08

               dfdu0 = (uperp(ni0) * dfduper0 + upara(mi0) * dfdupar0) &
                      / u_0
               dfdth0 = upara(mi0) * dfduper0 - uperp(ni0) * dfdupar0

               wdot_int(ni0, mi0) = (bqlavg(ni0, mi0, n) * dfdu0 &
                                   + cqlavg(ni0, mi0, n) * dfdth0) / u_0

            end do
         end do



!        ---------------------------------------------------
!        Do velocity space integral over midplane velocities
!        ---------------------------------------------------

         wdot_check(n) = 0.0


         call ugrate(wdot_int, uperp, upara, nuper, nupar, ans, myid,xm)


         wdot_check(n) = - 4.0 * pi * e * enorm / dldbavg(n) * ans


      end do


      deallocate(wdot_int)

!DLG:   Deallocate arrays
      deallocate ( aorsa_dfduPer_rho, aorsa_dfduPar_rho )


      return

 1311 format(1p9e12.4)
 1312 format(i10, 1p9e12.4)
 1313 format(2i10, 1p9e12.4)
  100 format (1p8e12.4)
  101 format (2i10, 1p8e12.4)

      end subroutine


!
!***************************************************************************
!



      subroutine ugrate(f, uperp, upara, nuper, nupar, fint, myid, xm)

      implicit none

      integer nuper, nupar, ni, mi, myid

      real uperp(nuper), upara(nupar), f(nuper, nupar)
      real favg, fint, duperp, dupara, xm

      duperp = (uperp(nuper) - uperp(1)) / (nuper - 1)
      dupara = (upara(nupar) - upara(1)) / (nupar - 1)

      fint = 0.0

      do ni = 1, nuper - 1
         do mi = 1, nupar - 1
            favg = (f(ni, mi)   + f(ni+1, mi) &
                  + f(ni, mi+1) + f(ni+1, mi+1)) / 4.0

            fint = fint + favg * uperp(ni) * duperp * dupara

!            if (myid.eq.0 .and. ni .eq. 32 .and. mi .eq. 95) then
!               write(6 ,1313)ni, mi, xm, favg, fint
!               write(15,1313)ni, mi, xm, favg, fint
!          end if

         end do
      end do

 1313 format(2i10, 1p9e12.4)

      return
      end subroutine  ugrate

!
!***************************************************************************
!

       end module ql_myra_mod
