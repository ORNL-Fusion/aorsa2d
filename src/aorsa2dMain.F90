
!-------------------------------------
!     AORSA2D: version 15 (02-29-2008)
!-------------------------------------
      program aorsa2dMain

      use size_mod
      use ql_myra_mod
      use profile_mod
      use aorsa2din_mod
      use eqdsk_dlg
      use gc_terms
      use interp
      use write_pf
      use plot_aorsa2dps
      use netcdf
      use read_particle_f
      use set_edge_density

!------------------------------------------------------------------------------
!     This version (4/01/03: newlab3) does not solve equations for the points
!     outside psi = 1.
!     Stix2 frame using SCALAPACK with new faster loading by E. DAzevedo
!------------------------------------------------------------------------

      implicit none

      integer n_theta_max, n_theta, n_psi_max, nt, k, nbessj

      integer :: nphi1 = -50
      integer :: nphi2 =  50
      real, parameter :: e = 1.6e-19
      integer :: what, ivalue, pnc_f_id

      real ydiff, ydiff_min, xkprl, xkprl_max, xkprl_min
      real dupara, duperp, dvperp, dvpara, vpara_in
      real vperp_in, xmi_slo
      real vce_mks_cql3d, vc1_mks_cql3d, vc2_mks_cql3d
      real vce_cgs_cql3d, vc1_cgs_cql3d, vc2_cgs_cql3d, vc3_cgs_cql3d
      real vc4_cgs_cql3d
      real bqlavg_i2_out, cqlavg_i2_out, eqlavg_i2_out, fqlavg_i2_out
      real xkalp, xkbet
      real xkprl_neg_max, xkprl_pos_min
      real sgn_kprl, akprl, gammab, y0, yb, fgam, xkprl_eff, descrim
      real xk_cutoff, term2, rmaxis, zmaxis, z0, ptot_wdot, xj

      integer jmid, nxdim, nydim, l, lmaxdim
      integer nnodex_loc, nx_overlap, nnodex_eff
      integer nnodey_loc, ny_overlap, nnodey_eff
      integer istart, ifinish, jstart, jfinish, &
        i_seg, j_seg, iseg1, iseg2, jseg1, jseg2, i0, j0

      integer nxmx, nymx, ieq
      integer incX, global_row, global_col
      integer ndfmax, &
          ninteg, nd, izoom1, izoom2, jzoom1, jzoom2, ndf, nmaxe, irnc
      integer nrow, ncol, norder

      integer mkdim1, mkdim2

      integer nkdim1, nkdim2, nkx1, nkx2, nldim, nldim3, &
         nky1, nky2, iant

      real :: kperp_max, valfven, kalfven, xky_ant
      real  kperp_max_actual
      real kperprhoi1, kperprhoi2, kperprhoi3, kperprhoi4, kperprhoi5

      real eslow, qi_slo, xmi_sl, yant_max, &
          xn_slo, a, epsa, epsa_potato, rmaxa, rmina

      real flimiter, parabola, dfdx, dfdy, caprmax, caprmin, rmax

      real tmem, tsys, tio, ttotal, time0, time, cpu, dummy, second1, dlgTime

      real tmin, gflops, gflopsp, ops, teev, t1ev, t2ev, t3ev
      real t4ev, t5ev, t6ev

      real sqx, signb, gausspsi, dpsiant, theta_antr

      real dxsqx, dxxsqx, dxysqx, dysqx, dyysqx, UminPara, UmaxPara

      real vce_mks, vc1_mks, vc2_mks, vc3_mks, vc4_mks, vc5_mks, vc6_mks
      real vce_cgs, vc1_cgs, vc2_cgs, vc3_cgs, vc4_cgs, vc5_cgs, vc6_cgs

      integer nrhomax, nthetamax, ndim

      parameter (n_theta_max = 200)
      parameter (n_psi_max = 200)

      parameter (nxmx = nmodesmax)
      parameter (nymx = mmodesmax)
      parameter (nrhomax = nmodesmax * 2)
      parameter (nthetamax = nmodesmax * 2)

      parameter (nkdim1 = - nmodesmax / 2)
      parameter (nkdim2 =   nmodesmax / 2)

      parameter (mkdim1 = - mmodesmax / 2)
      parameter (mkdim2 =   mmodesmax / 2)

      parameter (nldim  = nxmx * nymx)
      parameter (nldim3 = 3 * nldim)
      parameter (ndfmax = nldim3)

      parameter (lmaxdim = 25)

      parameter (ndim = 256)

      integer nsum

      real rhod(ndim), xneavg(ndim), tekev(ndim), tikev(ndim), &
           omega(ndim), zeffavg(ndim), xn_beam(ndim), e_beam(ndim), &
           drho_data, rhot, ti2kev(ndim), xn_maj(ndim), xn_min(ndim), &
           xn_min_cgs(ndim), xn_maj_cgs(ndim), xneavg_cgs(ndim), &
           xn_beam_cgs(ndim)

      real rhod2(ndim), xneavg2(ndim), tekev2(ndim), tikev2(ndim), &
          xn_beam2(ndim), e_beam2(ndim), zeffavg2(ndim)

      integer, parameter :: n_psi_dim = 64

      integer :: n_psi

      real :: dens(n_psi_dim), rho_a(n_psi_dim)

      real, dimension(:),   allocatable :: UPERP, UPARA
      real, dimension(:),   allocatable :: VPERP, VPARA
      real, dimension(:),   allocatable :: UPERP_work, UPARA_work

      real, dimension(:,:), allocatable :: DFDUPERe, DFDUPARe
      real, dimension(:,:), allocatable :: DFDUPER1, DFDUPAR1
      real, dimension(:,:), allocatable :: DFDUPER2, DFDUPAR2
      real, dimension(:,:), allocatable :: DFDUPER3, DFDUPAR3
      real, dimension(:,:), allocatable :: DFDUPER4, DFDUPAR4
      real, dimension(:,:), allocatable :: DFDUPER5, DFDUPAR5
      real, dimension(:,:), allocatable :: DFDUPER6, DFDUPAR6


      real, dimension(:,:), allocatable :: f

      real, dimension(:,:,:), allocatable :: fe_cql_cart
      real, dimension(:,:,:), allocatable :: dfe_cql_uprp
      real, dimension(:,:,:), allocatable :: dfe_cql_uprl

      real, dimension(:,:,:), allocatable :: f1_cql_cart
      real, dimension(:,:,:), allocatable :: df1_cql_uprp
      real, dimension(:,:,:), allocatable :: df1_cql_uprl

      real, dimension(:,:,:), allocatable :: f2_cql_cart
      real, dimension(:,:,:), allocatable :: df2_cql_uprp
      real, dimension(:,:,:), allocatable :: df2_cql_uprl

      real, dimension(:,:,:), allocatable :: f3_cql_cart
      real, dimension(:,:,:), allocatable :: df3_cql_uprp
      real, dimension(:,:,:), allocatable :: df3_cql_uprl

      real, dimension(:,:,:), allocatable :: f4_cql_cart
      real, dimension(:,:,:), allocatable :: df4_cql_uprp
      real, dimension(:,:,:), allocatable :: df4_cql_uprl

      real, dimension(:,:,:), allocatable :: f5_cql_cart
      real, dimension(:,:,:), allocatable :: df5_cql_uprp
      real, dimension(:,:,:), allocatable :: df5_cql_uprl

      real, dimension(:,:,:), allocatable :: f6_cql_cart
      real, dimension(:,:,:), allocatable :: df6_cql_uprp
      real, dimension(:,:,:), allocatable :: df6_cql_uprl

      real, dimension(:,:,:), allocatable :: bqlavg_e
      real, dimension(:,:,:), allocatable :: cqlavg_e
      real, dimension(:,:,:), allocatable :: eqlavg_e
      real, dimension(:,:,:), allocatable :: fqlavg_e

      real, dimension(:,:,:), allocatable :: bqlavg_i1
      real, dimension(:,:,:), allocatable :: cqlavg_i1
      real, dimension(:,:,:), allocatable :: eqlavg_i1
      real, dimension(:,:,:), allocatable :: fqlavg_i1

      real, dimension(:,:,:), allocatable :: bqlavg_i2
      real, dimension(:,:,:), allocatable :: cqlavg_i2
      real, dimension(:,:,:), allocatable :: eqlavg_i2
      real, dimension(:,:,:), allocatable :: fqlavg_i2

      real, dimension(:,:,:), allocatable :: bqlavg_i3
      real, dimension(:,:,:), allocatable :: cqlavg_i3
      real, dimension(:,:,:), allocatable :: eqlavg_i3
      real, dimension(:,:,:), allocatable :: fqlavg_i3

      real, dimension(:,:,:), allocatable :: bqlavg_i4
      real, dimension(:,:,:), allocatable :: cqlavg_i4
      real, dimension(:,:,:), allocatable :: eqlavg_i4
      real, dimension(:,:,:), allocatable :: fqlavg_i4

      real, dimension(:,:,:), allocatable :: bqlavg_i5
      real, dimension(:,:,:), allocatable :: cqlavg_i5
      real, dimension(:,:,:), allocatable :: eqlavg_i5
      real, dimension(:,:,:), allocatable :: fqlavg_i5

      real, dimension(:,:,:), allocatable :: bqlavg_i6
      real, dimension(:,:,:), allocatable :: cqlavg_i6
      real, dimension(:,:,:), allocatable :: eqlavg_i6
      real, dimension(:,:,:), allocatable :: fqlavg_i6

      real, dimension(:,:,:), allocatable :: bqlavg_work
      real, dimension(:,:,:), allocatable :: cqlavg_work
      real, dimension(:,:,:), allocatable :: eqlavg_work
      real, dimension(:,:,:), allocatable :: fqlavg_work


      integer :: i_uperp, i_upara, i_psi, i_psi_eq

      real factl_sav(0:lmaxdim), factrl

      character*1 trans

      integer info
      complex b, zi, cexpkxky
      complex alpha
      real norm2,norm2_max

      complex fdk, fek, ffk, fgk, fak, fpk, frk, fqk, fsk

      real xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2), &
           xkperp, xketa, &
           ptot, pcito2, pcrto2, powtot, pscale, &
           cosalp, sinalp, t1, gaussian, psimax, &
           rhomax, xnprl, rhoant, gaussantx, shapey, &
           dthetant, gaussantth, gaussiant, gaussthetap, &
           constantth, abstheta0, xnustar, thetamax, dtheta

      real gaussian_ne,  gaussian_te,  gaussian_ti1, gaussian_ti2, &
           gaussian_ti3, gaussian_ti4, gaussian_ti5, gaussian_ti6


      real shape, shapen,  shapen2,  shapen3, shapen_slo, &
         shapen4, shapen5, shapen6, &
         shapete, shapeti, shapeti2, shapeti3, &
         shapeti4, shapeti5, shapeti6

      real frho, dxfrho, dxxfrho, dxyfrho, dyfrho, dyyfrho


      complex xx(nkdim1 : nkdim2, 1 : nxmx), &
              yy(mkdim1 : mkdim2, 1 : nymx)

      real dldb_tot12(nxmx, nymx)

      real theta_(n_theta_max, n_psi_max)
      integer n_theta_(n_psi_max)

      integer ilayer,nlayer
      integer mask(nxmx, nymx), mask2(nxmx,nymx)
      integer mask_bbbs(nxmx,nymx)

      complex ealphak(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              ebetak(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              ebk(nkdim1 : nkdim2, mkdim1 : mkdim2)

      real ealphakmod(nkdim1 : nkdim2, mkdim1 : mkdim2), &
           ebetakmod(nkdim1 : nkdim2, mkdim1 : mkdim2), &
           ebkmod(nkdim1 : nkdim2, mkdim1 : mkdim2)

      complex fdksav(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              feksav(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              ffksav(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              fgksav(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              faksav(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              fpksav(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              frksav(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              fqksav(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              fsksav(nkdim1 : nkdim2, mkdim1 : mkdim2)

      complex fdksav2d(nxmx, nymx), &
              feksav2d(nxmx, nymx), &
              ffksav2d(nxmx, nymx), &
              fgksav2d(nxmx, nymx), &
              faksav2d(nxmx, nymx), &
              fpksav2d(nxmx, nymx), &
              frksav2d(nxmx, nymx), &
              fqksav2d(nxmx, nymx), &
              fsksav2d(nxmx, nymx)

      integer nnodex, nnodey, i, j, &
          jequat, iflag, liw, lw, nrhs, icenter
      integer nnodex_p, nnodey_p

      integer nnoderho, mnodetheta

      integer n, m, nphi, ndata, ndata2, nsmooth

      real xant, delta, xnurf, psio, psimag, psi_tor_max

      real telimj, tilimj, ti2limj, ti3limj, &
                 ti4limj, ti5limj, ti6limj


      real rhoplasm, reomg1, reomg2, reomg3, &
          reomg4, reomg5, reomg6, &
          r0, xiota0, rholim, &
          deltay, abs_deltay, delta_theta, &
          xwleft, xwright, psi_lim, psi1

      real xkphi(nxmx), xktau, xkrho, xkphi0, xnphi
      real rlim

      real xkthrho, wphase, vsound, domgk, xkthdx, omgestar, rhoi10, &
         v0i, vthi10, vthe, vthi, vphase, rhoi1overl, rnz, eta2, &
         xk0, shearedge, eta1, xn1, xmi1, xmh, xme, qi1, xmi2, &
         xmi3, t0i, t0i2, t0i3, t0e, q, teedge, clight, xmu0, xlnlam, &
         omgci10, omgrf, xmax, qe, qi2, pi, twopi, eps0, qi3, &
         costh, sinth, radius, rnx,  rny, rnphi

      real rhoi20, vthi20, omgci20
      real rhoi30, vthi30, omgci30
      real rhoi40, vthi40, omgci40
      real rhoi50, vthi50, omgci50



      real xmi4, t0i4, qi4
      real xmi5, t0i5, qi5
      real xmi6, t0i6, qi6



      real xjantx, xjanty, xjantz, xjant
      real xjx(nxmx, nymx), xjy(nxmx, nymx), xjz(nxmx, nymx)
      real djxdx, djydy
      complex rho_ant(nxmx, nymx), rho_pla(nxmx, nymx)

      complex xb(nxmx, nymx), xc(nxmx, nymx), xd(nxmx, nymx)


      real xprimec(nxmx), caprc(nxmx), xcourse(nxmx), capr(nxmx), &
         xprime(nxmx), x(nxmx), dx, dxc

      real capr_flux(nrhomax, nthetamax), capz_flux(nrhomax, nthetamax), &
         xgiv, ygiv
      complex eplus_flux(nrhomax, nthetamax), fout
      complex eminus_flux(nrhomax, nthetamax)
      real bmod_flux(nrhomax, nthetamax)

      real rhon(nrhomax), thetam(nthetamax), rho_in, theta_in, &
         capr_out, capz_out, &
         wdoti1avg(nrhomax), wdoti2avg(nrhomax), &
         wdoti3avg(nrhomax), wdoti4avg(nrhomax), &
         dvol(nrhomax), darea(nrhomax), rn(nrhomax), rh(nrhomax), &
         dn(nrhomax), dh(nrhomax), dr, an(nrhomax), bn(nrhomax), &
         cn(nrhomax), yn(nrhomax), epsn(nrhomax), &
         wdoteavg(nrhomax), wdotavg(nrhomax), drho, dpsi, &
         capr_bpol_mid(nrhomax), capr_bpol_midavg(nrhomax), &
         bmod_midavg(nrhomax), fvol(nrhomax), dvolp(nrhomax), &
         volume(nrhomax)
      real zeffavg1(nrhomax)

      real ank(nrhomax), bnk(nrhomax), cnk(nrhomax), ynk(nrhomax)


      real  wdote_ql(nrhomax),  wdoti1_ql(nrhomax), wdoti2_ql(nrhomax)
      real  wdoti3_ql(nrhomax), wdoti4_ql(nrhomax), wdoti5_ql(nrhomax)
      real  wdoti6_ql(nrhomax)

      real wdoti5avg(nrhomax), wdoti6avg(nrhomax)

      real xnavg(nrhomax), xn1avg(nrhomax), xn2avg(nrhomax), &
           xn3avg(nrhomax), xna_sloavg(nrhomax)
      real xn4avg(nrhomax), xn5avg(nrhomax), xn6avg(nrhomax)


      real xkteavg(nrhomax), xktiavg(nrhomax),  xkti2avg(nrhomax), &
           xkti3avg(nrhomax)
      real xkti4avg(nrhomax), xkti5avg(nrhomax), xkti6avg(nrhomax)
      real dldbavg(nrhomax)


      real fyavg(nrhomax), qhat, omgte, omgti, &
           vthe0, vthi0, xnuee, xnuii, xnu7omg

      real redotj1avg(nrhomax), redotj2avg(nrhomax), &
           redotjeavg(nrhomax), redotj3avg(nrhomax), &
           redotj4avg(nrhomax), &
           redotj5avg(nrhomax), redotj6avg(nrhomax), &
           redotjsavg(nrhomax), &
           redotjtavg(nrhomax), xjprlavg(nrhomax)

      real xjprl_int(nrhomax), fyp_int(nrhomax), fz0_int(nrhomax)

      real wdote_ql_int(nrhomax), wdoti1_ql_int(nrhomax), &
                                  wdoti2_ql_int(nrhomax), &
                                  wdoti3_ql_int(nrhomax), &
                                  wdoti4_ql_int(nrhomax), &
                                  wdoti5_ql_int(nrhomax), &
                                  wdoti6_ql_int(nrhomax)

      real wdoteavg_int(nrhomax), wdoti1avg_int(nrhomax), &
                                  wdoti2avg_int(nrhomax), &
                                  wdoti3avg_int(nrhomax), &
                                  wdoti4avg_int(nrhomax), &
                                  wdoti5avg_int(nrhomax), &
                                  wdoti6avg_int(nrhomax)


      real redotjeavg_int(nrhomax), redotj1avg_int(nrhomax), &
                                    redotj2avg_int(nrhomax), &
                                    redotj3avg_int(nrhomax), &
                                    redotj4avg_int(nrhomax), &
                                    redotj5avg_int(nrhomax), &
                                    redotj6avg_int(nrhomax)

      real fypavg(nrhomax), fypi1avg(nrhomax), fypi2avg(nrhomax), &
           fypeavg(nrhomax), fypi3avg(nrhomax)

      real fypi4avg(nrhomax), fypi5avg(nrhomax), fypi6avg(nrhomax)


      real fz0avg(nrhomax), fz0i1avg(nrhomax), fz0i2avg(nrhomax), &
           fz0eavg(nrhomax), fz0i3avg(nrhomax)

      real fz0i4avg(nrhomax), fz0i5avg(nrhomax), fz0i6avg(nrhomax)


      real capr_fzeta_avg(nrhomax), capr2_avg(nrhomax)
      real gpsi_avg(nrhomax), xjhat(nrhomax), dgdpsi_avg(nrhomax)
      real kpsi_avg(nrhomax), xkhat(nrhomax)

      real dpdpsi(nrhomax)
      real epsi_avg(nrhomax), dedpsi_avg(nrhomax)
      real epsig_avg(nrhomax), epsip_avg(nrhomax), epsir_avg(nrhomax)

      real dgdpsi(nxmx, nymx), dedpsi(nxmx, nymx)
      real vyavg(nrhomax),  vyi1avg(nrhomax),  vyi2avg(nrhomax)
      real muhat_avg(nrhomax), nu_star_avg(nrhomax), &
           pressiavg(nrhomax), &
           rhom1avg(nrhomax), qsafety_avg(nrhomax), ipsi_avg(nrhomax), &
           gamma_avg(nrhomax), beta_avg(nrhomax)

      real taup_b, taup_pl, taup_ps, taup, tauii, taup_inv


      real dvydrho(nrhomax), psi_dim_avg(nrhomax), &
           fpsi1_avg(nrhomax), ftheta1_avg(nrhomax), &
           gradprlb_avg(nrhomax), &
           bdotf_avg(nrhomax), bmod2_avg(nrhomax)

      real drhodx, drhodxx, drhody, drhodyy, gradrho

      real dxpsi, dxxpsi, dxypsi, dypsi, dyypsi

      real drhodr(nxmx, nymx), drhodz(nxmx, nymx), &
           dthedr(nxmx, nymx),  dthedz(nxmx, nymx)

      real reomg1a(nxmx, nymx), reomg2a(nxmx, nymx), reomg3a(nxmx, nymx)

      real psi(nxmx, nymx), rho(nxmx, nymx), theta(nxmx, nymx), &
           rhohatx(nxmx, nymx), rhohaty(nxmx, nymx), &
           theta0(nxmx, nymx), psi_dim(nxmx, nymx), &
           bx(nxmx, nymx), by(nxmx, nymx), bz(nxmx, nymx), &
           btau(nxmx, nymx), bzeta(nxmx, nymx), &
           dxdth(nxmx, nymx), dzdth(nxmx, nymx), xntau(nxmx, nymx), &
           xkte(nxmx, nymx), xkti(nxmx, nymx), &
           xkti2(nxmx, nymx), xkti3(nxmx, nymx), &
           xkti4(nxmx, nymx), xkti5(nxmx, nymx), xkti6(nxmx, nymx), &
           xn1a(nxmx, nymx), xnea(nxmx, nymx), xn2a(nxmx, nymx), &
           xn3a(nxmx, nymx), bmod(nxmx, nymx), omgce(nxmx, nymx), &
           xn4a(nxmx, nymx), xn5a(nxmx, nymx), xn6a(nxmx, nymx), &
           omgci1(nxmx, nymx), omgci2(nxmx, nymx), omgci3(nxmx, nymx), &
           omgci4(nxmx, nymx), omgci5(nxmx, nymx), omgci6(nxmx, nymx), &
           omgpe2(nxmx, nymx), bpol(nxmx, nymx), capr_bpol(nxmx, nymx), &
           omgp12(nxmx, nymx), omgp22(nxmx, nymx), omgp32(nxmx, nymx), &
           omgp42(nxmx, nymx), omgp52(nxmx, nymx), omgp62(nxmx, nymx), &
           xiota(nxmx, nymx), qsafety(nxmx, nymx), bmod_mid(nxmx, nymx), &
           capr_bpol_mid2(nxmx, nymx), rho_tor2d(nxmx, nymx), &
           psi_tor2d(nxmx, nymx)

      real zeff(nxmx, nymx)

      real rho_pol2d(nxmx, nymx), psi_pol2d(nxmx, nymx)

      real dxbx, dxxbx, dxybx, dybx, dyybx
      real dxby, dxxby, dxyby, dyby, dyyby
      real dxbz, dxxbz, dxybz, dybz, dyybz

      real dxbmod, dxxbmod, dxybmod, dybmod, dyybmod

      real xna_slo(nxmx, nymx), &
           omgp2_slo(nxmx, nymx), omgci_slo(nxmx, nymx)

      real rhome, rhomi1, rhomi2, rhomi3, rhomslo
      real rhomi4, rhomi5, rhomi6, rhom1(nxmx, nymx)

      real pressi(nxmx, nymx), presse, pressi1, pressi2, pressi3, &
                                       pressi4, pressi5, pressi6

      real muhat(nxmx, nymx), prod, nu_star(nxmx, nymx), &
         ipsi(nxmx, nymx)
      real dbdx, dbdy, d2bdx2, d2bdy2, gradprlb(nxmx, nymx)
      real bmod2(nxmx, nymx)

      real bxn(nxmx, nymx), byn(nxmx, nymx), bzn(nxmx, nymx)

      real dxbxn, dybxn, dxxbxn, dyybxn, dxybxn
      real dxbyn, dybyn, dxxbyn, dyybyn, dxybyn
      real dxbzn, dybzn, dxxbzn, dyybzn, dxybzn

      real uxx(nxmx, nymx), uxy(nxmx, nymx), uxz(nxmx, nymx), &
           uyx(nxmx, nymx), uyy(nxmx, nymx), uyz(nxmx, nymx), &
           uzx(nxmx, nymx), uzy(nxmx, nymx), uzz(nxmx, nymx)


      complex uxxk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              uxyk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              uxzk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              uyxk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              uyyk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              uyzk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              uzxk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              uzyk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              uzzk(nkdim1 : nkdim2, mkdim1 : mkdim2)

      real dxuxx(nxmx, nymx), dxuxy(nxmx, nymx), dxuxz(nxmx, nymx), &
           dxuyx(nxmx, nymx), dxuyy(nxmx, nymx), dxuyz(nxmx, nymx), &
           dxuzx(nxmx, nymx), dxuzy(nxmx, nymx), dxuzz(nxmx, nymx)

      real dyuxx(nxmx, nymx), dyuxy(nxmx, nymx), dyuxz(nxmx, nymx), &
           dyuyx(nxmx, nymx), dyuyy(nxmx, nymx), dyuyz(nxmx, nymx), &
           dyuzx(nxmx, nymx), dyuzy(nxmx, nymx), dyuzz(nxmx, nymx)


      real dyyuxx(nxmx, nymx), dyyuxy(nxmx, nymx), dyyuxz(nxmx, nymx), &
           dyyuyx(nxmx, nymx), dyyuyy(nxmx, nymx), dyyuyz(nxmx, nymx), &
           dyyuzx(nxmx, nymx), dyyuzy(nxmx, nymx), dyyuzz(nxmx, nymx)

      real dxyuxx(nxmx, nymx), dxyuxy(nxmx, nymx), dxyuxz(nxmx, nymx), &
           dxyuyx(nxmx, nymx), dxyuyy(nxmx, nymx), dxyuyz(nxmx, nymx), &
           dxyuzx(nxmx, nymx), dxyuzy(nxmx, nymx), dxyuzz(nxmx, nymx)

      real dxxuxx(nxmx, nymx), dxxuxy(nxmx, nymx), dxxuxz(nxmx, nymx), &
           dxxuyx(nxmx, nymx), dxxuyy(nxmx, nymx), dxxuyz(nxmx, nymx), &
           dxxuzx(nxmx, nymx), dxxuzy(nxmx, nymx), dxxuzz(nxmx, nymx)



      complex ex(nxmx, nymx), ey(nxmx, nymx), ez(nxmx, nymx)
      complex dezdx, dezdy, deydx, dexdy
      complex acold(nxmx, nymx), bcold(nxmx, nymx), ccold(nxmx, nymx)
      complex xkperp_cold(nxmx, nymx), xkperp_cold2
      complex xkperp_flux(nrhomax, nthetamax)
      real re_acold, rho_giv, the_giv

      complex bxwave(nxmx, nymx), bywave(nxmx, nymx), bzwave(nxmx, nymx)
      real spx(nxmx, nymx), spy(nxmx, nymx), spz(nxmx, nymx), isq2

      complex eplus_flux_plot(nxmx, nymx)
      complex eminus_flux_plot(nxmx, nymx)
      complex xkperp_flux_plot(nxmx, nymx)

      complex ealpha(nxmx, nymx), ebeta(nxmx, nymx), eb(nxmx, nymx)
      complex eplus(nxmx, nymx), eminus(nxmx, nymx)
      complex ealphakx(nxmx, nymx), ebetakx(nxmx, nymx)
      complex ealphaky(nxmx, nymx), ebetaky(nxmx, nymx)
      complex ebkx(nxmx, nymx), ebky(nxmx, nymx)

      complex xjpxe(nxmx, nymx), xjpye(nxmx, nymx), xjpze(nxmx, nymx)
      complex xjpx1(nxmx, nymx), xjpy1(nxmx, nymx), xjpz1(nxmx, nymx)
      complex xjpx2(nxmx, nymx), xjpy2(nxmx, nymx), xjpz2(nxmx, nymx)
      complex xjpx3(nxmx, nymx), xjpy3(nxmx, nymx), xjpz3(nxmx, nymx)
      complex xjpx4(nxmx, nymx), xjpy4(nxmx, nymx), xjpz4(nxmx, nymx)
      complex xjpx5(nxmx, nymx), xjpy5(nxmx, nymx), xjpz5(nxmx, nymx)
      complex xjpx6(nxmx, nymx), xjpy6(nxmx, nymx), xjpz6(nxmx, nymx)

      complex xjpx_ehst(nxmx, nymx), &
              xjpy_ehst(nxmx, nymx), &
              xjpz_ehst(nxmx, nymx)

      complex xj_slox(nxmx, nymx), &
              xj_sloy(nxmx, nymx), &
              xj_sloz(nxmx, nymx)

      complex xjpx(nxmx, nymx),  xjpy(nxmx, nymx),  xjpz(nxmx, nymx)
      complex xjpx_lab(nxmx, nymx), &
              xjpy_lab(nxmx, nymx), &
              xjpz_lab(nxmx, nymx)

      complex, allocatable :: xjpxe_lab(:,:), &
                                xjpye_lab(:,:), &
                                xjpze_lab(:,:)

      real, allocatable :: ntilda_e_real(:,:)
      complex, allocatable :: ntilda_e(:,:)

      complex pc(nxmx, nymx)

      real wdote(nxmx, nymx), wdoti1(nxmx, nymx), &
           wdoti2(nxmx, nymx), wdoti3(nxmx, nymx), wdot(nxmx, nymx)
      real wdoti4(nxmx, nymx), wdoti5(nxmx, nymx), wdoti6(nxmx, nymx)

      real fye, fyi1, fyi2, fyi3, fy
      real fyi4, fyi5, fyi6


      real fype(nxmx, nymx), fypi1(nxmx, nymx), &
           fypi2(nxmx, nymx), fypi3(nxmx, nymx), fyp(nxmx, nymx)
      real fypi4(nxmx, nymx), fypi5(nxmx, nymx), fypi6(nxmx, nymx)

      real fpol0e (nxmx, nymx), fpol0i1(nxmx, nymx), &
           fpol0i2(nxmx, nymx), fpol0i3(nxmx, nymx)
      real fpol0i4(nxmx, nymx), fpol0i5(nxmx, nymx), fpol0i6(nxmx,nymx)


      real fx0e (nxmx, nymx), fx0i1(nxmx, nymx), &
           fx0i2(nxmx, nymx), fx0i3(nxmx, nymx), &
           fx0i4(nxmx, nymx), fx0i5(nxmx, nymx), fx0i6(nxmx, nymx)

      real fy0e (nxmx, nymx), fy0i1(nxmx, nymx), &
           fy0i2(nxmx, nymx), fy0i3(nxmx, nymx), &
           fy0i4(nxmx, nymx), fy0i5(nxmx, nymx), fy0i6(nxmx, nymx)

      real fz0e (nxmx, nymx), fz0i1(nxmx, nymx), &
           fz0i2(nxmx, nymx), fz0i3(nxmx, nymx), &
           fz0i4(nxmx, nymx), fz0i5(nxmx, nymx), fz0i6(nxmx, nymx)

      real bdotf(nxmx, nymx)


      real capr_fzeta(nxmx, nymx), capr2(nxmx, nymx), &
           uzeta(nxmx, nymx), utheta(nxmx, nymx)
      real gpsi(nxmx, nymx), omgexb(nxmx, nymx)
      real jhat(nxmx, nymx)
      real kpsi(nxmx, nymx), epsi(nxmx, nymx)

      real fprl0e (nxmx, nymx), fprl0i1(nxmx, nymx), &
           fprl0i2(nxmx, nymx), fprl0i3(nxmx, nymx)
      real fprl0i4(nxmx, nymx), fprl0i5(nxmx, nymx), fprl0i6(nxmx, nymx)


      real fx0(nxmx, nymx), fy0(nxmx, nymx), fz0(nxmx, nymx), &
           fpsi0(nxmx, nymx), fpsi1(nxmx, nymx), &
           ftheta0(nxmx, nymx), ftheta1(nxmx, nymx)

      real fpol1e (nxmx, nymx), fpol1i1(nxmx, nymx), &
           fpol1i2(nxmx, nymx), fpol1i3(nxmx, nymx)
      real fpol1i4(nxmx, nymx), fpol1i5(nxmx, nymx), fpol1i6(nxmx, nymx)


      real dfpol1e (nxmx, nymx), dfpol1i1(nxmx, nymx), &
           dfpol1i2(nxmx, nymx), dfpol1i3(nxmx, nymx)
      real dfpol1i4(nxmx,nymx), dfpol1i5(nxmx,nymx), dfpol1i6(nxmx,nymx)


      real denom
      real redotj(nxmx, nymx), pcre(nxmx, nymx), pcim(nxmx, nymx)



      real fx_ant(nxmx, nymx), fy_ant(nxmx, nymx), fz_ant(nxmx, nymx)
      real fx_pla(nxmx, nymx), fy_pla(nxmx, nymx), fz_pla(nxmx, nymx)
      real fxtot_ant, fytot_ant, fztot_ant

      real rekxdotj(nxmx, nymx), rekydotj(nxmx, nymx), &
           rekzdotj(nxmx, nymx)
      real fxtot_plasma, fytot_plasma, fztot_plasma



      real xjprl(nxmx, nymx), xjtot, fyptot, fz0tot
      real redotj1(nxmx, nymx), redotj2(nxmx, nymx), redotj3(nxmx, nymx)
      real redotj4(nxmx, nymx), redotj5(nxmx, nymx), redotj6(nxmx, nymx)

      real redotje(nxmx, nymx), redotjt(nxmx, nymx), redotji(nxmx, nymx)
      real redotjs(nxmx, nymx)
      real redotj_ehst(nxmx, nymx)

      real pcedotj1, pcedotje, pcedotj2, pcedotjt, pcedotj3, pcedotjs
      real pcedotj4, pcedotj5, pcedotj6

      real pedotj1, pedotje, pedotj2, pedotjt, pedotj3, pedotjs
      real pedotj4, pedotj5, pedotj6

      real p, pi1, pi2, pit, pi3, pe, pt, pi4, pi5, pi6

      real p_ql, pi1_ql, pi2_ql, pit_ql, pi3_ql, pe_ql, &
                  pt_ql, pi4_ql, pi5_ql, pi6_ql

      real pcti1, pcte, pcti2, pctt, pcti3, pcti4, pcti5, pcti6

      real pcti1_ql, pcte_ql, pcti2_ql, pctt_ql, pcti3_ql, pcti4_ql, &
           pcti5_ql, pcti6_ql

      complex sigexx, sigexy, sigexz, &
              sigeyx, sigeyy, sigeyz, &
              sigezx, sigezy, sigezz

      complex sig1xx, sig1xy, sig1xz, &
              sig1yx, sig1yy, sig1yz, &
              sig1zx, sig1zy, sig1zz

      complex sig2xx, sig2xy, sig2xz, &
              sig2yx, sig2yy, sig2yz, &
              sig2zx, sig2zy, sig2zz

      complex sig3xx, sig3xy, sig3xz, &
              sig3yx, sig3yy, sig3yz, &
              sig3zx, sig3zy, sig3zz

      complex sig4xx, sig4xy, sig4xz, &
              sig4yx, sig4yy, sig4yz, &
              sig4zx, sig4zy, sig4zz

      complex sig5xx, sig5xy, sig5xz, &
              sig5yx, sig5yy, sig5yz, &
              sig5zx, sig5zy, sig5zz

      complex sig6xx, sig6xy, sig6xz, &
              sig6yx, sig6yy, sig6yz, &
              sig6zx, sig6zy, sig6zz




      complex sigsloxx, sigsloxy, sigsloxz, &
              sigsloyx, sigsloyy, sigsloyz, &
              sigslozx, sigslozy, sigslozz


      complex &
           sigxx, sigxy, sigxz, &
           sigyx, sigyy, sigyz, &
           sigzx, sigzy, sigzz

      complex &
           xkxx, xkxy, xkxz, &
           xkyx, xkyy, xkyz, &
           xkzx, xkzy, xkzz
      complex scap, capk_perp, capk_x

      complex xk1xx, xk1xy, xk1xz, &
              xk1yx, xk1yy, xk1yz, &
              xk1zx, xk1zy, xk1zz

      complex xk2xx, xk2xy, xk2xz, &
              xk2yx, xk2yy, xk2yz, &
              xk2zx, xk2zy, xk2zz


      complex &
           dxx, dxy, dxz, &
           dyx, dyy, dyz, &
           dzx, dzy, dzz

      complex &
           bxx, bxy, bxz, &
           byx, byy, byz, &
           bzx, bzy, bzz


      real yprimec(nymx), ycourse(nymx), &
           yprime(nymx), y(nymx), dy, dyc



! ------------------------------
! storage for parallel scalapack
! ------------------------------
      INTEGER            DLEN_
      PARAMETER          ( DLEN_ = 9 )
      INTEGER            CTXT_, M_, N_, MB_, NB_
      PARAMETER          ( CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6 )
      INTEGER            RSRC_, CSRC_, LLD_
      PARAMETER          ( RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )

      integer            p_amat_dim,p_brhs_dim,p_ipiv_dim
      integer            p_amat_size



      logical use_cur_mod
      parameter(use_cur_mod=.true.)

      complex,        allocatable::p_amat(:)


      integer, allocatable ::         p_ipiv(:)
      complex, allocatable ::         p_brhs(:)

      real err_max,err,b_nrm2,err_nrm2
      complex, dimension(ndfmax) ::brhs2,brhsk,brhs_tmp
      integer icnc
      logical ismine

      integer numroc
      external numroc
      integer lld,nrow_local,ncol_local

      character*4 suffix

      integer mb,nb,myid,nproc,wantnproc,myrow,mycol
      integer icontxt, wantnprocs
      integer lrindx,lcindx,rsrc,csrc,ipos,ii,jj
      integer desc_amat(dlen_), desc_brhs(dlen_)

      integer rsrc1, csrc1, irnc1, icnc1

!efd-begin
      integer idebug
      parameter(idebug=0)
      integer ia,ja, iastart,iaend, jastart,jaend
      integer ia_first,ja_first,   iproc
      integer, allocatable, dimension(:) :: &
               itable,jtable,ntable,mtable
      integer lrindx_base,lcindx_base
      integer lrindx1,lcindx1
      logical isok, found
      integer undefined
      parameter(undefined=-987654321)

!      ---------------------------------------------------
!      local variables for transformation to real space
!      ---------------------------------------------------
      integer :: ib,jb
      integer :: Locp, Locq, irow, icol, nia,nja
      integer :: niastart, niaend
      integer :: num_org, num_new, niasize
      integer :: org_nrow, org_ncol, mmb, nnb
      integer :: new_nrow, new_ncol
      integer :: num_inside, mm, ni, nn, iia, jja, niu, miu

      logical, dimension(:,:), allocatable :: is_inside
      integer, dimension(:), allocatable :: new_to_org
      integer, dimension(:,:), allocatable :: niabegin_all
      integer, dimension(:,:), allocatable :: isize_all
      integer, dimension(:,:,:), allocatable :: descBtmp_all
      integer, dimension(:,:,:), allocatable :: descbrhs_all
      integer, dimension(DLEN_) :: descBtmp, descbrhs

      real :: rho_thresh

        complex xx_inv(nkdim1 : nkdim2, 1 : nxmx), &
              yy_inv(mkdim1 : mkdim2, 1 : nymx)

      complex :: aij, beta
      complex, dimension(:,:), allocatable :: brhs

      complex, dimension(:,:), allocatable :: Btmp

        complex, dimension(ndfmax) :: row,rowk

      complex, dimension(:), allocatable :: work
      logical, parameter :: use_fft = .true.


!       ---------------------------------------------------------
!       further permutation to avoid load imbalance due to
!       assigning cells in resonance region to the same processor
!       ---------------------------------------------------------
        logical, parameter :: need_shuffle = .true.
        integer :: ip, nrow3,ilen
        real*8, allocatable, dimension(:) :: dtemp2
        integer, allocatable, dimension(:) :: itemp2, iperm


#ifdef  USE_HPL
        integer :: hpl_lld, hpl_ineed 
#endif




      integer indxg2p,indxg2l,indxl2g
      external indxg2p,indxg2l,indxl2g

!DLG:   Define variables
      character(len=100) :: eqdsk_fileName
      character(len=100) :: ncFileName
      integer :: nc_id, nuper_id, nupar_id, scalar_id, &
       nrho_id, f_vvp_id, uper_id, upar_id, rho_id, &
       vc_mks_id, pScale_id

    logical :: sane
         integer :: nR_id, nz_id, &
             ePlus_id, eMinu_id, &
             ePlus_img_id, eMinu_img_id, &
             kPer_cold_id, kPer_img_cold_id, &
             R_id, z_id



      allocate (xjpxe_lab(nxmx, nymx), &
                xjpye_lab(nxmx, nymx), &
                xjpze_lab(nxmx, nymx) )
      allocate ( ntilda_e_real ( nxmx, nymx ), &
                 ntilda_e ( nxmx, nymx ) )

!     --------------------------
!     setup parallel environment
!     --------------------------

!     -----------------------------------------
!     some environment may not require mpi_init
!     -----------------------------------------
      call mpi_init( info )



      open(unit=63, file='aorsa2d.in', &
                status='old', form='formatted')
      rewind(63)

      if (myid.eq.0) then
       open(unit=963, file='namelist_input', &
                                     status='unknown', form='formatted')
      end if


!     ------------------------
!     Read namelist input data
!     ------------------------
      read (63, aorsa2din)

!     ------------------------
!     Write namelist input data
!     ------------------------
      write (963, aorsa2din)

      if(abs(signbz).lt.1.0e-05)signbz = +1.0

!     ------------------------------------------- 
!     Set the remainder of the nphi_array to zero
!     -------------------------------------------  
      do i = nphi_number + 1, nphimx
            nphi_array(i) = 0.0
      enddo    

!     ------------------------------------------------------
!     if (nphi_number .gt. 1), set prfin=0.0 so that pscale=1.0
!     ------------------------------------------------------

      if(nphi_number .gt. 1) prfin=0.0

        if (myid.eq.0) then
            open(unit=15, file='out15', status='unknown', form='formatted')
            open(unit=166, file='bharvey_3d', status='unknown', form='formatted')
            if (.not. nPhi_sum_only) &
            open(unit=34, file='fpm',   status='unknown', form='formatted')
            open(unit = 167, file = 'dlg_profiling', status = 'unknown', form = 'formatted' )
        end if

!     -----------------------
!     write fpm (34) file
!     ------------------------
      if (myid.eq.0 .and. (.not. nPhi_sum_only)) then
         write(34,309) nphi1, nphi2
      end if


      idiag = nmodesx / 2
      jdiag = nmodesy / 2

      nbessj = lmax + 2



      call blacs_pinfo( myid, nproc )
      if (nproc .lt. 1) then
         write(6,*) '** blacs_pinfo returns:myid,nproc',myid,nproc
!         write(15,*) '** blacs_pinfo returns:myid,nproc',myid,nproc
       stop '** blacs not setup '
      endif
      call blacs_get(-1,0,icontxt)




!     --------------------
!     setup processor grid
!     --------------------
!     -------------------------------------------------------
!     1 x n processor grid would simplify pivoting
!     n x 1 processor grid would simplify eqnx,eqny,eqnz
!     sqrt(n) x sqrt(n) normally would give best performance
!     -------------------------------------------------------


      call blacs_gridinit( icontxt, 'column-major', nprow,npcol)
      call blacs_gridinfo( icontxt, nprow, npcol, myrow, mycol)


#ifdef  USE_HPL
!     --------------
!     initialize hpl
!     --------------
      call blacs_barrier( icontxt,'A')
      call HPL_blacsinit( icontxt )
#endif

      if (nphi_sum_only) go to 9001

!     ----------------
!     loop over nphi's
!     ----------------
      nt = 0

!     -------------------------------------------------------------
!     Leave all modes out of the sum except the pre-determined ones
!     -------------------------------------------------------------
      do 9000 nsum = 1, nphi_number

      nphi  = nphi_array(nsum)
      nt = nt + 1

      time0=second1(dummy)

!     ----------------------------------------------------------
!     simplifying assumption that all processors can open
!     same files and perform read
!     no need to re-broadcast input data
!     ----------------------------------------------------------
      if (myid .eq. 0) then
          write(6,*) 'blacs started: nprow,npcol,nproc ', &
                  nprow,npcol,nproc
!          write(15,*) 'blacs started: nprow,npcol,nproc ',
!     &            nprow,npcol,nproc
      endif

!      write(6, *) 'myid = ', myid


!     --------------------------------------------------
!     may need to open different files like 'out38.001'
!     for processor '001'
!     --------------------------------------------------
      suffix = '.000'
      suffix(4:4) = char(ichar('0') + mod( myid,    10))
      suffix(3:3) = char(ichar('0') + mod( myid/10, 10))
      suffix(2:2) = char(ichar('0') + mod( myid/100,10))

      dlgTime = second1(dummy)

      if (myid.eq.0) then
         open(unit=38,file='out38',status='unknown',form='formatted')
!         open(unit=29,file='out29',status='unknown',form='formatted')
!         open(unit=30,file='out30',status='unknown',form='formatted')

         open(unit=59,file='movie_eb',status='unknown',form='formatted')
         open(unit=60,file='movie_ealpha',status='unknown', &
            form='formatted')
         open(unit=51,file='rho',status='unknown',form='formatted')
         open(unit=69,file='movie_wdot',status='unknown', &
              form='formatted')


!         open(unit=32,file='diagnostic',status='unknown',
!     .      form='formatted')
         open(unit=53, file='out_fields', status='unknown', &
            form='formatted')
!         open(unit=54, file='fields_fourier', status='unknown',
!     .      form='formatted')
         open(unit=28, file='fields_local',   status='unknown', &
            form='formatted')
!DLG:   Adjust if statement for ndist?? >= 1
         if (ndiste .eq. 1)  open(unit=40, file='out_cql3d.coefe', &
                                  status='unknown', form='formatted')
         if (ndisti1 .ge. 1) open(unit=41, file='out_cql3d.coef1', &
                                  status='unknown', form='formatted')
         if (ndisti2 .ge. 1) open(unit=42, file='out_cql3d.coef2', &
                                  status='unknown', form='formatted')
         if (ndisti3 .ge. 1) open(unit=43, file='out_cql3d.coef3', &
                                  status='unknown', form='formatted')
         if (ndisti4 .ge. 1) open(unit=44, file='out_cql3d.coef4', &
                                  status='unknown', form='formatted')


      endif

      if (myId==0) write(167,12398) 'open various files:', &
        (dlgTime - second1(dummy))/60.0
12398 format(a30,2x,f6.2)
 
      nrhs = 1
      nmaxe=1



      nkx2 = nmodesx / 2
      nkx1 = - nmodesx / 2 + 1
      nnodex = nmodesx
      nnodey = nmodesy

      nnoderho = nnodex / 2
      mnodetheta = nnodey / 2

      nky2 = nmodesy / 2
      nky1 = - nmodesy / 2 + 1


!      jequat  = nnodey / 2
      icenter = nnodex / 2


      if (qavg0 .ne. 0.0) xiota0 = 1./qavg0

      rholim = sqrt(psilim)

      q = 1.6e-19
      if(te0 .eq. 0.0)te0 = ti0

      t0e = te0
      t0i = ti0
      t0i2 = ti02
      t0i3 = ti03
      t0i4 = ti0
      t0i5 = ti02
      t0i6 = ti03

      teedge = 400.0

      qhat = qavg0

      eslow = eslowev * q
      t0e = t0e   * q
      t0i = t0i   * q
      t0i2 = t0i2 * q
      t0i3 = t0i3 * q
      t0i4 = t0i4 * q
      t0i5 = t0i5 * q
      t0i6 = t0i6 * q

      teedge = teedge * q

      xme = 9.11e-31
      xmh = 1.67e-27
      xmi1 = amu1 * xmh
      xmi2 = amu2 * xmh
      xmi3 = amu3 * xmh
      xmi4 = amu4 * xmh
      xmi5 = amu5 * xmh
      xmi6 = amu6 * xmh
      xmi_slo = amu_slo * xmh

      qi1 = z1 * q
      qi2 = z2 * q
      qi3 = z3 * q
      qi4 = z4 * q
      qi5 = z5 * q
      qi6 = z6 * q
      qi_slo = z_slo * q

      qe = -q
      zi = cmplx(0.0,1.0)
      eps0 = 8.85e-12

      pi = 3.141592654
      twopi = 2.0 * pi

      xlnlam = 20.0
      xmu0 = 1.26e-06
      clight = 1.0 / sqrt(eps0 * xmu0)


      xkphi0 = nphi / rt
      if(xkphi0 .eq. 0.0)xkphi0 = 1.0e-05

      xnurf = freqcy

      omgrf = 2.0 * pi * xnurf

      vphase = omgrf / xkphi0
      vthe = sqrt(t0e / xme)
      vsound = sqrt(teedge / xmi1)
      wphase = 0.0
      if(vthe .ne. 0.0)wphase = vphase / vthe


      xkthrho = 0.2

      xkthdx = 1.0



      xk0 = omgrf / clight
      rnz = xkphi0 / xk0


      if (myid .eq. 0) then

         write(15, *) "AORSA2D: version_number = ", version_number
             write(15, *) "nphi_number = ", nphi_number
         write(15, *) "nphi_array = ", nphi_array


         write(6, *) "AORSA2D: version_number = ", version_number
         write(6, *) "nphi_number = ", nphi_number
         write(6, *) "nphi_array = ", nphi_array


         write (6, *)
         write (6, *) "xkperp_cutoff = ", xkperp_cutoff
         write (6, *) "damping       = ", damping
       write (6, *) "upshift = ", upshift
       write (6, *)
       write (6, *) "n_prof_flux = ", n_prof_flux
         write (6, *) "eqdsk       = ", eqdsk
       write (6, *) "netcdf_file1 = ", netcdf_file1
       write (6, *) "netcdf_file2 = ", netcdf_file2
       write (6, *)
         write (6, *) "nnodex  = ", nnodex
         write (6, *) "nnodey  = ", nnodey
         write (6, *)

       write (6, *)
         write (6, *) "nprow   = ", nprow
         write (6, *) "npcol   = ", npcol
       write (6, *)
       write (6, *) "i_antenna = ", i_antenna
       write (6, *) "prfin   = ", prfin
         write (6, *) "freqcy  = ", freqcy
         write (6, *) "nphi    = ", nphi
         write (6, *) "antlen  = ", antlen
         write (6, *) "antlc   = ", antlc
         write (6, *) "rant    = ", rant
         write (6, *) "yant    = ", yant
       write (6, *)
       write (6, *) "delta0  = ", delta0
         write (6, *) "xnuomg  = ", xnuomg
       write (6, *)
         write (6, *) "lmax    = ", lmax
         write (6, *) "zeffcd  = ", zeffcd
         write (6, *)
         write (6, *) "xn0     = ", xn0
         write (6, *) "xnlim   = ", xnlim
         write (6, *) "alphan  = ", alphan
         write (6, *) "betan   = ", betan
         write (6, *) "te0     = ", te0
       write (6, *) "telim   = ", telim
         write (6, *) "alphate = ", alphate
         write (6, *) "betate  = ", betate
       write (6, *)
         write (6, *) "ndisti1 = ", ndisti1
         write (6, *) "amu1    = ", amu1
         write (6, *) "z1      = ", z1
         write (6, *) "ti0     = ", ti0
         write (6, *) "tilim   = ", tilim
         write (6, *) "alphati = ", alphati
         write (6, *) "betati  = ", betati
       write (6, *)
         write (6, *) "ndisti2 = ", ndisti2
         write (6, *) "amu2    = ", amu2
         write (6, *) "z2      = ", z2
         write (6, *) "xn2     = ", xn2
         write (6, *) "xn2lim  = ", xn2lim
         write (6, *) "alphan2 = ", alphan2
         write (6, *) "betan2  = ", betan2
         write (6, *) "ti02    = ", ti02
         write (6, *) "ti2lim  = ", ti2lim
         write (6, *) "alphati2= ", alphati2
         write (6, *) "betati2 = ", betati2
         write (6, *)
         write (6, *) "ndisti3 = ", ndisti3
         write (6, *) "amu3    = ", amu3
         write (6, *) "z3      = ", z3
         write (6, *) "xn3     = ", xn3
         write (6, *) "xn3lim  = ", xn3lim
         write (6, *) "alphan3 = ", alphan3
         write (6, *) "betan3  = ", betan3
         write (6, *) "ti03    = ", ti03
         write (6, *) "ti3lim  = ", ti3lim
         write (6, *) "alphati3= ", alphati3
         write (6, *) "betati3 = ", betati3
         write (6, *)
         write (6, *) "ndisti4 = ", ndisti4
         write (6, *) "amu4    = ", amu4
         write (6, *) "z4      = ", z4
         write (6, *) "xn4     = ", xn4
         write (6, *) "xn4lim  = ", xn4lim
         write (6, *) "alphan4 = ", alphan4
         write (6, *) "betan4  = ", betan4
         write (6, *) "ti04    = ", ti04
         write (6, *) "ti4lim  = ", ti4lim
         write (6, *) "alphati4= ", alphati4
         write (6, *) "betati4 = ", betati4
         write (6, *)
         write (6, *) "ndisti5 = ", ndisti5
         write (6, *) "amu5    = ", amu5
         write (6, *) "z5      = ", z5
         write (6, *) "xn5     = ", xn5
         write (6, *) "xn5lim  = ", xn5lim
         write (6, *) "alphan5 = ", alphan5
         write (6, *) "betan5  = ", betan5
         write (6, *) "ti05    = ", ti05
         write (6, *) "ti5lim  = ", ti5lim
         write (6, *) "alphati5= ", alphati5
         write (6, *) "betati5 = ", betati5
         write (6, *)


         write (15, *)
         write (15, *) "xkperp_cutoff = ", xkperp_cutoff
         write (15, *) "damping       = ", damping
       write (15, *) "upshift = ", upshift
       write (15, *)
       write (15, *) "n_prof_flux = ", n_prof_flux
         write (15, *) "eqdsk       = ", eqdsk
       write (15, *) "netcdf_file1 = ", netcdf_file1
       write (15, *) "netcdf_file2 = ", netcdf_file2
       write (15, *)
         write (15, *) "nnodex  = ", nnodex
         write (15, *) "nnodey  = ", nnodey

       write (15, *)
       write (15, *) "nprow   = ", nprow
         write (15, *) "npcol   = ", npcol
       write (15, *)
         write (15, *) "i_antenna = ", i_antenna
       write (15, *) "prfin     = ", prfin
         write (15, *) "freqcy    = ", freqcy
         write (15, *) "nphi      = ", nphi
         write (15, *) "antlen    = ", antlen
         write (15, *) "rant      = ", rant
         write (15, *) "yant      = ", yant
       write (15, *)
       write (15, *) "delta0  = ", delta0
         write (15, *) "xnuomg  = ", xnuomg
       write (15, *)
         write (15, *) "lmax    = ", lmax
         write (15, *) "zeffcd  = ", zeffcd
         write (15, *)
         write (15, *) "xn0     = ", xn0
         write (15, *) "xnlim   = ", xnlim
         write (15, *) "alphan  = ", alphan
         write (15, *) "betan   = ", betan
         write (15, *) "te0     = ", te0
       write (15, *) "telim   = ", telim
         write (15, *) "alphate = ", alphate
         write (15, *) "betate  = ", betate
       write (15, *)
         write (15, *) "ndisti1 = ", ndisti1
         write (15, *) "amu1    = ", amu1
         write (15, *) "z1      = ", z1
         write (15, *) "ti0     = ", ti0
         write (15, *) "tilim   = ", tilim
         write (15, *) "alphati = ", alphati
         write (15, *) "betati  = ", betati
       write (15, *)
         write (15, *) "ndisti2 = ", ndisti2
         write (15, *) "amu2    = ", amu2
         write (15, *) "z2      = ", z2
         write (15, *) "xn2     = ", xn2
         write (15, *) "xn2lim  = ", xn2lim
         write (15, *) "alphan2 = ", alphan2
         write (15, *) "betan2  = ", betan2
         write (15, *) "ti02    = ", ti02
         write (15, *) "ti2lim  = ", ti2lim
         write (15, *) "alphati2= ", alphati2
         write (15, *) "betati2 = ", betati2
         write (15, *)
         write (15, *) "ndisti3 = ", ndisti3
         write (15, *) "amu3    = ", amu3
         write (15, *) "z3      = ", z3
         write (15, *) "xn3     = ", xn3
         write (15, *) "xn3lim  = ", xn3lim
         write (15, *) "alphan3 = ", alphan3
         write (15, *) "betan3  = ", betan3
         write (15, *) "ti03    = ", ti03
         write (15, *) "ti3lim  = ", ti3lim
         write (15, *) "alphati3= ", alphati3
         write (15, *) "betati3 = ", betati3
         write (15, *)
         write (15, *) "ndisti4 = ", ndisti4
         write (15, *) "amu4    = ", amu4
         write (15, *) "z4      = ", z4
         write (15, *) "xn4     = ", xn4
         write (15, *) "xn4lim  = ", xn4lim
         write (15, *) "alphan4 = ", alphan4
         write (15, *) "betan4  = ", betan4
         write (15, *) "ti04    = ", ti04
         write (15, *) "ti4lim  = ", ti4lim
         write (15, *) "alphati4= ", alphati4
         write (15, *) "betati4 = ", betati4
         write (15, *)
         write (15, *) "ndisti5 = ", ndisti5
         write (15, *) "amu5    = ", amu5
         write (15, *) "z5      = ", z5
         write (15, *) "xn5     = ", xn5
         write (15, *) "xn5lim  = ", xn5lim
         write (15, *) "alphan5 = ", alphan5
         write (15, *) "betan5  = ", betan5
         write (15, *) "ti05    = ", ti05
         write (15, *) "ti5lim  = ", ti5lim
         write (15, *) "alphati5= ", alphati5
         write (15, *) "betati5 = ", betati5
         write (15, *)



!         write(6, 7013)nmodesx, nmodesy
!         write(6, 7113)nwdot
!         write(6, 7213)nnodecx, nnodecy

!       write(6, *)"nprow = ", nprow
!        write(6, *)"npcol = ", npcol
!       write(6, *)"n_bin = ", n_bin

!         write(6, *) "i_write = ", i_write
!       write(15, *)"i_write = ", i_write


         write(6, 1815)te0
         write(6, 1821)ti0

         write(6, 7217)qavg0
         write(6, 7014)lmax
         write(6, 7015)ibessel
         write(6, 7115)nzfun
         write(6, 71160)xnuomg
!         write(6, 7016)rhoi1overl
!         write(6, 7017)rhoi10
         write(6, 7217)qavg0
         write(6, 1013)nphi
         write(6, 1321)wphase
         write(6, 1323)vthe
         write(6, 1021)rnz
         write(6, 1812)rt
         write(6, 1822)aplasm
         write(6, 1823)rant
         write(6, 1809)b0
         write(6, 6813)xn0
         write(6, 1813)xn1
         write(6, 1814)xn2
         write(6, 1834)xn3

         write(6, 1012)omgrf
         write(6, 1009)xnurf
         write(6, 1013)nphi
         write(6, 1321)wphase
         write(6, 1322)vphase
         write(6, 1323)vthe
         write(6, 1714)xk0
         write(6, 1016)xnuead
         write(6, 1017)xnu1ad
         write(6, 1018)xnu2ad
         write(6, 1020)nnodex, nnodey

       write(6, *) "iprofile = ", iprofile
         write(6, *) "psilim = ", psilim


         write(15, 162)
         write(15, 7013)nmodesx, nmodesy
         write(15, 7113)nwdot
         write(15, 7213)nnodecx, nnodecy

!         write(15,*)"nprow = ", nprow
!         write(15,*)"npcol = ", npcol

       write(15, 1815)t0e
         write(15, 1821)t0i

         write(15, 7217)qavg0
         write(15, 7014)lmax
         write(15, 7015)ibessel
         write(15, 7115)nzfun
         write(15, 71160)xnuomg
!         write(15, 7016)rhoi1overl
!         write(15, 7017)rhoi10
         write(15, 7217)qavg0
         write(15, 1013)nphi
         write(15, 1321)wphase
         write(15, 1323)vthe
         write(15, 1021)rnz
         write(15, 1812)rt
         write(15, 1822)aplasm
         write(15, 1823)rant
         write(15, 1809)b0
         write(15, 6813)xn0
         write(15, 1813)xn1
         write(15, 1814)xn2
         write(15, 1834)xn3

         write(15, 1012)omgrf
         write(15, 1009)xnurf
         write(15, 1013)nphi
         write(15, 1321)wphase
         write(15, 1322)vphase
         write(15, 1323)vthe
         write(15, 1714)xk0
         write(15, 1016)xnuead
         write(15, 1017)xnu1ad
         write(15, 1018)xnu2ad
         write(15, 1020)nnodex, nnodey

       write(15, *) "iprofile = ", iprofile
         write(15, *) "psilim = ", psilim
      endif

      telimj   = telim  * q
      tilimj   = tilim  * q
      ti2limj  = ti2lim * q
      ti3limj  = ti3lim * q
      ti4limj  = ti4lim * q
      ti5limj  = ti5lim * q
      ti6limj  = ti6lim * q


!------------------------------------
!     Store factorials from 0 to lmax
!------------------------------------
      do l = 0, lmax
         factl_sav(l) = factrl(l)
      end do


      t1 = second1(dummy)

!------------------------------------------
!     Flux surfaces from EQDSK (ga version)
!------------------------------------------

      if(igeom .eq. 5) then

      dlgTime = second1(dummy)

       call eqdsk_setup(myid, eqdsk, nmodesx, nmodesy, &
                rwleft, rwright, ytop, ybottom, &
            rmaxis, zmaxis, b0, psio, psimag, psi_tor_max, &
            bxn, byn, bzn, bmod, psi_pol2d, rho_pol2d, qsafety, &
            bmod_mid,  capr_bpol_mid2,  capr_bpol_mid, rho_tor2d, &
            i_psi_eq, dldb_tot12, dldbavg, n_prof_flux,dlg_yRange)

      if (myId==0) write(167,12398) 'eqdsk_setup:', &
        (dlgTime - second1(dummy))/60.0

       call blacs_barrier(icontxt, 'All')


 9318    format(a128)
       rt = rmaxis

       r0 = rmaxis
       z0 = zmaxis

       if (myid .eq. 0)then
          write (6, *) "rwleft  = ", rwleft
            write (6, *) "rwright = ", rwright
            write (6, *) "ybottom = ", ybottom
          write (6, *) "ytop    = ", ytop

          write (15, *)
            write (15, *) "rwleft  = ", rwleft
            write (15, *) "rwright = ", rwright
            write (15, *) "ybottom = ", ybottom
          write (15, *) "ytop    = ", ytop

          write (6, *) "eqdsk       = ", eqdsk
          write (15, *) "eqdsk       = ", eqdsk
       end if






!        ----------------------------------------------------
!        Default:  if n_prof_flux equals 0, use poloidal flux
!        ----------------------------------------------------

            do i = 1, nnodex
               do j = 1, nnodey
                rho(i,j) = rho_pol2d(i,j)
                psi(i,j) = psi_pol2d(i,j)
                  psi_dim(i,j) = psi(i,j) * psio
               end do
            end do


!        --------------------------------------------------------------
!        Alternate: if n_prof_flux does not equal 0,  use toroidal flux
!        Note:  Here we use toroidal flux  inside rho_pol = 1.0
!                       and poloidal flux outside rho_pol = 1.0!!
!        --------------------------------------------------------------
         if(n_prof_flux .ne. 0)then
            do i = 1, nnodex
               do j = 1, nnodey
                psi_tor2d(i,j) = rho_tor2d(i,j)**2

              if(rho(i,j) .lt. 1.0)then
                   rho(i,j) = rho_tor2d(i,j)
                   psi(i,j) = psi_tor2d(i,j)
                   psi_dim(i,j) = psi(i,j) * psi_tor_max
                     end if

               end do
            end do
         end if


         if(myid .eq. 0)then
            if(n_prof_flux .eq. 0)write(6, *)'use sqrt(poloidal flux)'
          if(n_prof_flux .ne. 0)write(6, *)'use sqrt(toroidal flux)'
            write(6, 1812)rt
            write(6, 1809)b0
          write(6, *) "rmaxis = ", rmaxis
          write(6, *) "zmaxis = ", zmaxis

            if(n_prof_flux .eq. 0)write(15,*)'use sqrt(poloidal flux)'
          if(n_prof_flux .ne. 0)write(15,*)'use sqrt(toroidal flux)'
            write(15, 1812)rt
            write(15, 1809)b0
         end if


       omgci10 = qi1 * b0 / xmi1
         vthi10 = sqrt(2.0 * t0i / xmi1)
         rhoi10 = vthi10 / omgci10
       if(myid .eq. 0)write(15, *)"rhoi10 = ", rhoi10

       omgci20 = qi2 * b0 / xmi2
         vthi20 = sqrt(2.0 * t0i2 / xmi2)
         rhoi20 = vthi20 / omgci20
       if(myid .eq. 0)write(15, *)"rhoi20 = ", rhoi20

       omgci30 = qi3 * b0 / xmi3
         vthi30 = sqrt(2.0 * t0i3 / xmi3)
         rhoi30 = vthi30 / omgci30
       if(myid .eq. 0)write(15, *)"rhoi30 = ", rhoi30

       omgci40 = qi4 * b0 / xmi4
         vthi40 = sqrt(2.0 * t0i4 / xmi4)
         rhoi40 = vthi40 / omgci40
       if(myid .eq. 0)write(15, *)"rhoi40 = ", rhoi40

       omgci50 = qi5 * b0 / xmi5
         vthi50 = sqrt(2.0 * t0i5 / xmi5)
         rhoi50 = vthi50 / omgci50
       if(myid .eq. 0)write(15, *)"rhoi50 = ", rhoi50


         do i = 1, nnodex
            do j = 1, nnodey
               btau(i,j) = sqrt(bxn(i,j)**2 + byn(i,j)**2)
             bpol(i,j) = btau(i,j) * bmod(i,j)
               bzeta(i,j) = bzn(i,j)


!              ---------------------------
!              Calculate rotation matrix U
!              ---------------------------

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

            end do
         end do

      end if


      xwleft = rwleft - rt
      xwright = rwright - rt



!     --------------------------------
!     Distribution function from CQL3D
!     --------------------------------

!     -----------------
!     allocate arrays
!     -----------------

      if(ndisti1   .eq. 0 .and. &
         ndisti2   .eq. 0 .and. &
         ndisti3   .eq. 0 .and. &
         ndisti4   .eq. 0 .and. &
         ndisti5   .eq. 0 .and. &
         ndisti6   .eq. 0 .and. &
         ndiste    .eq. 0)   then

         allocate( UPERP(nuper) )
         allocate( UPARA(nupar) )
       allocate( VPERP(nuper) )
         allocate( VPARA(nupar) )

         allocate( UPERP_work(nuper) )
         allocate( UPARA_work(nupar) )

       allocate( f(nuper, nupar) )

         allocate( dfdupere(nuper, nupar) )
         allocate( dfdupare(nuper, nupar) )

         allocate( dfduper1(nuper, nupar) )
         allocate( dfdupar1(nuper, nupar) )

       allocate( dfduper2(nuper, nupar) )
         allocate( dfdupar2(nuper, nupar) )

       allocate( dfduper3(nuper, nupar) )
         allocate( dfdupar3(nuper, nupar) )

       allocate( dfduper4(nuper, nupar) )
         allocate( dfdupar4(nuper, nupar) )

         allocate( dfduper5(nuper, nupar) )
         allocate( dfdupar5(nuper, nupar) )

       allocate( dfduper6(nuper, nupar) )
         allocate( dfdupar6(nuper, nupar) )




       UPERP = 0.0
       UPARA = 0.0
       UPERP_work = 0.0
       UPARA_work = 0.0

       f = 0.0

       dfdupere = 0.0
       dfdupare = 0.0

       dfduper1 = 0.0
       dfdupar1 = 0.0

       dfduper2 = 0.0
       dfdupar2 = 0.0

       dfduper3 = 0.0
       dfdupar3 = 0.0

       dfduper4 = 0.0
       dfdupar4 = 0.0

       dfduper5 = 0.0
       dfdupar5 = 0.0

       dfduper6 = 0.0
       dfdupar6 = 0.0

      end if

      if(ndisti1   .ne. 0 .or. &
         ndisti2   .ne. 0 .or. &
         ndisti3   .ne. 0 .or. &
         ndisti4   .ne. 0 .or. &
         ndisti5   .ne. 0 .or. &
         ndisti6   .ne. 0 .or. &
         ndiste    .ne. 0)   then


!        -----------------
!        allocate arrays
!        -----------------
         allocate( UPERP(nuper) )
         allocate( UPARA(nupar) )
       allocate( VPERP(nuper) )
         allocate( VPARA(nupar) )
         allocate( UPERP_work(nuper) )
         allocate( UPARA_work(nupar) )


       allocate( f(nuper, nupar) )

         allocate( dfdupere(nuper, nupar) )
         allocate( dfdupare(nuper, nupar) )

         allocate( dfduper1(nuper, nupar) )
         allocate( dfdupar1(nuper, nupar) )

       allocate( dfduper2(nuper, nupar) )
         allocate( dfdupar2(nuper, nupar) )

       allocate( dfduper3(nuper, nupar) )
         allocate( dfdupar3(nuper, nupar) )

       allocate( dfduper4(nuper, nupar) )
         allocate( dfdupar4(nuper, nupar) )

         allocate( dfduper5(nuper, nupar) )
         allocate( dfdupar5(nuper, nupar) )

       allocate( dfduper6(nuper, nupar) )
         allocate( dfdupar6(nuper, nupar) )


       UPERP = 0.0
       UPARA = 0.0
       UPERP_work = 0.0
       UPARA_work = 0.0

       f = 0.0

       dfdupere = 0.0
       dfdupare = 0.0

       dfduper1 = 0.0
       dfdupar1 = 0.0

       dfduper2 = 0.0
       dfdupar2 = 0.0

       dfduper3 = 0.0
       dfdupar3 = 0.0

       dfduper4 = 0.0
       dfdupar4 = 0.0

       dfduper5 = 0.0
       dfdupar5 = 0.0

       dfduper6 = 0.0
       dfdupar6 = 0.0

         allocate( fe_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( dfe_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( dfe_cql_uprl(nuper, nupar, n_psi_dim) )

         allocate( f1_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( df1_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( df1_cql_uprl(nuper, nupar, n_psi_dim) )

         allocate( f2_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( df2_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( df2_cql_uprl(nuper, nupar, n_psi_dim) )

       allocate( f3_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( df3_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( df3_cql_uprl(nuper, nupar, n_psi_dim) )

       allocate( f4_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( df4_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( df4_cql_uprl(nuper, nupar, n_psi_dim) )

       allocate( f5_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( df5_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( df5_cql_uprl(nuper, nupar, n_psi_dim) )

       allocate( f6_cql_cart(nuper, nupar, n_psi_dim) )
         allocate( df6_cql_uprp(nuper, nupar, n_psi_dim) )
         allocate( df6_cql_uprl(nuper, nupar, n_psi_dim) )



       fe_cql_cart = 0.0
       dfe_cql_uprp = 0.0
       dfe_cql_uprl = 0.0

       f1_cql_cart = 0.0
       df1_cql_uprp = 0.0
       df1_cql_uprl = 0.0

       f2_cql_cart = 0.0
       df2_cql_uprp = 0.0
       df2_cql_uprl = 0.0

       f3_cql_cart = 0.0
       df3_cql_uprp = 0.0
       df3_cql_uprl = 0.0

       f4_cql_cart = 0.0
       df4_cql_uprp = 0.0
       df4_cql_uprl = 0.0

       f5_cql_cart = 0.0
       df5_cql_uprp = 0.0
       df5_cql_uprl = 0.0

         f6_cql_cart = 0.0
       df6_cql_uprp = 0.0
       df6_cql_uprl = 0.0

       call blacs_barrier(icontxt, 'All')

          if(myid.eq.0)write(*,*) 'dlg: vce_mks pre cql: ', vce_mks
        if(myid.eq.0)write(*,*) 'dlg: vc1_mks pre cql: ', vc1_mks
        if(myid.eq.0)write(*,*) 'dlg: vc2_mks pre cql: ', vc2_mks


!DLG:   Adjust if statement for ndisti1 >= 1
         if (ndisti1 .ge. 1) then

    sanity_check1: &
    if ( enorm_factor_i1 .lt. 0 .and. ndisti1 .eq. 2) then 
        stop 'FAILED SANITY CHECK: plese set enorm_factor_i1'
    endif sanity_check1


          dlgTime = second1(dummy)

          if(myid .eq. 0) call cql3d_setup(netcdf_file1, nuper, nupar, &
                             xmi1, enorm_factor_i1, vc1_mks_cql3d)

          if (myId==0) write(167,12398) 'cql3d_setup 1:', &
            (dlgTime - second1(dummy))/60.0


            vc1_cgs_cql3d = vc1_mks_cql3d * 100.

            call blacs_barrier(icontxt, 'All')

            open(unit=50, file='cql3d.out', status='unknown', &
                                                       form='formatted')

!           ---------------------------------------------
!              Read data from CQL3D in u_perp and u_parallel
!           ---------------------------------------------

            rewind (50)

            read (50, 309) nuper
            read (50, 309) nupar
            read (50, 309) n_psi

            read (50, 3310) vc1_mks
            read (50, 3310) UminPara, UmaxPara

            read (50, 3310) (rho_a(i_psi), i_psi = 1, n_psi)
            read (50, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            read (50, 3310) (upara(i_upara), i_upara = 1, nupar)

            read (50, 3310) (((f1_cql_cart(i_uperp, i_upara, i_psi), &
             i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

            read (50, 3310) (((df1_cql_uprp(i_uperp, i_upara, i_psi), &
             i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

            read (50, 3310) (((df1_cql_uprl(i_uperp, i_upara, i_psi), &
             i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

             close(50)

       end if

       call blacs_barrier(icontxt, 'All')
!DLG:   Adjust if statement for ndisti2 >= 1

         if (ndisti2 .ge. 1) then

    sanity_check2: &
    if ( enorm_factor_i2 .lt. 0 .and. ndisti2 .eq. 2) then 
        stop 'FAILED SANITY CHECK: plese set enorm_factor_i2'
    endif sanity_check2

            dlgTime = second1(dummy)

          if(myid .eq. 0) call cql3d_setup(netcdf_file2, nuper, nupar, &
                              xmi2, eNorm_factor_i2, vc2_mks_cql3d)
            vc2_cgs_cql3d = vc2_mks_cql3d * 100.

            if (myId==0) write(167,12398) 'cql3d_setup 2:', &
            (dlgTime - second1(dummy))/60.0


            call blacs_barrier(icontxt, 'All')

            open(unit=50,file='cql3d.out', status='unknown', &
                                                       form='formatted')

!           ---------------------------------------------
!              Read data from CQL3D in u_perp and u_parallel
!           ---------------------------------------------

            rewind (50)

            read (50, 309) nuper
            read (50, 309) nupar
            read (50, 309) n_psi

            read (50, 3310) vc2_mks
            read (50, 3310) UminPara, UmaxPara

            read (50, 3310) (rho_a(i_psi), i_psi = 1, n_psi)
            read (50, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            read (50, 3310) (upara(i_upara), i_upara = 1, nupar)

          read (50, 3310) (((f2_cql_cart(i_uperp, i_upara, i_psi), &
             i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

            read (50, 3310) (((df2_cql_uprp(i_uperp, i_upara, i_psi), &
             i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

            read (50, 3310) (((df2_cql_uprl(i_uperp, i_upara, i_psi), &
             i_uperp = 1, nuper), i_upara = 1, nupar), i_psi = 1, n_psi)

            close(50)

        dlgTime = second1(dummy)
        if ( myId .eq. 0 ) then
!DLG:   Dump f(vPer,vPar,rho) for plotting
        write (*,*) 'WRITING output/f_vvp.nc ...'
        ncFileName = 'output/f_vvp.nc'

            call check ( &
       nf90_create ( ncFileName, nf90_clobber, nc_id ) )
            call check ( &
       nf90_def_dim ( nc_id, "nuper", nuper, nuper_id ) )
            call check ( &
       nf90_def_dim ( nc_id, "nupar", nupar, nupar_id ) )
            call check ( &
       nf90_def_dim ( nc_id, "scalar", 1, scalar_id ) )
            call check ( &
       nf90_def_dim ( nc_id, "nrho", n_psi, nrho_id ) )

            call check ( &
       nf90_def_var ( nc_id, "f_vvp", NF90_REAL, &
       (/ nuper_id, nupar_id, nrho_id /), f_vvp_id ) )
            call check ( &
       nf90_def_var ( nc_id, "uper", NF90_REAL, &
       (/ nuper_id /), uper_id ) )
            call check ( &
       nf90_def_var ( nc_id, "upar", NF90_REAL, &
       (/ nupar_id /), upar_id ) )
            call check ( &
       nf90_def_var ( nc_id, "rho", NF90_REAL, &
       (/ nrho_id /), rho_id ) )
            call check ( &
       nf90_def_var ( nc_id, "vc_mks", NF90_REAL, &
       scalar_id, vc_mks_id ) )
            call check ( nf90_enddef ( nc_id ) )

            call check ( &
       nf90_put_var ( nc_id, f_vvp_id, f2_cql_cart(:,:,1:n_psi) ) )
            call check ( &
       nf90_put_var ( nc_id, uper_id, uperp ) )
            call check ( &
       nf90_put_var ( nc_id, upar_id, upara ) )
            call check ( &
       nf90_put_var ( nc_id, rho_id, rho_a(1:n_psi) ) )
            call check ( &
       nf90_put_var ( nc_id, vc_mks_id, vc2_mks ) )

            call check ( nf90_close ( nc_id ) )
            write (*,*) 'DONE'
        endif
       end if

       if (myId==0) write(167,12398) 'write f_vvp.nc:', &
            (dlgTime - second1(dummy))/60.0




         if(myid .eq. 0)then
            WRITE (6,*)
            WRITE (6,*) "nuper = ",  nuper
            WRITE (6,*) "nupar = ",  nupar
            WRITE (6,*) "n_psi = ",  n_psi

            WRITE (6,*)
            WRITE (6,*) "vce_mks = ", vce_mks
          WRITE (6,*) "vc1_mks = ", vc1_mks
            WRITE (6,*) "vc2_mks = ", vc2_mks
            WRITE (6,*) "vc3_mks = ", vc3_mks
          WRITE (6,*) "vc4_mks = ", vc4_mks
            WRITE (6,*) "vc5_mks = ", vc5_mks
            WRITE (6,*) "vc6_mks = ", vc6_mks

            write(6, *)
            write(6, *) "rho/a(i_psi)"
            write(6, 310) (rho_a(i_psi), i_psi = 1, n_psi)

            write(6, *)
            write(6, *) "uperp ="
            write(6, 310) (uperp(i_uperp), i_uperp = 1, nuper)

            write(6, *)
            write(6, *) "upara ="
            write(6, 310) (upara(i_upara), i_upara = 1, nupar)



            WRITE (15,*)
            WRITE (15,*) "nuper = ",  nuper
            WRITE (15,*) "nupar = ",  nupar
            WRITE (15,*) "n_psi = ",  n_psi

            WRITE (15,*)
            WRITE (15,*) "vce_mks = ", vce_mks
          WRITE (15,*) "vc1_mks = ", vc1_mks
            WRITE (15,*) "vc2_mks = ", vc2_mks
            WRITE (15,*) "vc3_mks = ", vc3_mks
          WRITE (15,*) "vc4_mks = ", vc4_mks
            WRITE (15,*) "vc5_mks = ", vc5_mks
            WRITE (15,*) "vc6_mks = ", vc6_mks

            write(15, *)
            write(15, *) "rho/a(i_psi)"
            write(15, 310) (rho_a(i_psi), i_psi = 1, n_psi)

            write(15, *)
            write(15, *) "uperp ="
            write(15, 310) (uperp(i_uperp), i_uperp = 1, nuper)

            write(15, *)
            write(15, *) "upara ="
            write(15, 310) (upara(i_upara), i_upara = 1, nupar)
         end if

      end if

        if(myid.eq.0)write(*,*) 'dlg: vce_mks post cql: ', vce_mks
        if(myid.eq.0)write(*,*) 'dlg: vc1_mks post cql: ', vc1_mks
        if(myid.eq.0)write(*,*) 'dlg: vc2_mks post cql: ', vc2_mks

        if( eNorm_factor_e>0) &
        vce_mks = sqrt(2d0*1d3*eNorm_factor_e*e/xme)
        if( eNorm_factor_i1>0) &
        vc1_mks = sqrt(2d0*1d3*eNorm_factor_i1*e/xmi1)
        if( eNorm_factor_i2>0) &
        vc2_mks = sqrt(2d0*1d3*eNorm_factor_i2*e/xmi2)
        if( eNorm_factor_i3>0) &
        vc3_mks = sqrt(2d0*1d3*eNorm_factor_i3*e/xmi3)
        if( eNorm_factor_i4>0) &
        vc4_mks = sqrt(2d0*1d3*eNorm_factor_i4*e/xmi4)
        if( eNorm_factor_i5>0) &
        vc5_mks = sqrt(2d0*1d3*eNorm_factor_i5*e/xmi5)
        if( eNorm_factor_i6>0) &
        vc6_mks = sqrt(2d0*1d3*eNorm_factor_i6*e/xmi6)

        if(myid.eq.0)write(*,*) 'dlg: vce_mks post load1: ', vce_mks
        if(myid.eq.0)write(*,*) 'dlg: vc1_mks post load1: ', vc1_mks
        if(myid.eq.0)write(*,*) 'dlg: vc2_mks post load1: ', vc2_mks


      xmax = rwright - rwleft
      ymax = ytop - ybottom

      if(ndisti1   .ne. 0 .or. &
         ndisti2   .ne. 0 .or. &
         ndisti3   .ne. 0 .or. &
         ndisti4   .ne. 0 .or. &
         ndisti5   .ne. 0 .or. &
         ndisti6   .ne. 0 .or. &
         ndiste    .ne. 0)   then

         deallocate( fe_cql_cart )
         deallocate( f1_cql_cart )
         deallocate( f2_cql_cart )
         deallocate( f3_cql_cart )
         deallocate( f4_cql_cart )
         deallocate( f5_cql_cart )
         deallocate( f6_cql_cart )
      end if

!--------------------------------------------
!--   Define x mesh: x(i), xprime(i), capr(i)
!--------------------------------------------
!--   xprime: 0 to xmax
!--   x(i) : -xmax / 2.0   to   xmax / 2.0
      dx = xmax / nnodex



      do i = 1, nnodex
         xprime(i) = (i-1) * dx &
           + dx / 2.0
!--   Note: the code gives slightly smoother results with dx/2.0 added
         x(i) = xprime(i) + xwleft
         capr(i) = rt + x(i)
!         write(6, 1312)i, x(i), xprime(i), capr(i)

         xkphi(i) = nphi / capr(i)

      end do
      !dxc = xmax / nnodecx
      !do i = 1, nnodecx
      !   xprimec(i) = (i-1) * dxc + dxc / 2.0
      !   xcourse(i) = xprimec(i) + xwleft
      !   caprc(i) = rt + xcourse(i)

      !end do

      if(rzoom1 .eq. 0.0)rzoom1 = capr(1)
      if(rzoom2 .eq. 0.0)rzoom2 = capr(nnodex)

      izoom1 = int((rzoom1 - rwleft - dx / 2.0) / dx) + 1
      izoom2 = int((rzoom2 - rwleft - dx / 2.0) / dx) + 1
!      write(6, 1312) izoom1
!      write(6, 1312) izoom2

!-----------------------------------
!     Define y mesh: y(j), yprime(j)
!-----------------------------------
!--   yprime: 0 to ymax
!--   y(j) : -ymax / 2.0   to   ymax / 2.0
      dy = ymax / nnodey


      do j = 1, nnodey
         yprime(j) = (j-1) * dy &
            + dy / 2.0
!--      Note: the code gives slightly smoother results with dy/2.0 added
         y(j) = yprime(j) + ybottom

!         write(6, 1312)j, y(j), yprime(j)
      end do

      !dyc = ymax / nnodecy
      !do j = 1, nnodecy
      !   yprimec(j) = (j-1) * dyc + dyc / 2.0
      !   ycourse(j) = yprimec(j) + ybottom
!     !    write(6, 1312)j, caprc(j)
      !end do


      if(yzoom1 .eq. 0.0)yzoom1 = ybottom
      if(yzoom2 .eq. 0.0)yzoom2 = ytop

      jzoom1 = int((yzoom1 - ybottom - dy / 2.0) / dy) + 1
      jzoom2 = int((yzoom2 - ybottom - dy / 2.0) / dy) + 1
!      write(6, 1312) jzoom1
!      write(6, 1312) jzoom2


!--------------------------------------
!     find jmid = j in equatorial plane
!--------------------------------------
      ydiff_min = 100.0
      do j = 1, nnodey
         ydiff = abs(y(j) - zmaxis)
       if (ydiff .lt. ydiff_min) then
          jmid = j
          ydiff_min = ydiff
       end if
      end do

      jequat = jmid



!-------------------------------------------
!     define theta(i,j)
!        theta0(i,j) : -pi / 2.0  to  pi / 2.0
!        theta(i,j)  :       0.0  to  2.0 * pi
!-------------------------------------------
      do i = 1, nnodex
         do j = 1, nnodey
            if(y(j) .ne. 0.0 .or. x(i) .ne. 0.0) then
               theta0(i,j) = atan2(y(j), x(i))
               if(theta0(i,j) .ge. 0.0) theta(i,j) = theta0(i,j)
               if(theta0(i,j) .lt. 0.0) theta(i,j) = &
                                             theta0(i,j) + 2.0 * pi
            end if
         end do
      end do

!----------------------------------------------------
!     Calculate derivatives of rho(i,j) and theta(i,j)
!-----------------------------------------------------
      drhodr = 0.
      drhodz = 0.
      dthedr = 0.
      dthedz = 0.

      do i = 2, nnodex - 1
         do j = 2, nnodey - 1
            drhodr(i,j) = (rho(i+1, j) - rho(i-1, j)) / (2.0 * dx)
          drhodz(i,j) = (rho(i, j+1) - rho(i, j-1)) / (2.0 * dy)
          dthedr(i,j) = (theta(i+1, j) - theta(i-1, j)) / (2.0 * dx)
          dthedz(i,j) = (theta(i, j+1) - theta(i, j-1)) / (2.0 * dy)
         end do
      end do


!-----------------------------
!     Define rho mesh: rhon(n)
!------------------------------
!--   rhon: 0 to rhomax
      rhomax = 1.0
      drho = rhomax / (nnoderho - 1)
      do n = 1, nnoderho
         rhon(n) = (n-1) * drho
      end do

!---------------------------------
!     Define theta mesh: thetam(m)
!---------------------------------
!--   thetam: 0 to 2pi
      thetamax = 2.0 * pi
      dtheta = thetamax / (mnodetheta - 1)
      do m = 1, mnodetheta
         thetam(m) = (m-1) * dtheta
       if(myid .eq. 0)write(6, *)m, thetam(m)
      end do

!---------------------------------------------------
!     Invert mesh to get capr(nnoderho, nnodetheta)
!                    and capz(nnoderho, nnodetheta)
!---------------------------------------------------
      do n = 1, nnoderho
         do m = 1, mnodetheta
          rho_in = rhon(n)
          theta_in = thetam(m)

          call invert(rho_in, theta_in, capr_out, capz_out, &
               nnodex, nnodey, nxmx, nymx, rho, theta, capr, y, myid, &
               drhodr, drhodz, dthedr, dthedz, xj)


            capr_flux(n,m) = capr_out
          capz_flux(n,m) = capz_out
!          if(myid .eq. 0)write(6,*)n, m, capr_flux(n,m), xj
         end do
      end do


      delta = sqrt(dx**2 + dy**2)

      xant = rant - rt
      iant=int((rant - rwleft) / dx) + 1
      if(rant .ne. 0.0)psiant = psi(iant, jdiag)


!--note curden is in Amps per meter of toroidal length (2.*pi*rt).
      xjantx = curdnx / dx
      xjanty = curdny / dx
      xjantz = curdnz / dx
      xjant=sqrt(xjantx**2 + xjanty**2 + xjantz**2)

      rlim   = rt + alim




!----------------------------------------------------
!     Solovev flux surfaces with analytic derivatives:
!----------------------------------------------------
      if(igeom .eq. 2)then

         psi_lim = xiota0 * b0 / 2.0 * &
             ( (rlim**2 - rt**2)**2 / 4. / rt**2 )


         do i = 1, nnodex
            do j = 1, nnodey
               psi1 = b0 * xiota0 / 2. * &
                   ((capr(i) * y(j) / rt / ekappa)**2 &
                   + (capr(i)**2 - rt**2)**2 / 4. / rt**2 )

               psi(i,j) = psi1 / psi_lim
               if(psi(i,j) .eq. 0.0)psi(i,j) = 1.0e-10


               rho(i,j) = sqrt(psi(i,j))


               dxpsi = b0 * xiota0 / (2. * psi_lim) * &
                   (2. * capr(i) * (y(j) / rt / ekappa)**2 &
                   + (capr(i)**2 - rt**2) * capr(i) / rt**2 )

               dxxpsi = b0 * xiota0 / (2. * psi_lim) * &
                   (2. * (y(j) / rt / ekappa)**2 &
                   + 3.0 * capr(i)**2 / rt**2 - 1.0 )

               dxypsi = b0 * xiota0 / psi_lim * &
                   2.0 * capr(i) * y(j) / (rt * ekappa)**2

               dypsi = b0 * xiota0 / psi_lim * &
                   capr(i)**2 * y(j) / (rt * ekappa)**2

               dyypsi = b0 * xiota0 / psi_lim * &
                             (capr(i) / rt / ekappa)**2




               denom = ekappa / 2. / capr(i) *(capr(i)**2 - rt**2)
               if(y(j) .ne. 0.0 .or. denom .ne. 0.0) then
                  theta0(i,j) = atan2(y(j), denom)
                  if(theta0(i,j) .ge. 0.0) theta(i,j) = theta0(i,j)
                  if(theta0(i,j) .lt. 0.0) theta(i,j) = &
                                             theta0(i,j) + 2.0 * pi
               end if




               bx(i,j) = -b0 * xiota0 * capr(i) * y(j) &
                   / rt**2 / ekappa**2
               by(i,j) = xiota0 * b0 / 2.0 * &
                      ( 2.0 * y(j)**2 / rt**2 / ekappa**2 &
                       + capr(i)**2 / rt**2 - 1.0)



               gaussian =  exp(-psi(i,j) / psipne)

!--            use for regular runs (default):
               if (iqprof .eq. 1)then
                  frho = gaussian
                  dxfrho = - frho / psipne * dxpsi
                  dyfrho = - frho / psipne * dypsi
                  dxxfrho = dxfrho**2 / frho - frho / psipne * dxxpsi
                  dyyfrho = dyfrho**2 / frho - frho / psipne * dyypsi
                  dxyfrho = dxfrho * dyfrho / frho &
                                              - frho / psipne * dxypsi
               end if

!--            use for TAE modes:
               if (iqprof .eq. 2)then
                  frho = gaussian**(0.5)
               end if

               bx(i,j) = bx(i,j) * frho
               by(i,j) = by(i,j) * frho
               bz(i,j) = b0 * rt / capr(i)

               dxbx = bx(i,j) * (1. / capr(i) + dxfrho / frho)

               dybx = - b0 * xiota0 * capr(i) &
                   / rt**2 / ekappa**2 * frho + bx(i,j) / frho * dyfrho

               dxxbx = dxbx**2 / bx(i,j) + bx(i,j) * (-1.0 / capr(i)**2 &
                  - (dxfrho / frho)**2 + dxxfrho / frho)

               dyybx = - b0 * xiota0 * capr(i) &
                   / rt**2 / ekappa**2 * dyfrho + dybx * dyfrho / frho &
                   + bx(i,j) * (- dyfrho**2 / frho**2 + dyyfrho / frho)

               dxybx = dxbx * dybx / bx(i,j) + bx(i,j) * ( &
                  - dxfrho * dyfrho / frho**2 + dxyfrho / frho)



               dxby = xiota0 * b0 * capr(i) / rt**2 * frho &
                                             + by(i,j) / frho * dxfrho
               dyby = xiota0 * b0 * 2.0 * y(j) / rt**2 &
                          / ekappa**2 * frho + by(i,j) / frho * dyfrho

               dxxby =  xiota0 * b0 / rt**2 * (frho + capr(i) * dxfrho) &
                  + dxby * dxfrho / frho  + by(i,j) * &
                              ( - (dxfrho / frho)**2  + dxxfrho / frho)
               dyyby =  xiota0 * b0 * 2.0 / rt**2 / ekappa**2 &
                      * (frho + y(j) * dyfrho) + dyby * dyfrho / frho &
                  + by(i,j) * (- (dyfrho / frho)**2  + dyyfrho / frho)
               dxyby =  xiota0 * b0 * capr(i) / rt**2 * dyfrho + &
                   dyby * dxfrho / frho + by(i,j) * &
                      (- dxfrho * dyfrho / frho**2  + dxyfrho / frho)


               dxbz = - bz(i,j) / capr(i)
               dxxbz = - dxbz / capr(i) + bz(i,j) / capr(i)**2
               dybz = 0.0
               dyybz = 0.0
               dxybz = 0.0



               bmod(i, j) = sqrt(bx(i, j)**2 + by(i, j)**2 &
                            + bz(i, j)**2)


               dxdth(i, j) = - y(j) / ekappa
               dzdth(i, j) = ekappa * (capr(i)**2 - rt**2) / &
                  (2. *capr(i)) + y(j)**2 / (capr(i) * ekappa)
               xntau(i, j) = sqrt(dxdth(i, j)**2 + dzdth(i, j)**2)
               if(xntau(i,j) .eq. 0.0)xntau(i,j) = 1.0e-08


               btau(i,j) = 1.0 / xntau(i, j) * &
                     (dxdth(i, j) * bx(i, j) + dzdth(i, j) * by(i, j)) &
                     / bmod(i,j)
               bpol(i,j) = btau(i,j) * bmod(i,j)
               bzeta(i,j) = bz(i,j) / bmod(i,j)

               bxn(i,j) = bx(i,j) / bmod(i,j)
               byn(i,j) = by(i,j) / bmod(i,j)
               bzn(i,j) = bz(i,j) / bmod(i,j)

               dxbmod = bxn(i,j) * dxbx &
                      + byn(i,j) * dxby &
                      + bzn(i,j) * dxbz
               dybmod = bxn(i,j) * dybx &
                      + byn(i,j) * dyby &
                      + bzn(i,j) * dybz

               dxbxn = bxn(i,j) * (dxbx / bx(i,j) - dxbmod / bmod(i,j))
               dybxn = bxn(i,j) * (dybx / bx(i,j) - dybmod / bmod(i,j))

               dxbyn = byn(i,j) * (dxby / by(i,j) - dxbmod / bmod(i,j))
               dybyn = byn(i,j) * (dyby / by(i,j) - dybmod / bmod(i,j))

               dxbzn = bzn(i,j) * (dxbz / bz(i,j) - dxbmod / bmod(i,j))
               dybzn = bzn(i,j) * (dybz / bz(i,j) - dybmod / bmod(i,j))



               dxxbmod = dxbxn * dxbx &
                       + dxbyn * dxby &
                       + dxbzn * dxbz &
                       + bxn(i,j) * dxxbx &
                       + byn(i,j) * dxxby &
                       + bzn(i,j) * dxxbz
               dyybmod = dybxn * dybx &
                       + dybyn * dyby &
                       + dybzn * dybz &
                       + bxn(i,j) * dyybx &
                       + byn(i,j) * dyyby &
                       + bzn(i,j) * dyybz
               dxybmod = dybxn * dxbx &
                       + dybyn * dxby &
                       + dybzn * dxbz &
                       + bxn(i,j) * dxybx &
                       + byn(i,j) * dxyby &
                       + bzn(i,j) * dxybz




               dxxbxn = dxbxn**2 / bxn(i,j) &
                   + bxn(i,j) * (-dxbx**2 / bx(i,j)**2 + dxxbx / bx(i,j) &
                   + (dxbmod / bmod(i,j))**2  - dxxbmod / bmod(i,j))
               dyybxn = dybxn**2 / bxn(i,j) &
                   + bxn(i,j) * (-dybx**2 / bx(i,j)**2 + dyybx / bx(i,j) &
                   + (dybmod / bmod(i,j))**2  - dyybmod / bmod(i,j))
               dxybxn = dxbxn * dybxn / bxn(i,j) &
                   + bxn(i,j) * (- dxbx * dybx / bx(i,j)**2 &
                   + dxybx / bx(i,j) &
                   + dxbmod * dybmod / bmod(i,j)**2 &
                   - dxybmod / bmod(i,j))


               dxxbyn = dxbyn**2 / byn(i,j) &
                   + byn(i,j) * (-dxby**2 / by(i,j)**2 + dxxby / by(i,j) &
                   + (dxbmod / bmod(i,j))**2  - dxxbmod / bmod(i,j))
               dyybyn = dybyn**2 / byn(i,j) &
                   + byn(i,j) * (-dyby**2 / by(i,j)**2 + dyyby / by(i,j) &
                   + (dybmod / bmod(i,j))**2  - dyybmod / bmod(i,j))
               dxybyn = dxbyn * dybyn / byn(i,j) &
                   + byn(i,j) * (- dxby * dyby / by(i,j)**2 &
                   + dxyby / by(i,j) &
                   + dxbmod * dybmod / bmod(i,j)**2 &
                   - dxybmod / bmod(i,j))


               dxxbzn = dxbzn**2 / bzn(i,j) &
                   + bzn(i,j) * (-dxbz**2 / bz(i,j)**2 + dxxbz / bz(i,j) &
                   + (dxbmod / bmod(i,j))**2  - dxxbmod / bmod(i,j))
               dyybzn = dybzn**2 / bzn(i,j) &
                   + bzn(i,j) * (-dybz**2 / bz(i,j)**2 + dyybz / bz(i,j) &
                   + (dybmod / bmod(i,j))**2  - dyybmod / bmod(i,j))
               dxybzn = dxbzn * dybzn / bzn(i,j) &
                   + bzn(i,j) * (- dxbz * dybz / bz(i,j)**2 &
                   + dxybz / bz(i,j) &
                   + dxbmod * dybmod / bmod(i,j)**2 &
                   - dxybmod / bmod(i,j))


               xiota(i,j) = btau(i,j)/ bzeta(i,j) * capr(i) / xntau(i,j)
               if(xiota(i,j) .eq. 0.0) xiota(i,j) = 1.0e-06
               qsafety(i,j) = 1.0 / xiota(i,j)

!               if(j .eq. jequat)write(6, 1312) i, capr(i), btau(i,j),
!     1             bzeta(i,j)

               sqx = sqrt(1.0 - bxn(i,j)**2)

               dxsqx = - bxn(i,j) / sqx * dxbxn
               dysqx = - bxn(i,j) / sqx * dybxn
               dxxsqx =  dxsqx * (dxbxn / bxn(i,j) - dxsqx / sqx) &
                        - bxn(i,j) / sqx * dxxbxn
               dyysqx =  dysqx * (dybxn / bxn(i,j) - dysqx / sqx) &
                        - bxn(i,j) / sqx * dyybxn
               dxysqx =  dxsqx * (dybxn / bxn(i,j) - dysqx / sqx) &
                        - bxn(i,j) / sqx * dxybxn


               uxx(i, j) =   sqx
               uxy(i, j) = - bxn(i, j) * byn(i, j) / sqx
               uxz(i, j) = - bxn(i, j) * bzn(i, j) / sqx
               uyx(i, j) =   0.0
               uyy(i, j) =   bzn(i, j) / sqx
               uyz(i, j) = - byn(i, j) / sqx
               uzx(i, j) =   bxn(i, j)
               uzy(i, j) =   byn(i, j)
               uzz(i, j) =   bzn(i, j)

               dxuxx(i,j)  = dxsqx
               dxxuxx(i,j) = dxxsqx
               dyuxx(i,j)  = dysqx
               dyyuxx(i,j) = dyysqx
               dxyuxx(i,j) = dxysqx


               dxuxy(i,j)  = uxy(i,j) * (dxbxn / bxn(i,j) &
                                       + dxbyn / byn(i,j) - dxsqx / sqx)
               dyuxy(i,j)  = uxy(i,j) * (dybxn / bxn(i,j) &
                                       + dybyn / byn(i,j) - dysqx / sqx)
               dxxuxy(i,j) = dxuxy(i,j)**2 / uxy(i,j) &
                     + uxy(i,j) * (dxxbxn / bxn(i,j) &
                                 + dxxbyn / byn(i,j) &
                                 - dxxsqx / sqx  ) &
                     - uxy(i,j) * (dxbxn**2 / bxn(i,j)**2 &
                                 + dxbyn**2 / byn(i,j)**2 &
                                 - dxsqx**2 / sqx**2 )
               dyyuxy(i,j) = dyuxy(i,j)**2 / uxy(i,j) &
                     + uxy(i,j) * (dyybxn / bxn(i,j) &
                                 + dyybyn / byn(i,j) &
                                 - dyysqx / sqx  ) &
                     - uxy(i,j) * (dybxn**2 / bxn(i,j)**2 &
                                 + dybyn**2 / byn(i,j)**2 &
                                 - dysqx**2 / sqx**2 )
               dxyuxy(i,j) = dxuxy(i,j) * dyuxy(i,j) / uxy(i,j) &
                     + uxy(i,j) * (dxybxn / bxn(i,j) &
                                 + dxybyn / byn(i,j) &
                                 - dxysqx / sqx  ) &
                     - uxy(i,j) * (dxbxn * dybxn / bxn(i,j)**2 &
                                +  dxbyn * dybyn / byn(i,j)**2 &
                                -  dxsqx * dysqx / sqx**2   )


               dxuxz(i,j)  = uxz(i,j) * (dxbxn / bxn(i,j) &
                                       + dxbzn / bzn(i,j) - dxsqx / sqx)
               dyuxz(i,j)  = uxz(i,j) * (dybxn / bxn(i,j) &
                                       + dybzn / bzn(i,j) - dysqx / sqx)
               dxxuxz(i,j) = dxuxz(i,j)**2 / uxz(i,j) &
                     + uxz(i,j) * (dxxbxn / bxn(i,j) &
                                 + dxxbzn / bzn(i,j) &
                                 - dxxsqx / sqx  ) &
                     - uxz(i,j) * (dxbxn**2 / bxn(i,j)**2 &
                                 + dxbzn**2 / bzn(i,j)**2 &
                                 - dxsqx**2 / sqx**2 )
               dyyuxz(i,j) = dyuxz(i,j)**2 / uxz(i,j) &
                     + uxz(i,j) * (dyybxn / bxn(i,j) &
                                 + dyybzn / bzn(i,j) &
                                 - dyysqx / sqx  ) &
                     - uxz(i,j) * (dybxn**2 / bxn(i,j)**2 &
                                 + dybzn**2 / bzn(i,j)**2 &
                                 - dysqx**2 / sqx**2 )
               dxyuxz(i,j) = dxuxz(i,j) * dyuxz(i,j) / uxz(i,j) &
                     + uxz(i,j) * (dxybxn / bxn(i,j) &
                                 + dxybzn / bzn(i,j) &
                                 - dxysqx / sqx  ) &
                     - uxz(i,j) * (dxbxn * dybxn / bxn(i,j)**2 &
                                +  dxbzn * dybzn / bzn(i,j)**2 &
                                -  dxsqx * dysqx / sqx**2   )


               dxuyx(i,j)  = 0.0
               dxxuyx(i,j) = 0.0
               dyuyx(i,j)  = 0.0
               dyyuyx(i,j) = 0.0
               dxyuyx(i,j) = 0.0

               dxuyy(i,j) = uyy(i,j) * (dxbzn / bzn(i,j) - dxsqx / sqx)
               dyuyy(i,j) = uyy(i,j) * (dybzn / bzn(i,j) - dysqx / sqx)
               dxxuyy(i,j) = dxuyy(i,j)**2 / uyy(i,j) &
                  + uyy(i,j) * (dxxbzn / bzn(i,j) - dxxsqx / sqx) &
                  - uyy(i,j) * (dxbzn**2 / bzn(i,j)**2 &
                              - dxsqx**2 / sqx**2        )
               dyyuyy(i,j) = dyuyy(i,j)**2 / uyy(i,j) &
                  + uyy(i,j) * (dyybzn / bzn(i,j) - dyysqx / sqx) &
                  - uyy(i,j) * (dybzn**2 / bzn(i,j)**2 &
                              - dysqx**2 / sqx**2        )
               dxyuyy(i,j) = dyuyy(i,j) * dxuyy(i,j) / uyy(i,j) &
                  + uyy(i,j) * (dxybzn / bzn(i,j) - dxysqx / sqx) &
                  - uyy(i,j) * (dxbzn * dybzn / bzn(i,j)**2 &
                              - dxsqx * dysqx / sqx**2        )


               dxuyz(i,j) = uyz(i,j) * (dxbyn / byn(i,j) - dxsqx / sqx)
               dyuyz(i,j) = uyz(i,j) * (dybyn / byn(i,j) - dysqx / sqx)
               dxxuyz(i,j) = dxuyz(i,j)**2 / uyz(i,j) &
                  + uyz(i,j) * (dxxbyn / byn(i,j) - dxxsqx / sqx) &
                  - uyz(i,j) * (dxbyn**2 / byn(i,j)**2 &
                              - dxsqx**2 / sqx**2        )
               dyyuyz(i,j) = dyuyz(i,j)**2 / uyz(i,j) &
                  + uyz(i,j) * (dyybyn / byn(i,j) - dyysqx / sqx) &
                  - uyz(i,j) * (dybyn**2 / byn(i,j)**2 &
                              - dysqx**2 / sqx**2        )
               dxyuyz(i,j) = dxuyz(i,j) * dyuyz(i,j) / uyz(i,j) &
                  + uyz(i,j) * (dxybyn / byn(i,j) - dxysqx / sqx) &
                  - uyz(i,j) * (dxbyn * dybyn / byn(i,j)**2 &
                              - dxsqx * dysqx / sqx**2        )


               dxuzx(i,j)  = dxbxn
               dxxuzx(i,j) = dxxbxn
               dyuzx(i,j)  = dybxn
               dyyuzx(i,j) = dyybxn
               dxyuzx(i,j) = dxybxn


               dxuzy(i,j)  = dxbyn
               dxxuzy(i,j) = dxxbyn
               dyuzy(i,j)  = dybyn
               dyyuzy(i,j) = dyybyn
               dxyuzy(i,j) = dxybyn


               dxuzz(i,j)  = dxbzn
               dxxuzz(i,j) = dxxbzn
               dyuzz(i,j)  = dybzn
               dyyuzz(i,j) = dyybzn
               dxyuzz(i,j) = dxybzn

!              -----
!              Exact:
!              -----
!               gradprlb(i,j) = uzx(i,j) * dxbmod + uzy(i,j) * dybmod

!              -----------------------
!              Brambilla approximation:
!              -----------------------
             sinth = y(j) / sqrt(x(i)**2 + y(j)**2)
               gradprlb(i,j) = bmod(i,j) / capr(i) &
                                             * abs(btau(i,j) * sinth)

!              -------------------------------
!              Constant gradient scale length:
!              -------------------------------
!               gradprlb(i,j) = bmod(i,j) / capr(i)



               if (nzfun .eq. 0)gradprlb(i,j) = 1.0e-10


               drhodx = 0.5 / sqrt(psi(i,j)) * dxpsi
               drhody = 0.5 / sqrt(psi(i,j)) * dypsi

               gradrho = sqrt(drhodx**2 + drhody**2)
             if (gradrho .eq. 0.0) gradrho = 1.0e-08

               rhohatx(i,j) = drhodx / gradrho
               rhohaty(i,j) = drhody / gradrho

            end do
         end do


      end if


!     ----------------------------------------
!     Set up r mesh for circular flux surfaces
!     ----------------------------------------
      j = jequat
      do i = 2, nnodex -1
!         if (myid .eq. 0)write(6, 1312)i, rho(i,j), x(i), xprime(i),
!     .       capr(i)
         if(rho(i,j) .lt. 1. .and. rho(i-1,j) .ge. 1.)caprmin = capr(i)
         if(rho(i,j) .gt. 1. .and. rho(i-1,j) .le. 1.)caprmax =capr(i-1)
      end do

      rmax = (caprmax - caprmin) / 2.0

      if(myid .eq. 0)write(6, 311)caprmax, caprmin, rmax

      do n = 1, nnoderho
         rn(n) = rhon(n) * rmax
       epsn(n) = rn(n) / rt
         dn(n) = rmax**2 / taue
!         if (myid .eq. 0)write(6, 1312)n, rn(n), rhon(n), dn(n)
!         if (myid .eq. 0)write(15, 1312)n, rn(n), rhon(n), dn(n)
      end do

      do i = 1, nnodex
         do j = 1, nnodey
            n = int(rho(i,j) / drho) + 1
            if (n .eq. 1) then
               i0 = i
               j0 = j
            end if
         end do
      end do

      if(i0 .eq. 0 .and. j0 .eq. 0) then
         do i = 1, nnodex
            do j = 1, nnodey
               n = int(rho(i,j) / drho) + 1
               if (n .eq. 2) then
                  i0 = i
                  j0 = j
               end if
            end do
         end do

      end if


      dr = drho * rmax

      do n = 1, nnoderho - 1
         rh(n) = (rn(n) + rn(n+1)) / 2.0
         dh(n) = (dn(n) + dn(n+1)) / 2.0
!         if (myid .eq. 0)write(6, 1312)n, rh(n), dh(n)
      end do
      rh(nnoderho) = 0.0
      dh(nnoderho) = 0.0



!     -------------------------------------
!     Calculate capr * bpol in the midplane
!     -------------------------------------
      do i = 1, nnodex
         do j = 1, nnodey
          capr_bpol(i,j) = capr(i) * bpol(i,j)
          ipsi(i,j)= capr(i) * bzeta(i,j) * bmod(i,j)
       end do
      end do

      do i = 1, nnodex
         do j = 1, nnodey
            call midplane(i, j, capr_bpol, capr_bpol_mid2(i,j), &
                rho, nxmx, nymx, nnodex, nnodey, capr, rt, 0.0, jmid)
         end do
      end do

      call polavg(capr_bpol_mid2, capr_bpol_mid, rho, nxmx, nymx, &
         nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvolp, fvol)


      if (myid .eq. 0)then

         write(6, *) "jmid = ", jmid
       write(15, *) "jmid = ", jmid

         write(6, *)  "psio = ", psio
       write(15, *) "psio = ", psio
       write(6, *)  "psi_tor_max = ", psi_tor_max
       write(15, *) "psi_tor_max = ", psi_tor_max

         write(15, *)"rt = ", rt
       write(15, *)"b0 = ", b0
       write(15, *)"dx = ", dx
       write(15, *)"dy = ", dy
       write(15, *)"drho = ", drho

         write(15, *)
         write(15, *) "     i    capr  rhoij  capr_bpol  capr_bpol_mid2"
       write(15, *)

         write(6, *)
         write(6, *)  "     i    capr  rhoij  capr_bpol  capr_bpol_mid2"
       write(6, *)

         do i = 1, nnodex
            write(6,  1312)i, capr(i), rho(i, jequat), &
                         capr_bpol(i, jequat), capr_bpol_mid2(i, jequat)
            write(15, 1312)i, capr(i), rho(i, jequat), &
                         capr_bpol(i, jequat), capr_bpol_mid2(i, jequat)
         end do

         write(15, *)
         write(15, *) "     n   capr_bpol_mid"

         write(6, *)
         write(6, *)  "     n   capr_bpol_mid"

         do n = 1, nnoderho
            write(6,  1312)n, capr_bpol_mid(n)
            write(15, 1312)n, capr_bpol_mid(n)
         end do

      end if




!     -----------------------
!     Calculate bmod_mid(i,j)
!     -----------------------
      do i = 1, nnodex
         do j = 1, nnodey
            call midplane(i, j, bmod, bmod_mid(i,j), rho, &
                      nxmx, nymx, nnodex, nnodey, capr, rt, b0, jmid)
         end do
      end do

      call polavg(bmod_mid, bmod_midavg, rho, nxmx, nymx, &
         nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvolp, fvol)


!      if (myid .eq. 0)then

!         write(15, *)
!         write(15, *) "     n        bmod_midavg"

!         write(6, *)
!         write(6, *)  "     n        bmod_midavg"


!         do n = 1, nnoderho
!            write(6,  1312)n, bmod_midavg(n)
!            write(15, 1312)n, bmod_midavg(n)
!         end do

!      end if



!      if (myid .eq. 0)then

!         write(6, *)
!         write(15, *)

!         do i = 1, nnodex
!            do j = 1, nnodey
!              if(j .eq. jequat) then
!                 write(6,  1312)i, x(i), bmod(i, j), bmod_mid(i, j)
!                 write(15, 1312)i, x(i), bmod(i, j), bmod_mid(i, j)
!              end if
!            end do
!         end do

!         do i = 1, nnodex
!            do j = 1, nnodey
!               if (i .eq. 5 .and. j .eq. 5)then
!                  write(6, *) "bmod_mid(i, j) = ", bmod_mid(i, j)
!                  write(6, *) "bxn(i, j) = ", bxn(i, j)
!                  write(6, *) "byn(i, j) = ", byn(i, j)
!                  write(6, *) "bzn(i, j) = ", bzn(i, j)
!               end if
!            end do
!         end do

!      end if



!     ---------------
!     Antenna current
!     ---------------

      theta_antr = theta_ant / 180. * pi

      dpsiant = dpsiant0
      dthetant = dthetant0 / 360. * 2.0 * pi
      yant_max = antlen / 2.0

      do i = 1, nnodex
         do j = 1, nnodey

            delta_theta = theta0(i,j) - theta_antr
            gaussantth = exp(-delta_theta**2 / dthetant**2)
            gausspsi = exp(-(psi(i,j) - psiant)**2 / dpsiant**2)

          shapey = 0.0
          deltay = y(j) - yant

          if(capr(i) .gt. rt) then

!              ------------------------------------------------
!              if(i_antenna .eq. 0) antenna current is Gaussian
!              ------------------------------------------------

             if (i_antenna .eq. 0) then
                shapey = exp(- deltay**2 / yant_max**2)
             end if


!              ---------------------------------------------------------------
!              if(i_antenna .eq. 1) antenna current is cos(ky * y)  (DEFAULT)
!              ------------------------------------------------ --------------
             if (i_antenna .eq. 1 .and. abs(deltay) .lt. yant_max)then
                shapey = cos(xk0 * antlc * deltay)
             end if


          end if

            xjx(i,j) = 0.0
            xjy(i,j) = xjanty * shapey  * gausspsi
            xjz(i,j) = 0.0


         end do
      end do

!DLG: read eqdsk dlg style

    dlgTime = second1(dummy)

    if ( nt == 1 ) then 

        call read_geqdsk ( eqdsk, plot = .false. )
        call bCurvature ()
        call bGradient ()
        call init_interp ()

        if ( ndisti2 .eq. 2 ) call init_particleFile ( myId )

    endif

    if (myId==0) write(167,12398) 'read_geqdsk:', &
        (dlgTime - second1(dummy))/60.0

    
!     -----------
!     Mask array:
!     -----------

      dlgTime = second1 ( dummy )

      do i = 1, nnodex
         do j = 1, nnodey

            mask(i,j) = 1
            if (psi(i,j) .gt. psilim) mask(i,j) = 0

            if (limiter_boundary) then

                eqdsk_box_mask: &
                if ( is_inside_lim ( capR(i), y(j) ) ) then
                    mask(i,j) = 1
                else
                    mask(i,j) = 0
                endif eqdsk_box_mask

                eqdsk_bbbs_mask: &
                if ( is_inside_bbbs ( capR(i), y(j) ) ) then
                    mask_bbbs(i,j) = 1
                else
                    mask_bbbs(i,j) = 0
                endif eqdsk_bbbs_mask
                
                if ( bbbsMask ) mask = mask_bbbs

            endif

         end do
      end do


      if (myId==0) write(167,12398) 'calculate mask:', &
        (dlgTime - second1(dummy))/60.0

!    ----------------
!    plasma profiles:
!    ----------------

dlgTime = second1 ( dummy )

if(psimol .eq. 1.0) flimiter = 1.0
if(psimol .ne. 1.0) flimiter = 1.0 / (1.0 + (psi(i,j) / psimol)**16)

iprofile_eq_3: & 
if (iprofile .eq. 3) then
 
    do i = 1, nnodex
        do j = 1, nnodey

        shapen = 0.0
        shapen2 = 0.0
        shapen3 = 0.0
        shapen4 = 0.0
        shapen5 = 0.0
        shapen6 = 0.0
        shapen_slo = 0.0
        
        shapete = 0.0
        shapeti = 0.0
        shapeti2 = 0.0
        shapeti3 = 0.0
        shapeti4 = 0.0
        shapeti5 = 0.0
        shapeti6 = 0.0
        
        if ( rho(i,j) .le. 1.0 ) then
        
            shapen  = 1.0 - rho(i,j)**betan
            shapen2 = 1.0 - rho(i,j)**betan2
            shapen3 = 1.0 - rho(i,j)**betan3
            shapen4 = 1.0 - rho(i,j)**betan4
            shapen5 = 1.0 - rho(i,j)**betan5
            shapen6 = 1.0 - rho(i,j)**betan6
            shapen_slo =  1.0 - rho(i,j)**betan_slo
            
            shapete  = 1.0 - rho(i,j)**betate
            shapeti  = 1.0 - rho(i,j)**betati
            shapeti2 = 1.0 - rho(i,j)**betati2
            shapeti3 = 1.0 - rho(i,j)**betati3
            shapeti4 = 1.0 - rho(i,j)**betati4
            shapeti5 = 1.0 - rho(i,j)**betati5
            shapeti6 = 1.0 - rho(i,j)**betati6
            
        endif 

        xnea(i,j) = xnlim  + (xn0 - xnlim)  * shapen**alphan   * flimiter
        xn2a(i,j) = xn2lim + (xn2 - xn2lim) * shapen2**alphan2 * flimiter
        xn3a(i,j) = xn3lim + (xn3 - xn3lim) * shapen3**alphan3 * flimiter
        xn4a(i,j) = xn4lim + (xn4 - xn4lim) * shapen4**alphan4 * flimiter
        xn5a(i,j) = xn5lim + (xn5 - xn5lim) * shapen5**alphan5 * flimiter
        xn6a(i,j) = xn6lim + (xn6 - xn6lim) * shapen6**alphan6 * flimiter
        xna_slo(i,j) = xnslolim + (xnslo - xnslolim) * shapen_slo**alphan_slo 
   
        if (limiter_boundary) then   
            if (mask(i,j) .lt. 1) then 

                xnea(i,j) = 1e-10 
                xn2a(i,j) = 1e-10
                xn3a(i,j) = 1e-10 
                xn4a(i,j) = 1e-10 
                xn5a(i,j) = 1e-10 
                xn6a(i,j) = 1e-10 

            else
                if ( mask_bbbs(i,j) .lt. 1 ) then 
                    xnea(i,j) = xn_rho2lim!density_by_gradient ( capR(i), y(j), xn_rho2lim, gradient, 1e-10) 
                    xn2a(i,j) = xn2_rho2lim
                    xn3a(i,j) = xn3_rho2lim
                    xn4a(i,j) = xn4_rho2lim
                    xn5a(i,j) = xn5_rho2lim
                    xn6a(i,j) = xn6_rho2lim
                endif
           endif

        endif

        xkte(i,j) = telimj + (t0e - telimj) * shapete**alphate * flimiter
        xkti(i,j) = tilimj + (t0i - tilimj) * shapeti**alphati * flimiter
        xkti2(i,j) = ti2limj + (t0i2 - ti2limj) * shapeti2**alphati2 * flimiter
        xkti3(i,j) = ti3limj + (t0i3 - ti3limj) * shapeti3**alphati3 * flimiter
        xkti4(i,j) = ti4limj + (t0i4 - ti4limj) * shapeti4**alphati4 * flimiter
        xkti5(i,j) = ti5limj + (t0i5 - ti5limj) * shapeti5**alphati5 * flimiter
        xkti6(i,j) = ti6limj + (t0i6 - ti6limj) * shapeti6**alphati6 * flimiter

        if (limiter_boundary) then

            if (mask(i,j) .lt. 1) then 

                xkte(i,j)   =    te_rho2lim*q
                xkti(i,j)   =    ti_rho2lim*q
                xkti2(i,j)  =   ti2_rho2lim*q
                xkti3(i,j)  =   ti3_rho2lim*q
                xkti4(i,j)  =   ti4_rho2lim*q
                xkti5(i,j)  =   ti5_rho2lim*q
                xkti6(i,j)  =   ti6_rho2lim*q
 
            else
                if ( mask_bbbs(i,j) .lt. 1 ) then 

                    xkte(i,j)   =    te_rho2lim*q
                    xkti(i,j)   =    ti_rho2lim*q
                    xkti2(i,j)  =   ti2_rho2lim*q
                    xkti3(i,j)  =   ti3_rho2lim*q
                    xkti4(i,j)  =   ti4_rho2lim*q
                    xkti5(i,j)  =   ti5_rho2lim*q
                    xkti6(i,j)  =   ti6_rho2lim*q
 
                endif
            endif

      endif

        enddo
    enddo

    dlg_particle_density: &
    if ( nDisti2 .eq. 2 .and. particleDensity ) then
    
        xn2a = dlg_getDensity ( capR, y, &
                    rMAxis = rmaxis, zMAxis = zmaxis, &
                    densityAxis = xn2 )
    
    endif dlg_particle_density
    
    dlg_limiter_boundary: &
    if ( limiter_boundary ) then 
    
            xnea   = dlg_getDensity ( capR, y, &
                           ncFileNameIn = 'dlg_profiles.nc', &
                           ncVariableNameIn = 'ne', &
                           rMAxis = rmaxis, zMAxis = zmaxis, &
                           densityAxis = xn0 )

            where ( xnea < 1.0e-10 ) xnea = 1.0e-10

            xkte   = dlg_getDensity ( capR, y, &
                           ncFileNameIn = 'dlg_profiles.nc', &
                           ncVariableNameIn = 'te' ) * q
            xkti   = dlg_getDensity ( capR, y, &
                           ncFileNameIn = 'dlg_profiles.nc', &
                           ncVariableNameIn = 'ti' ) * q
    
    endif dlg_limiter_boundary
    
    xn1a    = (xnea - z2 * xn2a &
                - z3 * xn3a &
                - z4 * xn4a &
                - z5 * xn5a &
                - z6 * xn6a &
                - z_slo * xna_slo ) / z1

    where ( xn1a < 1.0e-10 ) xn1a = 1.0e-10
  
    eta2 = xn2 / xn0
    eta3 = xn3 / xn0
    eta4 = xn4 / xn0
    eta5 = xn5 / xn0
    eta6 = xn6 / xn0
    eta_slo = xnslo / xn0
    
    eta1 = 1.0 / z1 * (1.0 - z2 * eta2 - z3 * eta3 &
               - z4 * eta4 - z5 * eta5 - z6 * eta6 &
                                       - z_slo * eta_slo)
    xn1 = eta1 * xn0

    
end if iprofile_eq_3



do i = 1, nnodex
   do j = 1, nnodey

    if (iprofile .eq. 1 .or. iprofile .eq. 2) then
    
        xn1 = xn0 / z1 * (1.0 - z2 * eta  - z3 * eta3 &
                  - z4 * eta4 - z5 * eta5 - z6 * eta6 &
                                            - z_slo * eta_slo)
        eta1 = xn1 / xn0
        xn2 = xn0 * eta
        eta2 = xn2 / xn0
        xn3 = xn0 * eta3
        xn4 = xn0 * eta4
        xn5 = xn0 * eta5
        xn6 = xn0 * eta6
        
        xn_slo = xn0 * eta_slo
        
        gaussian_ne  = exp(-psi(i,j) / psipne)
        gaussian_te  = exp(-psi(i,j) / psipte)
        gaussian_ti1 = exp(-psi(i,j) / psipti1)
        gaussian_ti2 = exp(-psi(i,j) / psipti2)
        gaussian_ti3 = exp(-psi(i,j) / psipti3)
        gaussian_ti4 = exp(-psi(i,j) / psipti4)
        gaussian_ti5 = exp(-psi(i,j) / psipti5)
        gaussian_ti6 = exp(-psi(i,j) / psipti6)
        
        parabola = (1. - psi(i,j) / psilim)
        if (parabola .le. 0.0) parabola = 0.0
        
        
        if(iprofile .eq. 1) shape = gaussian_ne
        if(iprofile .eq. 2) shape = parabola * flimiter
        
        
        xnea(i,j)  = xnlim +  (xn0 - xnlim) * shape**alphan
        xn2a(i, j) = xnea(i, j) * eta2
        xn3a(i, j) = xnea(i, j) * eta3
        xn4a(i, j) = xnea(i, j) * eta4
        xn5a(i, j) = xnea(i, j) * eta5
        xn6a(i, j) = xnea(i, j) * eta6
        xna_slo(i, j) = xnea(i,j) * eta_slo
        
        xn1a(i, j) = (xnea(i, j) - z2 * xn2a(i, j) &
                                 - z3 * xn3a(i, j) &
                                 - z4 * xn4a(i, j) &
                                 - z5 * xn5a(i, j) &
                                 - z6 * xn6a(i, j) &
                              - z_slo * xna_slo(i, j)   ) / z1
        
        if (xn1a(i, j) .le. 0.0) xn1a(i, j) = 1.0e-10
        
        
        xkte(i,j)  = telimj +  (t0e - telimj) * shape**alphate
        xkti(i,j)  = tilimj +  (t0i - tilimj) * shape**alphati
        xkti2(i,j) = ti2limj + (t0i2 - ti2limj) * shape**alphati
        xkti3(i,j) = ti3limj + (t0i3 - ti3limj) * shape**alphati
        xkti4(i,j) = ti4limj + (t0i4 - ti4limj) * shape**alphati
        xkti5(i,j) = ti5limj + (t0i5 - ti5limj) * shape**alphati
        xkti6(i,j) = ti6limj + (t0i6 - ti6limj) * shape**alphati
    
    end if

   end do
end do


!     ---------------------------------------
!     read profile data from the plasma state:
!     ---------------------------------------
      if (iprofile .eq. 5) then
!        ------------------
!        read profile data:
!        ------------------
         rewind (63)
         read(63, state)



         if(myid .eq. 0)then
          write(6, *)
            write(6, *) "ndata = ", s_nrho_n
          write(6, 66164)
66164       format(1h ,5x,5h   n  ,3x,6h rhod , &
                       6x,6h  ne  ,6x,6h  te  ,6x,6h  ti  ,6x,6hn_beam, &
                       6x,6he_beam,6x,6h zeff ,6x,6h  n1  ,6x,6h  n2  , &
                       6x,6h  n3  )
            write(6, *)


               write(15, *)
          write(15, *) "ndata = ", s_nrho_n
          write(15, 66164)
            write(15, *)

            do n = 1, s_nrho_n
               write(6, 1312)n, s_rho_n_grid(n), s_n_s(n,0), &
                   s_t_s(n,0), s_t_s(n,1), &
                   s_n_s(n,3), s_t_s(n,3), s_zeff(n), s_n_s(n,1), &
                   s_n_s(n,2), s_n_s(n,3)
               write(15, 1312)n, s_rho_n_grid(n), s_n_s(n,0), &
                   s_t_s(n,0), s_t_s(n,1), &
                   s_n_s(n,3), s_t_s(n,3), s_zeff(n), s_n_s(n,1), &
                   s_n_s(n,2), s_n_s(n,3)
            end do

          ndata = s_nrho_n

          do n = 1, ndata

             rhod(n) = s_rho_n_grid(n)

               tekev(n)  = s_t_s(n,0)
             tikev(n)  = s_t_s(n,1)
             ti2kev(n) = s_t_s(n,2)
             e_beam(n) = s_t_s(n,3)

             xneavg(n) = s_n_s(n,0)
             xn_min(n) = s_n_s(n,2)
             xn_maj(n) =  s_n_s(n,1)
             xn_beam(n)=  s_n_s(n,3)

             if(eta .gt. 0.0)then
                xn_min(n) = eta * xneavg(n)
                xn_maj(n) = (xneavg(n) - z2 * xn_min(n) &
                                              - z3 * xn_beam(n)) / z1
!                  tikev(n) = tekev(n)/3.0
               end if


               xneavg_cgs(n) = xneavg(n) * 1.0e-06
             xn_min_cgs(n) = xn_min(n) * 1.0e-06
               xn_maj_cgs(n) = xn_maj(n) * 1.0e-06
             xn_beam_cgs(n) = xn_beam(n) * 1.0e-06

             zeffavg(n) = s_zeff(n)

!             rhod(n) = sqrt(rhod(n))
          end do

          write(6, *)
            write(6, *) "rhod = "
          write(6, 100) (rhod(n), n = 1, ndata)
          write(6, *)

          write(6, *) "xneavg = "
          write(6, 100) (xneavg(n), n = 1, ndata)
          write(6, *)

          write(6, *) "xn_majority = "
          write(6, 100) (xn_maj(n), n = 1, ndata)
          write(6, *)

          write(6, *) "xn_minority = "
          write(6, 100) (xn_min(n), n = 1, ndata)
          write(6, *)

          write(6, *) "xn_beam = "
          write(6, 100) (xn_beam(n), n = 1, ndata)
          write(6, *)

          write(6, *) "zeffavg = "
          write(6, 100) (zeffavg(n), n = 1, ndata)
          write(6, *)

          write(6, *) "tekev = "
          write(6, 100) (tekev(n), n = 1, ndata)
          write(6, *)

          write(6, *) "tikev = "
          write(6, 100) (tikev(n), n = 1, ndata)
          write(6, *)

          write(6, *) "e_beam = "
          write(6, 100) (e_beam(n), n = 1, ndata)
          write(6, *)

            write(6, *) "xn_majority_cgs = "
            write(6, 100) (xn_maj_cgs(n), n = 1, ndata)
          write(6, *)

          write(6, *) "xn_minority_cgs = "
            write(6, 100) (xn_min_cgs(n), n = 1, ndata)
          write(6, *)

          write(6, *) "xneavg_cgs = "
            write(6, 100) (xneavg_cgs(n), n = 1, ndata)
          write(6, *)


          write(6, *) "ti2kev = "
          write(6, 100) (ti2kev(n), n = 1, ndata)
          write(6, *)

          write(6, *) "xme = ", xme
          write(6, *) "xmi1 = ", xmi1
          write(6, *) "xmi2 = ", xmi2
          write(6, *) "xmi3 = ", xmi3
          write(6, *) "xmi4 = ", xmi4
          write(6, *) "xmi5 = ", xmi5
          write(6, *) "xmi6 = ", xmi6

          write(6, *) "qe = ", qe
          write(6, *) "qi1 = ", qi1
          write(6, *) "qi2 = ", qi2
          write(6, *) "qi3 = ", qi3
          write(6, *) "qi4 = ", qi4
          write(6, *) "qi5 = ", qi5
          write(6, *) "qi6 = ", qi6




          write(15, *)
            write(15, *) "rhod = "
          write(15, 100) (rhod(n), n = 1, ndata)
          write(15, *)

          write(15, *) "xneavg = "
          write(15, 100) (xneavg(n), n = 1, ndata)
          write(15, *)

          write(15, *) "xn_majority = "
          write(15, 100) (xn_maj(n), n = 1, ndata)
          write(15, *)

          write(15, *) "xn_minority = "
          write(15, 100) (xn_min(n), n = 1, ndata)
          write(15, *)

          write(15, *) "xn_beam = "
          write(15, 100) (xn_beam(n), n = 1, ndata)
          write(15, *)

          write(15, *) "zeffavg = "
          write(15, 100) (zeffavg(n), n = 1, ndata)
          write(15, *)

          write(15, *) "tekev = "
          write(15, 100) (tekev(n), n = 1, ndata)
          write(15, *)

          write(15, *) "tikev = "
          write(15, 100) (tikev(n), n = 1, ndata)
          write(15, *)

          write(15, *) "e_beam = "
          write(15, 100) (e_beam(n), n = 1, ndata)
          write(15, *)

            write(15, *) "xn_majority_cgs = "
            write(15, 100) (xn_maj_cgs(n), n = 1, ndata)
          write(15, *)

          write(15, *) "xn_minority_cgs = "
            write(15, 100) (xn_min_cgs(n), n = 1, ndata)
          write(15, *)

          write(15, *) "xneavg_cgs = "
            write(15, 100) (xneavg_cgs(n), n = 1, ndata)
          write(15, *)


          write(15, *) "ti2kev = "
          write(15, 100) (ti2kev(n), n = 1, ndata)
          write(15, *)

          write(15, *) "xme = ", xme
          write(15, *) "xmi1 = ", xmi1
          write(15, *) "xmi2 = ", xmi2
          write(15, *) "xmi3 = ", xmi3
          write(15, *) "xmi4 = ", xmi4
          write(15, *) "xmi5 = ", xmi5
          write(15, *) "xmi6 = ", xmi6

          write(15, *) "qe = ", qe
          write(15, *) "qi1 = ", qi1
          write(15, *) "qi2 = ", qi2
          write(15, *) "qi3 = ", qi3
          write(15, *) "qi4 = ", qi4
          write(15, *) "qi5 = ", qi5
          write(15, *) "qi6 = ", qi6



       end if
  100    format(1p6e12.3)

       s_t_s = q * s_t_s * 1.0e3    !input is kev TEMPS IN JOULES



!        -----------------------------
!        deposit data on 2D flux grid:
!        -----------------------------

       xnea = 1.0e-10

         call flux_to_rz(nnodex,nnodey,s_n_s(1:s_nrho_n,0), &
              xnea(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n), &
              s_nrho_n, rho(1:nnodex,1:nnodey))


         call flux_to_rz(nnodex,nnodey,s_t_s(1:s_nrho_n,0), &
              xkte(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n), &
              s_nrho_n, rho(1:nnodex,1:nnodey))

             xn1a = 1.0e-10
         call flux_to_rz(nnodex,nnodey,s_n_s(1:s_nrho_n,1), &
              xn1a(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n), &
              s_nrho_n, rho(1:nnodex,1:nnodey))

             xkti = 1.0e-10 *q
         call flux_to_rz(nnodex,nnodey,s_t_s(1:s_nrho_n,1), &
              xkti(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n), &
              s_nrho_n, rho(1:nnodex,1:nnodey))

             xn2a = 1.0e-10
         call flux_to_rz(nnodex,nnodey,s_n_s(1:s_nrho_n,2), &
              xn2a(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n), &
              s_nrho_n, rho(1:nnodex,1:nnodey))

             xkti2 = 1.0e-10 *q
         call flux_to_rz(nnodex,nnodey,s_t_s(1:s_nrho_n,2), &
              xkti2(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n), &
              s_nrho_n, rho(1:nnodex,1:nnodey))

             xn3a = 1.0e-10
         call flux_to_rz(nnodex,nnodey,s_n_s(1:s_nrho_n,3), &
              xn3a(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n), &
              s_nrho_n, rho(1:nnodex,1:nnodey))

             xkti3 = 1.0e-10 *q
         call flux_to_rz(nnodex,nnodey,s_t_s(1:s_nrho_n,3), &
              xkti3(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n), &
              s_nrho_n, rho(1:nnodex,1:nnodey))

             xn4a = 1.0e-10
         call flux_to_rz(nnodex,nnodey,s_n_s(1:s_nrho_n,4), &
              xn4a(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n), &
              s_nrho_n, rho(1:nnodex,1:nnodey))

             xkti4 = 1.0e-10 *q
         call flux_to_rz(nnodex,nnodey,s_t_s(1:s_nrho_n,4), &
              xkti4(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n), &
              s_nrho_n, rho(1:nnodex,1:nnodey))

             xn5a = 1.0e-10
         call flux_to_rz(nnodex,nnodey,s_n_s(1:s_nrho_n,5), &
              xn5a(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n), &
              s_nrho_n, rho(1:nnodex,1:nnodey))

             xkti5 = 1.0e-10 *q
         call flux_to_rz(nnodex,nnodey,s_t_s(1:s_nrho_n,5), &
              xkti5(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n), &
              s_nrho_n, rho(1:nnodex,1:nnodey))

             xn6a = 1.0e-10
         call flux_to_rz(nnodex,nnodey,s_n_s(1:s_nrho_n,6), &
              xn6a(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n), &
              s_nrho_n, rho(1:nnodex,1:nnodey))

             xkti6 = 1.0e-10 *q
         call flux_to_rz(nnodex,nnodey,s_t_s(1:s_nrho_n,6), &
              xkti6(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n), &
              s_nrho_n, rho(1:nnodex,1:nnodey))

         call flux_to_rz(nnodex,nnodey,s_zeff(1:s_nrho_n), &
              zeff(1:nnodex,1:nnodey),s_rho_n_grid(1:s_nrho_n), &
              s_nrho_n, rho(1:nnodex,1:nnodey))

         do i = 1, nnodex
            do j = 1, nnodey
               if(rho(i,j) .gt. rholim) then
                  xnea(i,j) = 1.0e-10
                zeff(i,j) = 1.0e-10
                  xn3a(i,j) = 1.0e-10
                  xna_slo(i,j) = 0.0e-10
                  xkte(i,j) = 1.0e-10 * q
                  xkti(i,j) = 1.0e-10 * q
                  xkti3(i,j) = 1.0e-10 * q
               end if
            end do
         end do

!        ---------------------------------
!        Minority ion (species 2) profiles:
!        ---------------------------------
!        for now leave the species 2 minority option open
!        if eta comes in as zero or small, use namelist values
!        else create a minority

         if(eta .ge. 1.0e-6) then
            do i = 1, nnodex
               do j = 1, nnodey
                  xn2a(i,j) = eta * xnea(i,j)

!                 ----------------
!                 Cold minority H:
!                 ----------------
                  xkti2(i,j) = xkti(i,j)

!                 ----------------
!                 Hot minority H:
!                 ----------------
!                 xkti2(i,j) = xkti3(i,j)

               end do
            end do

          do i = 1, nnodex
               do j = 1, nnodey
                  xn1a(i, j) = (xnea(i, j) - z2 * xn2a(i, j) &
                                        - z3 * xn3a(i, j) &
                                        - z4 * xn4a(i, j) &
                                        - z5 * xn5a(i, j) &
                                        - z6 * xn6a(i, j) &
                                     - z_slo * xna_slo(i, j)   ) / z1

                  if (xn1a(i, j) .le. 0.0) xn1a(i, j) = 1.0e-10
                  if (xn2a(i, j) .le. 0.0) xn2a(i, j) = 1.0e-10

               end do
            end do

       end if


         xn0 = xnea(i0, j0)
         xn1 = xn1a(i0, j0)
         xn2 = xn2a(i0, j0)
         xn3 = xn3a(i0, j0)
         xn_slo = xna_slo(i0, j0)
         eta1 = xn1 / xn0
         eta2 = xn2 / xn0
         eta3 = xn3 / xn0
         eta_slo = xn_slo / xn0

       if(myid .eq. 0)then
            write(6, 6834)eta1
            write(6, 6835)eta2
            write(6, 6836)eta3

            write(15, 6834)eta1
            write(15, 6835)eta2
            write(15, 6836)eta3
       end if


         call fluxavg3(xnea, xnavg, rho, nxmx, nymx, nrhomax, &
            nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(zeff, zeffavg1, rho, nxmx, nymx, nrhomax, &
            nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xn1a, xn1avg, rho, nxmx, nymx, nrhomax, &
            nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xn2a, xn2avg, rho, nxmx, nymx, nrhomax, &
            nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xn3a, xn3avg, rho, nxmx, nymx, nrhomax, &
            nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xna_slo, xna_sloavg, rho, nxmx, nymx, nrhomax, &
            nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xkte, xkteavg, rho, nxmx, nymx, nrhomax, &
            nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xkti, xktiavg, rho, nxmx, nymx, nrhomax, &
            nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xkti2, xkti2avg, rho, nxmx, nymx, nrhomax, &
            nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         call fluxavg3(xkti3, xkti3avg, rho, nxmx, nymx, nrhomax, &
            nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol)

         xn0 = xnea(i0, j0)
         xn1 = xn1a(i0, j0)
         xn2 = xn2a(i0, j0)
         xn3 = xn3a(i0, j0)
         xn4 = xn4a(i0, j0)
         xn5 = xn5a(i0, j0)
         xn6 = xn6a(i0, j0)
         xn_slo = xna_slo(i0, j0)

         eta1 = xn1 / xn0
         eta2 = xn2 / xn0
         eta3 = xn3 / xn0
         eta4 = xn4 / xn0
         eta5 = xn5 / xn0
         eta6 = xn6 / xn0
         eta_slo = xn_slo / xn0

      end if

!     --------------------------
!     End of iprofile = 5 option
!     --------------------------

      if (myId==0) write(167,12398) 'setup profiles:', &
        (dlgTime - second1(dummy))/60.0


!     write (6, *) "b0 = ", b0
!     write (6, *) "xmu0 = ", xmu0
!     write (6, *) "xn1 = ", xn1
!     write (6, *) "xmi1 = ", xmi1

      valfven = sqrt(b0**2/(xmu0 * xn1 * xmi1))
      kalfven = omgrf / valfven

      kperprhoi1 = 0.0
      kperprhoi2 = 0.0
      kperprhoi3 = 0.0
      kperprhoi4 = 0.0
      kperprhoi5 = 0.0

      if (eta1 .ne. 0.0) kperprhoi1 = kalfven * rhoi10
      if (eta2 .ne. 0.0) kperprhoi2 = kalfven * rhoi20
      if (eta3 .ne. 0.0) kperprhoi3 = kalfven * rhoi30
      if (eta4 .ne. 0.0) kperprhoi4 = kalfven * rhoi40
      if (eta5 .ne. 0.0) kperprhoi5 = kalfven * rhoi50


!     ----------------------------------------------
!     Calculate the differential volume on half mesh:
!     ----------------------------------------------

      dvol  = 0.0
      darea = 1.0e-05

      do i = 1, nnodex
         do j = 1, nnodey
            n = int(rho(i,j) / drho) + 1
            if(n .le. nnoderho)then
               dvol(n)  =  dvol(n) + dx * dy * 2.0 * pi * capr(i)
             darea(n) = darea(n) + dx * dy * capr(i) / rt
            end if
         end do
      end do


!     --------------------------------------------
!     Calculate the integrated volume on even mesh:
!     --------------------------------------------

      volume = 0.0

      volume(1) = 0.0
      do n = 1, nnoderho - 1
         volume(n+1) = volume(n) + dvol(n)
      end do




!     ------------------------
!     Do flux surface averages:
!     ------------------------

      call fluxavg(xnea, xnavg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xn1a, xn1avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xn2a, xn2avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xn3a, xn3avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xn4a, xn4avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xn5a, xn5avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xn6a, xn6avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xna_slo, xna_sloavg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)


      call fluxavg(xkte, xkteavg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xkti, xktiavg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xkti2, xkti2avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xkti3, xkti3avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xkti4, xkti4avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xkti5, xkti5avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(xkti6, xkti6avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)


!      call polavg(dldb_tot12, dldbavg, rho, nxmx, nymx, nrhomax,
!     .   nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvolp, fvol)


      if(myid .eq. 0) then
         write(6, *)
         write(6, *)"        i    taup        tauii        nustar", &
                      "     epsa"
         write(6, *)

         write(15,*)
         write(15,*)"        i    taup        tauii        nustar", &
                      "     epsa"
       write(15,*)
      end if

      do i = 1, nnodex
         do j = 1, nnodey

            rhome  = xnea(i, j) * xme
            rhomi1 = xn1a(i, j) * xmi1
            rhomi2 = xn2a(i, j) * xmi2
            rhomi3 = xn3a(i, j) * xmi3
            rhomi4 = xn4a(i, j) * xmi4
            rhomi5 = xn5a(i, j) * xmi5
            rhomi6 = xn6a(i, j) * xmi6

            rhom1(i,j) = xn1a(i, j) * xmi1


            presse  = xnea(i, j) * xkte(i, j)
          pressi1 = xn1a(i, j) * xkti(i, j)
          pressi2 = xn2a(i, j) * xkti2(i, j)
          pressi3 = xn3a(i, j) * xkti3(i, j)
          pressi4 = xn4a(i, j) * xkti4(i, j)
          pressi5 = xn5a(i, j) * xkti5(i, j)
          pressi6 = xn6a(i, j) * xkti6(i, j)

!          pressi(i, j) = presse + pressi1 pressi2 + pressi3
!     .                             + pressi4 pressi5 + pressi6
                pressi(i, j) = pressi1

          xnuee = 2.859e-27 * xnea(i,j) * xlnlam / &
                            (sqrt(xme) * (xkte(i,j)/q)**1.5)

            xnuii = 2.859e-27 * xn1a(i,j) * xlnlam / &
                            (sqrt(xmi1) * (xkti(i,j)/q)**1.5)

            if(xnuii .ne. 0.0)then

               vthi0 = sqrt(2.0 * xkti(i,j) / xmi1)
               vthe0 = sqrt(2.0 * xkte(i,j) / xme)

               tauii = 1.0 / xnuii
               omgti =  vthi0 / (rt * qhat)
               omgte =  vthe0 / (rt * qhat)
               xnu7omg = xnuii / omgti

!              -----------------------------
!              Assume circular flux surfaces
!              -----------------------------

             a = sqrt((x(i)-r0)**2 + (y(j)-z0)**2)

               rmaxa = rt + a
               rmina = rt - a

               epsa = (rmaxa - rmina) / (rmaxa + rmina)
             if( epsa .eq. 0.0) epsa = 1.0e-08

             xnustar = xnuii / (omgti * epsa**1.5)
             nu_star(i, j) = xnustar

!              -----------------------------------------------
!              bananna (collisionless) regime:   xnuii/omgti << eps**1.5
!              -----------------------------------------------
!              taup_b = 1.29 * epsa**1.5 * tauii

!              -------------------------------------------------
!              plateau regime:    eps**1.5  << xnuii/omgti << 1
!              -------------------------------------------------
!              taup_pl = 2.0 / sqrt(pi) / omgti

!              --------------------------------------------
!              Pfirsch - Schluter regime:  xnuii/omgti >> 1
!              --------------------------------------------
!              taup_ps = 1.0 / ( omgti * tauii) * (1.0 / omgti)


!              ----------------------------------------
!              Mult-regime model:  arbitray xnuii/omgti
!              ----------------------------------------
             taup_inv = 1.5 * omgti * xnustar / (1.0 + xnustar) &
                                   * 1.0 / (1.0 + xnustar * epsa**1.5)
               taup = 1.0 / taup_inv

             if(myid .eq. 0) then
                if(j .eq. 16) then
                     write(6, 1312) i, taup, tauii, xnustar, epsa
                   write(15,1312) i, taup, tauii, xnustar, epsa
              end if
               end if

!              ---------------------
!              poloidal viscosities:
!              ---------------------
               muhat(i, j) = epsa**2 / taup

          end if





            signb = b0 / abs(b0)

            omgce(i, j)  = qe  * bmod(i, j) / xme  * signb
            omgci1(i, j) = qi1 * bmod(i, j) / xmi1 * signb
            omgci2(i, j) = qi2 * bmod(i, j) / xmi2 * signb
            omgci3(i, j) = qi3 * bmod(i, j) / xmi3 * signb
            omgci4(i, j) = qi4 * bmod(i, j) / xmi4 * signb
            omgci5(i, j) = qi5 * bmod(i, j) / xmi5 * signb
            omgci6(i, j) = qi6 * bmod(i, j) / xmi6 * signb
            omgci_slo(i, j) = qi_slo * bmod(i, j) / xmi_slo * signb

            omgpe2(i, j) = xnea(i ,j) * qe**2  / (eps0 * xme)
            omgp12(i, j) = xn1a(i, j) * qi1**2 / (eps0 * xmi1)
            omgp22(i, j) = xn2a(i, j) * qi2**2 / (eps0 * xmi2)
            omgp32(i, j) = xn3a(i, j) * qi3**2 / (eps0 * xmi3)
            omgp42(i, j) = xn4a(i, j) * qi4**2 / (eps0 * xmi4)
            omgp52(i, j) = xn5a(i, j) * qi5**2 / (eps0 * xmi5)
            omgp62(i, j) = xn6a(i, j) * qi6**2 / (eps0 * xmi6)
            omgp2_slo(i, j) = xna_slo(i, j) * qi_slo**2 &
                                               / (eps0 * xmi_slo)


            xb(i,j) = -zi / omgrf / eps0 * xjx(i,j)
            xc(i,j) = -zi / omgrf / eps0 * xjy(i,j)
            xd(i,j) = -zi / omgrf / eps0 * xjz(i,j)

            capk_perp = 1.0 - omgpe2(i,j) / (omgrf**2 - omgce(i,j)**2) &
                            - omgp12(i,j) / (omgrf**2 - omgci1(i,j)**2) &
                            - omgp22(i,j) / (omgrf**2 - omgci2(i,j)**2) &
                            - omgp32(i,j) / (omgrf**2 - omgci3(i,j)**2) &
                            - omgp42(i,j) / (omgrf**2 - omgci4(i,j)**2) &
                            - omgp52(i,j) / (omgrf**2 - omgci5(i,j)**2) &
                            - omgp62(i,j) / (omgrf**2 - omgci6(i,j)**2)

            capk_x =(omgce(i,j) * omgpe2(i,j)/(omgrf**2-omgce(i,j)**2) &
                  + omgci1(i,j) * omgp12(i,j)/(omgrf**2-omgci1(i,j)**2) &
                  + omgci2(i,j) * omgp22(i,j)/(omgrf**2-omgci2(i,j)**2) &
                  + omgci3(i,j) * omgp32(i,j)/(omgrf**2-omgci3(i,j)**2) &
                  + omgci4(i,j) * omgp42(i,j)/(omgrf**2-omgci4(i,j)**2) &
                  + omgci5(i,j) * omgp52(i,j)/(omgrf**2-omgci5(i,j)**2) &
                  + omgci6(i,j) * omgp62(i,j)/(omgrf**2-omgci6(i,j)**2)) &
                   / omgrf


            xnphi = bzn(i,j) * xkphi(i) / xk0
!            acold(i,j) = xnphi**2 - capk_perp

          acold(i,j) = xk0**2 * capk_perp - xkphi(i)**2
          bcold(i,j) = xk0**2 * (capk_perp - capk_x) - xkphi(i)**2
          ccold(i,j) = xk0**2 * (capk_perp + capk_x) - xkphi(i)**2

          xkperp_cold2 = bcold(i,j) * ccold(i,j) / acold(i,j)
          xkperp_cold(i,j) = csqrt(xkperp_cold2)

         end do
      end do




      if (myid .eq. 0) then

         write (15, *) 'valfven = ', valfven, 'm/s'
       write (6,  *) 'valfven = ', valfven, 'm/s'

       write (15, *) 'kalfven = ', kalfven, 'm-1'
       write (6,  *) 'kalfven = ', kalfven, 'm-1'

       write (15, *) 'kperprhoi1 = ', kperprhoi1
       write (15, *) 'kperprhoi2 = ', kperprhoi2
       write (15, *) 'kperprhoi3 = ', kperprhoi3
       write (15, *) 'kperprhoi4 = ', kperprhoi4
       write (15, *) 'kperprhoi5 = ', kperprhoi5

       write (6, *) 'kperprhoi1 = ', kperprhoi1
       write (6, *) 'kperprhoi2 = ', kperprhoi2
       write (6, *) 'kperprhoi3 = ', kperprhoi3
       write (6, *) 'kperprhoi4 = ', kperprhoi4
       write (6, *) 'kperprhoi5 = ', kperprhoi5

         write (15, *) 'eta1 = ', eta1
         write (15, *) 'eta2 = ', eta2
         write (15, *) 'eta3 = ', eta3
         write (15, *) 'eta4 = ', eta4
         write (15, *) 'eta5 = ', eta5
         write (15, *) 'eta6 = ', eta6

         write (6, *) 'eta1 = ', eta1
         write (6, *) 'eta2 = ', eta2
         write (6, *) 'eta3 = ', eta3
         write (6, *) 'eta4 = ', eta4
         write (6, *) 'eta5 = ', eta5
         write (6, *) 'eta6 = ', eta6

         write(15, 920)
         write(15, 164)
         write(15, 920)

         write(6, 920)
         write(6, 164)
         write(6, 920)

      endif


      do i = 1, nnodex
       do j = 1, nnodey
            reomg1a(i,j) = omgrf / omgci1(i, j)
            reomg2a(i,j) = omgrf / omgci2(i, j)
            reomg3a(i,j) = omgrf / omgci3(i, j)
         end do
      end do

      j = jequat
      do i = 1, nnodex
         reomg1 = omgrf / omgci1(i, j)
         reomg2 = omgrf / omgci2(i, j)
         reomg3 = omgrf / omgci3(i, j)
         reomg4 = omgrf / omgci4(i, j)
         reomg5 = omgrf / omgci5(i, j)
         reomg6 = omgrf / omgci6(i, j)

         teev = xkte(i,j) / q
       re_acold = real(acold(i,j))

         if (myid .eq. 0) then
            write(6, 2163)i, capr(i), x(i), btau(i, j), &
               bmod(i,j), reomg1, reomg2, reomg3, reomg4, reomg5, &
               xnea(i,j), re_acold
            write(15, 2163)i, capr(i), x(i), btau(i, j), &
               bmod(i,j), reomg1, reomg2, reomg3, reomg4, reomg5, &
               xnea(i,j), re_acold
         endif
      end do

      if (myid .eq. 0) then
         write(15, 920)
         write(15, 9164)
         write(15, 920)

         write(6, 920)
         write(6, 9164)
         write(6, 920)
      endif

      i = icenter
      do j = 1, nnodey
         reomg1 = omgrf / omgci1(i, j)
         reomg2 = omgrf / omgci2(i, j)
         reomg3 = omgrf / omgci3(i, j)
         reomg4 = omgrf / omgci4(i, j)
         reomg5 = omgrf / omgci5(i, j)
         reomg6 = omgrf / omgci6(i, j)

         teev = xkte(i,j) / q

         if (myid .eq. 0) then
            write(6, 2163)j, capr(i), y(j), btau(i, j), &
               bmod(i,j), reomg1, reomg2, reomg3, reomg4, reomg5, &
               xnea(i,j), re_acold
            write(15, 2163)j, capr(i), y(j), btau(i, j), &
               bmod(i,j), reomg1, reomg2, reomg3, reomg4, reomg5, &
               xnea(i,j), re_acold
         endif
      end do



      if (myid .eq. 0) then
         write(15, 920)
         write(15, 8164)
         write(15, 920)

         write(6, 920)
         write(6, 8164)
         write(6, 920)
      endif

 8164 format(1h ,1x,4h i  ,  3x,6h R(x) , &
                 6x,6h  x   ,6x,6h bmod ,6x,6h   n1 ,6x,6h  T1  , &
                 6x,6h  n2  ,6x,6h  T2  ,6x,6h   n3 ,6x,6h  T3  , &
                 6x,6h      )

      j = jequat
      do i = 1, nnodex
         t1ev = xkti(i,j)  / q
         t2ev = xkti2(i,j) / q
         t3ev = xkti3(i,j) / q
         t4ev = xkti4(i,j) / q
         t5ev = xkti5(i,j) / q
         t6ev = xkti6(i,j) / q

         if (myid .eq. 0) then
            write(6, 2163) i, capr(i), x(i), bmod(i,j), &
                 xn1a(i,j), t1ev, xn2a(i,j), t2ev, xn3a(i,j), t3ev
            write(15, 2163)i, capr(i), x(i), bmod(i,j), &
                 xn1a(i,j), t1ev, xn2a(i,j), t2ev, xn3a(i,j), t3ev
         endif
      end do



      do n = nkx1, nkx2
         xkxsav(n) = 2.0 * pi * n / xmax
      end do

      do m = nky1, nky2
         xkysav(m) = 2.0 * pi * m / ymax
      end do

      kperp_max_actual = sqrt(xkxsav(nkx2)**2 + xkysav(nky2)**2)

      kperp_max = kperp_max_actual * 2.0

      xk_cutoff = kperp_max_actual * xkperp_cutoff


!     -------------------------------------
!     setup Z function table when nzfun = 3
!     -------------------------------------
      if (nzfun .eq. 3) then

         iflag = 0

         call ztable_setup(myid, iflag)

         call blacs_barrier(icontxt, 'All')

         if(iflag .eq. 1)call ztable_setup(myid, nproc)

      end if


!     ---------------------------
!     precompute xx(n,i), yy(m,j)
!     ---------------------------

      do i = 1, nnodex
         do n = nkx1, nkx2
            xx(n, i) = exp(zi * xkxsav(n) * xprime(i))
            xx_inv(n,i) = 1.0/xx(n,i)
         end do
      end do

      do j = 1, nnodey
         do m = nky1, nky2
            yy(m, j) = exp(zi * xkysav(m) * yprime(j))
            yy_inv(m,j) = 1.0/yy(m,j)
         end do
      end do

      if (igeom .eq. 5) then


!     --------------------------
!     Take numerical derivatives
!     --------------------------
      do i = 1, nnodex
         do j = 1, nnodey

            call deriv_x(bmod, nxmx, nymx, i, j, nnodex, nnodey, dx, &
               dbdx, d2bdx2)
            call deriv_y(bmod, nxmx, nymx, i, j, nnodex, nnodey, dy, &
               dbdy, d2bdy2)


!           -----
!           Exact:
!           -----
!            gradprlb(i,j) = abs(uzx(i,j) * dbdx + uzy(i,j) * dbdy)

!           -----------------------
!           Brambilla approximation:
!           -----------------------
          sinth = y(j) / sqrt(x(i)**2 + y(j)**2)
            gradprlb(i,j) = bmod(i,j) / capr(i) * abs(btau(i,j) * sinth)

!           -------------------------------
!           Constant gradient scale length:
!           -------------------------------
!            gradprlb(i,j) = bmod(i,j) / capr(i)


            if (nzfun .eq. 0)gradprlb(i,j) = 1.0e-10



            call deriv_x(rho, nxmx, nymx, i, j, nnodex, nnodey, dx, &
               drhodx, drhodxx)
            call deriv_y(rho, nxmx, nymx, i, j, nnodex, nnodey, dy, &
               drhody, drhodyy)

            gradrho = sqrt(drhodx**2 + drhody**2)
          if (gradrho .eq. 0.0) gradrho = 1.0e-08

            rhohatx(i,j) = drhodx / gradrho
            rhohaty(i,j) = drhody / gradrho


            call deriv_x(uxx, nxmx, nymx, i, j, nnodex, nnodey, dx, &
               dxuxx(i,j), dxxuxx(i,j))
            call deriv_x(uxy, nxmx, nymx, i, j, nnodex, nnodey, dx, &
               dxuxy(i,j), dxxuxy(i,j))
            call deriv_x(uxz, nxmx, nymx, i, j, nnodex, nnodey, dx, &
               dxuxz(i,j), dxxuxz(i,j))

            call deriv_x(uyx, nxmx, nymx, i, j, nnodex, nnodey, dx, &
               dxuyx(i,j), dxxuyx(i,j))
            call deriv_x(uyy, nxmx, nymx, i, j, nnodex, nnodey, dx, &
               dxuyy(i,j), dxxuyy(i,j))
            call deriv_x(uyz, nxmx, nymx, i, j, nnodex, nnodey, dx, &
               dxuyz(i,j), dxxuyz(i,j))

            call deriv_x(uzx, nxmx, nymx, i, j, nnodex, nnodey, dx, &
               dxuzx(i,j), dxxuzx(i,j))
            call deriv_x(uzy, nxmx, nymx, i, j, nnodex, nnodey, dx, &
               dxuzy(i,j), dxxuzy(i,j))
            call deriv_x(uzz, nxmx, nymx, i, j, nnodex, nnodey, dx, &
               dxuzz(i,j), dxxuzz(i,j))



            call deriv_y(uxx, nxmx, nymx, i, j, nnodex, nnodey, dy, &
               dyuxx(i,j), dyyuxx(i,j))
            call deriv_y(uxy, nxmx, nymx, i, j, nnodex, nnodey, dy, &
               dyuxy(i,j), dyyuxy(i,j))
            call deriv_y(uxz, nxmx, nymx, i, j, nnodex, nnodey, dy, &
               dyuxz(i,j), dyyuxz(i,j))

            call deriv_y(uyx, nxmx, nymx, i, j, nnodex, nnodey, dy, &
               dyuyx(i,j), dyyuyx(i,j))
            call deriv_y(uyy, nxmx, nymx, i, j, nnodex, nnodey, dy, &
               dyuyy(i,j), dyyuyy(i,j))
            call deriv_y(uyz, nxmx, nymx, i, j, nnodex, nnodey, dy, &
               dyuyz(i,j), dyyuyz(i,j))

            call deriv_y(uzx, nxmx, nymx, i, j, nnodex, nnodey, dy, &
               dyuzx(i,j), dyyuzx(i,j))
            call deriv_y(uzy, nxmx, nymx, i, j, nnodex, nnodey, dy, &
               dyuzy(i,j), dyyuzy(i,j))
            call deriv_y(uzz, nxmx, nymx, i, j, nnodex, nnodey, dy, &
               dyuzz(i,j), dyyuzz(i,j))




            call deriv_xy(uxx, nxmx, nymx, i, j, nnodex, nnodey, dx, dy, &
               dxyuxx(i,j))
            call deriv_xy(uxy, nxmx, nymx, i, j, nnodex, nnodey, dx, dy, &
               dxyuxy(i,j))
            call deriv_xy(uxz, nxmx, nymx, i, j, nnodex, nnodey, dx, dy, &
               dxyuxz(i,j))

            call deriv_xy(uyx, nxmx, nymx, i, j, nnodex, nnodey, dx, dy, &
               dxyuyx(i,j))
            call deriv_xy(uyy, nxmx, nymx, i, j, nnodex, nnodey, dx, dy, &
               dxyuyy(i,j))
            call deriv_xy(uyz, nxmx, nymx, i, j, nnodex, nnodey, dx, dy, &
               dxyuyz(i,j))

            call deriv_xy(uzx, nxmx, nymx, i, j, nnodex, nnodey, dx, dy, &
               dxyuzx(i,j))
            call deriv_xy(uzy, nxmx, nymx, i, j, nnodex, nnodey, dx, dy, &
               dxyuzy(i,j))
            call deriv_xy(uzz, nxmx, nymx, i, j, nnodex, nnodey, dx, dy, &
               dxyuzz(i,j))


         end do
      end do

      end if

      time = second1(dummy) - t1
      tmin = time / 60.
      if(myid .eq. 0)then
         write(6, 386)tmin
         write(15, 386) tmin
!       call flush
      end if
  386 format('time to call eqdsk_setup = ',f9.3, 4h min)


!       ------------------------------------
!       determine reordering and permutation
!       ------------------------------------

!      -------------------
!      add extra layer of points
!      -------------------
      nlayer = 0
      do ilayer=1,nlayer
      do j=1,nnodey
      do i=1,nnodex
          if (mask(i,j).eq.1) then
            mask2(i,j) = mask(i,j)

          if (i-1.ge.1)       mask2(i-1,j) = 1
          if (i+1.le.nnodex)  mask2(i+1,j) = 1
          if (j-1.ge.1)       mask2(i,j-1) = 1
          if (j+1.le.nnodey)  mask2(i,j+1) = 1
          endif
      enddo
      enddo
      mask(1:nnodex,1:nnodey) = mask2(1:nnodex,1:nnodey)
      enddo




        num_inside = 0
        do j=1,nnodey
        do i=1,nnodex
          if (mask(i,j).eq.1) then
             num_inside = num_inside + 1
          endif
        enddo
        enddo


      allocate( new_to_org(3*nnodex*nnodey) )

      num_new = 0
      num_org = 0
!      ------------------------------------
!      order the points in the plasma first
!      ------------------------------------
      do j=1,nnodey
      do i=1,nnodex
         if (mask(i,j).eq.1) then
                irnc = (j-1)*3 + (i-1)*nnodey*3+1
              num_org = irnc


            num_new = num_new + 1
            new_to_org( num_new ) = num_org
            num_new = num_new + 1
            new_to_org( num_new) = num_org + 1
            num_new = num_new + 1
            new_to_org( num_new) = num_org + 2
         endif
      enddo
      enddo

!      ------------------------------------
!      order the points in the plasma last
!      ------------------------------------
        do j=1,nnodey
        do i=1,nnodex
           if (mask(i,j).eq.0) then
                irnc = (j-1)*3 + (i-1)*nnodey*3+1
                num_org = irnc

                num_new = num_new + 1
                new_to_org( num_new ) = num_org
                num_new = num_new + 1
                new_to_org( num_new) = num_org + 1
                num_new = num_new + 1
                new_to_org( num_new) = num_org + 2
           endif
        enddo
        enddo





!      ----------------------------------------------
!      note we wish to be solving a smaller
!      matrix and ignoring points outside the plasma
!      ----------------------------------------------

      org_nrow = 3*nnodex*nnodey
      org_ncol = org_nrow

      new_nrow = 3*num_inside
        new_ncol = new_nrow

      nrow = new_nrow
      ncol = nrow

        if ((myrow.eq.0).and.(mycol.eq.0)) then
          write(*,*) 'new_nrow,org_nrow ', new_nrow,org_nrow
          write(15,*) 'new_nrow,org_nrow ', new_nrow,org_nrow
      endif


!     --------------------------------------------
!     perform reordering for better load balancing
!     --------------------------------------------
      if (need_shuffle) then
      if ((myrow.eq.0).and.(mycol.eq.0)) then
          nrow3 = new_nrow/3
          allocate( iperm( nrow3 ) )
          allocate( dtemp2(nrow3) )

          call random_number( dtemp2 )
          call dshell( nrow3, dtemp2, iperm )

          deallocate( dtemp2)
          allocate( itemp2(nrow3) )

          do i=1,nrow3
              itemp2(i) = new_to_org( (i-1)*3 + 1 )
          enddo
          do i=1,nrow3
              ip = (i-1)*3 + 1
              new_to_org( ip) = itemp2( iperm(i) )
              new_to_org(ip+1) = new_to_org(ip) + 1
              new_to_org(ip+2) = new_to_org(ip) + 2
          enddo

!         ---------------------------------------
!         broadcast to all processors so everyone
!         will agree with the new mapping
!         ---------------------------------------
          ilen = size(new_to_org)
          call igebs2d(icontxt,'A',' ',ilen,1,new_to_org, &
                              ilen)

          deallocate( itemp2 )
          deallocate( iperm )
       else
          ilen = size(new_to_org)
          call igebr2d(icontxt,'A',' ',ilen,1,new_to_org, &
                              ilen,0,0)

       endif
      endif


      t1=second1(dummy)


!     ----------------------------------
!     define scalapack array descriptors
!     ----------------------------------
      mb = 33 * 3
      nb = mb


!       --------------------------------------------------------
!       setting mmb = mb is convenient but may take up too much
!       memory, especially for large problems
!
!       mmb and mb should be divisible by 3, mod(mmb,3).eq.0
!       --------------------------------------------------------
        mmb = mb
        isok = (mmb.ge.3).and.(mod(mmb,3).eq.0)
        if (.not.isok) then
          mmb = 3
        endif

!       --------------------------------------------------------
!       4/23/2007:  set mmb = 3 to save memory in large problems
!       --------------------------------------------------------
      mmb = 3



!      -------------------------------
!      setup datastructures for Bmat_all(0:(nprow-1),0:(npcol-1))
!      each can be considered a trivial scalapack array of size mb by ncol
!      will all data for Bmat_all(irow,icol) residing
!      on processor (irow,icol)
!      -------------------------------

      allocate(niabegin_all(0:(nprow-1),0:(npcol-1)) )
      allocate( isize_all(0:(nprow-1),0:(npcol-1)) )
      niabegin_all(:,:) = 0
      isize_all(:,:) = 0

      allocate(Btmp(mmb,org_ncol))
      allocate(brhs(mmb,1))

      Btmp(:,:) = 0.0
      brhs(:,:) = 0.0

      allocate( descBtmp_all(DLEN_,0:(nprow-1),0:(npcol-1)) )
      allocate( descbrhs_all(DLEN_,0:(nprow-1),0:(npcol-1)) )
      descBtmp_all(:,:,:) = 0
      descbrhs_all(:,:,:) = 0

         Locp = mmb
         Locq = ncol
         lld = max(1,Locp)

      do icol=0,npcol-1
      do irow=0,nprow-1
         call descinit( descBtmp,mmb, ncol, mb, ncol, &
                        irow,icol,icontxt, lld, info)
         if (info.ne.0) then
          if ((myrow.eq.0).and.(mycol.eq.0)) then
            write(*,9710) irow,icol,info
 9710            format('descinit for descBtmp: irow,icol ', &
                       2(1x,i7),' info ',i7 )
          endif
         endif

         call descinit( descbrhs,mmb,1, mmb,1, &
                        irow,icol,icontxt, lld, info )
         descBtmp_all(:,irow,icol) = descBtmp(:)
         descbrhs_all(:,irow,icol) = descbrhs(:)
      enddo
      enddo

!     ------------------------------------------
!     define scalapack array descriptor for amat
!     ------------------------------------------


      nrow_local = max(1,numroc( nrow,mb,myrow,0,nprow ))
      ncol_local = max(1,numroc( ncol,nb,mycol,0,npcol ))
      lld = max(1,nrow_local)

      p_amat_dim = nrow_local * ncol_local + 1
      p_amat_size = p_amat_dim * 16 / 1.0e+06

      if (isolve .ne. 0) then


      if (myid.eq.0) then
         write (6, 6839) p_amat_dim, p_amat_size
         write (15,6839) p_amat_dim, p_amat_size
      end if
 6839 format('p_amat_dim = ', i10, ' words,   ', &
             'p_amat_size = ', i5, ' MBytes')


!      ----------------------------------------
!      allocate storage for scalapack matrices
!      ----------------------------------------
#ifdef USE_HPL
!      -----------------------------------
!      let HPL determine the storage needs
!      -----------------------------------
       call HPL_matinit( nrow, mb, hpl_lld, hpl_ineed )
       lld = max(lld,hpl_lld)
       p_amat_dim = max(p_amat_dim,hpl_ineed) 
#endif

      allocate(p_amat(p_amat_dim))
      p_amat(:) = cmplx(0.0,0.0)




      p_brhs_dim = nrow
      p_ipiv_dim = nrow
      allocate( p_brhs(p_brhs_dim), p_ipiv(p_ipiv_dim) )
      p_brhs(:) = cmplx(0.0,0.0)
      p_ipiv(:) = 0



!      if (nrow_local*ncol_local .gt. p_amat_dim) then
!      write(6,*) 'increase p_amat_dim to ',
!     &              nrow_local*ncol_local+1
!      write(15,*) 'increase p_amat_dim to ',
!     &              nrow_local*ncol_local+1
!        stop '** error ** '
!      endif

 
#ifdef  USE_HPL
      call descinit( desc_amat, nrow,ncol+1, mb,nb, 0,0, icontxt, &
               lld, info )
#else
      call descinit( desc_amat, nrow,ncol, mb,nb, 0,0, icontxt, &
               lld, info ) 
#endif
      if (info.ne.0) then
      write(6,*) '** descinit for amat returns info = ', info
      write(6,*) 'lld = ', lld
      write(15,*) '** descinit for amat returns info = ', info
      write(15,*) 'lld = ', lld
        stop '** error ** '
      endif


      do irnc=1,nrow_local*ncol_local
      p_amat(irnc) = cmplx(0., 0.)
      enddo

      if (lld .gt. p_brhs_dim) then
      write(6,*) 'increase p_brhs_dim to ', lld+1
      write(15,*) 'increase p_brhs_dim to ', lld+1
        stop '** error ** '
      endif


      call descinit( desc_brhs, nrow, 1, mb,nb, 0,0, icontxt, &
               lld, info )
      if (info.ne.0) then
        write(6,*) '** descinit for brhs returns info ', info
        write(15,*) '** descinit for brhs returns info ', info
        stop '** error ** '
      endif





! ------------------------------
! precompute mapping from ia->(i,j)
! precompute mapping from ja->(n,m)
! ------------------------------

       allocate( itable(3*nnodex*nnodey) )
       allocate( jtable(3*nnodex*nnodey) )
       allocate( mtable(3*nnodex*nnodey) )
       allocate( ntable(3*nnodex*nnodey) )

       do i=1,3*nnodex*nnodey
          itable(i) = undefined
          jtable(i) = undefined
          mtable(i) = undefined
          ntable(i) = undefined
       enddo

       do i=1,nnodex
       do j=1,nnodey
         irnc = (j-1)*3 + (i-1)*nnodey*3+1
         itable(irnc)   = i
       itable(irnc+1) = i
       itable(irnc+2) = i

         jtable(irnc)   = j
         jtable(irnc+1) = j
         jtable(irnc+2) = j
       enddo
       enddo

       do n=nkx1,nkx2
       do m=nky1,nky2
           icnc = (n-nkx1)*3*nnodey + (m-nky1)*3+1
           ntable(icnc)   = n
           ntable(icnc+1) = n
           ntable(icnc+2) = n

           mtable(icnc)   = m
           mtable(icnc+1) = m
           mtable(icnc+2) = m
        enddo
        enddo


!DLG: Initialise pf file

        dlgTime = second1 ( dummy )

        if ( write_f_file ) then

            if ( myId .eq. 0 ) &
                write(*,*) 'dlg: initialising full f file ...'
            what    = 10
            call blacs_get ( icontxt, what, ivalue )

            call init_pf_dlg ( iValue, myId, &
                nNodeX, nNodeY, nuPer, nuPar, &
                vc2_mks, capR, y, uPerp, uPara )

        endif

        if (myId==0) write(167,12398) 'initialise pf file:', &
           (dlgTime - second1(dummy))/60.0



!     ------------------------------------------------------------------------
!     Load x, y and z equations for spatial point (i,j) and mode number (n,m)
!     ------------------------------------------------------------------------




      dlgTime = second1 ( dummy )
    
      nnb = nprow*npcol*mmb
      nia_loop: &
      do niastart=1,desc_amat(M_),nnb
         niaend = min(desc_amat(M_), niastart+nnb-1)
         niasize = niaend-niastart+1
         if (niasize.le.0) exit


!      -------------------------------------------
!      determine what rows need to be computed for
!      (irow,icol) processor
!      -------------------------------------------
         nia = niastart
         do icol=0,(npcol-1)
         do irow=0,(nprow-1)
           niabegin_all(irow,icol) = nia
           isize_all(irow,icol) = max(0,min(mmb, niasize))

           nia = nia + isize_all(irow,icol)
           niasize = niasize - isize_all(irow,icol)
         enddo
         enddo


         ni_loop: &
         do ni=1,isize_all(myrow,mycol),3

!          -------------------------------------------
!          construct the entire row of original matrix
!          -------------------------------------------
          nia = (ni-1) + niabegin_all(myrow,mycol)
          ia = new_to_org( nia )

          ja_loop: &
          do ja=1,org_ncol,3


!  -----------------------------------------------
!  (ia,ja) are all local indices in a single block
!
!  note each ia and ja loop has increment by 3
!  to match the original code
!  -----------------------------------------------

! ------------------------------
! obtain (i,j,m,n) from (ia,ja)
! ------------------------------

                 i = itable(ia)
                 j = jtable(ia)
                 m = mtable(ja)
                 n = ntable(ja)


            isok = (1.le.i).and.(i.le.nnodex)  .and. &
                   (1.le.j).and.(j.le.nnodey)  .and. &
                   (nkx1.le.n).and.(n.le.nkx2) .and. &
                   (nky1.le.m).and.(m.le.nky2)
            if (.not.isok) then
                write(*,*) 'i,j ',i,j
                write(*,*) 'm,n ',m,n
                write(*,*) 'ia,ja ', ia,ja
            write(*,*) 'nkx1,nkx2 ', nkx1,nkx2
            write(*,*) 'nky1,nky2 ', nky1,nky2
                write(*,*) 'nnodex, nnodey ', nnodex, nnodey
                stop 'error in main '
            endif


! ------------
! double check
! ------------
            irnc = (j-1) * 3 + (i-1) * nnodey * 3 + 1
            icnc = (n - nkx1) * 3 * nnodey + (m - nky1) * 3 + 1


            isok = (irnc.eq.ia).and.(icnc.eq.ja)
            if (.not.isok) then
                write(*,*) 'main: ia,ja ', ia,ja
                write(*,*) 'i,j ',i,j
                write(*,*) 'n,m ',n,m
                write(*,*) 'main: irnc,icnc ',irnc,icnc
                stop 'error in current '
             endif


!            ------------------------------------------
!            local matrix entry is at
!            ipos = lrindx + (lcindx-1)*desc_amat(LLD_)
!            ------------------------------------------

             lrindx = (ia-iastart) + lrindx_base
             lcindx = (ja-jastart) + lcindx_base

             if (idebug.ge.2) then
             call infog2l(ia,ja, desc_amat,nprow,npcol, &
                myrow,mycol, lrindx1,lcindx1, rsrc,csrc)

             isok = (lrindx.eq.lrindx1).and.(lcindx.eq.lcindx1) &
                   .and. (rsrc.eq.myrow).and.(csrc.eq.mycol)
             if (.not.isok) then
                 write(*,*) 'lrindx,lcindx ', lrindx,lcindx
                 write(*,*) 'lrindx1,lcindx1 ', lrindx1,lcindx1
                 write(*,*) 'ia,ja ', ia,ja
                 write(*,*) 'rsrc,csrc ', rsrc,csrc
                 stop '** error ** '
             endif
             endif





!           -------------------------------------------------
!           copy variable to be compatible with previous code
!           -------------------------------------------------
            irnc1 = irnc
            icnc1 = icnc
            rsrc1 = myrow
            csrc1 = mycol

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

!                 ----------------------
!                 interior plasma region:
!                 ----------------------
                  if( mask(i,j) .eq. 1 .and. nboundary .ge. 1) then

                     if (isigma .eq. 1)then

                        call sigmad_cql3d(i, j, n, m, rho(i,j), rho_a, &
                        gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
                        xme, qe, xnea(i,j), xnuomg, &
                        xkte(i,j), omgce(i,j), omgpe2(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sigexx, sigexy, sigexz, &
                        sigeyx, sigeyy, sigeyz, &
                        sigezx, sigezy, sigezz, &
                        delta0, 0, nupar, nuper, n_psi, &
                        n_psi_dim, dfdupere, dfdupare, &
                        UminPara, UmaxPara, UPERP, UPARA, &
                        vce_mks, dfe_cql_uprp, dfe_cql_uprl, nbessj, &
                        nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max, &
                        i_sav, j_sav, upshift, damping, xk_cutoff,y(j), &
                        eNorm_factor_e)

                        call sigmad_cql3d(i, j, n, m, rho(i,j), rho_a, &
                        gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
                        xmi1, qi1, xn1a(i,j), xnuomg, &
                        xkti(i,j), omgci1(i,j), omgp12(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sig1xx, sig1xy, sig1xz, &
                        sig1yx, sig1yy, sig1yz, &
                        sig1zx, sig1zy, sig1zz, &
                        delta0, ndisti1, nupar, nuper, n_psi, &
                        n_psi_dim, dfduper1, dfdupar1, &
                        UminPara, UmaxPara, UPERP, UPARA, &
                        vc1_mks, df1_cql_uprp, df1_cql_uprl, nbessj, &
                        nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max, &
                        i_sav1,j_sav1,upshift,damping,xk_cutoff,y(j), &
                        eNorm_factor_i1)


                        if ( write_f_file ) then

                        if(eta2 .ne. 0.0) &
                        call sigmad_cql3d(i, j, n, m, rho(i,j), rho_a, &
                        gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
                        xmi2, qi2, xn2a(i,j), xnuomg, &
                        xkti2(i,j), omgci2(i,j), omgp22(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sig2xx, sig2xy, sig2xz, &
                        sig2yx, sig2yy, sig2yz, &
                        sig2zx, sig2zy, sig2zz, &
                        delta0, ndisti2, nupar, nuper, n_psi, &
                        n_psi_dim, dfduper2, dfdupar2, &
                        UminPara, UmaxPara, UPERP, UPARA, &
                        vc2_mks, df2_cql_uprp, df2_cql_uprl, nbessj, &
                        nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max, &
                        i_sav2,j_sav2,upshift,damping,xk_cutoff,y(j), &
                        eNorm_factor_i2, 1, iValue)

                        else

                        if(eta2 .ne. 0.0) &
                        call sigmad_cql3d(i, j, n, m, rho(i,j), rho_a, &
                        gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
                        xmi2, qi2, xn2a(i,j), xnuomg, &
                        xkti2(i,j), omgci2(i,j), omgp22(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sig2xx, sig2xy, sig2xz, &
                        sig2yx, sig2yy, sig2yz, &
                        sig2zx, sig2zy, sig2zz, &
                        delta0, ndisti2, nupar, nuper, n_psi, &
                        n_psi_dim, dfduper2, dfdupar2, &
                        UminPara, UmaxPara, UPERP, UPARA, &
                        vc2_mks, df2_cql_uprp, df2_cql_uprl, nbessj, &
                        nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max, &
                        i_sav2,j_sav2,upshift,damping,xk_cutoff,y(j), &
                        eNorm_factor_i2 )

                        endif


                        if(eta3 .ne. 0.0) &
                        call sigmad_cql3d(i, j, n, m, rho(i,j), rho_a, &
                        gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
                        xmi3, qi3, xn3a(i,j), xnuomg, &
                        xkti3(i,j), omgci3(i,j), omgp32(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sig3xx, sig3xy, sig3xz, &
                        sig3yx, sig3yy, sig3yz, &
                        sig3zx, sig3zy, sig3zz, &
                        delta0, ndisti3, nupar, nuper, n_psi, &
                        n_psi_dim, dfduper3, dfdupar3, &
                        UminPara, UmaxPara, UPERP, UPARA, &
                        vc3_mks, df3_cql_uprp, df3_cql_uprl, nbessj, &
                        nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max, &
                        i_sav,j_sav,upshift,damping,xk_cutoff,y(j), &
                        eNorm_factor_i3)

                        if(eta4 .ne. 0.0) &
                        call sigmad_cql3d(i, j, n, m, rho(i,j), rho_a, &
                        gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
                        xmi4, qi4, xn4a(i,j), xnuomg, &
                        xkti4(i,j), omgci4(i,j), omgp42(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sig4xx, sig4xy, sig4xz, &
                        sig4yx, sig4yy, sig4yz, &
                        sig4zx, sig4zy, sig4zz, &
                        delta0, ndisti4, nupar, nuper, n_psi, &
                        n_psi_dim, dfduper4, dfdupar4, &
                        UminPara, UmaxPara, UPERP, UPARA, &
                        vc4_mks, df4_cql_uprp, df4_cql_uprl, nbessj, &
                        nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max, &
                        i_sav,j_sav,upshift,damping,xk_cutoff,y(j), &
                        eNorm_factor_i4)

                        if(eta5 .ne. 0.0) &
                        call sigmad_cql3d(i, j, n, m, rho(i,j), rho_a, &
                        gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
                        xmi5, qi5, xn5a(i,j), xnuomg, &
                        xkti5(i,j), omgci5(i,j), omgp52(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sig5xx, sig5xy, sig5xz, &
                        sig5yx, sig5yy, sig5yz, &
                        sig5zx, sig5zy, sig5zz, &
                        delta0, ndisti5, nupar, nuper, n_psi, &
                        n_psi_dim, dfduper5, dfdupar5, &
                        UminPara, UmaxPara, UPERP, UPARA, &
                        vc5_mks, df5_cql_uprp, df5_cql_uprl, nbessj, &
                        nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max, &
                        i_sav,j_sav,upshift,damping,xk_cutoff,y(j), &
                        eNorm_factor_i5)

                        if(eta6 .ne. 0.0) &
                        call sigmad_cql3d(i, j, n, m, rho(i,j), rho_a, &
                        gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
                        xmi6, qi6, xn6a(i,j), xnuomg, &
                        xkti6(i,j), omgci6(i,j), omgp62(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sig6xx, sig6xy, sig6xz, &
                        sig6yx, sig6yy, sig6yz, &
                        sig6zx, sig6zy, sig6zz, &
                        delta0, ndisti6, nupar, nuper, n_psi, &
                        n_psi_dim, dfduper6, dfdupar6, &
                        UminPara, UmaxPara, UPERP, UPARA, &
                        vc6_mks, df6_cql_uprp, df6_cql_uprl, nbessj, &
                        nkperp, zi, eps0, v0i, omgrf, xk0, kperp_max, &
                        i_sav,j_sav,upshift,damping,xk_cutoff,y(j), &
                        eNorm_factor_i6)

                        if(eta_slo .ne. 0.0) &
                        call sigmah_slow(i, j, n, m, &
                        xmi_slo, qi_slo, xna_slo(i,j), xnuomg, &
                        eslow, omgci_slo(i,j), omgp2_slo(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sigsloxx, sigsloxy, sigsloxz, &
                        sigsloyx, sigsloyy, sigsloyz, &
                        sigslozx, sigslozy, sigslozz, &
                        xkte(i,j), zeffcd, zi, eps0, v0i, omgrf, xk0, &
                        kperp_max, i_sav, j_sav)


                     end if





                     if (isigma .eq. 0)then ! Cold plasma

                        call sigmac_stix(i, j, n, m, &
                        xme, qe, xnea(i,j), xnuomg, &
                        xkte(i,j), omgce(i,j), omgpe2(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sigexx, sigexy, sigexz, &
                        sigeyx, sigeyy, sigeyz, &
                        sigezx, sigezy, sigezz, &
                        delta0, zi, eps0, v0i, omgrf, xk0, kperp_max, &
                        i_sav, j_sav)

                        call sigmac_stix(i, j, n, m, &
                        xmi1, qi1, xn1a(i,j), xnuomg, &
                        xkti(i,j), omgci1(i,j), omgp12(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sig1xx, sig1xy, sig1xz, &
                        sig1yx, sig1yy, sig1yz, &
                        sig1zx, sig1zy, sig1zz, &
                        delta0, zi, eps0, v0i, omgrf, xk0, kperp_max, &
                        i_sav, j_sav)

                        if(eta2 .ne. 0.0) &
                        call sigmac_stix(i, j, n, m, &
                        xmi2, qi2, xn2a(i,j), xnuomg, &
                        xkti2(i,j), omgci2(i,j), omgp22(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sig2xx, sig2xy, sig2xz, &
                        sig2yx, sig2yy, sig2yz, &
                        sig2zx, sig2zy, sig2zz, &
                        delta0, zi, eps0, v0i, omgrf, xk0, kperp_max, &
                        i_sav, j_sav)

                        if(eta3 .ne. 0.0) &
                        call sigmac_stix(i, j, n, m, &
                        xmi3, qi3, xn3a(i,j), xnuomg, &
                        xkti3(i,j), omgci3(i,j), omgp32(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sig3xx, sig3xy, sig3xz, &
                        sig3yx, sig3yy, sig3yz, &
                        sig3zx, sig3zy, sig3zz, &
                        delta0, zi, eps0, v0i, omgrf, xk0, kperp_max, &
                        i_sav, j_sav)


                        if(eta4 .ne. 0.0) &
                        call sigmac_stix(i, j, n, m, &
                        xmi4, qi4, xn4a(i,j), xnuomg, &
                        xkti4(i,j), omgci4(i,j), omgp42(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sig4xx, sig4xy, sig4xz, &
                        sig4yx, sig4yy, sig4yz, &
                        sig4zx, sig4zy, sig4zz, &
                        delta0, zi, eps0, v0i, omgrf, xk0, kperp_max, &
                        i_sav, j_sav)

                        if(eta5 .ne. 0.0) &
                        call sigmac_stix(i, j, n, m, &
                        xmi5, qi5, xn5a(i,j), xnuomg, &
                        xkti5(i,j), omgci5(i,j), omgp52(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sig5xx, sig5xy, sig5xz, &
                        sig5yx, sig5yy, sig5yz, &
                        sig5zx, sig5zy, sig5zz, &
                        delta0, zi, eps0, v0i, omgrf, xk0, kperp_max, &
                        i_sav, j_sav)

                        if(eta6 .ne. 0.0) &
                        call sigmac_stix(i, j, n, m, &
                        xmi6, qi6, xn6a(i,j), xnuomg, &
                        xkti6(i,j), omgci6(i,j), omgp62(i,j), &
                        -lmax, lmax, nzfun, ibessel, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sig6xx, sig6xy, sig6xz, &
                        sig6yx, sig6yy, sig6yz, &
                        sig6zx, sig6zy, sig6zz, &
                        delta0, zi, eps0, v0i, omgrf, xk0, kperp_max, &
                        i_sav, j_sav)

                     end if


                  sigxx = sigexx + sig1xx + sig2xx + sig3xx &
                                 + sig4xx + sig5xx + sig6xx + sigsloxx
                  sigxy = sigexy + sig1xy + sig2xy + sig3xy &
                                 + sig4xy + sig5xy + sig6xy + sigsloxy
                  sigxz = sigexz + sig1xz + sig2xz + sig3xz &
                                 + sig4xz + sig5xz + sig6xz + sigsloxz

                  sigyx = sigeyx + sig1yx + sig2yx + sig3yx &
                                 + sig4yx + sig5yx + sig6yx + sigsloyx
                  sigyy = sigeyy + sig1yy + sig2yy + sig3yy &
                                 + sig4yy + sig5yy + sig6yy + sigsloyy
                  sigyz = sigeyz + sig1yz + sig2yz + sig3yz &
                                 + sig4yz + sig5yz + sig6yz + sigsloyz

                  sigzx = sigezx + sig1zx + sig2zx + sig3zx &
                                 + sig4zx + sig5zx + sig6zx + sigslozx
                  sigzy = sigezy + sig1zy + sig2zy + sig3zy &
                                 + sig4zy + sig5zy + sig6zy + sigslozy
                  sigzz = sigezz + sig1zz + sig2zz + sig3zz &
                                 + sig4zz + sig5zz + sig6zz + sigslozz


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
                     rnphi = xkphi(i) / xk0

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

                  end if

!                 --------------------------
!                 metal boundary edge region:
!                 --------------------------
                  if( mask(i,j) .eq. 0 .and. nboundary .ge. 1) then
                     fdk = cexpkxky
                     fek = 0.0
                     ffk = 0.0

                     fgk = 0.0
                     fak = cexpkxky
                     fpk = 0.0

                     frk = 0.0
                     fqk = 0.0
                     fsk = cexpkxky

                     xb(i,j) = 0.0
                     xc(i,j) = 0.0
                     xd(i,j) = 0.0
                  end if




                  if(i .eq. idiag .and. j .eq. jdiag)then
                     fdksav(n, m) = fdk
                     feksav(n, m) = fek
                     ffksav(n, m) = ffk

                     fgksav(n, m) = fgk
                     faksav(n, m) = fak
                     fpksav(n, m) = fpk

                     frksav(n, m) = frk
                     fqksav(n, m) = fqk
                     fsksav(n, m) = fsk
                  end if

!                 ------------------------------------------
!                 Load the collocation solution matrix, amat
!                 ------------------------------------------
!efd                  ipos = lrindx + (lcindx - 1) * desc_amat(LLD_)
!efd                  p_amat(ipos)     = fdk
!efd                  p_amat(ipos + 1) = fgk
!efd                  p_amat(ipos + 2) = frk
!efd
!efd                  if (icnc1 .eq. 1) then
!efd                     p_brhs(ipos)     = xb(i,j)
!efd                     p_brhs(ipos + 1) = xc(i,j)
!efd                     p_brhs(ipos + 2) = xd(i,j)
!efd                  end if
!efd
!efd                  ipos = lrindx + (lcindx + 1 - 1) * desc_amat(LLD_)
!efd                  p_amat(ipos)     = fek
!efd                  p_amat(ipos + 1) = fak
!efd                  p_amat(ipos + 2) = fqk
!efd
!efd                  ipos = lrindx + (lcindx + 2 - 1) * desc_amat(LLD_)
!efd                  p_amat(ipos)     = ffk
!efd                  p_amat(ipos + 1) = fpk
!efd                  p_amat(ipos + 2) = fsk

!       -------------------------------
!       copy the value into local array
!       -------------------------------
                        Btmp(ni,  ja) = fdk
                        Btmp(ni+1,ja) = fgk
                        Btmp(ni+2,ja) = frk

                        Btmp(ni,  ja+1) = fek
                        Btmp(ni+1,ja+1) = fak
                        Btmp(ni+2,ja+1) = fqk

                        Btmp(ni,  ja+2) = ffk
                        Btmp(ni+1,ja+2) = fpk
                        Btmp(ni+2,ja+2) = fsk



                        if (ja.eq.1) then
                                brhs(ni,1)   = xb(i,j)
                                brhs(ni+1,1) = xc(i,j)
                                brhs(ni+2,1) = xd(i,j)
                        endif

                   enddo ja_loop 
              enddo ni_loop


!      -------------------------------------
!      perform transformation to real space
!      -------------------------------------

      do ni=1,isize_all(myrow,mycol)
        row(1:org_ncol) = Btmp(ni,1:org_ncol)

        call convert2d_row(row, rowk, &
         xmax, ymax, nnodex, nnodey, &
         nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2, &
         nkx1, nkx2, nky1, nky2, xx, yy, dx, dy, ndfmax, &
         nmodesx, nmodesy,use_fft)



!      -----------------------------------------
!      select only subset of variables in plasma
!      -----------------------------------------
      do nja=1,ncol
        ja = new_to_org(nja)
          Btmp(ni,nja) = rowk(ja)
      enddo

      enddo



!      -------------------------
!      ready to copy into p_amat
!      -------------------------

      do icol=0,(npcol-1)
      do irow=0,(nprow-1)
        mm = isize_all(irow,icol)
        if (mm.le.0) cycle


        nn = ncol
         nia = niabegin_all(irow,icol)
         nja = 1
         descBtmp(:) = descBtmp_all(:,irow,icol)
         call pzgecopy( mm,nn, &
                  Btmp, 1,1, descBtmp, &
                p_amat, nia,nja,desc_amat )

         descbrhs(:) = descbrhs_all(:,irow,icol)
         call pzgecopy( mm,1, &
                brhs, 1,1, descbrhs, &
                  p_brhs, nia,1, desc_brhs )
      enddo
      enddo

      enddo nia_loop 

      if (myId==0) write(167,12398) 'nia_loop:', &
           (dlgTime - second1(dummy))/60.0

      time = second1(dummy) - t1
      tmin = time / 60.
      if (myid.eq.0) then
         write(6 ,835) tmin
         write(15,835) tmin
!       call flush
      endif

  835 format('time to load matrix =',f9.3,4h min)

      call blacs_barrier(icontxt, 'All')


      do n = nkx1, nkx2
         do m = nky1, nky2
            fdksav2d(n - nkx1 + 1, m - nky1 + 1) = fdksav(n,m)
            feksav2d(n - nkx1 + 1, m - nky1 + 1) = feksav(n,m)
            ffksav2d(n - nkx1 + 1, m - nky1 + 1) = ffksav(n,m)

            fgksav2d(n - nkx1 + 1, m - nky1 + 1) = fgksav(n,m)
            faksav2d(n - nkx1 + 1, m - nky1 + 1) = faksav(n,m)
            fpksav2d(n - nkx1 + 1, m - nky1 + 1) = fpksav(n,m)

            frksav2d(n - nkx1 + 1, m - nky1 + 1) = frksav(n,m)
            fqksav2d(n - nkx1 + 1, m - nky1 + 1) = fqksav(n,m)
            fsksav2d(n - nkx1 + 1, m - nky1 + 1) = fsksav(n,m)
         end do
      end do

      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fdksav2d, &
         nxmx, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, feksav2d, &
         nxmx, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, ffksav2d, &
         nxmx, -1, -1)

      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fgksav2d, &
         nxmx, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, faksav2d, &
         nxmx, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fpksav2d, &
         nxmx, -1, -1)

      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, frksav2d, &
         nxmx, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fqksav2d, &
         nxmx, -1, -1)
      call zgsum2d(icontxt, 'All', ' ', nnodex, nnodey, fsksav2d, &
         nxmx, -1, -1)


      do n = nkx1, nkx2
         do m = nky1, nky2
            fdksav(n,m) = fdksav2d(n - nkx1 + 1, m - nky1 + 1)
            feksav(n,m) = feksav2d(n - nkx1 + 1, m - nky1 + 1)
            ffksav(n,m) = ffksav2d(n - nkx1 + 1, m - nky1 + 1)

            fgksav(n,m) = fgksav2d(n - nkx1 + 1, m - nky1 + 1)
            faksav(n,m) = faksav2d(n - nkx1 + 1, m - nky1 + 1)
            fpksav(n,m) = fpksav2d(n - nkx1 + 1, m - nky1 + 1)

            frksav(n,m) = frksav2d(n - nkx1 + 1, m - nky1 + 1)
            fqksav(n,m) = fqksav2d(n - nkx1 + 1, m - nky1 + 1)
            fsksav(n,m) = fsksav2d(n - nkx1 + 1, m - nky1 + 1)
         end do
      end do






!      --------------------------
!      scale each row by its norm
!      --------------------------
      global_col = 1
      incX = desc_amat(M_)

      do global_row = 1, nrow
        call pdznrm2( ncol, norm2, p_amat, global_row, global_col, &
             desc_amat, incX )

        alpha = 1.0
        if (norm2 .ne. 0.0) alpha = 1.0/norm2

        call pzscal ( ncol, alpha, p_amat, global_row, global_col, &
             desc_amat, incX )

        call pzscal ( 1, alpha, p_brhs, global_row, 1, &
             desc_brhs, desc_brhs(M_))


      enddo




      dlgTime = second1 ( dummy )

      call blacs_barrier(icontxt, 'All')
      t1 = second1(dummy) 
#ifdef  USE_HPL
!      ------------------------------------------------
!      HPL uses A(:,ncol+1) as storage for rhs vector and solution
!      so we need extra vector copies before and after to p_brhs
!      ------------------------------------------------
       call pzcopy(nrow,p_brhs,1,1,desc_brhs,1, &
                        p_amat,1,ncol+1,desc_amat,1)

       call HPL_pzgesv( nrow, p_amat, desc_amat, info)

       call pzcopy(nrow,p_amat,1,ncol+1,desc_amat,1, &
                        p_brhs,1,1,desc_brhs,1)

       if (info.ne.0) then
         write(6,*) 'hpl_pzgesv returns info = ',info
         write(15,*) 'hpl_pzgesv returns info = ',info
         stop '** error **'
       endif
#else

!      ---------------------
!      Call scalapack solver
!      ---------------------
      call pzgetrf( nrow,ncol, p_amat, 1, 1, desc_amat, p_ipiv, info )
        if (info.ne.0) then
         write(6,*) 'pzgetrf returns info = ', info
         write(15,*) 'pzgetrf returns info = ', info
           stop '** error ** '
      endif

      call blacs_barrier(icontxt, 'All')

      call pzgetrs( 'notrans', nrow, 1, p_amat, 1, 1, desc_amat, &
                      p_ipiv, p_brhs, 1,1,desc_brhs, info )


      if (info.ne.0) then
         write(6,*) 'pzgetrs returns info = ', info
         write(15,*) 'pzgetrs returns info = ', info
         stop '** error ** '
      endif 
#endif
        call blacs_barrier(icontxt, 'All')
      time = second1(dummy) - t1
        time = max(1.0, time)

      if (myId==0) write(167,12398) 'field solve:', &
           (dlgTime - second1(dummy))/60.0


      deallocate(p_amat)
      deallocate(Btmp)

!--   Operations for complex matrix factor and solve:
      ops = 8. / 3. * (real(nrow))**3 + 7. * (real(nrow))**2

      gflops  = ops / time / 1.0e+09
      gflopsp = gflops / nproc
      tmin = time / 60.
      if (myid.eq.0) then
        write(6 ,839) nrow, nproc
        write(15,839) nrow, nproc
        write(6 ,833) tmin
        write(15,833) tmin
        write(6 ,837) gflops
        write(15,837) gflops
        write(6 ,838) gflopsp
        write(15,838) gflopsp
      endif

  833 format('time taken by ScaLAPACK =',f9.3,4h min)
  834 format('info =', i10)
  839 format('nrow =', i10, 5x, 'nproc =', i10)
  837 format('operations per sec by pzgetrf & pzgetrs =' &
         ,f11.3,11h Gflops/sec)
  838 format('operations per sec per processor =' &
         ,f9.3,21h Gflops/sec/processor)

      call blacs_barrier(icontxt, 'All')


!     ------------------
!     broadcast solution
!     ------------------

      brhs2(:) = 0.0
      brhsk(:) = 0.0
      brhs_tmp(:) = 0.0


      icnc = 1
      do irnc=1,nrow
         call infog2l( irnc, icnc, desc_brhs, nprow,npcol, &
                          myrow,mycol, lrindx, lcindx, rsrc,csrc )
         ismine = (rsrc .eq. myrow) .and. (csrc .eq. mycol)
         if (ismine) then
             ipos = lrindx + (lcindx-1)*desc_amat(LLD_)
             brhs_tmp(irnc) =  p_brhs(ipos)
         endif
      enddo
      call zgsum2d(icontxt, 'All', ' ', nrow, 1, brhs_tmp, nrow, -1, -1)


!     --------------------------------
!     expand solution to original size
!     --------------------------------
      do nia=1,nrow
      ia = new_to_org(nia)
        brhs2(ia) = brhs_tmp(nia)
      enddo


!     ---------------------------------------
!     adjust problem size
!     ---------------------------------------
      nrow = org_nrow
      ncol = nrow


      do i = 1, nnodex
         do j = 1, nnodey

            irnc = (i - 1) * 3 * nnodey + (j - 1) * 3 + 1
            ealpha(i, j) = brhs2(irnc)

            irnc = (i - 1) * 3 * nnodey + (j - 1) * 3 + 2
            ebeta(i, j) = brhs2(irnc)

            irnc = (i - 1) * 3 * nnodey + (j - 1) * 3 + 3
            eb(i, j) = brhs2(irnc)

         end do
      end do


      end if

      call blacs_barrier(icontxt, 'All')

!     --------------------------------
!     write out solution in real space
!     --------------------------------

      if (myid .eq.0) then
         write(53, 1314) nnodex, nnodey
       write(53, 310) (capr(i), i = 1, nnodex)
       write(53, 310) (y(j), j = 1, nnodey)
         write(53, 8310)((ealpha(i, j), i = 1, nnodex), j = 1, nnodey)
         write(53, 8310)((ebeta(i, j), i = 1, nnodex), j = 1, nnodey)
         write(53, 8310)((eb(i, j), i = 1, nnodex), j = 1, nnodey)
         close(53)
      end if




!      call checkx(myid, idiag, jdiag, nkx1, nkx2, nky1, nky2,
!     1   fdksav, feksav, ffksav, ealpha, ebeta, eb,
!     1   xb, xc, xd, nkdim1, nkdim2, mkdim1, mkdim2,
!     1   nnodex, nnodey, nxmx, nymx, xprime, yprime,
!     1   idiag, jdiag, 1)

!      call checky(myid, idiag, jdiag, nkx1, nkx2, nky1, nky2,
!     1   fgksav, faksav, fpksav, ealpha, ebeta, eb,
!     1   xb, xc, xd, nkdim1, nkdim2, mkdim1, mkdim2,
!     1   nnodex, nnodey, nxmx, nymx, xprime, yprime,
!     1   idiag, jdiag, 1)

!      call checkz(myid, idiag, jdiag, nkx1, nkx2, nky1, nky2,
!     1   frksav, fqksav, fsksav, ealpha, ebeta, eb,
!     1   xb, xc, xd, nkdim1, nkdim2, mkdim1, mkdim2,
!     1   nnodex, nnodey, nxmx, nymx, xprime, yprime,
!     1   idiag, jdiag, 1)



      t1 = second1(dummy)

!     ------------------------------------------
!     Fourier transform the Stix electric fields:
!     ------------------------------------------

      call sft2d_parallel(ealpha, xmax, ymax, nnodex, nnodey, &
         nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2, &
         nkx1, nkx2, nky1, nky2, xx_inv, yy_inv, ealphak, dx, dy, &
         myid, nproc, icontxt)

      call sft2d_parallel(ebeta, xmax, ymax, nnodex, nnodey, &
         nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2, &
         nkx1, nkx2, nky1, nky2, xx_inv, yy_inv, ebetak, dx, dy, &
         myid, nproc, icontxt)

      call sft2d_parallel(eb, xmax, ymax, nnodex, nnodey, &
         nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2, &
         nkx1, nkx2, nky1, nky2, xx_inv, yy_inv, ebk, dx, dy, &
         myid, nproc, icontxt)



!      -----------------------------------
!      write out solution in Fourier modes
!      -----------------------------------

!      if (myid .eq.0) then
!         write(54, 309) nkx1, nkx2, nky1, nky2
!         write(54, 8310)((ealphak(n, m), n = nkx1, nkx2), m = nky1,nky2)
!         write(54, 8310)((ebetak(n, m), n = nkx1, nkx2), m = nky1, nky2)
!         write(54, 8310)((ebk(n, m), n = nkx1, nkx2), m = nky1, nky2)
!         write(54, 8310)(xkxsav(n), n = nkx1, nkx2)
!         write(54, 8310)(xkysav(m), m = nky1, nky2)

!         close(54)

!      end if



      do n = nkx1, nkx2
         do m = nky1, nky2

            ealphakmod(n, m) = sqrt(conjg(ealphak(n, m))* ealphak(n, m))
            ebetakmod(n, m) = sqrt(conjg(ebetak(n, m)) * ebetak(n, m))
            ebkmod(n, m) = sqrt(conjg(ebk(n, m)) * ebk(n, m))

         end do
      end do



      time = second1(dummy) - t1
      tmin = time / 60.
      if (myid.eq.0) then
         write(6,  1837) tmin
          write(15, 1837) tmin
      endif

 1837 format('time to Fourier transform the fields =',f9.3,4h min)

      if (myid .eq. 0) then
         write(15, 910)
         write(15, 930)
         write(15, 920)
         write(15, 940)
         write(15, 920)

         write(6, 910)
         write(6, 930)
         write(6, 920)
         write(6, 940)
         write(6, 920)


         i = idiag
         do j = 1, nnodey
            write(6,  900)i, j, capr(i), y(j), &
                                       ealpha(i,j), ebeta(i,j), eb(i,j)
            write(15, 900)i, j, capr(i), y(j), &
                                       ealpha(i,j), ebeta(i,j), eb(i,j)
         end do

         write(15, 910)
         write(15, 930)
         write(15, 920)
         write(15, 940)
         write(15, 920)

         write(6, 910)
         write(6, 930)
         write(6, 920)
         write(6, 940)
         write(6, 920)

         j = jdiag
         do i = 1, nnodex
            write(6,  900)i, j, capr(i), y(j), &
                                       ealpha(i,j), ebeta(i,j), eb(i,j)
            write(15, 900)i, j, capr(i), y(j), &
                                       ealpha(i,j), ebeta(i,j), eb(i,j)
         end do
      endif

  900 format(i5, i5, 1p9e12.3)
  940 format(1h , 4h   i, 5h    j, &
                          5x, 7h  R(i) , 5x, 7h  Z(j) , &
                          3x, 7hre ealp, 5x, 7him ealp, &
                          5x, 7hre ebet, 5x, 7him ebet, &
                          5x, 7h re eb , 5x, 7h im eb )
  910 format (1h1)
  920 format (1h0)
  930 format(3x, 'rf electric field in Stix frame')


!     ----------------------------------------------
!     Calculate E in the Lab frame and eplus, eminus
!     ----------------------------------------------
      isq2 = SQRT(0.5)
      do i = 1, nnodex
         do j = 1, nnodey

            ex(i,j)   = uxx(i,j) * ealpha(i,j) &
                      + uyx(i,j) * ebeta(i,j) &
                      + uzx(i,j) * eb(i,j)

            ey(i,j)   = uxy(i,j) * ealpha(i,j) &
                      + uyy(i,j) * ebeta(i,j) &
                      + uzy(i,j) * eb(i,j)

            ez(i,j)   = uxz(i,j) * ealpha(i,j) &
                      + uyz(i,j) * ebeta(i,j) &
                      + uzz(i,j) * eb(i,j)

            eplus(i,j)  = isq2 * (ealpha(i,j) + zi * ebeta(i,j))
          eminus(i,j) = isq2 * (ealpha(i,j) - zi * ebeta(i,j))

         end do
      end do

!-------------------------------------------------------
!     Interpolate Eplus, Eminus and mod B onto (rho, theta) grid
!-------------------------------------------------------
      do n = 1, nnoderho
         do m = 1, mnodetheta

          xgiv = capr_flux(n, m) - rt - xwleft
          ygiv = capz_flux(n, m) - ybottom

          call intplt(xgiv, ygiv, fout, nnodex, nnodey, eplus, &
               nxmx, nymx, dx, dy)
            eplus_flux(n, m) = fout

          call intplt(xgiv, ygiv, fout, nnodex, nnodey, eminus, &
               nxmx, nymx, dx, dy)
            eminus_flux(n, m) = fout

          call intplt(xgiv, ygiv, fout, nnodex, nnodey, xkperp_cold, &
               nxmx, nymx, dx, dy)
            xkperp_flux(n, m) = fout

          call intplt_re(xgiv, ygiv, fout, nnodex, nnodey, bmod, &
               nxmx, nymx, dx, dy)
            bmod_flux(n, m) = fout



!          if(myid .eq. 0)write(6,*)n, m, eplus_flux(n,m)

         end do
      end do

!---------------------------------------------------
!     Interpolate back onto (R, Z) grid for plotting
!---------------------------------------------------
      do i = 1, nnodex
         do j = 1, nnodey

          rho_giv = rho(i, j)
          the_giv = theta(i,j)

          call intplt(rho_giv, the_giv, fout, nnoderho, mnodetheta, &
               eplus_flux, nrhomax, nthetamax, drho, dtheta)
            eplus_flux_plot(i, j) = fout

          call intplt(rho_giv, the_giv, fout, nnoderho, mnodetheta, &
               eminus_flux, nrhomax, nthetamax, drho, dtheta)
            eminus_flux_plot(i, j) = fout

          call intplt(rho_giv, the_giv, fout, nnoderho, mnodetheta, &
               xkperp_flux, nrhomax, nthetamax, drho, dtheta)
            xkperp_flux_plot(i, j) = fout

!          if(myid .eq. 0)write(6,*)i, j, eplus_flux_plot(i, j)

         end do
      end do


      if (myid .eq. 0) then
         write(15, 910)
         write(15, 1930)
         write(15, 920)
         write(15, 1940)
         write(15, 920)

         write(6, 910)
         write(6, 1930)
         write(6, 920)
         write(6, 1940)
         write(6, 920)

         i = idiag
         do j = 1, nnodey
            write(6,  900)i, j, capr(i), y(j), ex(i,j), ey(i,j), ez(i,j)
            write(15, 900)i, j, capr(i), y(j), ex(i,j), ey(i,j), ez(i,j)
         end do

         write(15, 910)
         write(15, 1930)
         write(15, 920)
         write(15, 1940)
         write(15, 920)

         write(6, 910)
         write(6, 1930)
         write(6, 920)
         write(6, 1940)
         write(6, 920)

         j = jdiag
         do i = 1, nnodex
            write(6,  900)i, j, capr(i), y(j), ex(i,j), ey(i,j), ez(i,j)
            write(15, 900)i, j, capr(i), y(j), ex(i,j), ey(i,j), ez(i,j)
         end do
      endif


 1940 format(1h , 4h   i, 5h    j, &
                          5x, 7h  R(i) , 5x, 7h  Z(j) , &
                          4x, 6h re ex, 6x, 6h im ex, &
                          6x, 6h re ey, 6x, 6h im ey, &
                          6x, 6h re ez, 6x, 6h im ez)

 1930 format(3x, 'rf electric field in Lab frame')

        
!     -----------------------------------
!     Calculate minimum and maximum xkprl
!     -----------------------------------

      dlgTime = second1 ( dummy )

      xkprl_min = 10000.0
      xkprl_max = -10000.0
      xkprl_pos_min = 10000.0
      xkprl_neg_max = -10000.0

      do i = 1, nnodex
         do j = 1, nnodey

            if( mask(i,j) .eq. 1 .and. nboundary .ge. 1)then

               do n = nkx1, nkx2
                  do m = nky1, nky2

                 xkalp = uxx(i, j) * xkxsav(n) &
                           + uxy(i, j) * xkysav(m) &
                           + uxz(i, j) * xkphi(i)
                     xkbet = uyx(i, j) * xkxsav(n) &
                           + uyy(i, j) * xkysav(m) &
                           + uyz(i, j) * xkphi(i)
                     xkprl = uzx(i, j) * xkxsav(n) &
                           + uzy(i, j) * xkysav(m) &
                           + uzz(i, j) * xkphi(i)
                     xkperp = sqrt(xkalp**2 + xkbet**2)

!                    ------------------------------------
!                    Optional: leave out upshift in xkprl
!                    --------------------------------- --
                     if (upshift .eq. 0) xkprl = uzz(i,j) * xkphi(i)

                     if (upshift .eq. -1)then
                    if (xkperp  .gt. xk_cutoff) &
                                         xkprl = uzz(i,j) * xkphi(i)
                     end if

                     if (xkprl  .eq. 0.0) xkprl  = 1.0e-08
                 if (xkperp .eq. 0.0) xkperp = 1.0e-08

                 sgn_kprl = sign(1.0, xkprl)
                     akprl = abs(xkprl)


                 gammab = abs(omgci2(i,j)/ (2.0 * vthi20 * xkprl**2) &
                                           * gradprlb(i,j) / bmod(i,j))

!                    ------------------------------------------------
!                    Calculate Brambillas xkrpl_eff using l = 1 only
!                    ------------------------------------------------
                     y0 = 1.5
                     yb = y0

                     if(sgn_kprl .ge. 0.0)then
                        fgam = 1.0
                        if(gammab .gt. 1.0e-05)then
                           yb = y0
                           fgam = (sqrt(1. +  4. * gammab * yb) - 1.) &
                              / (2. * gammab * yb)
                        endif
                        xkprl_eff = xkprl / fgam
                     end if

                     if(sgn_kprl .lt. 0.0)then
                        fgam = 1.0
                        if(gammab .gt. 1.0e-05)then
                           descrim = 1. - 4. * gammab * y0
                           if (descrim .ge. 0.0) yb =   y0
                           if (descrim .lt. 0.0) yb = - y0
                           fgam = (1. - sqrt(1. -  4. * gammab * yb) ) &
                                / (2. * gammab * yb)
                        endif
                        xkprl_eff = xkprl / fgam
                     end if


                 xkprl = xkprl_eff


                     if(xkprl .gt. xkprl_max) xkprl_max = xkprl
                 if(xkprl .lt. xkprl_min) xkprl_min = xkprl

                 if(xkprl .gt. 0.0)then
                   if(xkprl .lt. xkprl_pos_min)xkprl_pos_min = xkprl
                 end if

                 if(xkprl .lt. 0.0)then
                   if(xkprl .gt. xkprl_neg_max)xkprl_neg_max = xkprl
                 end if

                  end do
             end do

          end if

         end do
      end do

      if (myId==0) write(167,12398) 'calculate min and max kprl:', &
           (dlgTime - second1(dummy))/60.0


      if (myid .eq. 0) then
             write(15, *) "xkprl_max = ", xkprl_max
         write(15, *) "xkprl_pos_min = ", xkprl_pos_min
       write(15, *) "xkprl_neg_max = ", xkprl_neg_max
       write(15, *) "xkprl_min = ", xkprl_min

             write(6, *) "xkprl_max = ", xkprl_max
         write(6, *) "xkprl_pos_min = ", xkprl_pos_min
       write(6, *) "xkprl_neg_max = ", xkprl_neg_max
       write(6, *) "xkprl_min = ", xkprl_min

      end if

        if(myid.eq.0)write(*,*) 'dlg: vce_mks post sigma: ', vce_mks
        if(myid.eq.0)write(*,*) 'dlg: vc1_mks post sigma: ', vc1_mks
        if(myid.eq.0)write(*,*) 'dlg: vc2_mks post sigma: ', vc2_mks



!     --------------------------------------
!     calculate individual species currents:
!     --------------------------------------

      dlgTime = second1 ( dummy )

      t1 = second1(dummy)


      call current_cd(xjpxe, xjpye, xjpze, &
         nxmx, nymx, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xme, qe, xnea, xkte, omgce, omgpe2, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         ealphak, ebetak, ebk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, 0, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfdupere, dfdupare, &
         UminPara, UmaxPara, UPERP, UPARA, &
         vce_mks, dfe_cql_uprp, dfe_cql_uprl, rho, rho_a, nbessj, &
         zeffcd, clight, x, y, rt, b0, ftrap, omgrf, &
         xjpx_ehst, xjpy_ehst, xjpz_ehst, nkperp, zi, eps0, &
         v0i, xk0, kperp_max, i_sav, j_sav, upshift, damping, xk_cutoff, &
         y, eNorm_factor_e, mask )




      call current_1(xjpx1, xjpy1, xjpz1, nxmx, nymx, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xmi1, qi1, xn1a, xkti, omgci1, omgp12, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         ealphak, ebetak, ebk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndisti1, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper1, dfdupar1, &
         UminPara, UmaxPara, UPERP, UPARA, &
         vc1_mks, df1_cql_uprp, df1_cql_uprl, rho, rho_a, nbessj,nkperp, &
         zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav1, j_sav1, upshift, &
         damping, xk_cutoff,y, eNorm_factor_i1, mask)



      if(eta2 .ne. 0.0) then

         call current_2(xjpx2, xjpy2, xjpz2, nxmx, nymx, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xmi2, qi2, xn2a, xkti2, omgci2, omgp22, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         ealphak, ebetak, ebk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndisti2, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper2, dfdupar2, &
         UminPara, UmaxPara, UPERP, UPARA, &
         vc2_mks, df2_cql_uprp, df2_cql_uprl, rho, rho_a, nbessj,nkperp, &
         zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav2, j_sav2, upshift, &
         damping, xk_cutoff,y, eNorm_factor_i2, mask)

      endif

      if(eta3 .ne. 0.0) then

         call current(xjpx3, xjpy3, xjpz3, nxmx, nymx, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xmi3, qi3, xn3a, xkti3, omgci3, omgp32, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         ealphak, ebetak, ebk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndisti3, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper3, dfdupar3, &
         UminPara, UmaxPara, UPERP, UPARA, &
         vc3_mks, df3_cql_uprp, df3_cql_uprl, rho, rho_a, nbessj,nkperp, &
         zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift, &
         damping, xk_cutoff,y, eNorm_factor_i3, mask)

      endif

      if(eta4 .ne. 0.0) then

         call current(xjpx4, xjpy4, xjpz4, nxmx, nymx, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xmi4, qi4, xn4a, xkti4, omgci4, omgp42, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         ealphak, ebetak, ebk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndisti4, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper4, dfdupar4, &
         UminPara, UmaxPara, UPERP, UPARA, &
         vc4_mks, df4_cql_uprp, df4_cql_uprl, rho, rho_a, nbessj,nkperp, &
         zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift, &
         damping, xk_cutoff,y, eNorm_factor_i4, mask)

      endif

      if(eta5 .ne. 0.0) then

         call current(xjpx5, xjpy5, xjpz5, nxmx, nymx, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xmi5, qi5, xn5a, xkti5, omgci5, omgp52, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         ealphak, ebetak, ebk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndisti5, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper5, dfdupar5, &
         UminPara, UmaxPara, UPERP, UPARA, &
         vc5_mks, df5_cql_uprp, df5_cql_uprl, rho, rho_a, nbessj,nkperp, &
         zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift, &
         damping, xk_cutoff,y, eNorm_factor_i5, mask)

      endif

      if(eta6 .ne. 0.0) then

         call current(xjpx6, xjpy6, xjpz6, nxmx, nymx, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xmi6, qi6, xn6a, xkti6, omgci6, omgp62, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         ealphak, ebetak, ebk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndisti6, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper6, dfdupar6, &
         UminPara, UmaxPara, UPERP, UPARA, &
         vc6_mks, df6_cql_uprp, df6_cql_uprl, rho, rho_a, nbessj,nkperp, &
         zi, eps0, v0i, omgrf, xk0, kperp_max, i_sav, j_sav, upshift, &
         damping, xk_cutoff,y, eNorm_factor_i6, mask)

      endif


      if(eta_slo .ne. 0.0) &
         call cur_slo(xj_slox, xj_sloy, xj_sloz, nxmx, nymx, &
         nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xmi_slo, qi_slo, xna_slo, eslow, omgci_slo, omgp2_slo, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         ealphak, ebetak, ebk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         xkte, zeffcd, myid, nproc)



      time = second1(dummy) - t1
      tmin = time / 60.
      if (myid .eq. 0) then
         write(6 ,836) tmin
       write(15,836) tmin
      endif

      if (myId==0) write(167,12398) 'currents:', &
           (dlgTime - second1(dummy))/60.0


  836 format('time to calculate currents =',f9.3,4h min)

        if(myid.eq.0)write(*,*) 'dlg: vce_mks post current: ', vce_mks
        if(myid.eq.0)write(*,*) 'dlg: vc1_mks post current: ', vc1_mks
        if(myid.eq.0)write(*,*) 'dlg: vc2_mks post current: ', vc2_mks


!     --------------------------------------------
!     calculate total plasma current in Stix frame:
!     --------------------------------------------
      do i = 1, nnodex
         do j = 1, nnodey
            xjpx(i,j) = xjpxe(i,j) + xjpx1(i,j) + xjpx2(i,j)+ xjpx3(i,j) &
                      + xjpx4(i,j) + xjpx5(i,j) + xjpx6(i,j) &
                      + xj_slox(i, j)
            xjpy(i,j) = xjpye(i,j) + xjpy1(i,j) + xjpy2(i,j)+ xjpy3(i,j) &
                      + xjpy4(i,j) + xjpy5(i,j) + xjpy6(i,j) &
                      + xj_sloy(i, j)
            xjpz(i,j) = xjpze(i,j) + xjpz1(i,j) + xjpz2(i,j)+ xjpz3(i,j) &
                      + xjpz4(i,j) + xjpz5(i,j) + xjpz6(i,j) &
                      + xj_sloz(i, j)

         end do
      end do

!     -------------------------------------------
!     Calculate total plasma current in Lab frame
!     -------------------------------------------
      do i = 1, nnodex
         do j = 1, nnodey

            xjpx_lab(i,j)   = uxx(i,j) * xjpx(i,j) &
                            + uyx(i,j) * xjpy(i,j) &
                            + uzx(i,j) * xjpz(i,j)

            xjpy_lab(i,j)   = uxy(i,j) * xjpx(i,j) &
                            + uyy(i,j) * xjpy(i,j) &
                            + uzy(i,j) * xjpz(i,j)

            xjpz_lab(i,j)   = uxz(i,j) * xjpx(i,j) &
                            + uyz(i,j) * xjpy(i,j) &
                            + uzz(i,j) * xjpz(i,j)

            xjpxe_lab(i,j)   = uxx(i,j) * xjpxe(i,j) + uyx(i,j) * xjpye(i,j) + uzx(i,j) * xjpze(i,j) 
            xjpye_lab(i,j)   = uxy(i,j) * xjpxe(i,j) + uyy(i,j) * xjpye(i,j) + uzy(i,j) * xjpze(i,j)
            xjpze_lab(i,j)   = uxz(i,j) * xjpxe(i,j) + uyz(i,j) * xjpye(i,j) + uzz(i,j) * xjpze(i,j)

         end do
      end do

!     ----------------------------------------------------------------
!     Calculate fluctuating charge density (rho_plasma) for electrons:
!     ----------------------------------------------------------------
      do i = 2, nnodex - 1
         do j = 2, nnodey - 1

            djxdx = (xjpxe_lab(i+1, j) - xjpxe_lab(i-1, j)) / (2.0 * dx)
            djydy = (xjpye_lab(i, j+1) - xjpye_lab(i, j-1)) / (2.0 * dy)

            rho_pla(i,j) = (djxdx + djydy + xjpxe_lab(i,j) / capr(i) &
                          + zi * nphi / capr(i) * xjpze_lab(i,j) ) &
                       / (zi * omgrf)

!           ntilda_e(i,j) = rho_pla(i,j) / qe

         end do
      end do

!     ---------------------------------------------------
!     Calculate charge density (rho_plasma) in the plasma:
!     ---------------------------------------------------
      do i = 2, nnodex - 1
         do j = 2, nnodey - 1

            djxdx = (xjpx_lab(i+1, j) - xjpx_lab(i-1, j)) / (2.0 * dx)
            djydy = (xjpy_lab(i, j+1) - xjpy_lab(i, j-1)) / (2.0 * dy)

            rho_pla(i,j) = (djxdx + djydy + xjpx_lab(i,j) / capr(i) &
                           + zi * nphi / capr(i) * xjpz_lab(i,j) ) &
                        / (zi * omgrf)

         end do
      end do

!     -------------------------------------------------
!     Calculate charge density (rho_ant) on the antenna:
!     -------------------------------------------------
      do i = 2, nnodex - 1
         do j = 2, nnodey - 1

            djxdx = (xjx(i+1, j) - xjx(i-1, j)) / (2.0 * dx)
            djydy = (xjy(i, j+1) - xjy(i, j-1)) / (2.0 * dy)

            rho_ant(i,j) = (djxdx + djydy + xjx(i,j)/capr(i) &
                           + zi * nphi / capr(i) * xjz(i,j) ) &
                        / (zi * omgrf)

         end do
      end do


      t1 = second1(dummy)

!     --------------------------------
!     Distribution function from CQL3D
!     --------------------------------

!     -----------------
!     allocate arrays
!     -----------------

      allocate( bqlavg_e(nuper, nupar, nnoderho) )
      allocate( cqlavg_e(nuper, nupar, nnoderho) )
      allocate( eqlavg_e(nuper, nupar, nnoderho) )
      allocate( fqlavg_e(nuper, nupar, nnoderho) )
      bqlavg_e = 0.0
      cqlavg_e = 0.0
      eqlavg_e = 0.0
      fqlavg_e = 0.0


      allocate( bqlavg_i1(nuper, nupar, nnoderho) )
      allocate( cqlavg_i1(nuper, nupar, nnoderho) )
      allocate( eqlavg_i1(nuper, nupar, nnoderho) )
      allocate( fqlavg_i1(nuper, nupar, nnoderho) )
      bqlavg_i1 = 0.0
      cqlavg_i1 = 0.0
      eqlavg_i1 = 0.0
      fqlavg_i1 = 0.0

      allocate( bqlavg_work(nuper, nupar, nnoderho) )
      allocate( cqlavg_work(nuper, nupar, nnoderho) )
      allocate( eqlavg_work(nuper, nupar, nnoderho) )
      allocate( fqlavg_work(nuper, nupar, nnoderho) )
      bqlavg_work = 0.0
      cqlavg_work = 0.0
      eqlavg_work = 0.0
      fqlavg_work = 0.0


      if(eta2 .ne. 0.0)then
         allocate( bqlavg_i2(nuper, nupar, nnoderho) )
         allocate( cqlavg_i2(nuper, nupar, nnoderho) )
         allocate( eqlavg_i2(nuper, nupar, nnoderho) )
         allocate( fqlavg_i2(nuper, nupar, nnoderho) )
         bqlavg_i2 = 0.0
         cqlavg_i2 = 0.0
         eqlavg_i2 = 0.0
         fqlavg_i2 = 0.0
      end if

      if(eta3 .ne. 0.0)then
         allocate( bqlavg_i3(nuper, nupar, nnoderho) )
         allocate( cqlavg_i3(nuper, nupar, nnoderho) )
         allocate( eqlavg_i3(nuper, nupar, nnoderho) )
         allocate( fqlavg_i3(nuper, nupar, nnoderho) )
         bqlavg_i3 = 0.0
         cqlavg_i3 = 0.0
         eqlavg_i3 = 0.0
         fqlavg_i3 = 0.0
      end if

      if(eta4 .ne. 0.0)then
         allocate( bqlavg_i4(nuper, nupar, nnoderho) )
         allocate( cqlavg_i4(nuper, nupar, nnoderho) )
         allocate( eqlavg_i4(nuper, nupar, nnoderho) )
         allocate( fqlavg_i4(nuper, nupar, nnoderho) )
         bqlavg_i4 = 0.0
         cqlavg_i4 = 0.0
         eqlavg_i4 = 0.0
         fqlavg_i4 = 0.0
      end if

      if(eta5 .ne. 0.0)then
         allocate( bqlavg_i5(nuper, nupar, nnoderho) )
         allocate( cqlavg_i5(nuper, nupar, nnoderho) )
         allocate( eqlavg_i5(nuper, nupar, nnoderho) )
         allocate( fqlavg_i5(nuper, nupar, nnoderho) )
         bqlavg_i5 = 0.0
         cqlavg_i5 = 0.0
         eqlavg_i5 = 0.0
         fqlavg_i5 = 0.0
      end if

      if(eta6 .ne. 0.0)then
         allocate( bqlavg_i6(nuper, nupar, nnoderho) )
         allocate( cqlavg_i6(nuper, nupar, nnoderho) )
         allocate( eqlavg_i6(nuper, nupar, nnoderho) )
         allocate( fqlavg_i6(nuper, nupar, nnoderho) )
         bqlavg_i6 = 0.0
         cqlavg_i6 = 0.0
         eqlavg_i6 = 0.0
         fqlavg_i6 = 0.0
      end if



    dlgTime = second1 ( dummy )

      if (nzeta_wdot .gt. 0) then
!DLG:   Setup particle tracing


!        if( eNorm_factor_e>0)
!     .  vce_mks = sqrt(2d0*1d3*eNorm_factor_e*e/xme)
!        if( eNorm_factor_i1>0)
!     .  vc1_mks = sqrt(2d0*1d3*eNorm_factor_i1*e/xmi1)
!        if( eNorm_factor_i2>0)
!     .  vc2_mks = sqrt(2d0*1d3*eNorm_factor_i2*e/xmi2)
!        if( eNorm_factor_i3>0)
!     .  vc3_mks = sqrt(2d0*1d3*eNorm_factor_i3*e/xmi3)
!        if( eNorm_factor_i4>0)
!     .  vc4_mks = sqrt(2d0*1d3*eNorm_factor_i4*e/xmi4)
!        if( eNorm_factor_i5>0)
!     .  vc5_mks = sqrt(2d0*1d3*eNorm_factor_i5*e/xmi5)
!        if( eNorm_factor_i6>0)
!     .  vc6_mks = sqrt(2d0*1d3*eNorm_factor_i6*e/xmi6)


!      if (myid .eq. 1) write(*,*)
!     . 'dlg: Initialising particle tracing variables'
      !eqdsk_fileName  = 'g122080.03100'


      ! ------------------ !
      ! --  Electrons   -- !
      ! ------------------ !

      call ql_myra_write(bqlavg_e, cqlavg_e, eqlavg_e, fqlavg_e, &
         wdote, fx0e, fy0e, fz0e, darea, nrhomax, nnoderho, drho, &
         dx, dy, rt, nxmx, nymx, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xme, qe, xnea, xkte, omgce, omgpe2, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         ealphak, ebetak, ebk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, 0, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfdupere, dfdupare, &
         UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work, &
         vce_mks, dfe_cql_uprp, dfe_cql_uprl, rho, rho_a, nbessj, &
         zeffcd, clight, x, y, rt, b0, ftrap, omgrf, &
         nkperp, lmaxdim, nzeta_wdot, theta_, &
         n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg, &
         n_bin, upshift, 0, xk_cutoff, y, eNorm_factor_e, mask )



      ! ------------------ !
      ! --Majority ions -- !
      ! ------------------ !

      call ql_myra_write(bqlavg_i1, cqlavg_i1, eqlavg_i1, fqlavg_i1, &
         wdoti1, fx0i1, fy0i1, fz0i1, darea, nrhomax, nnoderho, drho, &
         dx, dy, rt, nxmx, nymx, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xmi1, qi1, xn1a, xkti, omgci1, omgp12, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         ealphak, ebetak, ebk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndisti1, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper1, dfdupar1, &
         UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work, &
         vc1_mks, df1_cql_uprp, df1_cql_uprl, rho, rho_a, nbessj, &
         zeffcd, clight, x, y, rt, b0, ftrap, omgrf, &
         nkperp, lmaxdim, nzeta_wdot, theta_, &
         n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg, &
         n_bin, upshift, 0, xk_cutoff, y, eNorm_factor_i1, mask)

      ! ------------------ !
      ! --Minority ions -- !
      ! ------------------ !

      if(eta2 .ne. 0.0)then

      call ql_myra_write(bqlavg_i2, cqlavg_i2, eqlavg_i2, fqlavg_i2, &
         wdoti2, fx0i2, fy0i2, fz0i2, darea, nrhomax, nnoderho, drho, &
         dx, dy, rt, nxmx, nymx, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xmi2, qi2, xn2a, xkti2, omgci2, omgp22, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         ealphak, ebetak, ebk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndisti2, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper2, dfdupar2, &
         UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work, &
         vc2_mks, df2_cql_uprp, df2_cql_uprl, rho, rho_a, nbessj, &
         zeffcd, clight, x, y, rt, b0, ftrap, omgrf, &
         nkperp, lmaxdim, nzeta_wdot, theta_, &
         n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg, &
         n_bin, upshift, i_write, xk_cutoff, y, eNorm_factor_i2, mask)
      end if

      ! --------------------- !
      ! --Third ion species-- !
      ! --------------------- !

      if(eta3 .ne. 0.0)then
      call ql_myra_write(bqlavg_i3, cqlavg_i3, eqlavg_i3, fqlavg_i3, &
         wdoti3, fx0i3, fy0i3, fz0i3, darea, nrhomax, nnoderho, drho, &
         dx, dy, rt, nxmx, nymx, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xmi3, qi3, xn3a, xkti3, omgci3, omgp32, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         ealphak, ebetak, ebk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndisti3, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper3, dfdupar3, &
         UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work, &
         vc3_mks, df3_cql_uprp, df3_cql_uprl, rho, rho_a, nbessj, &
         zeffcd, clight, x, y, rt, b0, ftrap, omgrf, &
         nkperp, lmaxdim, nzeta_wdot, theta_, &
         n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg, &
         n_bin, upshift, 0, xk_cutoff, y, eNorm_factor_i3, mask)
      end if




      ! ---------------------- !
      ! --Fourth ion species-- !
      ! ---------------------- !

      if(eta4 .ne. 0.0)then
      call ql_myra_write(bqlavg_i4, cqlavg_i4, eqlavg_i4, fqlavg_i4, &
         wdoti4, fx0i4, fy0i4, fz0i4, darea, nrhomax, nnoderho, drho, &
         dx, dy, rt, nxmx, nymx, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xmi4, qi4, xn4a, xkti4, omgci4, omgp42, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         ealphak, ebetak, ebk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndisti4, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper4, dfdupar4, &
         UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work, &
         vc4_mks, df4_cql_uprp, df4_cql_uprl, rho, rho_a, nbessj, &
         zeffcd, clight, x, y, rt, b0, ftrap, omgrf, &
         nkperp, lmaxdim, nzeta_wdot, theta_, &
         n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg, &
         n_bin, upshift, 0, xk_cutoff, y, eNorm_factor_i4, mask)
      end if


      ! --------------------- !
      ! --Fifth ion species-- !
      ! --------------------- !

      if(eta5 .ne. 0.0)then
      call ql_myra_write(bqlavg_i5, cqlavg_i5, eqlavg_i5, fqlavg_i5, &
         wdoti5, fx0i5, fy0i5, fz0i5, darea, nrhomax, nnoderho, drho, &
         dx, dy, rt, nxmx, nymx, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xmi5, qi5, xn5a, xkti5, omgci5, omgp52, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         ealphak, ebetak, ebk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndisti5, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper5, dfdupar5, &
         UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work, &
         vc5_mks, df5_cql_uprp, df5_cql_uprl, rho, rho_a, nbessj, &
         zeffcd, clight, x, y, rt, b0, ftrap, omgrf, &
         nkperp, lmaxdim, nzeta_wdot, theta_, &
         n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg, &
         n_bin, upshift, 0, xk_cutoff, y, eNorm_factor_i5, mask)
      end if


      ! --------------------- !
      ! --Sixth ion species-- !
      ! --------------------- !

      if(eta6 .ne. 0.0)then
      call ql_myra_write(bqlavg_i6, cqlavg_i6, eqlavg_i6, fqlavg_i6, &
         wdoti6, fx0i6, fy0i6, fz0i6, darea, nrhomax, nnoderho, drho, &
         dx, dy, rt, nxmx, nymx, nnodex, nnodey, &
         nkx1, nkx2, nky1, nky2, nkdim1, nkdim2, mkdim1, mkdim2, &
         xmi6, qi6, xn6a, xkti6, omgci6, omgp62, lmax, &
         xkxsav, xkysav, nzfun, ibessel, &
         ealphak, ebetak, ebk, nphi, capr, &
         bxn, byn, bzn, &
         uxx, uxy, uxz, &
         uyx, uyy, uyz, &
         uzx, uzy, uzz, &
         myrow, mycol, nprow, npcol, icontxt, desc_amat, dlen_, &
         xx, yy, isigma, xnuomg, psi, psilim, nboundary, &
         myid, nproc, delta0, gradprlb, bmod, ndisti6, bmod_mid, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper6, dfdupar6, &
         UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work, &
         vc6_mks, df6_cql_uprp, df6_cql_uprl, rho, rho_a, nbessj, &
         zeffcd, clight, x, y, rt, b0, ftrap, omgrf, &
         nkperp, lmaxdim, nzeta_wdot, theta_, &
         n_theta_max, n_psi_max, i_psi_eq, n_theta_, dldbavg, &
         n_bin, upshift, 0, xk_cutoff, y, eNorm_factor_i6, mask)
      end if


      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(6 , 2838) tmin
         write(15, 2838) tmin
      endif

 2838 format('time to calculate quasilinear operator =', f9.3, ' min')

        if(myid.eq.0)write(*,*) 'dlg: vce_mks post qlWrite: ', &
         vce_mks
        if(myid.eq.0)write(*,*) 'dlg: vc1_mks post qlWrite: ', &
        vc1_mks
        if(myid.eq.0)write(*,*) 'dlg: vc2_mks post qlWrite: ', &
        vc2_mks

        if (myId==0) write(167,12398) 'ql:', &
           (dlgTime - second1(dummy))/60.0


!        if( eNorm_factor_e>0)
!     .  vce_mks = sqrt(2d0*1d3*eNorm_factor_e*e/xme)
!        if( eNorm_factor_i1>0)
!     .  vc1_mks = sqrt(2d0*1d3*eNorm_factor_i1*e/xmi1)
!        if( eNorm_factor_i2>0)
!     .  vc2_mks = sqrt(2d0*1d3*eNorm_factor_i2*e/xmi2)
!        if( eNorm_factor_i3>0)
!     .  vc3_mks = sqrt(2d0*1d3*eNorm_factor_i3*e/xmi3)
!        if( eNorm_factor_i4>0)
!     .  vc4_mks = sqrt(2d0*1d3*eNorm_factor_i4*e/xmi4)
!        if( eNorm_factor_i5>0)
!     .  vc5_mks = sqrt(2d0*1d3*eNorm_factor_i5*e/xmi5)
!        if( eNorm_factor_i6>0)
!     .  vc6_mks = sqrt(2d0*1d3*eNorm_factor_i6*e/xmi6)


      t1 = second1(dummy)

      ! ------------------ !
      ! --  Electrons   -- !
      ! ------------------ !
      dlgTime = second1 ( dummy )

      call wdot_qlcheck(wdote_ql, &
         nnoderho, nrhomax, &
         bqlavg_e, cqlavg_e, xme, omgrf, xkteavg, 0, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfdupere, dfdupare, &
         UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work, &
         vce_mks, dfe_cql_uprp, dfe_cql_uprl, rhon, rho_a, myid, &
         dldbavg, eNorm_factor_e)



      ! ------------------ !
      ! --Majority ions -- !
      ! ------------------ !

      call wdot_qlcheck(wdoti1_ql, &
         nnoderho, nrhomax, &
         bqlavg_i1, cqlavg_i1, xmi1, omgrf, xktiavg, ndisti1, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper1, dfdupar1, &
         UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work, &
         vc1_mks, df1_cql_uprp, df1_cql_uprl, rhon, rho_a, myid, &
         dldbavg, eNorm_factor_i1)



      ! ------------------ !
      ! --Minority ions -- !
      ! ------------------ !

      if(eta2 .ne. 0.0)then

      call wdot_qlcheck(wdoti2_ql, &
         nnoderho, nrhomax, &
         bqlavg_i2, cqlavg_i2, xmi2, omgrf, xkti2avg, ndisti2, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper2, dfdupar2, &
         UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work, &
         vc2_mks, df2_cql_uprp, df2_cql_uprl, rhon, rho_a, myid, &
         dldbavg, eNorm_factor_i2)

      end if



      ! --------------------- !
      ! --Third ion species-- !
      ! --------------------- !

      if(eta3 .ne. 0.0)then

      call wdot_qlcheck(wdoti3_ql, &
         nnoderho, nrhomax, &
         bqlavg_i3, cqlavg_i3, xmi3, omgrf, xkti3avg, ndisti3, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper3, dfdupar3, &
         UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work, &
         vc3_mks, df3_cql_uprp, df3_cql_uprl, rhon, rho_a, myid, &
         dldbavg, eNorm_factor_i3)

      end if

      ! --------------------- !
      ! --Fourth ion species-- !
      ! --------------------- !

      if(eta4 .ne. 0.0)then

      call wdot_qlcheck(wdoti4_ql, &
         nnoderho, nrhomax, &
         bqlavg_i4, cqlavg_i4, xmi4, omgrf, xkti4avg, ndisti4, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper4, dfdupar4, &
         UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work, &
         vc4_mks, df4_cql_uprp, df4_cql_uprl, rhon, rho_a, myid, &
         dldbavg, eNorm_factor_i4)

      end if

      ! --------------------- !
      ! --Fifth ion species-- !
      ! --------------------- !

      if(eta5 .ne. 0.0)then

      call wdot_qlcheck(wdoti5_ql, &
         nnoderho, nrhomax, &
         bqlavg_i5, cqlavg_i5, xmi5, omgrf, xkti5avg, ndisti5, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper5, dfdupar5, &
         UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work, &
         vc5_mks, df5_cql_uprp, df5_cql_uprl, rhon, rho_a, myid, &
         dldbavg, eNorm_factor_i5)

      end if

      ! --------------------- !
      ! --Sixth ion species-- !
      ! --------------------- !

      if(eta6 .ne. 0.0)then

      call wdot_qlcheck(wdoti6_ql, &
         nnoderho, nrhomax, &
         bqlavg_i6, cqlavg_i6, xmi6, omgrf, xkti6avg, ndisti6, &
         nupar, nuper, n_psi, &
         n_psi_dim, dfduper6, dfdupar6, &
         UminPara, UmaxPara, UPERP, UPARA, UPERP_work, UPARA_work, &
         vc6_mks, df6_cql_uprp, df6_cql_uprl, rhon, rho_a, myid, &
         dldbavg, eNorm_factor_i6)

      end if
      if(myid.eq.0) write(*,*) 'dlg: vc_mks2', vc2_mks

      end if


      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(6 , 2839) tmin
         write(15, 2839) tmin
      endif

        if(myid.eq.0)write(*,*) 'dlg: vce_mks post qlCheck: ', &
         vce_mks
        if(myid.eq.0)write(*,*) 'dlg: vc1_mks post qlCheck: ', &
        vc1_mks
        if(myid.eq.0)write(*,*) 'dlg: vc2_mks post qlCheck: ', &
        vc2_mks


 2839 format('time to check the quasilinear operator =', f9.3, ' min')

    if (myId==0) write(167,12398) 'ql check:', &
           (dlgTime - second1(dummy))/60.0


      t1 = second1(dummy)

!     -------------------------------------------------------
!     Calculate E* dot J, Poynting flux (sp), and zero order
!     poloidal force (fz0) in real space:
!     -------------------------------------------------------

      do i = 1, nnodex
         do j = 1, nnodey

          fx0(i,j) = fx0e(i,j) + fx0i1(i,j) + fx0i2(i,j) + fx0i3(i,j) &
                                 + fx0i4(i,j) + fx0i5(i,j) + fx0i6(i,j)
          fy0(i,j) = fy0e(i,j) + fy0i1(i,j) + fy0i2(i,j) + fy0i3(i,j) &
                                 + fy0i4(i,j) + fy0i5(i,j) + fy0i6(i,j)
          fz0(i,j) = fz0e(i,j) + fz0i1(i,j) + fz0i2(i,j) + fz0i3(i,j) &
                                 + fz0i4(i,j) + fz0i5(i,j) + fz0i6(i,j)


          bdotf(i,j) = bmod(i,j) * (bxn(i,j) * fx0(i,j) &
                                   +  byn(i,j) * fy0(i,j) &
                                   +  bzn(i,j) * fz0(i,j))

          bmod2(i,j) = bmod(i,j)**2

          capr_fzeta(i,j) = capr(i) * fz0(i,j)
          capr2(i,j) = capr(i)**2

          wdot(i,j) = wdote(i,j) &
                      + wdoti1(i,j) + wdoti2(i,j) + wdoti3(i,j) &
                      + wdoti4(i,j) + wdoti5(i,j) + wdoti6(i,j)

          redotj(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpx(i, j) &
                + conjg(ebeta(i, j)) * xjpy(i, j) &
                + conjg(eb(i, j)) * xjpz(i, j) )

          rekxdotj(i, j) = 0.5 / omgrf * &
                real(conjg(ealphakx(i, j)) * xjpx(i, j) &
                + conjg(ebetakx(i, j)) * xjpy(i, j) &
                + conjg(ebkx(i, j)) * xjpz(i, j) )

          rekydotj(i, j) = 0.5 / omgrf * &
                real(conjg(ealphaky(i, j)) * xjpx(i,j) &
                + conjg(ebetaky(i, j)) * xjpy(i, j) &
                + conjg(ebky(i, j)) * xjpz(i, j) )

          rekzdotj(i, j) = xkphi(i) / omgrf * redotj(i, j)

          pc(i, j) = -0.5 * (conjg(ex(i, j)) * xjx(i, j) &
                             + conjg(ey(i, j)) * xjy(i, j) &
                             + conjg(ez(i, j)) * xjz(i, j) )

          pcre(i, j) = real(pc(i, j))
          pcim(i, j) = aimag(pc(i, j))


          bxwave(i, j) = 0.0
          bzwave(i, j) = 0.0
          bywave(i, j) = 0.0

!          rhohatx(i,j) = 0.0
!            rhohaty(i,j) = 0.0

            if (i .ne. 1 .and. i .ne. nnodex) then
               if (j .ne. 1 .and. j .ne. nnodey) then

!                  drhodx = (rho(i+1, j) - rho(i-1, j)) / (2.0 * dx)
!                   drhody = (rho(i, j+1) - rho(i, j-1)) / (2.0 * dy)

!                   gradrho = sqrt(drhodx**2 + drhody**2)
!                if (gradrho .eq. 0.0) gradrho = 1.0e-08

!                   rhohatx(i,j) = drhodx / gradrho
!                   rhohaty(i,j) = drhody / gradrho

              dezdx = (ez(i+1, j) - ez(i-1, j)) / (2.0 * dx)
              dezdy = (ez(i, j+1) - ez(i, j-1)) / (2.0 * dy)
              deydx = (ey(i+1, j) - ey(i-1, j)) / (2.0 * dx)
              dexdy = (ex(i, j+1) - ex(i, j-1)) / (2.0 * dy)

              bxwave(i, j) = 1. / (zi * omgrf) &
                 * (zi * xkphi(i) * ey(i,j) - dezdy)
              bzwave(i, j) = 1. / (zi * omgrf) &
                 * (dexdy - deydx)
              bywave(i, j) = 1. / (zi * omgrf) * (ez(i,j) / capr(i) &
                  + dezdx - zi * xkphi(i) * ex(i,j))

              spx(i, j) = 1. / (2. * xmu0) * real( &
                          conjg(ez(i,j)) * bywave(i,j) &
                        - conjg(ey(i,j)) * bzwave(i,j)  )
              spz(i, j) = 1. / (2. * xmu0) * real( &
                          conjg(ey(i,j)) * bxwave(i,j) &
                        - conjg(ex(i,j)) * bywave(i,j)  )
              spy(i, j) = 1. / (2. * xmu0) * real( &
                          conjg(ex(i,j)) * bzwave(i,j) &
                        - conjg(ez(i,j)) * bxwave(i,j)  )

!                  dspxdx(i) = 1. / (2. * xmu0 * omgrf) * imag(
!    1                conjg(deydx(i)) * (deydx(i) - zi * xky * ex(i))
!    1              + conjg(ey(i)) * (d2eydx(i) - zi * xky * dexdx(i))
!    1              + conjg(dezdx(i)) * (dezdx(i) - zi * xkz * ex(i))
!    1              + conjg(ez(i)) * (d2ezdx(i) - zi * xkz * dexdx(i)))

!                  divsp(i) = dspxdx(i)
!     1               + zi * xky * spy(i) + zi * xkz * spz(i)
               end if
            end if

!     --------------------------------
!     Calculate radial and theta force:
!     --------------------------------

          fpsi0(i, j)  = fx0i2(i, j) * rhohatx(i, j) &
                         + fy0i2(i, j) * rhohaty(i, j)

            ftheta0(i,j) = fy0i2(i, j) * rhohatx(i, j) &
                         - fx0i2(i, j) * rhohaty(i, j)


            if ( eta2 /= 0.0 ) & 
            fpsi1(i,j) = 1.0 / (qi2 * xn2a(i, j) ) &
                       * fpsi0(i, j) / (capr(i) * bpol(i, j))

            ftheta1(i,j) = ftheta0(i,j) * bpol(i,j)

         end do
      end do

      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(6 , 2840) tmin
         write(15, 2840) tmin
      endif

 2840 format('time to calculate E* dot J and Poynting flux =', f9.3, &
       ' min')


      t1 = second1(dummy)


!     -------------------------
!     Integrate antenna current
!     -------------------------
      call sgrator(x, y, pcre, 1, nnodex, 1, nnodey, pcrto2, &
         capr, rt,  nxmx, nymx)
      call sgrator(x, y, pcim, 1, nnodex, 1, nnodey, pcito2, &
         capr, rt, nxmx, nymx)


!     ---------------------------
!     Integrate 1/2 Re (E* dot J):
!     ---------------------------
      call sgrator(x, y, redotj, 1, nnodex, 1, nnodey, ptot, &
         capr, rt,  nxmx, nymx)

!     ---------------
!     Integrate Wdot:
!     ---------------
      call sgrator(x, y, wdot, 1, nnodex, 1, nnodey, ptot_wdot, &
         capr, rt,  nxmx, nymx)


!     ------------------
!     Normalize to jdote:
!     ------------------
      powtot = 2.0 * pi * rt * ptot

!     -----------------
!     Normalize to wdot:
!     -----------------
!      powtot = 2.0 * pi * rt * ptot_wdot

      if (powtot .ne. 0.0) pscale = prfin / powtot
      if (powtot .eq. 0.0) pscale = 1.0

      if (nphi_number .gt. 1) prfin = 0.0

      if (prfin .lt. 1e-5) pscale = 1.0
      if (pscale .le. 0.0) pscale = 1.0

      if (myid .eq. 0) then
         write(6, 1218) pscale
         write(15, 1218) pscale
      endif

!DLG:   Append ql file with pScale

    dlgTime = second1 ( dummy )

        if ( myId .eq. 0 .and. i_write .ne. 0 ) then

            write (*,*) 'APPENDING output/p_ql.nc ...'
            ncFileName = 'output/p_ql.nc'

            call check ( &
       nf90_open ( ncFileName, nf90_write, nc_id ) )
            call check ( nf90_inq_varId ( nc_id, "pscale", pScale_id ))
            call check ( &
       nf90_put_var ( nc_id, pScale_id, pscale ) )
            call check ( nf90_close ( nc_id ) )

        endif

    if (myId==0) write(167,12398) 'append p_ql.nc:', &
           (dlgTime - second1(dummy))/60.0


!     -----------------------
!     Scale powers to Prf,in:
!     -----------------------
      powtot = powtot * pscale
      ptot = ptot * pscale
      pcrto2 = pcrto2 * pscale
      pcito2 = pcito2 * pscale

      fxtot_ant = fxtot_ant * pscale
      fytot_ant = fytot_ant * pscale
      fztot_ant = fztot_ant * pscale

      fxtot_plasma = fxtot_plasma * pscale
      fytot_plasma = fytot_plasma * pscale
      fztot_plasma = fztot_plasma * pscale


!     ------------------------------------
!     Scale fields and current to Prf, in:
!     ------------------------------------

      do n = nkx1, nkx2
         do m = nky1, nky2
            ealphakmod(n, m) = ealphakmod(n, m) * sqrt(pscale)
            ebetakmod(n, m) = ebetakmod(n, m) * sqrt(pscale)
            ebkmod(n, m) = ebkmod(n, m) * sqrt(pscale)
         end do
      end do

      do n = 1, nnoderho
         do m = 1, mnodetheta
            eplus_flux(n, m) = eplus_flux(n, m) * sqrt(pscale)
          eminus_flux(n, m) = eminus_flux(n, m) * sqrt(pscale)
       end do
      end do

      do i = 1, nnodex
         do j = 1, nnodey

            eplus_flux_plot(i,j) = eplus_flux_plot(i,j) * sqrt(pscale)
            eminus_flux_plot(i,j) = eminus_flux_plot(i,j) * sqrt(pscale)

            ex(i,j) = ex(i,j) * sqrt(pscale)
            ey(i,j) = ey(i,j) * sqrt(pscale)
            ez(i,j) = ez(i,j) * sqrt(pscale)

            bxwave(i,j) = bxwave(i,j) * sqrt(pscale)
            bywave(i,j) = bywave(i,j) * sqrt(pscale)
            bzwave(i,j) = bzwave(i,j) * sqrt(pscale)

            ealpha(i,j) = ealpha(i,j) * sqrt(pscale)
            ebeta(i,j) = ebeta(i,j) * sqrt(pscale)
            eb(i,j) = eb(i,j) * sqrt(pscale)

            eplus(i,j) = eplus(i,j) * sqrt(pscale)
            eminus(i,j) = eplus(i,j) * sqrt(pscale)         
            
            ntilda_e(i,j) = ntilda_e(i,j) * sqrt(pscale)
            ntilda_e_real(i,j) = real(ntilda_e(i,j))

            xjpxe(i,j) = xjpxe(i,j) * sqrt(pscale)
            xjpye(i,j) = xjpye(i,j) * sqrt(pscale)
            xjpze(i,j) = xjpze(i,j) * sqrt(pscale)

            xjpxe_lab(i,j) = xjpxe_lab(i,j) * sqrt(pscale)
            xjpye_lab(i,j) = xjpye_lab(i,j) * sqrt(pscale)
            xjpze_lab(i,j) = xjpze_lab(i,j) * sqrt(pscale)

            xjpx_ehst(i,j) = xjpx_ehst(i,j) * sqrt(pscale)
            xjpy_ehst(i,j) = xjpy_ehst(i,j) * sqrt(pscale)
            xjpz_ehst(i,j) = xjpz_ehst(i,j) * sqrt(pscale)

            xjpx1(i,j) = xjpx1(i,j) * sqrt(pscale)
            xjpy1(i,j) = xjpy1(i,j) * sqrt(pscale)
            xjpz1(i,j) = xjpz1(i,j) * sqrt(pscale)

            xjpx2(i,j) = xjpx2(i,j) * sqrt(pscale)
            xjpy2(i,j) = xjpy2(i,j) * sqrt(pscale)
            xjpz2(i,j) = xjpz2(i,j) * sqrt(pscale)

            xjpx3(i,j) = xjpx3(i,j) * sqrt(pscale)
            xjpy3(i,j) = xjpy3(i,j) * sqrt(pscale)
            xjpz3(i,j) = xjpz3(i,j) * sqrt(pscale)

            xjpx4(i,j) = xjpx4(i,j) * sqrt(pscale)
            xjpy4(i,j) = xjpy4(i,j) * sqrt(pscale)
            xjpz4(i,j) = xjpz4(i,j) * sqrt(pscale)

            xjpx5(i,j) = xjpx5(i,j) * sqrt(pscale)
            xjpy5(i,j) = xjpy5(i,j) * sqrt(pscale)
            xjpz5(i,j) = xjpz5(i,j) * sqrt(pscale)

            xjpx6(i,j) = xjpx6(i,j) * sqrt(pscale)
            xjpy6(i,j) = xjpy6(i,j) * sqrt(pscale)
            xjpz6(i,j) = xjpz6(i,j) * sqrt(pscale)

            xj_slox(i,j) = xj_slox(i,j) * sqrt(pscale)
            xj_sloy(i,j) = xj_sloy(i,j) * sqrt(pscale)
            xj_sloz(i,j) = xj_sloz(i,j) * sqrt(pscale)


            redotj(i,j) = redotj(i,j) * pscale

            pc(i,j) = pc(i,j) * pscale
            pcre(i,j) = pcre(i,j) * pscale
            pcim(i,j) = pcim(i,j) * pscale

            wdote(i,j)  = wdote(i,j)  * pscale
            wdoti1(i,j) = wdoti1(i,j) * pscale
            wdoti2(i,j) = wdoti2(i,j) * pscale
            wdoti3(i,j) = wdoti3(i,j) * pscale
            wdoti4(i,j) = wdoti4(i,j) * pscale
            wdoti5(i,j) = wdoti5(i,j) * pscale
            wdoti6(i,j) = wdoti6(i,j) * pscale

            wdot(i,j)   = wdot(i,j)  * pscale

            xjprl(i,j)  = xjprl(i,j)  * pscale

            fype(i,j)  = fype(i,j)  * pscale
            fypi1(i,j) = fypi1(i,j) * pscale
            fypi2(i,j) = fypi2(i,j) * pscale
            fypi3(i,j) = fypi3(i,j) * pscale
            fypi4(i,j) = fypi4(i,j) * pscale
            fypi5(i,j) = fypi5(i,j) * pscale
            fypi6(i,j) = fypi6(i,j) * pscale
            fyp(i,j)   = fyp(i,j)   * pscale


            fz0e(i,j)  = fz0e(i,j)  * pscale
            fz0i1(i,j) = fz0i1(i,j) * pscale
            fz0i2(i,j) = fz0i2(i,j) * pscale
            fz0i3(i,j) = fz0i3(i,j) * pscale
            fz0i4(i,j) = fz0i4(i,j) * pscale
            fz0i5(i,j) = fz0i5(i,j) * pscale
            fz0i6(i,j) = fz0i6(i,j) * pscale
            fz0(i,j)   = fz0(i,j)   * pscale

            capr_fzeta(i,j) = capr_fzeta(i,j) * pscale

          fpsi0(i, j)  = fpsi0(i, j)  * pscale
            ftheta0(i,j) = ftheta0(i,j) * pscale
          bdotf(i,j) = bdotf(i,j) * pscale

            fpsi1(i,j)= fpsi1(i,j) * pscale
            ftheta1(i,j) = ftheta1(i,j) * pscale

         end do
      end do

      do n = 1, nnoderho
         wdote_ql(n)  = wdote_ql(n)  * pscale
         wdoti1_ql(n) = wdoti1_ql(n) * pscale
         wdoti2_ql(n) = wdoti2_ql(n) * pscale
         wdoti3_ql(n) = wdoti3_ql(n) * pscale
       wdoti4_ql(n) = wdoti4_ql(n) * pscale
       wdoti5_ql(n) = wdoti5_ql(n) * pscale
       wdoti6_ql(n) = wdoti6_ql(n) * pscale

         do niu = 1, nuper
            do miu = 1, nupar

               bqlavg_e(niu, miu, n) = bqlavg_e(niu, miu, n) * pscale
               cqlavg_e(niu, miu, n) = cqlavg_e(niu, miu, n) * pscale
               eqlavg_e(niu, miu, n) = eqlavg_e(niu, miu, n) * pscale
               fqlavg_e(niu, miu, n) = fqlavg_e(niu, miu, n) * pscale


               bqlavg_i1(niu, miu, n) = bqlavg_i1(niu, miu, n) * pscale
               cqlavg_i1(niu, miu, n) = cqlavg_i1(niu, miu, n) * pscale
               eqlavg_i1(niu, miu, n) = eqlavg_i1(niu, miu, n) * pscale
               fqlavg_i1(niu, miu, n) = fqlavg_i1(niu, miu, n) * pscale

               if(eta2 .ne. 0.0)then
               bqlavg_i2(niu, miu, n) = bqlavg_i2(niu, miu, n) * pscale
               cqlavg_i2(niu, miu, n) = cqlavg_i2(niu, miu, n) * pscale
               eqlavg_i2(niu, miu, n) = eqlavg_i2(niu, miu, n) * pscale
               fqlavg_i2(niu, miu, n) = fqlavg_i2(niu, miu, n) * pscale
             end if

             if(eta3 .ne. 0.0)then
             bqlavg_i3(niu, miu, n) = bqlavg_i3(niu, miu, n) * pscale
               cqlavg_i3(niu, miu, n) = cqlavg_i3(niu, miu, n) * pscale
               eqlavg_i3(niu, miu, n) = eqlavg_i3(niu, miu, n) * pscale
               fqlavg_i3(niu, miu, n) = fqlavg_i3(niu, miu, n) * pscale
             end if

             if(eta4 .ne. 0.0)then
             bqlavg_i4(niu, miu, n) = bqlavg_i4(niu, miu, n) * pscale
               cqlavg_i4(niu, miu, n) = cqlavg_i4(niu, miu, n) * pscale
               eqlavg_i4(niu, miu, n) = eqlavg_i4(niu, miu, n) * pscale
               fqlavg_i4(niu, miu, n) = fqlavg_i4(niu, miu, n) * pscale
             end if

             if(eta5 .ne. 0.0)then
             bqlavg_i5(niu, miu, n) = bqlavg_i5(niu, miu, n) * pscale
               cqlavg_i5(niu, miu, n) = cqlavg_i5(niu, miu, n) * pscale
               eqlavg_i5(niu, miu, n) = eqlavg_i5(niu, miu, n) * pscale
               fqlavg_i5(niu, miu, n) = fqlavg_i5(niu, miu, n) * pscale
             end if

             if(eta6 .ne. 0.0)then
             bqlavg_i6(niu, miu, n) = bqlavg_i6(niu, miu, n) * pscale
               cqlavg_i6(niu, miu, n) = cqlavg_i6(niu, miu, n) * pscale
               eqlavg_i6(niu, miu, n) = eqlavg_i6(niu, miu, n) * pscale
               fqlavg_i6(niu, miu, n) = fqlavg_i6(niu, miu, n) * pscale
             end if


            end do
         end do

      end do

      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(6 , 2841) tmin
         write(15, 2841) tmin
      endif

 2841 format('time to do power scaling =', f9.3, ' min')



      t1 = second1(dummy)

!DLG: Write the mchoi file contents into a netCDF file but
!     on the R,z grid instead of the psi/theta grid ;-)

        dlgTime = second1 ( dummy )

        if ( myId .eq. 0 ) then

        write (*,*) 'WRITING output/mchoi_dlg.nc ...'
        ncFileName = 'output/mchoi_dlg.nc'

            call check ( &
       nf90_create ( ncFileName, nf90_clobber, nc_id ) )
            call check ( &
       nf90_def_dim ( nc_id, "nR", nnodex, nR_id ) )
            call check ( &
       nf90_def_dim ( nc_id, "nZ", nnodey, nz_id ) )

            call check ( &
       nf90_def_var ( nc_id, "ePlus", NF90_REAL, &
       (/ nR_id, nz_id /), ePlus_id ) )
            call check ( &
       nf90_def_var ( nc_id, "eMinu", NF90_REAL, &
       (/ nR_id, nz_id /), eMinu_id ) )
            call check ( &
       nf90_def_var ( nc_id, "ePlus_img", NF90_REAL, &
       (/ nR_id, nz_id /), ePlus_img_id ) )
            call check ( &
       nf90_def_var ( nc_id, "eMinu_img", NF90_REAL, &
       (/ nR_id, nz_id /), eMinu_img_id ) )

            call check ( &
       nf90_def_var ( nc_id, "kPer_cold", NF90_REAL, &
       (/ nR_id, nz_id /), kPer_cold_id ) )
            call check ( &
       nf90_def_var ( nc_id, "kPer_img_cold", NF90_REAL, &
       (/ nR_id, nz_id /), kPer_img_cold_id ) )

            call check ( &
       nf90_def_var ( nc_id, "R", NF90_REAL, &
       (/ nR_id /), R_id ) )
            call check ( &
       nf90_def_var ( nc_id, "z", NF90_REAL, &
       (/ nz_id /), z_id ) )

            call check ( nf90_enddef ( nc_id ) )

            call check ( &
       nf90_put_var ( nc_id, R_id, capr(1:nnodex) ) )
            call check ( &
       nf90_put_var ( nc_id, z_id, y(1:nnodey) ) )
            call check ( &
       nf90_put_var ( nc_id, ePlus_id, &
       real ( eplus(1:nnodex,1:nnodey) * sqrt(pscale) ) ) )
            call check ( &
       nf90_put_var ( nc_id, eMinu_id, &
       real ( eminus(1:nnodex,1:nnodey) * sqrt(pscale) ) ) )
            call check ( &
       nf90_put_var ( nc_id, ePlus_img_id, &
       aimag ( eplus(1:nnodex,1:nnodey) * sqrt(pscale) ) ) )
            call check ( &
       nf90_put_var ( nc_id, eMinu_img_id, &
       aimag ( eminus(1:nnodex,1:nnodey) * sqrt(pscale) ) ) )

            call check ( &
       nf90_put_var ( nc_id, kPer_cold_id, &
       real ( xkperp_cold(1:nnodex,1:nnodey) ) ) )
            call check ( &
       nf90_put_var ( nc_id, kPer_img_cold_id, &
       aimag ( xkperp_cold(1:nnodex,1:nnodey) ) ) )

            call check ( nf90_close ( nc_id ) )
            write (*,*) 'DONE'

        endif

    if (myId==0) write(167,12398) 'write mchoi_dlg.nc:', &
           (dlgTime - second1(dummy))/60.0




!     -------------------------------------------------------------------
!     Convert quasilinear diffusion coefficients to CGS units for
!     Bob Harvey and CQL3D:
!
!     1. divide out the density since CQL includes density in f
!     2. multiply by 100 to change xlamda to cm
!     3. multiply by powers of vc_cgs to change to a dimensional velocity
!     -------------------------------------------------------------------

      vce_cgs = vce_mks * 100.
      vc1_cgs = vc1_mks * 100.
      vc2_cgs = vc2_mks * 100.
      vc3_cgs = vc3_mks * 100.
      vc4_cgs = vc4_mks * 100.
      vc5_cgs = vc5_mks * 100.
      vc6_cgs = vc6_mks * 100.

      do n = 1, nnoderho
         do niu = 1, nuper
            do miu = 1, nupar

               if(xnavg(n) .ne. 0.0)then
               bqlavg_e(niu, miu, n) = bqlavg_e(niu, miu, n) / xnavg(n) &
                  * vce_cgs**4 * 100.
               cqlavg_e(niu, miu, n) = cqlavg_e(niu, miu, n) / xnavg(n) &
                  * vce_cgs**3 * 100.
               eqlavg_e(niu, miu, n) = eqlavg_e(niu, miu, n) / xnavg(n) &
                  * vce_cgs**3 * 100.
               fqlavg_e(niu, miu, n) = fqlavg_e(niu, miu, n) / xnavg(n) &
                  * vce_cgs**2 * 100.

               else
                bqlavg_e(niu, miu, n) = 0.
                  cqlavg_e(niu, miu, n) = 0.
                  eqlavg_e(niu, miu, n) = 0.
                  fqlavg_e(niu, miu, n) = 0.
               end if

             if(xn1avg(n) .ne. 0.0)then
               bqlavg_i1(niu, miu, n)= bqlavg_i1(niu, miu, n)/ xn1avg(n) &
                  * vc1_cgs**4 * 100.
               cqlavg_i1(niu, miu, n)= cqlavg_i1(niu, miu, n)/ xn1avg(n) &
                  * vc1_cgs**3 * 100.
               eqlavg_i1(niu, miu, n)= eqlavg_i1(niu, miu, n)/ xn1avg(n) &
                  * vc1_cgs**3 * 100.
               fqlavg_i1(niu, miu, n)= fqlavg_i1(niu, miu, n)/ xn1avg(n) &
                  * vc1_cgs**2 * 100.

               else
                bqlavg_i1(niu, miu, n)= 0.
                  cqlavg_i1(niu, miu, n)= 0.
                  eqlavg_i1(niu, miu, n)= 0.
                  fqlavg_i1(niu, miu, n)= 0.
               end if

             if(eta2 .ne. 0.0)then
               if(xn2avg(n) .ne. 0.0)then
               bqlavg_i2(niu, miu, n)= bqlavg_i2(niu, miu, n)/ xn2avg(n) &
                  * vc2_cgs**4 * 100.
               cqlavg_i2(niu, miu, n)= cqlavg_i2(niu, miu, n)/ xn2avg(n) &
                  * vc2_cgs**3 * 100.
               eqlavg_i2(niu, miu, n)= eqlavg_i2(niu, miu, n)/ xn2avg(n) &
                  * vc2_cgs**3 * 100.
               fqlavg_i2(niu, miu, n)= fqlavg_i2(niu, miu, n)/ xn2avg(n) &
                  * vc2_cgs**2 * 100.

               else
                  bqlavg_i2(niu, miu, n)= 0.
                  cqlavg_i2(niu, miu, n)= 0.
                  eqlavg_i2(niu, miu, n)= 0.
                  fqlavg_i2(niu, miu, n)= 0.
               end if
             end if


             if(eta3 .ne. 0.0)then
               if(xn3avg(n) .ne. 0.0)then
               bqlavg_i3(niu, miu, n)= bqlavg_i3(niu, miu, n)/ xn3avg(n) &
                  * vc3_cgs**4 * 100.
               cqlavg_i3(niu, miu, n)= cqlavg_i3(niu, miu, n)/ xn3avg(n) &
                  * vc3_cgs**3 * 100.
               eqlavg_i3(niu, miu, n)= eqlavg_i3(niu, miu, n)/ xn3avg(n) &
                  * vc3_cgs**3 * 100.
               fqlavg_i3(niu, miu, n)= fqlavg_i3(niu, miu, n)/ xn3avg(n) &
                  * vc3_cgs**2 * 100.

               else
                  bqlavg_i3(niu, miu, n)= 0.
                  cqlavg_i3(niu, miu, n)= 0.
                  eqlavg_i3(niu, miu, n)= 0.
                  fqlavg_i3(niu, miu, n)= 0.
               end if
             end if


             if(eta4 .ne. 0.0)then
             if(xn4avg(n) .ne. 0.0)then
               bqlavg_i4(niu, miu, n)= bqlavg_i4(niu, miu, n)/ xn4avg(n) &
                  * vc4_cgs**4 * 100.
               cqlavg_i4(niu, miu, n)= cqlavg_i4(niu, miu, n)/ xn4avg(n) &
                  * vc4_cgs**3 * 100.
               eqlavg_i4(niu, miu, n)= eqlavg_i4(niu, miu, n)/ xn4avg(n) &
                  * vc4_cgs**3 * 100.
               fqlavg_i4(niu, miu, n)= fqlavg_i4(niu, miu, n)/ xn4avg(n) &
                  * vc4_cgs**2 * 100.

               else
                  bqlavg_i4(niu, miu, n)= 0.
                  cqlavg_i4(niu, miu, n)= 0.
                  eqlavg_i4(niu, miu, n)= 0.
                  fqlavg_i4(niu, miu, n)= 0.
               end if
             end if

             if(eta5 .ne. 0.0)then
             if(xn5avg(n) .ne. 0.0)then
               bqlavg_i5(niu, miu, n)= bqlavg_i5(niu, miu, n)/ xn5avg(n) &
                  * vc5_cgs**4 * 100.
               cqlavg_i5(niu, miu, n)= cqlavg_i5(niu, miu, n)/ xn5avg(n) &
                  * vc5_cgs**3 * 100.
               eqlavg_i5(niu, miu, n)= eqlavg_i5(niu, miu, n)/ xn5avg(n) &
                  * vc5_cgs**3 * 100.
               fqlavg_i5(niu, miu, n)= fqlavg_i5(niu, miu, n)/ xn5avg(n) &
                  * vc5_cgs**2 * 100.

               else
                  bqlavg_i5(niu, miu, n)= 0.
                  cqlavg_i5(niu, miu, n)= 0.
                  eqlavg_i5(niu, miu, n)= 0.
                  fqlavg_i5(niu, miu, n)= 0.
               end if
             end if

             if(eta6 .ne. 0.0)then
             if(xn6avg(n) .ne. 0.0)then
               bqlavg_i6(niu, miu, n)= bqlavg_i6(niu, miu, n)/ xn6avg(n) &
                  * vc6_cgs**4 * 100.
               cqlavg_i6(niu, miu, n)= cqlavg_i6(niu, miu, n)/ xn6avg(n) &
                  * vc6_cgs**3 * 100.
               eqlavg_i6(niu, miu, n)= eqlavg_i6(niu, miu, n)/ xn6avg(n) &
                  * vc6_cgs**3 * 100.
               fqlavg_i6(niu, miu, n)= fqlavg_i6(niu, miu, n)/ xn6avg(n) &
                  * vc6_cgs**2 * 100.

               else
                  bqlavg_i6(niu, miu, n)= 0.
                  cqlavg_i6(niu, miu, n)= 0.
                  eqlavg_i6(niu, miu, n)= 0.
                  fqlavg_i6(niu, miu, n)= 0.
               end if
             end if


            end do
         end do
      end do

!     --------------------------------------------------------
!     Interpolate quasilinear coefficients onto the CQL3D grid
!     --------------------------------------------------------

      if(enorm_factor_i2 .ne. 0.0) then

          if(eta1 .ne. 0.0 .and. ndisti1 .eq. 1)then

            vperp = uperp * vc1_mks
          vpara = upara * vc1_mks

            do n = 1, nnoderho
               do i_uperp = 1, nuper
                  do i_upara = 1, nupar
                   vperp_in = uperp(i_uperp) * vc1_mks_cql3d
                   vpara_in = upara(i_upara) * vc1_mks_cql3d

                   call intplt_to_cql3d(vperp_in, vpara_in, &
                              bqlavg_work(i_uperp, i_upara, n), &
                              nuper, nupar, nnoderho, bqlavg_i1, &
                              nuper, nupar, nnoderho, vperp, vpara, n)

                        call intplt_to_cql3d(vperp_in, vpara_in, &
                              cqlavg_work(i_uperp, i_upara, n), &
                              nuper, nupar, nnoderho, cqlavg_i1, &
                              nuper, nupar, nnoderho, vperp, vpara, n)

                        call intplt_to_cql3d(vperp_in, vpara_in, &
                              eqlavg_work(i_uperp, i_upara, n), &
                              nuper, nupar, nnoderho, eqlavg_i1, &
                              nuper, nupar, nnoderho, vperp, vpara, n)

                        call intplt_to_cql3d(vperp_in, vpara_in, &
                              fqlavg_work(i_uperp, i_upara, n), &
                              nuper, nupar, nnoderho, fqlavg_i1, &
                              nuper, nupar, nnoderho, vperp, vpara, n)
                end do
               end do
            end do

            bqlavg_i1 = bqlavg_work
          cqlavg_i1 = cqlavg_work
          eqlavg_i1 = eqlavg_work
          fqlavg_i1 = fqlavg_work

         end if



         if(eta2 .ne. 0.0 .and. ndisti2 .eq. 1) then

            vperp = uperp * vc2_mks
          vpara = upara * vc2_mks

            do n = 1, nnoderho
               do i_uperp = 1, nuper
                  do i_upara = 1, nupar
                   vperp_in = uperp(i_uperp) * vc2_mks_cql3d
                   vpara_in = upara(i_upara) * vc2_mks_cql3d

                   call intplt_to_cql3d(vperp_in, vpara_in, &
                              bqlavg_work(i_uperp, i_upara, n), &
                              nuper, nupar, nnoderho, bqlavg_i2, &
                              nuper, nupar, nnoderho, vperp, vpara, n)

                        call intplt_to_cql3d(vperp_in, vpara_in, &
                              cqlavg_work(i_uperp, i_upara, n), &
                              nuper, nupar, nnoderho, cqlavg_i2, &
                              nuper, nupar, nnoderho, vperp, vpara, n)

                        call intplt_to_cql3d(vperp_in, vpara_in, &
                              eqlavg_work(i_uperp, i_upara, n), &
                              nuper, nupar, nnoderho, eqlavg_i2, &
                              nuper, nupar, nnoderho, vperp, vpara, n)

                        call intplt_to_cql3d(vperp_in, vpara_in, &
                              fqlavg_work(i_uperp, i_upara, n), &
                              nuper, nupar, nnoderho, fqlavg_i2, &
                              nuper, nupar, nnoderho, vperp, vpara, n)
                end do
               end do
            end do

            bqlavg_i2 = bqlavg_work
          cqlavg_i2 = cqlavg_work
          eqlavg_i2 = eqlavg_work
          fqlavg_i2 = fqlavg_work

         end if


      end if





      if (myid .eq. 0) then

!        -------------------------------------------------------------------------
!        Write electron quasilinear diffusion coefficients to file out_cql3d.coefe
!        -------------------------------------------------------------------------
       if(ndiste .eq. 1)then

          write (40, 309) nuper
            write (40, 309) nupar
            write (40, 309) nnoderho

            write (40, 3310) vce_cgs
            write (40, 3310) UminPara, UmaxPara

            write (40, 3310) (rhon(n), n = 1, nnoderho)
            write (40, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            write (40, 3310) (upara(i_upara), i_upara = 1, nupar)


            write (40, 3310) (((bqlavg_e(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (40, 3310) (((cqlavg_e(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (40, 3310) (((eqlavg_e(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (40, 3310) (((fqlavg_e(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (40, 3310) xme

         end if

!        ----------------------------------------------------------------
!        Write quasilinear diffusion coefficients to file out_cql3d.coef1
!        ----------------------------------------------------------------
!DLG:   Adjust if statement for ndisti1 >= 1
       if(ndisti1 .ge. 1)then

          write (41, 309) nuper
            write (41, 309) nupar
            write (41, 309) nnoderho

            write (41, 3310) vc1_cgs_cql3d
            write (41, 3310) UminPara, UmaxPara

            write (41, 3310) (rhon(n), n = 1, nnoderho)
            write (41, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            write (41, 3310) (upara(i_upara), i_upara = 1, nupar)


            write (41, 3310) (((bqlavg_i1(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (41, 3310) (((cqlavg_i1(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (41, 3310) (((eqlavg_i1(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (41, 3310) (((fqlavg_i1(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (41, 3310) xmi1

         end if


!        ----------------------------------------------------------------
!        Write quasilinear diffusion coefficients to file out_cql3d.coef2
!        ----------------------------------------------------------------
!DLG:   Adjust if statement for ndisti2 >= 1
       if(eta2 .ne. 0.0 .and. ndisti2 .ge. 1)then

          write (42, 309) nuper
            write (42, 309) nupar
            write (42, 309) nnoderho

            write (42, 3310) vc2_cgs_cql3d
            write (42, 3310) UminPara, UmaxPara

            write (42, 3310) (rhon(n), n = 1, nnoderho)
            write (42, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            write (42, 3310) (upara(i_upara), i_upara = 1, nupar)

            write (42, 3310) (((bqlavg_i2(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (42, 3310) (((cqlavg_i2(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (42, 3310) (((eqlavg_i2(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (42, 3310) (((fqlavg_i2(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (42, 3310) xmi2

         end if

!        ----------------------------------------------------------------
!        Write quasilinear diffusion coefficients to file out_cql3d.coef3
!        ----------------------------------------------------------------
!DLG:   Adjust if statement for ndisti3 >= 1
       if(eta3 .ne. 0.0 .and. ndisti3 .ge. 1)then

          write (43, 309) nuper
            write (43, 309) nupar
            write (43, 309) nnoderho

            write (43, 3310) vc3_cgs_cql3d
            write (43, 3310) UminPara, UmaxPara

            write (43, 3310) (rhon(n), n = 1, nnoderho)
            write (43, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            write (43, 3310) (upara(i_upara), i_upara = 1, nupar)

            write (43, 3310) (((bqlavg_i3(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (43, 3310) (((cqlavg_i3(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (43, 3310) (((eqlavg_i3(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (43, 3310) (((fqlavg_i3(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (43, 3310) xmi3

         end if


!        ----------------------------------------------------------------
!        Write quasilinear diffusion coefficients to file out_cql3d.coef4
!        ----------------------------------------------------------------
!DLG:   Adjust if statement for ndist4 >= 1
       if(eta4 .ne. 0.0 .and. ndisti4 .ge. 1)then

          write (44, 309) nuper
            write (44, 309) nupar
            write (44, 309) nnoderho

            write (44, 3310) vc4_cgs_cql3d
            write (44, 3310) UminPara, UmaxPara

            write (44, 3310) (rhon(n), n = 1, nnoderho)
            write (44, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            write (44, 3310) (upara(i_upara), i_upara = 1, nupar)

            write (44, 3310) (((bqlavg_i4(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (44, 3310) (((cqlavg_i4(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (44, 3310) (((eqlavg_i4(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

            write (44, 3310) (((fqlavg_i4(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

          write (44, 3310) xmi4


         end if

       close (40)
         close (41)
       close (42)
       close (43)
       close (44)

      end if


      tmin = (second1(dummy) - t1) / 60.

      if (myid.eq.0) then
         write(6 , 2842) tmin
         write(15, 2842) tmin
      endif

 2842 format('time to write quasilinear operator in cgs units =', f9.3, &
          ' min')


      t1 = second1(dummy)
!     ------------------------------------------------------------------------
!     calculate individual species edotj's and driven current on solution mesh
!     ------------------------------------------------------------------------

      do i = 1, nnodex
         do j = 1, nnodey
            redotje(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpxe(i, j) &
                                     + conjg(ebeta(i, j)) * xjpye(i, j) &
                                     + conjg(eb(i, j)) * xjpze(i, j) )
            redotj_ehst(i, j) = &
                       0.5 * real(conjg(ealpha(i, j)) * xjpx_ehst(i, j) &
                                 + conjg(ebeta(i, j)) * xjpy_ehst(i, j) &
                                  + conjg(eb(i, j)) * xjpz_ehst(i, j) )

            redotj1(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpx1(i, j) &
                                     + conjg(ebeta(i, j)) * xjpy1(i, j) &
                                     + conjg(eb(i, j)) * xjpz1(i, j) )

            redotj2(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpx2(i, j) &
                                     + conjg(ebeta(i, j)) * xjpy2(i, j) &
                                     + conjg(eb(i, j)) * xjpz2(i, j) )

            redotj3(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpx3(i, j) &
                                     + conjg(ebeta(i, j)) * xjpy3(i, j) &
                                     + conjg(eb(i, j)) * xjpz3(i, j) )

            redotj4(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpx4(i, j) &
                                     + conjg(ebeta(i, j)) * xjpy4(i, j) &
                                     + conjg(eb(i, j)) * xjpz4(i, j) )

            redotj5(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpx5(i, j) &
                                     + conjg(ebeta(i, j)) * xjpy5(i, j) &
                                     + conjg(eb(i, j)) * xjpz5(i, j) )

            redotj6(i, j) = 0.5 * real(conjg(ealpha(i, j)) * xjpx6(i, j) &
                                     + conjg(ebeta(i, j)) * xjpy6(i, j) &
                                     + conjg(eb(i, j)) * xjpz6(i, j) )

            redotjs(i, j) = 0.5 * &
                                real(conjg(ealpha(i, j)) * xj_slox(i, j) &
                                    + conjg(ebeta(i, j)) * xj_sloy(i, j) &
                                    + conjg(eb(i, j)) * xj_sloz(i, j) )

            redotjt(i, j) = redotje(i, j) + redotj1(i, j) &
                          + redotj2(i, j) + redotj3(i, j) &
                          + redotj4(i, j) + redotj5(i, j) &
                          + redotj6(i, j) + redotjs(i, j)

            redotji(i, j) = redotj1(i, j) &
                          + redotj2(i, j) + redotj3(i, j) &
                          + redotj4(i, j) + redotj5(i, j) &
                          + redotj6(i, j)



           if(nnode_local .le. 0)xjprl(i,j) = redotj_ehst(i,j)

         end do
      end do


!     ----------------------------------------------------------
!     Integrate individual species Estar dot J on solution mesh:
!     ----------------------------------------------------------

      call sgrator(x, y, redotje, 1, nnodex, 1, nnodey, pedotje, &
         capr, rt, nxmx, nymx)
      call sgrator(x, y, redotj1, 1, nnodex, 1, nnodey, pedotj1, &
         capr, rt, nxmx, nymx)
      call sgrator(x, y, redotj2, 1, nnodex, 1, nnodey, pedotj2, &
         capr, rt, nxmx, nymx)
      call sgrator(x, y, redotj3, 1, nnodex, 1, nnodey, pedotj3, &
         capr, rt, nxmx, nymx)
      call sgrator(x, y, redotj4, 1, nnodex, 1, nnodey, pedotj4, &
         capr, rt, nxmx, nymx)
      call sgrator(x, y, redotj5, 1, nnodex, 1, nnodey, pedotj5, &
         capr, rt, nxmx, nymx)
      call sgrator(x, y, redotj6, 1, nnodex, 1, nnodey, pedotj6, &
         capr, rt, nxmx, nymx)


      call sgrator(x, y, redotjs, 1, nnodex, 1, nnodey, pedotjs, &
         capr, rt, nxmx, nymx)
      call sgrator(x, y, redotjt, 1, nnodex, 1, nnodey, pedotjt, &
         capr, rt, nxmx, nymx)

      call sgrator(x, y, xjprl, 1, nnodex, 1, nnodey, xjtot, &
         capr, rt,  nxmx, nymx)

      call sgrator(x, y, wdote, 1, nnodex, 1, nnodey, pe, &
         capr, rt, nxmx, nymx)
      call sgrator(x, y, wdoti1, 1, nnodex, 1, nnodey, pi1, &
         capr, rt, nxmx, nymx)
      call sgrator(x, y, wdoti2, 1, nnodex, 1, nnodey, pi2, &
         capr, rt, nxmx, nymx)
      call sgrator(x, y, wdoti3, 1, nnodex, 1, nnodey, pi3, &
         capr, rt, nxmx, nymx)
      call sgrator(x, y, wdoti4, 1, nnodex, 1, nnodey, pi4, &
         capr, rt, nxmx, nymx)
      call sgrator(x, y, wdoti5, 1, nnodex, 1, nnodey, pi5, &
         capr, rt, nxmx, nymx)
      call sgrator(x, y, wdoti6, 1, nnodex, 1, nnodey, pi6, &
         capr, rt, nxmx, nymx)



      call sgrator(x, y, wdot, 1, nnodex, 1, nnodey, pt, &
         capr, rt, nxmx, nymx)

      call sgrator(x, y, fyp, 1, nnodex, 1, nnodey, fyptot, &
         capr, rt, nxmx, nymx)

      call sgrator(x, y, fz0, 1, nnodex, 1, nnodey, fz0tot, &
         capr, rt, nxmx, nymx)

!     -------------------------
!     Do flux surface averages:
!     -------------------------

      call fluxavg(redotje, redotjeavg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotj1, redotj1avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotj2, redotj2avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotj3, redotj3avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotj4, redotj4avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotj5, redotj5avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotj6, redotj6avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotjs, redotjsavg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(redotjt, redotjtavg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(xjprl, xjprlavg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(wdote, wdoteavg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(wdoti1, wdoti1avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(wdoti2, wdoti2avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(wdoti3, wdoti3avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(wdoti4, wdoti4avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(wdoti5, wdoti5avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(wdoti6, wdoti6avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(wdot, wdotavg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)





!     --------------------------------------------
!     Integrate flux averaged quantities over rhon:
!     --------------------------------------------

      call rhograte(rhon, wdoteavg,  1, nnoderho, wdoteavg_int, &
         nrhomax, dvol)
      call rhograte(rhon, wdoti1avg, 1, nnoderho, wdoti1avg_int, &
         nrhomax, dvol)
      call rhograte(rhon, wdoti2avg, 1, nnoderho, wdoti2avg_int, &
         nrhomax, dvol)
      call rhograte(rhon, wdoti3avg, 1, nnoderho, wdoti3avg_int, &
         nrhomax, dvol)
      call rhograte(rhon, wdoti4avg, 1, nnoderho, wdoti4avg_int, &
         nrhomax, dvol)
      call rhograte(rhon, wdoti5avg, 1, nnoderho, wdoti5avg_int, &
         nrhomax, dvol)
      call rhograte(rhon, wdoti6avg, 1, nnoderho, wdoti6avg_int, &
         nrhomax, dvol)


      call rhograte(rhon, wdote_ql,  1, nnoderho, wdote_ql_int, &
         nrhomax, dvol)
      call rhograte(rhon, wdoti1_ql, 1, nnoderho, wdoti1_ql_int, &
         nrhomax, dvol)
      call rhograte(rhon, wdoti2_ql, 1, nnoderho, wdoti2_ql_int, &
         nrhomax, dvol)
      call rhograte(rhon, wdoti3_ql, 1, nnoderho, wdoti3_ql_int, &
         nrhomax, dvol)
      call rhograte(rhon, wdoti4_ql, 1, nnoderho, wdoti4_ql_int, &
         nrhomax, dvol)
      call rhograte(rhon, wdoti5_ql, 1, nnoderho, wdoti5_ql_int, &
         nrhomax, dvol)
      call rhograte(rhon, wdoti6_ql, 1, nnoderho, wdoti6_ql_int, &
         nrhomax, dvol)


      call rhograte(rhon, redotjeavg, 1, nnoderho, redotjeavg_int, &
         nrhomax, dvol)
      call rhograte(rhon, redotj1avg, 1, nnoderho, redotj1avg_int, &
         nrhomax, dvol)
      call rhograte(rhon, redotj2avg, 1, nnoderho, redotj2avg_int, &
         nrhomax, dvol)
      call rhograte(rhon, redotj3avg, 1, nnoderho, redotj3avg_int, &
         nrhomax, dvol)
      call rhograte(rhon, redotj4avg, 1, nnoderho, redotj4avg_int, &
         nrhomax, dvol)
      call rhograte(rhon, redotj5avg, 1, nnoderho, redotj5avg_int, &
         nrhomax, dvol)
      call rhograte(rhon, redotj6avg, 1, nnoderho, redotj6avg_int, &
         nrhomax, dvol)






31109    format(3x, "wdoteavg_int = ", 1pe12.5, 1x, "Watts")
31110    format(3x, "wdote_ql_int = ", 1pe12.5, 1x, "Watts")



41109    format(3x, "wdoti1avg_int = ", 1pe12.5, 1x, "Watts")
41110    format(3x, "wdoti1_ql_int = ", 1pe12.5, 1x, "Watts")


51109    format(3x, "wdoti2avg_int = ", 1pe12.5, 1x, "Watts")
51110    format(3x, "wdoti2_ql_int = ", 1pe12.5, 1x, "Watts")



61109    format(3x, "wdoti3avg_int = ", 1pe12.5, 1x, "Watts")
61110    format(3x, "wdoti3_ql_int = ", 1pe12.5, 1x, "Watts")



91109    format(3x, "wdoti4avg_int = ", 1pe12.5, 1x, "Watts")
91110    format(3x, "wdoti4_ql_int = ", 1pe12.5, 1x, "Watts")




      call fluxavg(fyp, fypavg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fype, fypeavg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fypi1, fypi1avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fypi2, fypi2avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fypi3, fypi3avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fypi4, fypi4avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fypi5, fypi5avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fypi6, fypi6avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)




      call fluxavg(fz0, fz0avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fz0e, fz0eavg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fz0i1, fz0i1avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fz0i2, fz0i2avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fz0i3, fz0i3avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fz0i4, fz0i4avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fz0i5, fz0i5avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(fz0i6, fz0i6avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)


      call fluxavg(capr_fzeta, capr_fzeta_avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)
      call fluxavg(capr2, capr2_avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)


      call fluxavg(rhom1, rhom1avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)


      call fluxavg(muhat, muhat_avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(nu_star, nu_star_avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(qsafety, qsafety_avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(ipsi, ipsi_avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)



      call fluxavg(pressi, pressiavg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(psi_dim, psi_dim_avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(fpsi1, fpsi1_avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(ftheta1, ftheta1_avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)


      call fluxavg(bdotf, bdotf_avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(gradprlb, gradprlb_avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)

      call fluxavg(bmod2, bmod2_avg, rho, nxmx, nymx, nrhomax, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, rt, dvol, fvol)




      dpdpsi(1) = 0.0
      do n = 2, nnoderho - 1
       dpdpsi(n) = (pressiavg(n+1) - pressiavg(n-1)) &
                    / (psi_dim_avg(n+1) - psi_dim_avg(n-1))
      end do

  300 format (1p8e11.3)

!     --------------------------------------------
!     Integrate flux averaged quantities over rhon:
!     --------------------------------------------

      call rhograte(rhon, xjprlavg, 1, nnoderho, xjprl_int, &
         nrhomax, darea)

      call rhograte(rhon, fypavg, 1, nnoderho, fyp_int, &
         nrhomax, darea)

      call rhograte(rhon, fz0avg, 1, nnoderho, fz0_int, &
         nrhomax, darea)

!      write(6, 300) (xjprl_int(n), n = 1, nnoderho)

      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(6 , 2843) tmin
         write(15, 2843) tmin
      endif

 2843 format('time to do flux averages =', f9.3, ' min')


      if(muhat_avg(1) .eq. 0.0) go to 9945


!     ----------------------------------------
!     Calculate force and flow velocities:
!     ----------------------------------------
      t1 = second1(dummy)


!     ----------------------------------------
!     Calculate coefficients for Khat equation
!     ----------------------------------------

      do n = 1, nnoderho
         if (rhom1avg(n) .gt. 1.0e-15)then
          gamma_avg(n) = muhat_avg(n) &
                         * bmod2_avg(n) / (ipsi_avg(n) / rt)**2
         endif
       beta_avg(n) = (epsn(n) / qsafety_avg(n))**2 &
                                    * (1. + 2.0 * qsafety_avg(n)**2)
      end do


!     -----------------------------
!     Solve for Khat with diffusion
!     -----------------------------

      ank(1) = 0.0
      bnk(1) = - gamma_avg(1) - 4.0 * dn(1) * beta_avg(1) / dr**2
      cnk(1) =   4.0 * dn(1) * beta_avg(2) / dr**2
      ynk(1) = capr_fzeta_avg(1) / rt - rt * bdotf_avg(1) / ipsi_avg(1)

      do n = 2, nnoderho - 1
       ank(n) = rh(n-1) * beta_avg(n-1) * dh(n-1) / (rn(n) * dr**2)
         cnk(n) = rh(n)   * beta_avg(n+1) * dh(n)   / (rn(n) * dr**2)
         bnk(n) = - gamma_avg(n) &
                  - rh(n-1) * beta_avg(n) * dh(n-1) / (rn(n) * dr**2) &
                  - rh(n)   * beta_avg(n) * dh(n)   / (rn(n) * dr**2)
       ynk(n) = capr_fzeta_avg(n) /rt - rt * bdotf_avg(n) /ipsi_avg(n)

      end do

      ank(nnoderho) = 0.0
      bnk(nnoderho) = 1.0
      cnk(nnoderho) = 0.0
      ynk(nnoderho) = 0.0

      call tridag(ank, bnk, cnk, ynk, xkhat, nnoderho)

      do n = 1, nnoderho
             kpsi_avg(n) = 0.0
         if (rhom1avg(n) .gt. 1.0e-10)then
          kpsi_avg(n) = xkhat(n) / (rhom1avg(n) * ipsi_avg(n) / rt)

         endif
      end do


      if(myid .eq. 0) write(6,  *)
      if(myid .eq. 0) write(15,  *)
      if(myid .eq. 0) write(6,  *)"        n   fzeta-fprl     beta"
      if(myid .eq. 0) write(15, *)"        n   fzeta-fprl     beta"

      do n = 1, nnoderho - 1
         if(myid .eq. 0)write(6, 1312)n, ynk(n), beta_avg(n), &
                                         qsafety_avg(n), epsn(n)
         if(myid .eq. 0)write(15,1312)n, ynk(n), beta_avg(n), &
                                         qsafety_avg(n), epsn(n)
      end do




!     --------------
!     Solve for Jhat
!     --------------
      an(1) = 0.0
      bn(1) = - 4.0 * dn(1) / dr**2
      cn(1) = - bn(1)
      yn(1) = - capr_fzeta_avg(1) / rt

      do n = 2, nnoderho - 1
         an(n) = rh(n-1) * dh(n-1) / (rn(n) * dr**2)
         cn(n) = rh(n) * dh(n)     / (rn(n) * dr**2)
         bn(n) = - an(n) - cn(n)
       yn(n) = - capr_fzeta_avg(n) / rt
      end do

      an(nnoderho) = 0.0
      bn(nnoderho) = 1.0
      cn(nnoderho) = 0.0
      yn(nnoderho) = 0.0

      call tridag(an, bn, cn, yn, xjhat, nnoderho)



!     -------------------------
!     Calculate G(psi) and Epsi
!     -------------------------
      if(myid .eq. 0) write(6,  *)
      if(myid .eq. 0) write(15,  *)
      if(myid .eq. 0) write(6,  *)  "        n    xjhat      yn"
      if(myid .eq. 0) write(15, *)  "        n    xjhat      yn"

      do n = 1, nnoderho
         gpsi_avg(n) = 0.0
         epsi_avg(n) = 0.0

         if (rhom1avg(n) .gt. 1.0e-10)then

               gpsi_avg(n) = (xjhat(n) - rhom1avg(n) * kpsi_avg(n) * b0) &
               / (rhom1avg(n) * rt * (1. + 1.5 * epsn(n)**2))

            term2 = rhom1avg(n) * kpsi_avg(n) * b0

          if(myid .eq. 0) then
               write(6,  1312)n, xjhat(n), yn(n)
               write(15, 1312)n, xjhat(n), yn(n)
          end if

            epsi_avg(n) = gpsi_avg(n)
!     .                 - fpsi1_avg(n)
!     .                 + 1. / (qi2 * xn2avg(n)) * dpdpsi(n)

         endif
      end do



!     --------------------
!     Calculate dE / dpsi
!     --------------------

      dedpsi_avg(1) = 0.0

      do n = 2, nnoderho - 1
             dedpsi_avg(n) = (epsi_avg(n+1) - epsi_avg(n-1)) &
                       / (psi_dim_avg(n+1) - psi_dim_avg(n-1))
      end do

      dedpsi_avg(nnoderho) = (epsi_avg(nnoderho) - epsi_avg(nnoderho-1)) &
                     / (psi_dim_avg(nnoderho) - psi_dim_avg(nnoderho-1))



!     -----------------------------------
!     Put G, K, E, Jhat and dE/dpsi on 2D grid
!     -----------------------------------

      do i = 1, nnodex
         do j = 1, nnodey
            n = int(rho(i,j) / drho) + 1

            gpsi(i,j) = 0.0
          kpsi(i,j) = 0.0
            epsi(i,j) = 0.0
            dedpsi(i,j) = 0.0
          jhat(i,j) = 0.0

            if(n .lt. nnoderho .and. rho(i,j) .le. rholim)then

               gpsi(i,j) = gpsi_avg(n) + (gpsi_avg(n + 1) - gpsi_avg(n)) &
                   / drho * (rho(i,j) - rhon(n))
               kpsi(i,j) = kpsi_avg(n) + (kpsi_avg(n + 1) - kpsi_avg(n)) &
                   / drho * (rho(i,j) - rhon(n))

               epsi(i,j) = epsi_avg(n) + (epsi_avg(n + 1) - epsi_avg(n)) &
                   / drho * (rho(i,j) - rhon(n))
               dedpsi(i,j) = dedpsi_avg(n) + (dedpsi_avg(n + 1) &
                   - dedpsi_avg(n)) / drho * (rho(i,j) - rhon(n))
               jhat(i,j) = xjhat(n) + (xjhat(n + 1) &
                   - xjhat(n)) / drho * (rho(i,j) - rhon(n))

            end if

            if(n .eq. nnoderho .and. rho(i,j) .le. rholim)then

               gpsi(i,j) = gpsi_avg(n)
             kpsi(i,j) = kpsi_avg(n)
               epsi(i,j) = epsi_avg(n)
                    dedpsi(i,j) = dedpsi_avg(n)
             jhat(i,j) = xjhat(n)

            end if
         end do
      end do


!     -----------------------------------------------
!     Calculate ExB shearing rate and flow velocities
!     -----------------------------------------------

      do i = 2, nnodex - 1
         do j = 1, nnodey

            omgexb(i,j) = (capr(i) * bpol(i,j))**2 / bmod(i,j) &
                                                 * dedpsi(i,j)
            uzeta(i,j) = gpsi(i,j) * capr(i) + &
                            kpsi(i,j) * bzeta(i,j) * bmod(i,j)

!            uzeta(i,j) = jhat(i,j) / rhom1(i,j)

            utheta(i,j) = kpsi(i,j) * bpol(i,j)


         end do
      end do


      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(6 , 2844) tmin
         write(15, 2844) tmin
      endif

 2844 format('time to calculate omg_ExB =', f9.3, ' min')

 9945 continue

      if (myid .eq. 0) then
         write(6,  1309) ptot
         write(15, 1309) ptot

         write(6,  1110) pcrto2, pcito2
         write(15, 1110) pcrto2, pcito2

         write(6,  1216)powtot
         write(15, 1216)powtot

         write(6,  1217) xjtot
         write(15, 1217) xjtot

         write(6, 169)
         write(6, 1120) fxtot_ant, fxtot_plasma
         write(6, 1121) fytot_ant, fytot_plasma
         write(6, 1122) fztot_ant, fztot_plasma

         write(6, 1123) fz0tot

         write(15, 169)
         write(15, 1120) fxtot_ant, fxtot_plasma
         write(15, 1121) fytot_ant, fytot_plasma
         write(15, 1122) fztot_ant, fztot_plasma

         write(15, 1123) fz0tot

      endif

!     ---------------------------------------------
!     Calculate species fractions from Estar dot J:
!     ---------------------------------------------
      pcedotje = pedotje / pedotjt * 100.
      pcedotj1 = pedotj1 / pedotjt * 100.
      pcedotj2 = pedotj2 / pedotjt * 100.
      pcedotj3 = pedotj3 / pedotjt * 100.
      pcedotj4 = pedotj4 / pedotjt * 100.
      pcedotj5 = pedotj5 / pedotjt * 100.
      pcedotj6 = pedotj6 / pedotjt * 100.
      pcedotjs = pedotjs / pedotjt * 100.

      pcedotjt = pcedotje + pcedotj1 + pcedotj2 + pcedotj3 &
               + pcedotj4 + pcedotj5 + pcedotj6 + pcedotjs


      pedotje = pedotje * twopi * rt
      pedotj1 = pedotj1 * twopi * rt
      pedotj2 = pedotj2 * twopi * rt
      pedotj3 = pedotj3 * twopi * rt
      pedotj4 = pedotj4 * twopi * rt
      pedotj5 = pedotj5 * twopi * rt
      pedotj6 = pedotj6 * twopi * rt
      pedotjs = pedotjs * twopi * rt
      pedotjt = pedotjt * twopi * rt


71667 format(' Species absorption from Estar dot J:')
      if (myid .eq. 0) then
         write(15, 169)
         write(15, 71667)
         write(15, 1109)  pedotje, pcedotje, &
                          pedotj1, pcedotj1, &
                          pedotj2, pcedotj2, &
                          pedotj3, pcedotj3, &
                          pedotj4, pcedotj4, &
                          pedotj5, pcedotj5, &
                          pedotj6, pcedotj6, &
                          pedotjs, pcedotjs, &
                          pedotjt, pcedotjt

 1109 format( &
         3x, 35h      power absorbed by electrons = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       / &
         3x, 35h  power absorbed by majority ions = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       / &
         3x, 35h  power absorbed by minority ions = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       / &
         3x, 35hpower absorbed by 3rd ion species = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       / &
         3x, 35hpower absorbed by 4th ion species = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       / &
         3x, 35hpower absorbed by 5th ion species = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       / &
         3x, 35hpower absorbed by 6th ion species = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       / &
         3x, 35hpower absorbed by slowing species = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       / &
         3x, 35h             total power absorbed = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       )

         write(6, 169)
         write(6, 71667)
         write(6, 1109)   pedotje, pcedotje, &
                          pedotj1, pcedotj1, &
                          pedotj2, pcedotj2, &
                          pedotj3, pcedotj3, &
                          pedotj4, pcedotj4, &
                          pedotj5, pcedotj5, &
                          pedotj6, pcedotj6, &
                          pedotjs, pcedotjs, &
                          pedotjt, pcedotjt


!--      Calculate species fractions from Wdot:
         if (pt .ne. 0.0) then
            pcte =   pe / pt  * 100.
            pcti1 =  pi1 / pt * 100.
            pcti2 =  pi2 / pt * 100.
            pcti3 =  pi3 / pt * 100.
            pcti4 =  pi4 / pt * 100.
            pcti5 =  pi5 / pt * 100.
            pcti6 =  pi6 / pt * 100.
         end if

         pctt = pcte + pcti1 + pcti2 + pcti3 &
                     + pcti4 + pcti5 + pcti6


81667 format(' Species absorption from Wdot:')

         pe  = pe  * twopi * rt
       pi1 = pi1 * twopi * rt
       pi2 = pi2 * twopi * rt
       pi3 = pi3 * twopi * rt
       pi4 = pi4 * twopi * rt
       pi5 = pi5 * twopi * rt
       pi6 = pi6 * twopi * rt
       pt  = pt  * twopi * rt

         write(15, 169)
         write(15, 81667)
         write(15, 71109)  pe, pcte, &
                          pi1, pcti1, &
                          pi2, pcti2, &
                          pi3, pcti3, &
                          pi4, pcti4, &
                          pi5, pcti5, &
                          pi6, pcti6, &
                          pt,  pctt
         write(6, 169)
         write(6, 81667)
         write(6, 71109)  pe,  pcte, &
                          pi1, pcti1, &
                          pi2, pcti2, &
                          pi3, pcti3, &
                          pi4, pcti4, &
                          pi5, pcti5, &
                          pi6, pcti6, &
                          pt,  pctt


81668 format('Species absorption from quasilinear Wdot:')

         pe_ql  =  wdote_ql_int(nnoderho)
       pi1_ql = wdoti1_ql_int(nnoderho)
       pi2_ql = wdoti2_ql_int(nnoderho)
       pi3_ql = wdoti3_ql_int(nnoderho)
       pi4_ql = wdoti4_ql_int(nnoderho)
       pi5_ql = wdoti5_ql_int(nnoderho)
       pi6_ql = wdoti6_ql_int(nnoderho)

       pt_ql  = pe_ql + pi1_ql + pi2_ql + pi3_ql &
                         + pi4_ql +pi5_ql + pi6_ql


         if (pt_ql .ne. 0.0) then
            pcte_ql =   pe_ql / pt_ql  * 100.
            pcti1_ql =  pi1_ql / pt_ql * 100.
            pcti2_ql =  pi2_ql / pt_ql * 100.
            pcti3_ql =  pi3_ql / pt_ql * 100.
            pcti4_ql =  pi4_ql / pt_ql * 100.
            pcti5_ql =  pi5_ql / pt_ql * 100.
            pcti6_ql =  pi6_ql / pt_ql * 100.
         end if

         pctt_ql = pcte_ql + pcti1_ql + pcti2_ql + pcti3_ql &
                           + pcti4_ql + pcti5_ql + pcti6_ql

         write(15, 169)
         write(15, 81668)
         write(15, 71109)  pe_ql, pcte_ql, &
                          pi1_ql, pcti1_ql, &
                          pi2_ql, pcti2_ql, &
                          pi3_ql, pcti3_ql, &
                          pi4_ql, pcti4_ql, &
                          pi5_ql, pcti5_ql, &
                          pi6_ql, pcti6_ql, &
                          pt_ql,  pctt_ql
         write(6, 169)
         write(6, 81668)
         write(6, 71109)  pe_ql,  pcte_ql, &
                          pi1_ql, pcti1_ql, &
                          pi2_ql, pcti2_ql, &
                          pi3_ql, pcti3_ql, &
                          pi4_ql, pcti4_ql, &
                          pi5_ql, pcti5_ql, &
                          pi6_ql, pcti6_ql, &
                          pt_ql,  pctt_ql


71109 format( &
         3x, 35h      power absorbed by electrons = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       / &
         3x, 35h  power absorbed by majority ions = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       / &
         3x, 35h  power absorbed by minority ions = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       / &
         3x, 35hpower absorbed by 3rd ion species = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       / &
         3x, 35hpower absorbed by 4th ion species = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       / &
         3x, 35hpower absorbed by 5th ion species = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       / &
         3x, 35hpower absorbed by 6th ion species = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       / &
         3x, 35h             total power absorbed = ,1e12.5, &
         7h Watts , 3h = , f10.4, 9h %       )

!        ---------------------------------------------------
!        Calculate the differential volume elements, dvol(n):
!        ---------------------------------------------------

!         dvol  = 0.0
!         darea = 0.0

!         do i = 1, nnodex
!            do j = 1, nnodey
!               n = int(rho(i,j) / drho) + 1
!               if(n .lt. nnoderho)then
!                  dvol(n)  =  dvol(n) + dx * dy * 2.0 * pi * capr(i)
!                darea(n) = darea(n) + dx * dy * capr(i) / rt
!               end if
!            end do
!         end do



         write(38, 309) nnodex, nnodey, jmid
         write(38, 310) psilim, prfin
         write(38, 310) (x(i), i = 1, nnodex)
         write(38, 310) (y(j), j = 1, nnodey)
         write(38, 310) (capr(i), i = 1, nnodex)

         write(38, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((theta(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((xnea(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((xkte(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((xkti(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((xkperp_cold(i,j), i = 1, nnodex),j = 1,nnodey)
         write(38, 310) ((acold(i, j), i = 1, nnodex), j = 1, nnodey)

       write (38, 310) ((reomg1a(i, j), i = 1, nnodex), j = 1, nnodey)
         write (38, 310) ((reomg2a(i, j), i = 1, nnodex), j = 1, nnodey)
         write (38, 310) ((reomg3a(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((xjy(i, j),  i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((xjx(i, j),  i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((ex(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((ey(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((ez(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 309) nkx1, nkx2
         write(38, 309) nky1, nky2

         write(38, 310) (xkxsav(n), n = nkx1, nkx2)
         write(38, 310) (xkysav(m), m = nky1, nky2)

         write(38, 310) ((ealphakmod(n, m), n = nkx1, nkx2), &
             m = nky1, nky2)
         write(38, 310) ((ebetakmod(n, m), n = nkx1, nkx2), &
             m = nky1, nky2)
         write(38, 310) ((ebkmod(n, m), n = nkx1, nkx2), &
             m = nky1, nky2)

         write(38, 310) ((redotje(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((redotj1(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((redotj2(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((redotj3(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((redotjt(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((bmod(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((qsafety(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((ealpha(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((ebeta(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((eb(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((btau(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((bzeta(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((spx(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((spy(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((spz(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((ealphak(n, m), n = nkx1, nkx2), m= nky1, nky2)
         write(38, 310) ((ebetak(n, m), n = nkx1, nkx2), m = nky1, nky2)
         write(38, 310) ((ebk(n, m), n = nkx1, nkx2), m = nky1, nky2)

         write(38, 310) ((wdote(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((wdoti1(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((wdoti2(i, j), i = 1, nnodex), j = 1, nnodey)
       write(38, 310) ((wdot(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 309) nnoderho
         write(38, 310) (rhon(n), n = 1, nnoderho)
       write(38, 310) (gradprlb_avg(n), n = 1, nnoderho)

       write(38, 309) mnodetheta
       write(38, 310) (thetam(m), m = 1, mnodetheta)

       write(38, 310) (dvol(n), n = 1, nnoderho)
       write(38, 310) (volume(n), n = 1, nnoderho)

         write(38, 310) (xnavg(n), n = 1, nnoderho)
         write(38, 310) (wdoteavg(n), n = 1, nnoderho)
         write(38, 310) (wdoti1avg(n), n = 1, nnoderho)
         write(38, 310) (wdoti2avg(n), n = 1, nnoderho)
         write(38, 310) (vyavg(n), n = 1, nnoderho)
         write(38, 310) (dvydrho(n), n = 1, nnoderho)

         write(38, 310) (redotjeavg(n), n = 1, nnoderho)
         write(38, 310) (redotj1avg(n), n = 1, nnoderho)
         write(38, 310) (redotj2avg(n), n = 1, nnoderho)
         write(38, 310) (redotj3avg(n), n = 1, nnoderho)
         write(38, 310) (redotjtavg(n), n = 1, nnoderho)

         write(38, 310) ((redotjs(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((wdoti3(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) (wdoti3avg(n), n = 1, nnoderho)


         write(38, 310) ((fz0e(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((fz0i1(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((fz0i2(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((fz0i3(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((fz0(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (fz0eavg(n), n = 1, nnoderho)
         write(38, 310) (fz0i1avg(n), n = 1, nnoderho)
         write(38, 310) (fz0i2avg(n), n = 1, nnoderho)
         write(38, 310) (fz0i3avg(n), n = 1, nnoderho)
         write(38, 310) (fz0avg(n), n = 1, nnoderho)

         write(38, 310) (xjprlavg(n), n = 1, nnoderho)

         write(38, 310) (xjprl_int(n), n = 1, nnoderho)
         write(38, 310) (fz0_int(n), n = 1, nnoderho)

         write(38, 310) (xjhat(n), n = 1, nnoderho)
       write(38, 310) (xkhat(n), n = 1, nnoderho)

       write(38, 310) (gpsi_avg(n), n = 1, nnoderho)
         write(38, 310) ((gpsi(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (kpsi_avg(n), n = 1, nnoderho)
         write(38, 310) ((kpsi(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (muhat_avg(n), n = 1, nnoderho)
         write(38, 310) ((muhat(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (nu_star_avg(n), n = 1, nnoderho)
         write(38, 310) ((nu_star(i, j), i = 1, nnodex), j = 1, nnodey)

       write(38, 310) (ipsi_avg(n), n = 1, nnoderho)
       write(38, 310) ((ipsi(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) ((omgexb(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((uzeta(i,j), i = 1, nnodex), j = 1, nnodey)
       write(38, 310) ((utheta(i,j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (xn1avg(n), n = 1, nnoderho)
         write(38, 310) (xn2avg(n), n = 1, nnoderho)
         write(38, 310) (xn3avg(n), n = 1, nnoderho)
         write(38, 310) (xna_sloavg(n), n = 1, nnoderho)
         write(38, 310) (xkteavg(n), n = 1, nnoderho)
         write(38, 310) (xktiavg(n), n = 1, nnoderho)
         write(38, 310) (xkti2avg(n), n = 1, nnoderho)
         write(38, 310) (xkti3avg(n), n = 1, nnoderho)

         write(38, 310) (redotjsavg(n), n = 1, nnoderho)

         write(38, 310) ((bmod_mid(i,j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) (bmod_midavg(n), n = 1, nnoderho)

         write(38, 310) ((capr_bpol_mid2(i, j), i = 1,nnodex), &
                                                          j = 1,nnodey)
         write(38, 310) (capr_bpol_mid(n), n = 1, nnoderho)




         write(38, 310) ((redotj4(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((redotj5(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((redotj6(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (redotj4avg(n), n = 1, nnoderho)
         write(38, 310) (redotj5avg(n), n = 1, nnoderho)
         write(38, 310) (redotj6avg(n), n = 1, nnoderho)

         write(38, 310) ((wdoti4(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((wdoti5(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((wdoti6(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (wdoti4avg(n), n = 1, nnoderho)
         write(38, 310) (wdoti5avg(n), n = 1, nnoderho)
         write(38, 310) (wdoti6avg(n), n = 1, nnoderho)

         write(38, 310) ((fz0i4(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((fz0i5(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((fz0i6(i, j), i = 1, nnodex), j = 1, nnodey)

         write(38, 310) (fz0i4avg(n), n = 1, nnoderho)
         write(38, 310) (fz0i5avg(n), n = 1, nnoderho)
         write(38, 310) (fz0i6avg(n), n = 1, nnoderho)

         write(38, 310) (xn4avg(n), n = 1, nnoderho)
         write(38, 310) (xn5avg(n), n = 1, nnoderho)
         write(38, 310) (xn6avg(n), n = 1, nnoderho)

         write(38, 310) (xkti4avg(n), n = 1, nnoderho)
         write(38, 310) (xkti5avg(n), n = 1, nnoderho)
         write(38, 310) (xkti6avg(n), n = 1, nnoderho)

       write(38, 310) (wdoteavg_int(n),  n = 1, nnoderho)
         write(38, 310) (wdoti1avg_int(n), n = 1, nnoderho)
       write(38, 310) (wdoti2avg_int(n), n = 1, nnoderho)
       write(38, 310) (wdoti3avg_int(n), n = 1, nnoderho)
       write(38, 310) (wdoti4avg_int(n), n = 1, nnoderho)
       write(38, 310) (wdoti5avg_int(n), n = 1, nnoderho)
       write(38, 310) (wdoti6avg_int(n), n = 1, nnoderho)

       write(38, 310) (wdote_ql_int(n),  n = 1, nnoderho)
         write(38, 310) (wdoti1_ql_int(n), n = 1, nnoderho)
       write(38, 310) (wdoti2_ql_int(n), n = 1, nnoderho)
       write(38, 310) (wdoti3_ql_int(n), n = 1, nnoderho)
       write(38, 310) (wdoti4_ql_int(n), n = 1, nnoderho)
       write(38, 310) (wdoti5_ql_int(n), n = 1, nnoderho)
       write(38, 310) (wdoti6_ql_int(n), n = 1, nnoderho)

       write(38, 310) (redotjeavg_int(n), n = 1, nnoderho)
         write(38, 310) (redotj1avg_int(n), n = 1, nnoderho)
       write(38, 310) (redotj2avg_int(n), n = 1, nnoderho)
       write(38, 310) (redotj3avg_int(n), n = 1, nnoderho)
       write(38, 310) (redotj4avg_int(n), n = 1, nnoderho)
       write(38, 310) (redotj5avg_int(n), n = 1, nnoderho)
       write(38, 310) (redotj6avg_int(n), n = 1, nnoderho)


         write(59,309) nnodex, nnodey
         write(59,310) ((eb(i, j), i = 1, nnodex), j = 1, nnodey)

         write(60,309) nnodex, nnodey
         write(60,310) ((ealpha(i, j), i = 1, nnodex), j = 1, nnodey)

         write(51,309) nnodex, nnodey
         write(51,310) psilim
         write(51,310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)

         write(69,309) nnodex, nnodey
         write(69,310) ((wdote(i, j), i = 1, nnodex), j = 1, nnodey)
         write(69,310) ((wdoti1(i, j), i = 1, nnodex), j = 1, nnodey)
         write(69,310) ((wdoti2(i, j), i = 1, nnodex), j = 1, nnodey)
         write(69,310) ((wdoti3(i, j), i = 1, nnodex), j = 1, nnodey)


!        -------------------------------
!        write edotj data to edotj.out
!        -------------------------------

         open(unit=39, file='out_edotj', status='unknown', &
            form='formatted')
         rewind(39)

         write(39, 1314) nnodex, nnodey
         write(39, 310) ((redotje(i, j), i = 1, nnodex), j = 1, nnodey)
         write(39, 310) ((redotji(i, j), i = 1, nnodex), j = 1, nnodey)

         write(39, 310) ((redotj1(i, j), i = 1, nnodex), j = 1, nnodey)
         write(39, 310) ((redotj2(i, j), i = 1, nnodex), j = 1, nnodey)
         write(39, 310) ((redotj3(i, j), i = 1, nnodex), j = 1, nnodey)
         write(39, 310) ((redotj4(i, j), i = 1, nnodex), j = 1, nnodey)
         write(39, 310) ((redotj5(i, j), i = 1, nnodex), j = 1, nnodey)
         write(39, 310) ((redotj6(i, j), i = 1, nnodex), j = 1, nnodey)

       close (39)

!       -----------------------------------------------
!       Write file 46, "naoto_2d" for synthetic diagnostic
!       -----------------------------------------------
        open(unit=46, file='naoto_2d', status='unknown', form= 'formatted')
        write (46, 309) nnodex, nnodey, nphi
        write (46, 310) (capr(i),  i = 1, nnodex)
        write (46, 310) (y(j),  j = 1, nnodey)
        write (46, 310) ((ealpha(i, j), i = 1, nnodex), j = 1, nnodey) 
        write (46, 310) ((ebeta(i, j), i = 1, nnodex), j = 1, nnodey)  
        write (46, 310) ((eb(i, j), i = 1, nnodex), j = 1, nnodey)
        write (46, 310) ((xjpxe(i,j) , i = 1, nnodex), j = 1, nnodey)
        write (46, 310) ((xjpye(i,j) , i = 1, nnodex), j = 1, nnodey)
        write (46, 310) ((xjpze(i,j) , i = 1, nnodex), j = 1, nnodey)
        write (46, 310) ((xjpxe_lab(i,j) , i = 1, nnodex), j = 1, nnodey) 
        write (46, 310) ((xjpye_lab(i,j) , i = 1, nnodex), j = 1, nnodey) 
        write (46, 310) ((xjpze_lab(i,j) , i = 1, nnodex), j = 1, nnodey)
        write (46, 310) ((ntilda_e(i, j), i = 1, nnodex), j = 1, nnodey)
        close (46)   

      endif


      t1 = second1(dummy)

      nxdim = nxmx
      nydim = nymx

      call run_rf2x(nmodesx, nmodesy, rwleft, rwright, &
         ytop, ybottom, myid, nxdim, nydim, &
         rt, b0, rho, &
         redotje, redotji, &
         redotj1, redotj2, &
         redotj3, redotj4, &
         redotj5, redotj6, &
         xjprl, wdote, &
         wdoti1, wdoti2, &
         wdoti3, wdoti4, &
         wdoti5, wdoti6)

      if (myid.eq.0) then

         write(38, 310) (wdote_ql(n),  n = 1, nnoderho)
         write(38, 310) (wdoti1_ql(n), n = 1, nnoderho)
         write(38, 310) (wdoti2_ql(n), n = 1, nnoderho)
         write(38, 310) (wdoti3_ql(n), n = 1, nnoderho)
       write(38, 310) (wdoti4_ql(n), n = 1, nnoderho)
       write(38, 310) (wdoti5_ql(n), n = 1, nnoderho)
       write(38, 310) (wdoti6_ql(n), n = 1, nnoderho)

         write(38, 310) (dldbavg(n), n = 1, nnoderho)



         write(38, 310) ((bxwave(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((bywave(i, j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((bzwave(i, j), i = 1, nnodex), j = 1, nnodey)

       write(38, 310) ((fpsi0(i, j), i = 1, nnodex), j = 1, nnodey)
       write(38, 310) ((ftheta0(i, j), i = 1, nnodex), j = 1, nnodey)
       write(38, 310) ((pressi(i,j), i = 1, nnodex), j = 1, nnodey)
         write(38, 310) ((dldb_tot12(i, j), i = 1, nnodex), &
                                            j = 1, nnodey)

       write (38, 309) ndisti1, ndisti2, ndisti3

       write (38, 310) ((eplus_flux_plot(i, j), i = 1, nnodex), &
                                                   j = 1, nnodey)

         write (38, 310) ((eminus_flux_plot(i, j), i = 1, nnodex), &
                                                   j = 1, nnodey)

         write (38, 310) ((xkperp_flux_plot(i, j), i = 1, nnodex), &
                                                   j = 1, nnodey)




       write (38, 310) ((eplus_flux(n, m), n = 1, nnoderho), &
                                             m = 1, mnodetheta)

         write (38, 310) ((eminus_flux(n, m), n = 1, nnoderho), &
                                              m = 1, mnodetheta)

         write (38, 310) ((xkperp_flux(n, m), n = 1, nnoderho), &
                                              m = 1, mnodetheta)

         write (38, 310) ((capr_flux(n,m), n = 1, nnoderho), &
                                             m = 1, mnodetheta)
         write (38, 310) ((capz_flux(n,m), n = 1, nnoderho), &
                                             m = 1, mnodetheta)
         write (38, 310) ((bmod_flux(n,m), n = 1, nnoderho), &
                                             m = 1, mnodetheta)






         close (38)


      end if

      tmin = (second1(dummy) - t1) / 60.
      if (myid.eq.0) then
         write(6 , 2845) tmin
         write(15, 2845) tmin
      endif

 2845 format('time to call run_rf2x =', f9.3, ' min')

      call blacs_barrier(icontxt, 'All')

!     -------------------
!     write fpm (34) file
!     -------------------
      if (myid.eq.0) then

         if(nphi .eq. nphi1 .or. nt .eq. 1)then
         
            write(34,309) nmodesx, nmodesy
            write(34,310) rwleft, rwright, ytop, ybottom

            write(34,309) nnodex, nnodey
            write(34,310) (capr(i), i = 1, nnodex)
            write(34,310) (y(j), j = 1, nnodey)

            write(166, 309) nnodex, nnodey
            write(166, 310) (capr(i), i = 1, nnodex)     
            write(166, 310) (y(j), j = 1, nnodey) 
            write(166, 310) psilim      
            write(166, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey) 
            write(166, 309) nphi_number
            write(166, 309) (nphi_array(n), n = 1, nphi_number)      
            
            write(34,310) psilim, rt, b0
            write(34,310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
            write(34,310) ((xnea(i, j), i = 1, nnodex), j = 1, nnodey)

            write (34, 309) nuper
            write (34, 309) nupar
            write (34, 309) nnoderho

!DLG:   A   djust if statement for ndisti? >= 1
            if(ndisti1 .ge. 1)write (34, 3310) vc1_cgs
            if(ndisti2 .ge. 1)write (34, 3310) vc2_cgs

            write (34, 3310) UminPara, UmaxPara

            write (34, 3310) (rhon(n), n = 1, nnoderho)
            write (34, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
            write (34, 3310) (upara(i_upara), i_upara = 1, nupar)

        end if

        write(34,309) nphi

        write(34,310) ((ealpha(i, j), i = 1, nnodex), j = 1, nnodey)
        write(34,310) ((ebeta(i, j),  i = 1, nnodex), j = 1, nnodey)
        write(34,310) ((eb(i, j),     i = 1, nnodex), j = 1, nnodey)

        write(166,310) ((eplus(i, j), i = 1, nnodex), j = 1, nnodey)
        write(166,310) ((eminus(i, j),  i = 1, nnodex), j = 1, nnodey)
        write(166,310) ((eb(i, j),     i = 1, nnodex), j = 1, nnodey)
 
        write(166,310) ((ex(i, j), i = 1, nnodex), j = 1, nnodey)
        write(166,310) ((ey(i, j),  i = 1, nnodex), j = 1, nnodey)
        write(166,310) ((ez(i, j),     i = 1, nnodex), j = 1, nnodey)

        write(166,310) ((bxwave(i, j), i = 1, nnodex), j = 1, nnodey)
        write(166,310) ((bywave(i, j), i = 1, nnodex), j = 1, nnodey)
        write(166,310) ((bzwave(i, j), i = 1, nnodex), j = 1, nnodey)

        write(34,310) ((xjpx(i, j), i = 1, nnodex), j = 1, nnodey)
        write(34,310) ((xjpy(i, j), i = 1, nnodex), j = 1, nnodey)
        write(34,310) ((xjpz(i, j), i = 1, nnodex), j = 1, nnodey)

        write(34,310) ((xjpxe_lab(i, j), i = 1, nnodex), j = 1, nnodey)
        write(34,310) ((xjpye_lab(i, j), i = 1, nnodex), j = 1, nnodey)
        write(34,310) ((xjpze_lab(i, j), i = 1, nnodex), j = 1, nnodey)
 
        write(34,310) ((ntilda_e(i, j), i = 1, nnodex), j = 1, nnodey)
        write(34,310) ptot, pcrto2, pcito2, xjtot

        write(34, 309) nnoderho
        write(34, 310) (rhon(n), n = 1, nnoderho)

        write(34, 310) (redotjeavg(n), n = 1, nnoderho)
        write(34, 310) (redotj1avg(n), n = 1, nnoderho)
        write(34, 310) (redotj2avg(n), n = 1, nnoderho)
        write(34, 310) (redotj3avg(n), n = 1, nnoderho)
        write(34, 310) (redotj4avg(n), n = 1, nnoderho)
        write(34, 310) (redotj5avg(n), n = 1, nnoderho)
        write(34, 310) (redotj6avg(n), n = 1, nnoderho)

        write(34, 310) (wdoteavg(n),  n = 1, nnoderho)
        write(34, 310) (wdoti1avg(n), n = 1, nnoderho)
        write(34, 310) (wdoti2avg(n), n = 1, nnoderho)
        write(34, 310) (wdoti3avg(n), n = 1, nnoderho)
        write(34, 310) (wdoti4avg(n), n = 1, nnoderho)
        write(34, 310) (wdoti5avg(n), n = 1, nnoderho)
        write(34, 310) (wdoti6avg(n), n = 1, nnoderho)

        write(34, 310) (wdote_ql(n),  n = 1, nnoderho)
        write(34, 310) (wdoti1_ql(n), n = 1, nnoderho)
        write(34, 310) (wdoti2_ql(n), n = 1, nnoderho)
        write(34, 310) (wdoti3_ql(n), n = 1, nnoderho)
        write(34, 310) (wdoti4_ql(n), n = 1, nnoderho)
        write(34, 310) (wdoti5_ql(n), n = 1, nnoderho)

        write(34, 310) (xjprlavg(n),  n = 1, nnoderho)

!DLG:   Adjust if statement for ndisti1 >= 1
        if(ndisti1 .ge. 1)then
            write (34, 3310) (((bqlavg_i1(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) (((cqlavg_i1(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) (((eqlavg_i1(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) (((fqlavg_i1(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) xmi1
       end if

!DLG:   Adjust if statement for ndisti2 >= 1
       if(ndisti2 .ge. 1)then
            write (34, 3310) (((bqlavg_i2(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) (((cqlavg_i2(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) (((eqlavg_i2(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) (((fqlavg_i2(i_uperp, i_upara, n), &
              i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)
            write (34, 3310) xmi2
       end if

       write(6,*) "finished writing fpm; starting plots"

      end if

      call blacs_barrier(icontxt, 'All')

!     -----------------------
!     do plotting with pgplot
!     -----------------------
      if (myid.eq.0) then
         t1 = second1(dummy)

       call fieldws(prfin,mask)

       tmin = (second1(dummy) - t1) / 60.
       write(6 , 2846) tmin
         write(15, 2846) tmin
      end if

 2846 format('time to do plots =', f9.3, ' min')
 2847 format('time to sum modes =', f9.3, ' min')

!      write(6,*) "finished plots; starting deallocation"
      call blacs_barrier(icontxt, 'All')


  310 format(1p6e12.4)
 3310 format(1p6e18.10)
 8310 format(1p6e14.6)
  309 format(10i10)
 9310 format(1p7e12.4)
  311 format(1p10e12.4)
 3117 format(1p11e12.4)
 1312 format(i10, 1p10e12.4)
 1314 format (2i10, 1p11e12.4)
13149 format (4i10, 1p8e12.4)
 1313 format(10i10)
 1311 format(1p9e12.4)
   10 format(i10,1p4e10.3,i10,1pe10.3)
 1010 format(1f4.0,4f8.3,3f7.3,1f8.3,1f9.3,2f8.3)
 1009 format(3x,21h        frequency  = ,1pe12.4,7h hertz )
 1012 format(3x,21h             omgrf = ,1pe12.4,7h hertz )
 2012 format(3x,51h2*pi*rt/Ly * real part of impedance (resistance) = , &
         1pe12.4,10h ohms     )
 2013 format(3x, &
         55h2*pi*rt/Ly * imaginary part of impedance (reactance) = , &
         1pe12.4,10h ohms     )
 1014 format(3x,21h               xkz = ,1pe12.4,7h m-1   )
 1321 format(3x,21h         vph / vth = ,1pe12.4,7h       )
 1391 format(3x,21h critical shear(0) = ,1pe12.4,7h s-1   )
 1392 format(3x,21h critical shear(a) = ,1pe12.4,7h s-1   )
 1393 format(3x,21h         mu neo(a) = ,1pe12.4,7h s-1   )
 1322 format(3x,21h               vph = ,1pe12.4,7h       )
 1323 format(3x,21h               vth = ,1pe12.4,7h       )
 1714 format(3x,21h               xk0 = ,1pe12.4,7h m-1   )
 1021 format(3x,21h        n parallel = ,1pe12.4,7h       )
 1812 format(3x,21h                rt = ,1pe12.4,7h m     )
 1822 format(3x,21h            aplasm = ,1pe12.4,7h m     )
 1823 format(3x,21h              rant = ,1pe12.4,7h m     )
 1809 format(3x,21h                b0 = ,1pe12.4,7h T     )
 1813 format(3x,21h              xn10 = ,1pe12.4,7h m-3   )
 6813 format(3x,21h              xne0 = ,1pe12.4,7h m-3   )
 1814 format(3x,21h              xn20 = ,1pe12.4,7h m-3   )
 1834 format(3x,21h              xn30 = ,1pe12.4,7h m-3   )

 1815 format(3x,21h               te0 = ,1pe12.4,7h eV    )
 1821 format(3x,21h               ti0 = ,1pe12.4,7h eV    )
 1016 format(3x,21h xnue/omgrf ad hoc = ,1pe12.4,7h       )
 1017 format(3x,21h xnu1/omgrf ad hoc = ,1pe12.4,7h       )
 1018 format(3x,21h xnu2/omgrf ad hoc = ,1pe12.4,7h       )
 1013 format(3x,21h              nphi = ,i12,7h       )
 7013 format(3x,21h           nmodesx = ,i12,7h       / &
             3x,21h           nmodesy = ,i12,7h       )

 7113 format(3x,21h             nwdot = ,i12,7h       )
 7213 format(3x,21h           nnodecx = ,i12,7h       / &
             3x,21h           nnodecy = ,i12,7h       )
 7014 format(3x,21h              lmax = ,i12,7h       )
 7015 format(3x,21h           ibessel = ,i12,7h       )
 7115 format(3x,21h             nzfun = ,i12,7h       )
 7016 format(3x,21h         rhoi1 / L = ,1pe12.4,7h       )
 7017 format(3x,21h             rhoi1 = ,1pe12.4,7h m     )
 7217 format(3x,21h             qavg0 = ,1pe12.4,7h       )
 1020 format(3x,21h            nnodex = ,i12, 7h       / &
             3x,21h            nnodey = ,i12, 7h       )

 3013 format(3x,21h                i0 = ,i12,7h       )
30131 format(3x,21h             ileft = ,i12,7h       )
30132 format(3x,21h            iright = ,i12,7h       )
 3014 format(3x,21h xnuii(0)/omgti(0) = ,1pe12.4,7h s-1   )
 3015 format(3x,21h           vthi(0) = ,1pe12.4,7h m/s   )
 3016 format(3x,21h          omgti(0) = ,1pe12.4,7h s-1   )
 3017 format(3x,21h           xnup(0) = ,1pe12.4,7h s-1   )
 3018 format(3x,21h          eps**1.5 = ,1pe12.4,7h       )
71160 format(3x,21h         xnu / omg = ,1pe12.4,7h       )

 1309 format(3x, &
         35h             total power absorbed = ,1pe12.4,9h watts/m )
11091 format(3x, &
         35h       total power absorbed (old) = ,1pe12.4,9h watts/m )





81109 format( &
         3x, 35h      power absorbed by electrons = ,1e12.5, &
         9h watts/m / &
         3x, 35h  power absorbed by majority ions = ,1e12.5, &
         9h watts/m / &
         3x, 35h  power absorbed by minority ions = ,1e12.5, &
         9h watts/m )

 1112 format( &
         3x,35h      power absorbed by electrons = ,1f12.4,9h %       / &
         3x,35h  power absorbed by majority ions = ,1f12.4,9h %       / &
         3x,35h  power absorbed by minority ions = ,1f12.4,9h %       / &
         3x,35hpower absorbed by 3rd ion species = ,1f12.4,9h %       / &
         3x,35h             total power absorbed = ,1f12.4,9h %       )
 1189 format( &
         3x,36h                             ref1 = ,f12.3,2h %/ &
         3x,36h                            trans = ,f12.3,2h %/ &
         3x,36h                           trans2 = ,f12.3,2h %/ &
         3x,36h                             conv = ,f12.3,2h %/ &
         3x,36h                            conv2 = ,f12.3,2h %/ &
         3x,36h  power absorbed by minority ions = ,f12.3,2h %/ &
         3x,36h  power absorbed by majority ions = ,f12.3,2h %/ &
         3x,36hpower absorbed by 3rd ion species = ,f12.3,2h %/ &
         3x,36h      power absorbed by electrons = ,f12.3,2h %/ &
         3x,36h             total power absorbed = ,f12.3,2h %/ &
         3x,36h                        total sum = ,f12.3,2h %)
 1289 format( &
         3x,36h                         x(iedge) = ,1pe12.4,2h m/ &
         3x,36h                               bz = ,1pe12.4,2h T/ &
         3x,36h              real(eps parallel)  = ,1pe12.4,2h  / &
         3x,36h    estimate: real(eps parallel)  = ,1pe12.4,2h  / &
         3x,36h                 real (eps left)  = ,1pe12.4,2h  / &
         3x,36h     estimate:   real (eps left)  = ,1pe12.4,2h  / &
         3x,36h                real (eps right)  = ,1pe12.4,2h  / &
         3x,36h    estimate:   real (eps right)  = ,1pe12.4,2h  / &
         3x,36h                               xn = ,1pe12.4,4h m-3/ &
         3x,36h                               LN = ,1pe12.4,2h m/ &
         3x,36h                        mod2 (ez) = ,1pe12.4,9h (V/m)**2/ &
         3x,36h                            L rf  = ,1pe12.4,2h m/ &
         3x,36h                x=omgrf/omgci(i)  = ,1pe12.4,4h    / &
         3x,36h                             xkti = ,1pe12.4,4h deg/ &
         3x,36h                             xkte = ,1pe12.4,4h deg)
 1389 format( &
         3x,36h                           dlnndr = ,1pe12.4,4h m-1/ &
         3x,36h                            deprl = ,1pe12.4, &
         11h (V/m)**2/m/ &
         3x,36h                           deleft = ,1pe12.4, &
         11h (V/m)**2/m/ &
         3x,36h                           derght = ,1pe12.4, &
         11h (V/m)**2/m/ &
         3x,36h                        real (ez) = ,1pe12.4,4h V/m/ &
         3x,36h                         imag(ez) = ,1pe12.4,4h V/m)
 1290 format( &
         3x,36h             alpha rf parallel    = ,1pe12.4,2h  / &
         3x,36h estimate of alpha rf parallel    = ,1pe12.4,2h  / &
         3x,36h                 alpha rf left    = ,1pe12.4,2h  / &
         3x,36h     estimate of alpha rf left    = ,1pe12.4,2h  / &
         3x,36h                 alpha rf right   = ,1pe12.4,2h  / &
         3x,36h     estimate of alpha rf right   = ,1pe12.4,2h  / &
         3x,36h                 alpha rf total   = ,1pe12.4,2h  / &
         3x,36h electron ponderomotive potential = ,1pe12.4,3h eV/ &
         3x,36h      ion ponderomotive potential = ,1pe12.4,3h eV)
 1110 format( &
         3x,35h               real power emitted = ,1pe12.4,9h watts/m / &
         3x,35h          imaginary power emitted = ,1pe12.4,9h watts/m )
 1113 format(3x,'driven current per meter = ',1pe12.4,' Amps/m')
 1114 format(3x,'current driven per watt = ',1pe12.4,' A/W')
 1215 format(3x,'current drive efficiency = ',1pe12.4,' A/W/m**2')
 1216 format(3x,'total RF power = ',1pe12.4,' Watts')
 1217 format(3x,'total driven current = ',1pe12.4,' Amps')

 1220 format(3x,'total x force on antenna = ',1pe12.4,' Nt/m'/ &
             3x,'total y force on antenna = ',1pe12.4,' Nt/m'/ &
             3x,'total z force on antenna = ',1pe12.4,' Nt/m')


 1120 format(3x,'total x force on antenna = ',1pe12.4,' Nt/m'/ &
             3x,'total x force on  plasma = ',1pe12.4,' Nt/m'/)

 1121 format(3x,'total y force on antenna = ',1pe12.4,' Nt/m'/ &
             3x,'total y force on  plasma = ',1pe12.4,' Nt/m'/)

 1122 format(3x,'total z force on antenna = ',1pe12.4,' Nt/m'/ &
             3x,'total z force on  plasma = ',1pe12.4,' Nt/m'/)

 1123 format(3x,'total toroidal force on plasma = ',1pe12.4,' Nt/m'/)


  163 format(1h0)


 1218 format(3x,'pscale = ',1pe12.4,' ')
 1219 format(3x,'gamma = ',1pe12.4,' ')
 1111 format(3x,8h At x = ,1pe12.4,2h m, &
             3x,9h ifail = ,i5)
  162 format(1h1)
  169 format(1h )

 2162 format(1p8e12.4)
 2163 format(i5, 1p11e12.3)
 2165 format(3i10,5e10.3)
 1002 format(2i10,7e10.3)
11313 format(2i10, 1p8e12.4)
 1000 format(1i10,7e10.3)

 1001 format(8e10.3)


  164 format(1h ,1x,4h i  ,  3x,6h R(x) , &
                 6x,6h  x   ,6x,6h btau ,6x,6h bmod ,6x,6hre om1, &
                 6x,6hre om2,6x,6hre om3,6x,6hre om4,6x,6hre om5, &
                 6x,6h xne  ,6x,6h Acold)

 9164 format(1h ,1x,4h j  ,  3x,6h R(x) , &
                 6x,6h  y   ,6x,6h btau ,6x,6h bmod ,6x,6hre om1, &
                 6x,6hre om2,6x,6hre om3,6x,6hre om4,6x,6hre om5, &
                 6x,6h xne  ,6x,6h Acold)



 1164 format(1h ,3x,6h R(x) , &
                 6x,6h  x   , 6x,6hvymean,6x,6hmu*vyx, &
                 5x,7hmu*vzxx,6x,6hmu*vzx, &
                 6x,6h      ,6x,6h      ,6x,6h      , &
                 6x,6h      )

 4000 format(1h ,3x,7h  x    ,5x,7hre k**2,5x,7hre k**2, &
                 5x,7hre k**2,5x,7hre k**2,5x,7hre k**2, &
                 5x,7hre k**2,5x,7hre k**2,5x,7hre k**2)

 7000 format(1h ,3x,6h  x   ,6x,6h amp1 ,6x,6h amp2 , &
                 6x,6h amp3 ,6x,6h amp4 ,6x,6h amp5 , &
                 6x,6h amp6 ,6x,6h amp7 ,6x,6h amp8 )

 7100 format(1h ,3x,6h  x   ,6x,6h  sx1 ,6x,6h  sx2 , &
                 6x,6h  sx3 ,6x,6h  sx4 ,6x,6h  sx5 , &
                 6x,6h  sx6 ,6x,6h  sx7 ,6x,6h  sx8 , &
                 6x,6hsx sum)
 4002 format(1h ,3x,6h  x   ,6x,6hre dt1,6x,6hre dt2, &
                 6x,6hre dt3,6x,6hre dt4,6x,6hre dt5, &
                 6x,6hre dt6,6x,6hre dt7,6x,6hre dt8)

 4001 format(1h ,3x,6h  x   , 5x,7him k**2,5x,7him k**2, &
                 5x,7him k**2,5x,7him k**2,5x,7him k**2, &
                 5x,7him k**2,5x,7him k**2,5x,7him k**2)

 4003 format(1h ,3x,6h  x   ,6x,6him dt1,6x,6him dt2, &
                 6x,6him dt3,6x,6him dt4,6x,6him dt5, &
                 6x,6him dt6,6x,6him dt7,6x,6him dt8)


 1650 format(3x,32hdispersion relation             )
 7650 format(3x,32hamplitude of modes              )
 7750 format(3x,32hflux of power carried by modes  )
 1651 format(3x,32hdispersion relation check       )
 1115 format(3x,13h        kh = ,1pe12.4,7h m-1   )
 1116 format(3x,13h        lp = ,1pe12.4,7h m     )
 1117 format(3x,13h      eps2 = ,1pe12.4,7h tesla )
 1118 format(3x,13havg iotabr = ,1pe12.4,7h       )
 1119 format(2x,14hellipticity = ,1pe12.4,7h       )

 5555 continue

      call blacs_barrier(icontxt, 'All')

      time=second1(dummy)-time0

      ttotal = time/60.

      if (myid .eq. 0) then
         write(15,162)
         write(6, 162)
         write(15,899) ttotal
         write(6, 899) ttotal
!       close (15)
      endif

  899 format('total cpu time used =',f9.3,4h min)




 9930 format(3x, 'electron power and flow')
 1213 format(9x, 'i', 6x, 'R', 9x, 'wdote', 7x, 'fye')
  933 format(3x, 'species #3 ion power and flow')
 1315 format(9x, 'i', 6x, 'R', 9x, 'wdoti3', 7x, 'fyi3')
  932 format(3x, 'minority ion power and flow')
 1214 format(9x, 'i', 6x, 'R', 9x, 'wdoti2', 7x, 'fyi2')
  931 format(3x, 'majority ion power and flow')
 1319 format(9x, 'i', 6x, 'R', 9x, 'wdoti1', 7x, 'fyi1')

 9321 format(3x, 'Electric field:')
 1394 format(9x, 'i', 6x, 'R', 9x, 'Re Ex', 7x, 'Im Ex', &
                               7x, 'Re Ey', 7x, 'Im Ey', &
                               7x, 'Re Ez', 7x, 'Im Ez')

 2313 format(1p8e12.4)
 6834 format(3x,21h              eta1 = ,1pe12.4,7h       )
 6835 format(3x,21h              eta2 = ,1pe12.4,7h       )
 6836 format(3x,21h              eta3 = ,1pe12.4,7h       )



!efd-begin
!      if (myid.eq.0) then
!      call profstat()
!      endif
!efd-end


    call blacs_barrier(icontxt, 'All')
      close(5)


      if (myid.eq.0) then

         close(29)
         close(30)

         close(59)
         close(60)
         close(51)
         close(69)

         close(28)

       close(40)
         close(41)
         close(42)
       close(43)
       close(44)
       close(167)

      endif

!     -----------------
!     deallocate arrays
!     -----------------

      if(ndisti1   .eq. 0 .and. &
         ndisti2   .eq. 0 .and. &
         ndisti3   .eq. 0 .and. &
         ndisti4   .eq. 0 .and. &
         ndisti5   .eq. 0 .and. &
         ndisti6   .eq. 0 .and. &
         ndiste    .eq. 0)   then

         deallocate( UPERP )
         deallocate( UPARA )
         deallocate ( vPerp, vPara )
         deallocate( UPERP_work )
         deallocate( UPARA_work )

         deallocate( f )

       deallocate( dfdupere )
         deallocate( dfdupare )

         deallocate( dfduper1 )
         deallocate( dfdupar1 )

       deallocate( dfduper2 )
         deallocate( dfdupar2 )

         deallocate( dfduper3 )
         deallocate( dfdupar3 )

       deallocate( dfduper4 )
         deallocate( dfdupar4 )

         deallocate( dfduper5 )
         deallocate( dfdupar5 )

       deallocate( dfduper6 )
         deallocate( dfdupar6 )

         deallocate( bqlavg_e )
         deallocate( cqlavg_e )
         deallocate( eqlavg_e )
         deallocate( fqlavg_e )

         deallocate( bqlavg_i1 )
         deallocate( cqlavg_i1 )
         deallocate( eqlavg_i1 )
         deallocate( fqlavg_i1 )

         deallocate ( bqlavg_work, cqlavg_work, eqlavg_work, fqlavg_work )

         if(eta2 .ne. 0.0) then
            deallocate( bqlavg_i2 )
            deallocate( cqlavg_i2 )
            deallocate( eqlavg_i2 )
            deallocate( fqlavg_i2 )
       end if

         if(eta3 .ne. 0.0) then
          deallocate( bqlavg_i3 )
            deallocate( cqlavg_i3 )
            deallocate( eqlavg_i3 )
            deallocate( fqlavg_i3 )
       end if

       if(eta4 .ne. 0.0) then
            deallocate( bqlavg_i4 )
            deallocate( cqlavg_i4 )
            deallocate( eqlavg_i4 )
            deallocate( fqlavg_i4 )
       end if

       if(eta5 .ne. 0.0) then
          deallocate( bqlavg_i5 )
            deallocate( cqlavg_i5 )
            deallocate( eqlavg_i5 )
            deallocate( fqlavg_i5 )
       end if

       if(eta6 .ne. 0.0) then
          deallocate( bqlavg_i6 )
            deallocate( cqlavg_i6 )
            deallocate( eqlavg_i6 )
            deallocate( fqlavg_i6 )
       end if

      end if

      if(ndisti1   .ne. 0 .or. &
         ndisti2   .ne. 0 .or. &
         ndisti3   .ne. 0 .or. &
         ndisti4   .ne. 0 .or. &
         ndisti5   .ne. 0 .or. &
         ndisti6   .ne. 0 .or. &
         ndiste    .ne. 0)   then

         deallocate( UPERP )
         deallocate( UPARA )
         deallocate( UPERP_work )
         deallocate( UPARA_work )

         deallocate( f )

       deallocate( dfdupere )
         deallocate( dfdupare )

         deallocate( dfduper1 )
         deallocate( dfdupar1 )

       deallocate( dfduper2 )
         deallocate( dfdupar2 )

         deallocate( dfduper3 )
         deallocate( dfdupar3 )

       deallocate( dfduper4 )
         deallocate( dfdupar4 )

         deallocate( dfduper5 )
         deallocate( dfdupar5 )

       deallocate( dfduper6 )
         deallocate( dfdupar6 )

         deallocate( dfe_cql_uprp )
         deallocate( dfe_cql_uprl )

         deallocate( df1_cql_uprp )
         deallocate( df1_cql_uprl )

         deallocate( df2_cql_uprp )
         deallocate( df2_cql_uprl )

         deallocate( df3_cql_uprp )
         deallocate( df3_cql_uprl )

         deallocate( df4_cql_uprp )
         deallocate( df4_cql_uprl )

         deallocate( df5_cql_uprp )
         deallocate( df5_cql_uprl )

         deallocate( df6_cql_uprp )
         deallocate( df6_cql_uprl )

         deallocate( bqlavg_e )
         deallocate( cqlavg_e )
         deallocate( eqlavg_e )
         deallocate( fqlavg_e )

         deallocate( bqlavg_i1 )
         deallocate( cqlavg_i1 )
         deallocate( eqlavg_i1 )
         deallocate( fqlavg_i1 )

         if(eta2 .ne. 0.0) then
            deallocate( bqlavg_i2 )
            deallocate( cqlavg_i2 )
            deallocate( eqlavg_i2 )
            deallocate( fqlavg_i2 )
       end if

         if(eta3 .ne. 0.0) then
          deallocate( bqlavg_i3 )
            deallocate( cqlavg_i3 )
            deallocate( eqlavg_i3 )
            deallocate( fqlavg_i3 )
       end if

       if(eta4 .ne. 0.0) then
            deallocate( bqlavg_i4 )
            deallocate( cqlavg_i4 )
            deallocate( eqlavg_i4 )
            deallocate( fqlavg_i4 )
       end if

       if(eta5 .ne. 0.0) then
          deallocate( bqlavg_i5 )
            deallocate( cqlavg_i5 )
            deallocate( eqlavg_i5 )
            deallocate( fqlavg_i5 )
       end if

       if(eta6 .ne. 0.0) then
          deallocate( bqlavg_i6 )
            deallocate( cqlavg_i6 )
            deallocate( eqlavg_i6 )
            deallocate( fqlavg_i6 )
       end if

      end if

      deallocate( new_to_org )

      deallocate(niabegin_all )
      deallocate( isize_all )

!     deallocate(Btmp)
      deallocate(brhs)

      deallocate( descBtmp_all )
      deallocate( descbrhs_all )

      deallocate( p_brhs )
      deallocate( p_ipiv )

      deallocate( itable )
      deallocate( jtable )
      deallocate( mtable )
      deallocate( ntable )

 !     end if

 9000 continue
!    -----------------------
!     End of loop over nphi's
!     ------------------------

 9001 continue

      deallocate (xjpxe_lab, xjpye_lab, xjpze_lab)
      deallocate ( ntilda_e, ntilda_e_real )

      if (myid .eq. 0) then
            close ( 963 )
            close (  34 )
            close ( 166 )
      endif

      close (63)

      toroidalModeSum: &
      if ( myid == 0 .and. nphi_number > 1 ) then
              t1    = second1 ( dummy )
              call aorsa2dSum ( myId )
              tMin = ( second1 ( dummy ) - t1 ) / 60.0
              write ( 6, 2847 ) tmin
              write ( 15, 2847 ) tmin
      endif toroidalModeSum

      if ( myId == 0 ) close (15)
      close(6)

!     -------------------------
!     stop parallel environment
!     -------------------------

      call blacs_barrier(icontxt, 'All')

      call blacs_gridexit( icontxt )
      call blacs_exit(0)


      end

!
!*************************************************************************
!

        subroutine pzgecopy( m,n, A,ia,ja,descA, B,ib,jb,descB)
        integer m,n,  ia,ja,descA(*), ib,jb,descB(*)
        complex A(*), B(*)

        complex alpha,beta

        alpha = 1.0
        beta = 0.0
!       ----------------------------------------------------
!       pzgeadd is new capability found in PBLAS V2 or later
!       otherwise, it can be implmented less efficiently by
!       repeated calls to pzcopy
!       ----------------------------------------------------
        call pzgeadd( 'N', m,n, alpha, A,ia,ja,descA, &
              beta, B,ib,jb,descB )
        return
        end

!
!***************************************************************************
!


      subroutine fluxavg(f, favg, rho, nxdim, nydim, nrhodim, &
         nnodex, nnodey, nnoderho, drho, dx, dy, capr, r0, vol, fvol)

      implicit none

      integer nxdim, nydim, nrhodim, nnodex, nnodey, nnoderho
      integer n, i, j

      real f(nxdim, nydim), favg(nrhodim), rho(nxdim, nydim), r0, &
         drho, dx, dy, fvol(nrhodim), vol(nrhodim), capr(nxdim), pi, &
         twopi

      pi = 3.141592654
      twopi = 2.0 * pi

      fvol = 0.0
      vol = 0.0
      favg = 0.0

      do i = 1, nnodex
         do j = 1, nnodey
            n = int(rho(i,j) / drho) + 1
            if(n .le. nnoderho)then
               fvol(n) = fvol(n) + dx * dy * twopi * capr(i) * f(i,j)
               vol(n) =  vol(n)  + dx * dy * twopi * capr(i)
            end if
         end do
      end do


      do n = 1, nnoderho
         favg(n) = 0.0
         if(vol(n) .ne. 0.0)favg(n) = fvol(n) / vol(n)
      end do

      do n = 2, nnoderho - 1
         if(favg(n) .eq. 0.0)then
            favg(n) = (favg(n-1) + favg(n+1)) / 2.0
         end if
      end do

      if (favg(1) .eq. 0.0) favg(1) = favg(2)
      if (favg(nnoderho) .eq. 0.0) favg(nnoderho) = favg(nnoderho - 1)


  100 format (1i10, 1p8e12.4)
  102 format (2i10)
      return
      end


!
!***************************************************************************
!

      subroutine rhograte(x, f, nx1, nx2, ans, nxdim, dvol)

      implicit none

      integer n, nx1, nx2, nxdim

      real f(nxdim), x(nxdim), ans(nxdim), dvol(nxdim)

      ans(nx1) = 0.0

      do n = nx1, nx2 - 1

       ans(n+1) = ans(n) + f(n) * dvol(n)

      end do


  100 format (1i10, 1p8e12.4)
  102 format (2i10)
      return
      end
!
!***************************************************************************
!



