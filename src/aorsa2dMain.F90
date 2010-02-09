program aorsa2dMain
    
    use constants
    use eqdsk_setup_mod
    use eqdsk_dlg
    use aorsasubs_mod
    use sigma_mod
    use fourier_mod

    !use ql_myra_mod
    !use profile_mod
    !use aorsa2din_mod
    !use interp
    !use write_pf
    !use plot_aorsa2dps
    !use netcdf
    !use read_particle_f
    !use set_edge_density
    !use dlg_ant
    !use current_module
    !use sigma_module
    !use parameters
    !use rf2x_mod
    !use cql3d_setup_mod
    !use nc_check 
        
    implicit none

!   Variable list
!   -------------

    integer  n_theta, nt, k

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
    real xk_cutoff, term2, ptot_wdot, xj
    real :: qsafety_avg(nRhoMax)
    integer jmid, l
    integer nnodex_loc, nx_overlap, nnodex_eff
    integer nnodey_loc, ny_overlap, nnodey_eff
    integer istart, ifinish, jstart, jfinish, &
      i_seg, j_seg, iseg1, iseg2, jseg1, jseg2, i0, j0

    integer ieq
    integer incX, global_row, global_col
    integer ninteg, nd, izoom1, izoom2, jzoom1, jzoom2, ndf, irnc
    integer nrow, ncol, norder

    integer iant

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

    integer :: i_uperp, i_upara, i_psi_eq

    character*1 trans

    integer info
    complex b, cexpkxky
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

    real theta_(n_theta_max, n_psi_max)
    integer n_theta_(n_psi_max)

    integer ilayer,nlayer
    integer mask(nxmx, nymx), mask2(nxmx,nymx)
    integer mask_bbbs(nxmx,nymx), mask_domain(nxmx,nymx), &
      mask_lim(nxmx,nymx)

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

    integer i, j, jequat, iflag, liw, lw
    integer nnodex_p, nnodey_p

    integer nnoderho, mnodetheta

    integer n, m,  ndata, ndata2, nsmooth

    real xant, delta

    real telimj, tilimj, ti2limj, ti3limj, &
               ti4limj, ti5limj, ti6limj


    real rhoplasm, reomg1, reomg2, reomg3, &
        reomg4, reomg5, reomg6, &
        xiota0, &
        deltay, abs_deltay, delta_theta, &
        xwright, psi_lim, psi1

    real xkphi(nxmx), xktau, xkrho, xkphi0, xnphi
    real rlim

    real xkthrho, wphase, vsound, domgk, xkthdx, omgestar, rhoi10, &
       v0i, vthi10, vthe, vthi, vphase, rhoi1overl, rnz, eta2, &
       xk0, shearedge, eta1, xn1, xmi1, qi1, xmi2, &
       xmi3, t0i, t0i2, t0i3, t0e, teedge, xlnlam, &
       omgci10, omgrf, xmax, qe, qi2, qi3, &
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
    real dlg_limSigma(nxmx, nymx), dlg_wallSig
    real dlg_xnuomg(nxmx,nymx)
    real djxdx, djydy
    complex rho_ant(nxmx, nymx), rho_pla(nxmx, nymx)

    complex xb(nxmx, nymx), xc(nxmx, nymx), xd(nxmx, nymx)

    real caprc(nxmx), xcourse(nxmx), capr(nxmx), &
       x(nxmx), dxc

    real capr_flux(nrhomax, nthetamax), capz_flux(nrhomax, nthetamax), &
       xgiv, ygiv
    complex eplus_flux(nrhomax, nthetamax), fout
    real :: fOut_re
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
       capr_bpol_midavg(nrhomax), &
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
         rhom1avg(nrhomax), ipsi_avg(nrhomax), &
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

    real rho(nxmx, nymx), theta(nxmx, nymx), &
         rhohatx(nxmx, nymx), rhohaty(nxmx, nymx), &
         theta0(nxmx, nymx), psi_dim(nxmx, nymx), &
         bx(nxmx, nymx), by(nxmx, nymx), bz(nxmx, nymx), &
         btau(nxmx, nymx), bzeta(nxmx, nymx), &
         dxdth(nxmx, nymx), dzdth(nxmx, nymx), xntau(nxmx, nymx), &
         xkte(nxmx, nymx), xkti(nxmx, nymx), &
         xkti2(nxmx, nymx), xkti3(nxmx, nymx), &
         xkti4(nxmx, nymx), xkti5(nxmx, nymx), xkti6(nxmx, nymx), &
         xn1a(nxmx, nymx), xnea(nxmx, nymx), xn2a(nxmx, nymx), &
         xn3a(nxmx, nymx), omgce(nxmx, nymx), &
         xn4a(nxmx, nymx), xn5a(nxmx, nymx), xn6a(nxmx, nymx), &
         omgci1(nxmx, nymx), omgci2(nxmx, nymx), omgci3(nxmx, nymx), &
         omgci4(nxmx, nymx), omgci5(nxmx, nymx), omgci6(nxmx, nymx), &
         omgpe2(nxmx, nymx), bpol(nxmx, nymx), capr_bpol(nxmx, nymx), &
         omgp12(nxmx, nymx), omgp22(nxmx, nymx), omgp32(nxmx, nymx), &
         omgp42(nxmx, nymx), omgp52(nxmx, nymx), omgp62(nxmx, nymx), &
         xiota(nxmx, nymx), &
         psi_tor2d(nxmx, nymx)

    real dlcfs(nxmx,nymx)

    real zeff(nxmx, nymx)

    real psi_pol2d(nxmx, nymx)

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
         y(nymx), dyc

    ! storage for parallel scalapack
    ! ------------------------------

    integer            dlen_
    parameter          ( dlen_ = 9 )
    integer            ctxt_, m_, n_, mb_, nb_
    parameter          ( ctxt_ = 2, m_ = 3, n_ = 4, mb_ = 5, nb_ = 6 )
    integer            rsrc_, csrc_, lld_
    parameter          ( rsrc_ = 7, csrc_ = 8, lld_ = 9 )

    integer            p_amat_dim,p_brhs_dim,p_ipiv_dim
    integer            p_amat_size

    logical use_cur_mod
    parameter(use_cur_mod=.true.)

    complex,        allocatable::p_amat(:)

    integer, allocatable ::         p_ipiv(:), cula_ipiv(:)
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

    !   local variables for transformation to real space
    !   ------------------------------------------------

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

    !   further permutation to avoid load imbalance due to
    !   assigning cells in resonance region to the same processor
    !   ---------------------------------------------------------

    logical, parameter :: need_shuffle = .true.
    integer :: ip, nrow3,ilen
    real*8, allocatable, dimension(:) :: dtemp2
    integer, allocatable, dimension(:) :: itemp2, iperm

    character(len=3) :: ntStr

    integer indxg2p,indxg2l,indxl2g
    external indxg2p,indxg2l,indxl2g

    character(len=100) :: eqdsk_fileName
    character(len=100) :: ncFileName
    integer :: nc_id, nuper_id, nupar_id, scalar_id, &
        nrho_id, f_vvp_id, uper_id, upar_id, rho_id, &
        vc_mks_id, pScale_id

    logical :: sane, gridMatch
    integer :: nR_id, nz_id, &
        ePlus_id, eMinu_id, &
        ePlus_img_id, eMinu_img_id, &
        kPer_cold_id, kPer_img_cold_id, &
        R_id, z_id, rho_pla_id, rho_ant_id, &
        xjx_id, xjy_id, xjpx_lab_id, xjpy_lab_id, &
        scalar_2_id, antOmega_id, nPhi_id, &
        xjz_id, xjpz_lab_id, &
        bxn_id, byn_id, bzn_id


!   Begin main program
!   ------------------

    allocate (xjpxe_lab(nxmx, nymx), &
        xjpye_lab(nxmx, nymx), &
        xjpze_lab(nxmx, nymx) )
    allocate ( ntilda_e_real ( nxmx, nymx ), &
        ntilda_e ( nxmx, nymx ) )


    real :: omgrf, xk0
    real :: xmi1, xmi2
    real :: qe, qi1, qi2
    real :: te, ti1, ti2
    real, dimension(:,:) :: omgci1, omgci2
    real, dimension(:,:) :: btau, sqx
    real, dimension(:,:) :: uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz
    real :: dx, dy
    real, dimension(:) :: xPrime, x, capR, xkphi
    real, dimension(:) :: yPrime, y
    integer :: i, j, m, n
    real, dimension(:) :: xkxsav, xkysav


!   read namelist input data
!   ------------------------

    call read_nameList ()


!   read g-eqdsk file
!   -----------------

    call eqdsk_setup ( myid, eqdsk )
    call read_geqdsk ( eqdsk, plot = .false. )


!   calculate the ion cyclotron freqs
!   ---------------------------------

    omgrf = 2.0 * pi * freqcy
    xk0 = omgrf / clight

    xmi1 = amu1 * xmh
    xmi2 = amu2 * xmh
   
    qe = -q
    qi1 = z1 * q
    qi2 = z2 * q
    
    te = 1e3 * q
    ti1 = 1e3 * q
    ti2 = 1e3 * q
   
    allocate ( &
        omgci1 ( nmodesx, nmodesy ), & 
        omgci1 ( nmodesx, nmodesy ) )

    omgci1 = qi1 * b0 / xmi1
    omgci2 = qi2 * b0 / xmi2

    
!   calculate rotation matrix U
!   ---------------------------

    allocate ( &
        btau ( nmodesx, nmodesy ) &
        sqx ( nmodesx, nmodesy ), )

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


!   define x mesh: x(i), xprime(i), capr(i)
!   --------------------------------------------

    !   xprime: 0 to xmax
    !   x(i) : -xmax / 2.0   to   xmax / 2.0


    allocate ( &
        xPrime ( nmodesx ), &
        x ( nmodesx ), &
        capR ( nmodesx ), &
        xkphi ( nmodesx ) )

    dx = xmax / nmodesx
    do i = 1, nmodesx

        xprime(i) = (i-1) * dx + dx / 2.0

        !   note: the code gives slightly smoother results with dx/2.0 added

        x(i) = xprime(i) + xwleft
        capr(i) = rmaxis + x(i)

        xkphi(i) = nphi / capr(i)

    enddo


!   define y mesh: y(j), yprime(j)
!---------------------------------

    !   yprime: 0 to ymax
    !   y(j) : -ymax / 2.0   to   ymax / 2.0

    allocate ( &
        yPrime ( nmodesy ), &
        y ( nmodesy ) )

    dy = ymax / nmodesy
    do j = 1, nmodesy

        yprime(j) = (j-1) * dy + dy / 2.0

        !   note: the code gives slightly smoother results with dy/2.0 added

        y(j) = yprime(j) + ybottom

    enddo


!   Calculate the integrated volume on even mesh:
!   --------------------------------------------


    allocate ( &
        xkxsav ( nmodesx + 1 ), &
        xkysav ( nmodesy ) )

    do m = -nmodesx/2, nmodesx/2
        xkxsav(m) = 2.0 * pi * m / xmax
    enddo

    do n = -nmodesy/2 + 1, nmodesy/2
        xkysav(n) = 2.0 * pi * n / ymax
    enddo

    kperp_max_actual = sqrt(xkxsav(nmodesx /2)**2 + xkysav(nmodesy /2)**2)
    kperp_max = kperp_max_actual * 2.0
    xk_cutoff = kperp_max_actual * xkperp_cutoff


!   Take numerical derivatives
!   --------------------------

    do i = 1, nmodesx
        do j = 1, nmodesy

            call deriv_x(bmod, nxmx, nymx, i, j, nmodesx, nmodesy, dx, &
                dbdx, d2bdx2)
            call deriv_y(bmod, nxmx, nymx, i, j, nmodesx, nmodesy, dy, &
                dbdy, d2bdy2)

            !   Brambilla approximation:
            !   -----------------------

            sinth = y(j) / sqrt(x(i)**2 + y(j)**2)
            gradprlb(i,j) = bmod(i,j) / capr(i) * abs(btau(i,j) * sinth)

            if (nzfun .eq. 0)gradprlb(i,j) = 1.0e-10

            call deriv_x(rho, nxmx, nymx, i, j, nmodesx, nmodesy, dx, &
                drhodx, drhodxx)
            call deriv_y(rho, nxmx, nymx, i, j, nmodesx, nmodesy, dy, &
                drhody, drhodyy)

            gradrho = sqrt(drhodx**2 + drhody**2)
            if (gradrho .eq. 0.0) gradrho = 1.0e-08

            rhohatx(i,j) = drhodx / gradrho
            rhohaty(i,j) = drhody / gradrho

            call deriv_x(uxx, nxmx, nymx, i, j, nmodesx, nmodesy, dx, &
                dxuxx(i,j), dxxuxx(i,j))
            call deriv_x(uxy, nxmx, nymx, i, j, nmodesx, nmodesy, dx, &
                dxuxy(i,j), dxxuxy(i,j))
            call deriv_x(uxz, nxmx, nymx, i, j, nmodesx, nmodesy, dx, &
                dxuxz(i,j), dxxuxz(i,j))

            call deriv_x(uyx, nxmx, nymx, i, j, nmodesx, nmodesy, dx, &
                dxuyx(i,j), dxxuyx(i,j))
            call deriv_x(uyy, nxmx, nymx, i, j, nmodesx, nmodesy, dx, &
                dxuyy(i,j), dxxuyy(i,j))
            call deriv_x(uyz, nxmx, nymx, i, j, nmodesx, nmodesy, dx, &
                dxuyz(i,j), dxxuyz(i,j))

            call deriv_x(uzx, nxmx, nymx, i, j, nmodesx, nmodesy, dx, &
                dxuzx(i,j), dxxuzx(i,j))
            call deriv_x(uzy, nxmx, nymx, i, j, nmodesx, nmodesy, dx, &
                dxuzy(i,j), dxxuzy(i,j))
            call deriv_x(uzz, nxmx, nymx, i, j, nmodesx, nmodesy, dx, &
                dxuzz(i,j), dxxuzz(i,j))

            call deriv_y(uxx, nxmx, nymx, i, j, nmodesx, nmodesy, dy, &
                dyuxx(i,j), dyyuxx(i,j))
            call deriv_y(uxy, nxmx, nymx, i, j, nmodesx, nmodesy, dy, &
                dyuxy(i,j), dyyuxy(i,j))
            call deriv_y(uxz, nxmx, nymx, i, j, nmodesx, nmodesy, dy, &
                dyuxz(i,j), dyyuxz(i,j))

            call deriv_y(uyx, nxmx, nymx, i, j, nmodesx, nmodesy, dy, &
                dyuyx(i,j), dyyuyx(i,j))
            call deriv_y(uyy, nxmx, nymx, i, j, nmodesx, nmodesy, dy, &
                dyuyy(i,j), dyyuyy(i,j))
            call deriv_y(uyz, nxmx, nymx, i, j, nmodesx, nmodesy, dy, &
                dyuyz(i,j), dyyuyz(i,j))

            call deriv_y(uzx, nxmx, nymx, i, j, nmodesx, nmodesy, dy, &
                dyuzx(i,j), dyyuzx(i,j))
            call deriv_y(uzy, nxmx, nymx, i, j, nmodesx, nmodesy, dy, &
                dyuzy(i,j), dyyuzy(i,j))
            call deriv_y(uzz, nxmx, nymx, i, j, nmodesx, nmodesy, dy, &
                dyuzz(i,j), dyyuzz(i,j))

            call deriv_xy(uxx, nxmx, nymx, i, j, nmodesx, nmodesy, dx, dy, &
                dxyuxx(i,j))
            call deriv_xy(uxy, nxmx, nymx, i, j, nmodesx, nmodesy, dx, dy, &
                dxyuxy(i,j))
            call deriv_xy(uxz, nxmx, nymx, i, j, nmodesx, nmodesy, dx, dy, &
                dxyuxz(i,j))

            call deriv_xy(uyx, nxmx, nymx, i, j, nmodesx, nmodesy, dx, dy, &
                dxyuyx(i,j))
            call deriv_xy(uyy, nxmx, nymx, i, j, nmodesx, nmodesy, dx, dy, &
                dxyuyy(i,j))
            call deriv_xy(uyz, nxmx, nymx, i, j, nmodesx, nmodesy, dx, dy, &
                dxyuyz(i,j))

            call deriv_xy(uzx, nxmx, nymx, i, j, nmodesx, nmodesy, dx, dy, &
                dxyuzx(i,j))
            call deriv_xy(uzy, nxmx, nymx, i, j, nmodesx, nmodesy, dx, dy, &
                dxyuzy(i,j))
            call deriv_xy(uzz, nxmx, nymx, i, j, nmodesx, nmodesy, dx, dy, &
                dxyuzz(i,j))

       enddo
    enddo


!   Load x, y and z equations for spatial point (i,j) and mode number (n,m)
!   ------------------------------------------------------------------------

    i_loop: &
    do i=1,nmodesx
        j_loop: &
        do j=1,nmodesy
            m_loop: &
            do m=1,nmodesx
                n_loop: &
                do n=1,nmodesy

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

                hot_plasma: &
                if (isigma .eq. 1)then

                    call sigma_maxwellian(i, j, n, m, &
                        gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
                        xme, qe, xnea(i,j), xnuomg, &
                        xkte(i,j), omgce(i,j), omgpe2(i,j), &
                        -lmax, lmax, nzfun, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sigexx, sigexy, sigexz, &
                        sigeyx, sigeyy, sigeyz, &
                        sigezx, sigezy, sigezz, &
                        delta0, 0, omgrf, xk0, &
                        upshift, damping, xk_cutoff )

                    call sigma_maxwellian(i, j, n, m,  &
                        gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
                        xmi1, qi1, xn1a(i,j), xnuomg, &
                        xkti(i,j), omgci1(i,j), omgp12(i,j), &
                        -lmax, lmax, nzfun, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sig1xx, sig1xy, sig1xz, &
                        sig1yx, sig1yy, sig1yz, &
                        sig1zx, sig1zy, sig1zz, &
                        delta0, 0, omgrf, xk0, &
                        upshift,damping,xk_cutoff )

                    call sigma_maxwellian(i, j, n, m, &
                        gradprlb(i,j), bmod(i,j), bmod_mid(i,j), &
                        xmi2, qi2, xn2a(i,j), dlg_xnuomg(i,j), &
                        xkti2(i,j), omgci2(i,j), omgp22(i,j), &
                        -lmax, lmax, nzfun, &
                        xkxsav(n), xkysav(m), nphi, capr(i), &
                        bxn(i,j), byn(i,j), bzn(i,j), &
                        uxx(i,j), uxy(i,j), uxz(i,j), &
                        uyx(i,j), uyy(i,j), uyz(i,j), &
                        uzx(i,j), uzy(i,j), uzz(i,j), &
                        sig2xx, sig2xy, sig2xz, &
                        sig2yx, sig2yy, sig2yz, &
                        sig2zx, sig2zy, sig2zz, &
                        delta0, 0, omgrf, xk0, &
                        upshift,damping,xk_cutoff )

                endif hot_plasma

                cold_plasma: &
                if (isigma .eq. 0)then 

                    call sigmac_stix(i, j, n, m, &
                        xme, qe, xnea(i,j), dlg_xnuomg(i,j), &
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
                        xmi1, qi1, xn1a(i,j), dlg_xnuomg(i,j), &
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
                        xmi2, qi2, xn2a(i,j), dlg_xnuomg(i,j), &
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

                endif cold_plasma

                sigxx = sigexx + sig1xx + sig2xx 
                sigxy = sigexy + sig1xy + sig2xy 
                sigxz = sigexz + sig1xz + sig2xz 

                sigyx = sigeyx + sig1yx + sig2yx 
                sigyy = sigeyy + sig1yy + sig2yy 
                sigyz = sigeyz + sig1yz + sig2yz 

                sigzx = sigezx + sig1zx + sig2zx 
                sigzy = sigezy + sig1zy + sig2zy 
                sigzz = sigezz + sig1zz + sig2zz 

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

                enddo n_loop
            enddo m_loop 
        enddo j_loop
    enddo i_loop 


!   Antenna current
!   ---------------

    xant    = rant - rmaxis
    iant    = int((rant - rwleft) / dx) + 1

    !   note curden is in Amps per meter of toroidal length (2.*pi*rt).

    xjantx = curdnx / dx
    xjanty = curdny / dx
    xjantz = curdnz / dx
    xjant=sqrt(xjantx**2 + xjanty**2 + xjantz**2)

    theta_antr = theta_ant / 180. * pi

    dpsiant = dpsiant0
    dthetant = dthetant0 / 360. * 2.0 * pi
    yant_max = antlen / 2.0

    do i = 1, nmodesx
        do j = 1, nmodesy

            delta_theta = theta0(i,j) - theta_antr
            gaussantth = exp(-delta_theta**2 / dthetant**2)
            gausspsi = exp(-(psi(i,j) - psiant)**2 / dpsiant**2)

            shapey = 0.0
            deltay = y(j) - yant

            if(capr(i) .gt. rmaxis) then

                !   if(i_antenna .eq. 1) antenna current is cos(ky * y)  (DEFAULT)
                !   ------------------------------------------------ --------------

                if (i_antenna .eq. 1 .and. abs(deltay) .lt. yant_max)then
                    shapey = cos(xk0 * antlc * deltay)
                endif

            endif

            xjx(i,j) = 0.0
            xjy(i,j) = xjanty * shapey  * gausspsi
            xjz(i,j) = 0.0

       enddo
    enddo


!   precompute xx(n,i), yy(m,j)
!   ---------------------------

    do i = 1, nmodesx
        do n = -nmodesx /2, nmodesx /2
            xx(n, i) = exp(zi * xkxsav(n) * xprime(i))
            xx_inv(n,i) = 1.0/xx(n,i)
        enddo
    enddo

    do j = 1, nmodesy
        do m = -nmodesy /2 + 1, nmodesy /2
            yy(m, j) = exp(zi * xkysav(m) * yprime(j))
            yy_inv(m,j) = 1.0/yy(m,j)
        enddo
    enddo


!     Fourier transform the Stix electric fields:
!     ------------------------------------------

      call sft2d_parallel(ealpha, xmax, ymax, nmodesx, nmodesy, &
         nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2, &
         -nmodesx /2, nmodesx /2, -nmodesy /2 + 1, nmodesy /2, xx_inv, yy_inv, ealphak, dx, dy, &
         myid, nproc, icontxt)

      call sft2d_parallel(ebeta, xmax, ymax, nmodesx, nmodesy, &
         nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2, &
         -nmodesx /2, nmodesx /2, -nmodesy /2 + 1, nmodesy /2, xx_inv, yy_inv, ebetak, dx, dy, &
         myid, nproc, icontxt)

      call sft2d_parallel(eb, xmax, ymax, nmodesx, nmodesy, &
         nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2, &
         -nmodesx /2, nmodesx /2, -nmodesy /2 + 1, nmodesy /2, xx_inv, yy_inv, ebk, dx, dy, &
         myid, nproc, icontxt)



      do n = -nmodesx /2, nmodesx /2
         do m = -nmodesy /2 + 1, nmodesy /2

            ealphakmod(n, m) = sqrt(conjg(ealphak(n, m))* ealphak(n, m))
            ebetakmod(n, m) = sqrt(conjg(ebetak(n, m)) * ebetak(n, m))
            ebkmod(n, m) = sqrt(conjg(ebk(n, m)) * ebk(n, m))

         end do
      end do



!     ----------------------------------------------
!     Calculate E in the Lab frame and eplus, eminus
!     ----------------------------------------------
      isq2 = SQRT(0.5)
      do i = 1, nmodesx
         do j = 1, nmodesy

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

end program aorsa2dMain


