module plot_aorsa2dps

contains

      subroutine fieldws(prfin, mask)
      use netcdf
      use size_mod

      implicit none

      integer, intent(in) :: mask(:,:)
      character*32 title
      character*32 titll
      character*32 titlr
      character*32 titlb
      character*32 titx
      character*32 tity
      character*32 titz

      !DLG
      character(len=100) :: ncFileName
      integer :: nc_id, nR_id, nz_id, scalar_id, &
       capR_id, y_id, nRho_id, rho_id, &
       wdote_id, wdoti1_id, wdoti2_id, &
       wdote_rz_id, wdoti1_rz_id, wdoti2_rz_id, &
       pscale_id, ePlus_real_id, ePlus_imag_id, &
       eMinu_real_id, eMinu_imag_id, &
       bx_wave_real_id, bx_wave_imag_id, &
       bz_wave_real_id, bz_wave_imag_id, &
       mask_id, density_id, janty_id, jantx_id, &
       jeDotE_id, eAlpha_id, eBeta_id, eParallel_id, &
       ex_id, ey_id, ez_id, xkperp_imag_id, xkperp_real_id, &
       bmod_id


      real logmax, ycut, dy, xmax, ymax, xmi, E_eV, vperp_mks, &
       vpara_mks, vperp_cgs, uperp_1kev, duperp, dz, dx
      real exkmin, exkmax 
      integer, parameter :: DBL = selected_real_kind ( 13, 300 )
      real(kind=DBL) :: prfin, prfinMOD

      integer pgopen, pgbeg, ier, nmid, mmid
      integer n_theta_max, n_u_max, n_psi_max

      integer ndisti1, ndisti2, ndisti3, number_points, k, nnodez
      integer nxmx, nymx, nkdim1, nkdim2, mkdim1, mkdim2, nlevmax
      integer ibackground, nkx1, nkx2, nky1, nky2, n, m, &
         nkpltdim, mkpltdim, nkxplt, nkyplt, ipage, n1, n2, n3, n4
      integer i_psi, i_psi1, i_psi2, i_psi3, i_psi4, i_psi5, i_psi6
      integer i_psi_array(6), i_psi_index

      integer :: nuper, nupar, i_uperp, i_upara, nz
      real uminpara, umaxpara, vc_cgs, vpara_cgs

      integer  ncolln10, ncolln9, nwheat, ngrey, naqua, &
         npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen, &
         nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5, &
         ncolln8, nwhite, ncolbox, nmagenta, nsalmon, ncolln2, &
         ncolln3, ncolbrd, ncolln1, ncollin, ncollab, ncolion, &
         ncolelec, norange

      integer nrhomax, nnoderho, nnoderho2, iflag, nnoderho_half, &
         nnoderho2_half, mnodetheta, nthetamax


      parameter (nlevmax = 101)

      parameter (nxmx = nmodesmax)
      parameter (nymx = mmodesmax)

      parameter (nrhomax = nxmx)
      parameter (nthetamax = nymx)

      parameter (nkdim1 = - nxmx / 2)
      parameter (nkdim2 =   nxmx / 2)

      parameter (mkdim1 = - nxmx / 2)
      parameter (mkdim2 =   nxmx / 2)

      parameter (nkpltdim = 2 * nkdim2)
      parameter (mkpltdim = 2 * mkdim2)

      parameter (n_theta_max = 150)
      parameter (n_u_max = 150)
      parameter (n_psi_max = 150)

      real u(n_u_max), theta_u(n_theta_max)
      real f_cql(n_theta_max, n_u_max, n_psi_max)
      real f_cql_2d(n_u_max, n_theta_max)
      real f_cql_1d_1(n_u_max), f_cql_1d_2(n_u_max)
      real f_cql_1d_3(n_u_max), f_cql_1d_4(n_u_max)
      real f_cql_1d_5(n_u_max), f_cql_1d_6(n_u_max)
      real f_cql_1d_7(n_u_max)


      integer n_theta_(n_psi_max)
      real theta_(n_theta_max, n_psi_max)

      integer n_theta, n_u, n_psi
      integer i_theta, i_u
      integer i_theta1, i_theta2, i_theta3, i_theta4, &
                        i_theta5, i_theta6, i_theta7
      real vc, r0

      real, dimension(:),   allocatable :: UPERP, UPARA

      real, dimension(:,:), allocatable :: bqlavg_i1_2d
      real, dimension(:,:), allocatable :: E_kick_2d
      real, dimension(:,:), allocatable :: f_cql_cart_2d
      real, dimension(:,:,:), allocatable :: f_cql_cart

      real, dimension(:), allocatable :: wperp_cql
      real, dimension(:), allocatable :: wpar_cql



      real*8, dimension(:,:,:), allocatable :: bqlavg_i1
      real*8, dimension(:,:,:), allocatable :: cqlavg_i1
      real*8, dimension(:,:,:), allocatable :: eqlavg_i1
      real*8, dimension(:,:,:), allocatable :: fqlavg_i1


      real capr_bpol_mid(nrhomax), capr_bpol_midavg(nrhomax)
      real capr_bpol_mid2(nxmx, nymx), bmod_midavg(nrhomax)
      real bmod_mid(nxmx, nymx), bratio(nxmx, nymx)

      real xkxsav(nkpltdim), xkysav(mkpltdim), pscale
      real wdoti1avg(nrhomax), wdoti2avg(nrhomax), wdoti3avg(nrhomax)
      real wdoti4avg(nrhomax), wdoti5avg(nrhomax), wdoti6avg(nrhomax)
      real zdummy(3)

      real xjprl_int(nrhomax)

      real redotje_int(nrhomax), &
           redotj1_int(nrhomax), redotj2_int(nrhomax), &
           redotj3_int(nrhomax), redotj4_int(nrhomax), &
           redotj5_int(nrhomax), redotj6_int(nrhomax)

      real wdote_int(nrhomax), &
           wdot1_int(nrhomax), wdot2_int(nrhomax), &
           wdot3_int(nrhomax), wdot4_int(nrhomax), &
           wdot5_int(nrhomax), wdot6_int(nrhomax)

      real wdoti1_dvol(nrhomax), wdoti2_dvol(nrhomax), &
           wdoti3_dvol(nrhomax), wdoti4_dvol(nrhomax), &
           wdoti5_dvol(nrhomax), wdoti6_dvol(nrhomax), &
           wdote_dvol(nrhomax)

      real redotj1_dvol(nrhomax), redotj2_dvol(nrhomax), &
           redotj3_dvol(nrhomax), redotj4_dvol(nrhomax), &
           redotj5_dvol(nrhomax), redotj6_dvol(nrhomax), &
           redotje_dvol(nrhomax)


      real wdoteavg_int(nrhomax), wdoti1avg_int(nrhomax), &
                                  wdoti2avg_int(nrhomax), &
                                  wdoti3avg_int(nrhomax), &
                                  wdoti4avg_int(nrhomax), &
                                  wdoti5avg_int(nrhomax), &
                                  wdoti6avg_int(nrhomax)

      real wdote_ql_int(nrhomax), wdoti1_ql_int(nrhomax), &
                                  wdoti2_ql_int(nrhomax), &
                                  wdoti3_ql_int(nrhomax), &
                                  wdoti4_ql_int(nrhomax), &
                                  wdoti5_ql_int(nrhomax), &
                                  wdoti6_ql_int(nrhomax)


      real redotjeavg_int(nrhomax), redotj1avg_int(nrhomax), &
                                    redotj2avg_int(nrhomax), &
                                    redotj3avg_int(nrhomax), &
                                    redotj4avg_int(nrhomax), &
                                    redotj5avg_int(nrhomax), &
                                    redotj6avg_int(nrhomax)

      real rhon(nrhomax), thetam(nthetamax)

      real xnavg(nrhomax), xn1avg(nrhomax), xkteavg(nrhomax), &
           xktiavg(nrhomax), &
           xkti1avg(nrhomax), xkti2avg(nrhomax), xkti3avg(nrhomax), &
           xn2avg(nrhomax), xn3avg(nrhomax), xk3avg(nrhomax), &
           xna_sloavg(nrhomax), &
           vyi1avg(nrhomax),  vyi2avg(nrhomax), &
           dvol(nrhomax), rhon_half(nrhomax), volume(nrhomax)
      real rhon_save(nrhomax), rhon_half_save(nrhomax)
      real xkti4avg(nrhomax),  xkti5avg(nrhomax), xkti6avg(nrhomax)
      real wdote_ql(nrhomax), wdoti1_ql(nrhomax), wdoti2_ql(nrhomax), &
           wdoti3_ql(nrhomax), wdoti4_ql(nrhomax), wdoti5_ql(nrhomax), &
           wdoti6_ql(nrhomax), dldbavg(nrhomax), gradprlb2_avg(nrhomax)
      real xn4avg(nrhomax), xn5avg(nrhomax), xn6avg(nrhomax)
      real vyavg(nrhomax), dvydrho(nrhomax), wdoteavg(nrhomax)
      real fz0i1avg(nrhomax), fz0i2avg(nrhomax), fz0i3avg(nrhomax), &
           fz0eavg(nrhomax), fz0avg(nrhomax)

      real fz0i4avg(nrhomax), fz0i5avg(nrhomax), fz0i6avg(nrhomax)

      real gpsi_avg(nrhomax), kpsi_avg(nrhomax), xjhat(nrhomax), &
           muhat_avg(nrhomax), nu_star_avg(nrhomax), xkhat(nrhomax)
      real gpsi(nxmx, nymx),  kpsi(nxmx, nymx), &
           omgexb(nxmx, nymx), uzeta(nxmx, nymx), utheta(nxmx, nymx), &
           fpsi0(nxmx, nymx), ftheta0(nxmx, nymx), &
           muhat(nxmx, nymx), nu_star(nxmx, nymx)


      real redotjeavg(nrhomax), redotj1avg(nrhomax), &
           redotj2avg(nrhomax), redotj3avg(nrhomax), &
           redotj4avg(nrhomax), &
           redotj5avg(nrhomax), redotj6avg(nrhomax), &
           redotjsavg(nrhomax), &
           redotjtavg(nrhomax), xjprlavg(nrhomax), &
           redotjiavg(nrhomax), ipsi_avg(nrhomax)

      real redotj4(nxmx, nymx), redotj5(nxmx, nymx), redotj6(nxmx, nymx)

      real fz0_int(nrhomax)

      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
         ncolion,ncolelec,ncollin,ncollab,ncolln2
      common/boundcom/rhoplasm
      common/zoom/ xmaxz, xminz, ymaxz, yminz

      real xmaxz, xminz, ymaxz, yminz

      real exkmod(nkpltdim, mkpltdim), &
           eykmod(nkpltdim, mkpltdim), &
           ezkmod(nkpltdim, mkpltdim)

      real exklog(nkpltdim, mkpltdim), &
           eyklog(nkpltdim, mkpltdim), &
           ezklog(nkpltdim, mkpltdim)

      real exklogmin, exklogmax, &
           eyklogmin, eyklogmax, &
           ezklogmin, ezklogmax

      real fmodm(mkpltdim)

      complex exk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              eyk(nkdim1 : nkdim2, mkdim1 : mkdim2), &
              ezk(nkdim1 : nkdim2, mkdim1 : mkdim2)


      real reomg1a(nxmx, nymx), reomg2a(nxmx, nymx), reomg3a(nxmx, nymx)

      real capr_bpol(nxmx, nymx), pressi(nxmx, nymx)
      real redotj1(nxmx, nymx), redotj2(nxmx, nymx), redotj3(nxmx, nymx)
      real redotje(nxmx, nymx), redotjt(nxmx, nymx), redotjs(nxmx, nymx)

      real wdoti1(nxmx, nymx), wdoti2(nxmx, nymx), wdoti3(nxmx, nymx)
      real wdoti4(nxmx, nymx), wdoti5(nxmx, nymx), wdoti6(nxmx, nymx)
      real wdote(nxmx, nymx), wdott(nxmx, nymx)

      real fz0e(nxmx, nymx), fz0i1(nxmx, nymx), &
           fz0i2(nxmx, nymx), fz0i3(nxmx, nymx), fz0(nxmx, nymx), &
           fz0i4(nxmx, nymx), fz0i5(nxmx, nymx), fz0i6(nxmx, nymx)


      real fmin, fmax, fminre, fmaxre, fminim, fmaxim, fmin1, &
         fmax1, fmax2, fmax3, fmaxs, fmaxt, fmin2, fmin3, fmins,fmint
      real x(nxmx), capr(nxmx), y(nymx), xn(nxmx, nymx)
      real xkti(nxmx, nymx), xkte(nxmx, nymx)
      real rho(nxmx, nymx), theta(nxmx, nymx), psi(nxmx, nymx)
      real xjy(nxmx, nymx), bmod(nxmx, nymx)
      real xjx(nxmx, nymx)
      real xiota(nxmx, nymx), qsafety(nxmx, nymx), ipsi(nxmx, nymx)
      real btau(nxmx, nymx), bzeta(nxmx, nymx), isq2


      complex xkperp_cold(nxmx, nymx), acold(nxmx, nymx)
      real freal(nxmx, nymx), fimag(nxmx, nymx), fmod(nxmx, nymx)
      real mod_Eplus(nxmx, nymx), mod_Eminus(nxmx, nymx)
      real mod_Eb(nxmx, nymx)

      complex ex(nxmx, nymx), ey(nxmx, nymx), ez(nxmx, nymx)
      complex bxwave(nxmx, nymx), bywave(nxmx, nymx), bzwave(nxmx, nymx)


      complex ealpha(nxmx, nymx), ebeta(nxmx, nymx), eb(nxmx, nymx)
      complex eplus(nxmx, nymx), eminus(nxmx, nymx)

      complex eplus_flux(nrhomax, nthetamax)
      complex eminus_flux(nrhomax, nthetamax)
      complex xkperp_flux(nrhomax, nthetamax)

      real bmod_flux(nrhomax, nthetamax)
      real bmod_plot(nrhomax)

      real capr_flux(nrhomax, nthetamax), capz_flux(nrhomax, nthetamax)
      complex eplus_flux_plot(nxmx, nymx)
      complex eminus_flux_plot(nxmx, nymx)
      complex xkperp_flux_plot(nxmx, nymx)


      real dldb_tot12(nxmx, nymx)

      real   reex_dx(nxmx, nymx)
      real   reey_dx(nxmx, nymx)
      real   reez_dx(nxmx, nymx)

      real  ximex_dx(nxmx, nymx)
      real  ximey_dx(nxmx, nymx)
      real  ximez_dx(nxmx, nymx)



      real spx(nxmx, nymx), spy(nxmx, nymx), spz(nxmx, nymx)

      real xnmid(nxmx), xktimid(nxmx),  xktemid(nxmx), qmid(nxmx)
      real bpmid(nxmx), xiotamid(nxmx)
      real fmidre(nxmx), fmidim(nxmx), fmid1(nxmx), &
           fmid2(nxmx), fmid3(nxmx), fmid4(nxmx), fmid5(nxmx), &
           fmids(nxmx), fmidt(nxmx), fmid6(nxmx)

      real mod_Eplus_mid(nxmx),  mod_Eminus_mid(nxmx),  mod_Eb_mid(nxmx)

      real ff(101)
      real q, omgrf, xk0, n0, clight, xmu0, eps0, rhoplasm
      real temax, temin, timin, tmin, tmax, timax, caprmaxp, &
         caprmin, caprminp, caprmax, xnmax, xnmin, qmin, qmax
      real bpmin, bpmax
      integer nnodex, j, i, nnodey, numb, jmid, jcut
      complex zi

      namelist/fieldin/ibackground, xminz, xmaxz, yminz, ymaxz, logmax, &
         ipage, &
         numb, ycut, i_psi1, i_psi2, i_psi3, i_psi4, i_psi5, i_psi6



!--set default values of input data:

      ibackground = 0

      xminz = .58
      xmaxz = .72
      ymaxz = .30
      yminz = -99.99

      logmax = 2.0
      ipage = 2
      numb = 20
      ycut = -0.0

      i_psi1 = 2
      i_psi2 = 10
      i_psi3 = 20
      i_psi4 = 30
      i_psi5 = 40
      i_psi6 = 50

      i_psi_array(1) = 2
      i_psi_array(2) = 10
      i_psi_array(3) = 20
      i_psi_array(4) = 30
      i_psi_array(5) = 40
      i_psi_array(6) = 50


      if(yminz .eq. -99.99)yminz = -ymaxz


      open(unit=38,file='out38',status='old',form='formatted')

      open(unit=140,file='Efield_2D.vtk',status='unknown', &
                                             form='formatted')
      open(unit=141,file='E_kicks_2D.vtk',status='unknown', &
                                             form='formatted')
      open(unit=142,file='Bql_avg_2D.vtk',status='unknown', &
                                             form='formatted')

!      open(unit=53,file='movie_ex',status='unknown',form='formatted')
!      open(unit=58,file='movie_ey',status='unknown',form='formatted')
!      open(unit=52,file='movie_ez',status='unknown',form='formatted')
!      open(unit=59,file='movie_eb',status='unknown',form='formatted')
!      open(unit=60,file='movie_ealpha',status='unknown',
!     .   form='formatted')

      open(unit=51,file='rho',status='unknown',form='formatted')
!      open(unit=61,file='acold',status='unknown',form='formatted')


      open(unit=66,file='bharvey',status='unknown',form='formatted')


      open(unit=57,file='swain',status='unknown',form='formatted')
      open(unit=62,file='murakami',status='unknown', form='formatted')

      open(unit=65,file='mchoi',status='unknown', form='formatted')
      open(unit=67,file='mchoi2',status='unknown', form='formatted')
      open(unit=69,file='mchoi3',status='unknown', form='formatted')
      open(unit=68,file='zowens',status='unknown', form='formatted')


      open(unit=64,file='E_lab_frame',status='unknown',form='formatted')


      zi = cmplx(0.0, 1.0)
      eps0 = 8.85e-12
      xmu0 = 1.26e-06
      clight = 1./sqrt(eps0 * xmu0)
      xk0 = omgrf / clight
      q = 1.6e-19

      fmid1 = 0.0
      fmid2 = 0.0
      fmid3 = 0.0
      fmid4 = 0.0
      fmid5 = 0.0
      fmid6 = 0.0

      read(38, 309) nnodex, nnodey, jmid
      read(38, 310) rhoplasm, prfin
      read(38, 310) (x(i), i = 1, nnodex)
      read(38, 310) (y(j), j = 1, nnodey)
      read(38, 310) (capr(i), i = 1, nnodex)

      xmax = x(nnodex) - x(1)
      ymax = y(nnodex) - y(1)


      read(38, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((theta(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((xn(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((xkte(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((xkti(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((xkperp_cold(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((acold(i, j), i = 1, nnodex), j = 1, nnodey)

      read (38, 310) ((reomg1a(i, j), i = 1, nnodex), j = 1, nnodey)
      read (38, 310) ((reomg2a(i, j), i = 1, nnodex), j = 1, nnodey)
      read (38, 310) ((reomg3a(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((xjy(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((xjx(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((ex(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((ey(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((ez(i, j), i = 1, nnodex), j = 1, nnodey)



      write(51,309) nnodex, nnodey
      write(51,310) rhoplasm
      write(51,310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)



      read(38, 309) nkx1, nkx2
      read(38, 309) nky1, nky2

      nkxplt = 2 * nkx2
      nkyplt = 2 * nky2

      read (38, 310) (xkxsav(n), n = 1, nkxplt)
      read (38, 310) (xkysav(m), m = 1, nkyplt)

      read(38, 310) ((exkmod(n, m), n = 1, nkxplt), m = 1, nkyplt)
      read(38, 310) ((eykmod(n, m), n = 1, nkxplt), m = 1, nkyplt)
      read(38, 310) ((ezkmod(n, m), n = 1, nkxplt), m = 1, nkyplt)


    ! catch for multiple nPhi setting prfin=0
    if ( prfin == 0 ) then
            prfinMOD = 1 
    else 
            prfinMOD = prfin
    endif


!
!--calculate log exkmod, eykmod, ezkmod
!
      do n = 1, nkxplt
         do m = 1, nkyplt


            exklog(n, m) = 0.0
            eyklog(n, m) = 0.0
            ezklog(n, m) = 0.0

            if(abs(exkmod(n,m)) .ne. 0.0) &
               exklog(n,m)=alog10(abs(exkmod(n,m)))
            if(abs(eykmod(n,m)) .ne. 0.0) &
               eyklog(n,m)=alog10(abs(eykmod(n,m)))
            if(abs(ezkmod(n,m)) .ne. 0.0) &
               ezklog(n,m)=alog10(abs(ezkmod(n,m)))
         end do
      end do




      call a2dmnmx_r4(exklog, nkpltdim, mkpltdim, nkxplt, nkyplt, &
          exklogmin, exklogmax)

      call a2dmnmx_r4(eyklog, nkpltdim, mkpltdim, nkxplt, nkyplt, &
          eyklogmin, eyklogmax)

      call a2dmnmx_r4(ezklog, nkpltdim, mkpltdim, nkxplt, nkyplt, &
          ezklogmin, ezklogmax)


!--   Normalize log's to a maximum value of logmax
      do n = 1, nkxplt
         do m = 1, nkyplt
            if(exklog(n,m) .ne. 0.0) &
               exklog(n,m) = exklog(n,m) + (logmax - exklogmax)
            if(eyklog(n,m) .ne. 0.0) &
               eyklog(n,m) = eyklog(n,m) + (logmax - eyklogmax)
            if(ezklog(n,m) .ne. 0.0) &
               ezklog(n,m) = ezklog(n,m) + (logmax - ezklogmax)
         end do
      end do


      read(38, 310) ((redotje(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((redotj1(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((redotj2(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((redotj3(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((redotjt(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((bmod(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((qsafety(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((ealpha(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((ebeta(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((eb(i, j), i = 1, nnodex), j = 1, nnodey)



      read(38, 310) ((btau(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((bzeta(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((spx(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((spy(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((spz(i, j), i = 1, nnodex), j = 1, nnodey)



      read(38, 310) ((exk(n, m), n = nkx1, nkx2), m = nky1, nky2)
      read(38, 310) ((eyk(n, m), n = nkx1, nkx2), m = nky1, nky2)
      read(38, 310) ((ezk(n, m), n = nkx1, nkx2), m = nky1, nky2)

!      write(6, *)"ealphak(32,32) = ", exk(32,32)
!      write(6, *)"ebetak(32,32)  = ", eyk(32,32)
!      write(6, *)"ebk(32,32)     = ", ezk(32,32)

      read(38, 310) ((wdote(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((wdoti1(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((wdoti2(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((wdott(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 309) nnoderho
      read(38, 310) (rhon(n), n = 1, nnoderho)
      read(38, 310) (gradprlb2_avg(n), n = 1, nnoderho)

      read(38, 309) mnodetheta
      read(38, 310) (thetam(m), m = 1, mnodetheta)

      read(38, 310) (dvol(n), n = 1, nnoderho)
      read(38, 310) (volume(n), n = 1, nnoderho)

      read(38, 310) (xnavg(n), n = 1, nnoderho)
      read(38, 310) (wdoteavg(n), n = 1, nnoderho)
      read(38, 310) (wdoti1avg(n), n = 1, nnoderho)
      read(38, 310) (wdoti2avg(n), n = 1, nnoderho)
      read(38, 310) (vyavg(n), n = 1, nnoderho)
      read(38, 310) (dvydrho(n), n = 1, nnoderho)

      read(38, 310) (redotjeavg(n), n = 1, nnoderho)
      read(38, 310) (redotj1avg(n), n = 1, nnoderho)
      read(38, 310) (redotj2avg(n), n = 1, nnoderho)
      read(38, 310) (redotj3avg(n), n = 1, nnoderho)
      read(38, 310) (redotjtavg(n), n = 1, nnoderho)

      read(38, 310) ((redotjs(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((wdoti3(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) (wdoti3avg(n), n = 1, nnoderho)



      read(38, 310) ((fz0e(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((fz0i1(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((fz0i2(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((fz0i3(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((fz0(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) (fz0eavg(n), n = 1, nnoderho)
      read(38, 310) (fz0i1avg(n), n = 1, nnoderho)
      read(38, 310) (fz0i2avg(n), n = 1, nnoderho)
      read(38, 310) (fz0i3avg(n), n = 1, nnoderho)
      read(38, 310) (fz0avg(n), n = 1, nnoderho)

      read(38, 310) (xjprlavg(n), n = 1, nnoderho)

      read(38, 310) (xjprl_int(n), n = 1, nnoderho)
      read(38, 310) (fz0_int(n), n = 1, nnoderho)



      read(38, 310) (xjhat(n), n = 1, nnoderho)
      read(38, 310) (xkhat(n), n = 1, nnoderho)

      read(38, 310) (gpsi_avg(n), n = 1, nnoderho)
      read(38, 310) ((gpsi(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) (kpsi_avg(n), n = 1, nnoderho)
      read(38, 310) ((kpsi(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) (muhat_avg(n), n = 1, nnoderho)
      read(38, 310) ((muhat(i, j), i = 1, nnodex), j = 1, nnodey)



      read(38, 310) (nu_star_avg(n), n = 1, nnoderho)
      read(38, 310) ((nu_star(i, j), i = 1, nnodex), j = 1, nnodey)



      read(38, 310) (ipsi_avg(n), n = 1, nnoderho)
      read(38, 310) ((ipsi(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((omgexb(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((uzeta(i,j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((utheta(i,j), i = 1, nnodex), j = 1, nnodey)


      read(38, 310) (xn1avg(n), n = 1, nnoderho)
      read(38, 310) (xn2avg(n), n = 1, nnoderho)
      read(38, 310) (xn3avg(n), n = 1, nnoderho)
      read(38, 310) (xna_sloavg(n), n = 1, nnoderho)
      read(38, 310) (xkteavg(n), n = 1, nnoderho)
      read(38, 310) (xktiavg(n), n = 1, nnoderho)
      read(38, 310) (xkti2avg(n), n = 1, nnoderho)
      read(38, 310) (xkti3avg(n), n = 1, nnoderho)


      read(38, 310) (redotjsavg(n), n = 1, nnoderho)

      read(38, 310) ((bmod_mid(i,j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) (bmod_midavg(n), n = 1, nnoderho)

      read(38, 310) ((capr_bpol_mid2(i, j), i = 1,nnodex), &
                                                          j = 1,nnodey)
      read(38, 310) (capr_bpol_midavg(n), n = 1, nnoderho)



      read(38, 310) ((redotj4(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((redotj5(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((redotj6(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) (redotj4avg(n), n = 1, nnoderho)
      read(38, 310) (redotj5avg(n), n = 1, nnoderho)
      read(38, 310) (redotj6avg(n), n = 1, nnoderho)

      read(38, 310) ((wdoti4(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((wdoti5(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((wdoti6(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) (wdoti4avg(n), n = 1, nnoderho)
      read(38, 310) (wdoti5avg(n), n = 1, nnoderho)
      read(38, 310) (wdoti6avg(n), n = 1, nnoderho)

      read(38, 310) ((fz0i4(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((fz0i5(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((fz0i6(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) (fz0i4avg(n), n = 1, nnoderho)
      read(38, 310) (fz0i5avg(n), n = 1, nnoderho)
      read(38, 310) (fz0i6avg(n), n = 1, nnoderho)

      read(38, 310) (xn4avg(n), n = 1, nnoderho)
      read(38, 310) (xn5avg(n), n = 1, nnoderho)
      read(38, 310) (xn6avg(n), n = 1, nnoderho)

      read(38, 310) (xkti4avg(n), n = 1, nnoderho)
      read(38, 310) (xkti5avg(n), n = 1, nnoderho)
      read(38, 310) (xkti6avg(n), n = 1, nnoderho)

      read(38, 310) (wdoteavg_int(n),  n = 1, nnoderho)
      read(38, 310) (wdoti1avg_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti2avg_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti3avg_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti4avg_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti5avg_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti6avg_int(n), n = 1, nnoderho)

      read(38, 310) (wdote_ql_int(n),  n = 1, nnoderho)
      read(38, 310) (wdoti1_ql_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti2_ql_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti3_ql_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti4_ql_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti5_ql_int(n), n = 1, nnoderho)
      read(38, 310) (wdoti6_ql_int(n), n = 1, nnoderho)

      read(38, 310) (redotjeavg_int(n), n = 1, nnoderho)
      read(38, 310) (redotj1avg_int(n), n = 1, nnoderho)
      read(38, 310) (redotj2avg_int(n), n = 1, nnoderho)
      read(38, 310) (redotj3avg_int(n), n = 1, nnoderho)
      read(38, 310) (redotj4avg_int(n), n = 1, nnoderho)
      read(38, 310) (redotj5avg_int(n), n = 1, nnoderho)
      read(38, 310) (redotj6avg_int(n), n = 1, nnoderho)




!     -----------------------
!     Calculate bratio(i,j)
!     -----------------------
      do i = 1, nnodex
         do j = 1, nnodey
          bratio(i,j) = bmod_mid(i,j) / bmod(i,j)
          if (bratio(i,j) .gt. 1.0) bratio(i,j) = 1.0
       end do
      end do



      nnoderho_half = nnoderho - 1

      do n = 1, nnoderho_half
       rhon_half(n) = (rhon(n) + rhon(n+1)) / 2.
      end do

      do n = 1, nnoderho
       rhon_save(n) = rhon(n)
       rhon_half_save(n) = rhon_half(n)
      end do


      write(57, 309) nnoderho
      do n = 1, nnoderho
         write(57,1312)n, rhon(n), redotjeavg(n), redotj1avg(n), &
                    redotj2avg(n), redotj3avg(n), redotj4avg(n), &
                    redotj5avg(n), xjprlavg(n)
      end do


      write(62, 309) nnoderho
      do n = 1, nnoderho
         write(62,1312)n, rhon(n), wdoteavg(n), wdoti1avg(n), &
                    wdoti2avg(n), wdoti3avg(n), xjprlavg(n)
      end do




      do i = 1, nnodex
         xnmid(i) = xn(i, jmid)
         xktemid(i) = xkte(i, jmid) / q
         xktimid(i) = xkti(i, jmid) / q
         qmid(i) = qsafety(i, jmid)
         xiotamid(i) = 1.0 / qmid(i)
         bpmid(i) = btau(i, jmid)
      end do



      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8


      ncollin=ncyan
      ncolln2=nblueviolet
      ncolln3=ngreen
      ncolln4=naqua
      ncolln5=nyellow
      ncolbrd=nwheat
      ncollab=ncyan

      ncolbox=nred
      ncolbrd=nwheat
!     ncolbox=nblack
!     ncolbrd=nblack

!     ncolion=nblue
!     ncolelec=nred
      ncolion=naqua
      ncolelec=nyellow

      if (ibackground.eq.1)then
         ncolbox=nred
         ncolbrd=nwheat
      end if

      if (ibackground.eq.0)then
         ncolbox=nblack
         ncolbrd=nblack
      end if

! Open graphics device

!      IER = PGBEG(0, 'aorsa2d.ps/cps', 2, 2)
      IER = PGBEG(0, 'aorsa2d.ps/vcps', 1, 1)

      IF (IER.NE.1) STOP

!      if (pgopen('aorsa2d.cps/cps') .lt. 1) stop

      titx = 'R (m)'
      tity = 'Z (m)'

       !--------------------------------------------
       !DLG - write nc file for nice plotting

        write (*,*) 'WRITING output/plotData.nc ...'
        ncFileName = 'output/plotData.nc'

            call check ( &
       nf90_create ( ncFileName, nf90_clobber, nc_id ) )

            call check ( &
       nf90_def_dim ( nc_id, "nR", nnodex, nR_id ) )
            call check ( &
       nf90_def_dim ( nc_id, "nz", nnodey, nz_id ) )
            call check ( &
       nf90_def_dim ( nc_id, "scalar", 1, scalar_id ) )
            call check ( &
       nf90_def_dim ( nc_id, "nRho", nnoderho_half, nRho_id ) )


            call check ( &
       nf90_def_var ( nc_id, "capR", NF90_REAL, &
       (/ nR_id /), capR_id ) )
            call check ( &
       nf90_def_var ( nc_id, "zLoc", NF90_REAL, &
       (/ nz_id /), y_id ) )
            call check ( &
       nf90_def_var ( nc_id, "rho", NF90_REAL, &
       (/ nRho_id /), rho_id ) )
            call check ( &
       nf90_def_var ( nc_id, "wdote", NF90_REAL, &
       (/ nRho_id /), wdote_id ) )
            call check ( &
       nf90_def_var ( nc_id, "wdoti1", NF90_REAL, &
       (/ nRho_id /), wdoti1_id ) )
            call check ( &
       nf90_def_var ( nc_id, "wdoti2", NF90_REAL, &
       (/ nRho_id /), wdoti2_id ) )
            call check ( &
       nf90_def_var ( nc_id, "wdote_rz", NF90_REAL, &
       (/ nR_id, nz_id /), wdote_rz_id ) )
            call check ( &
       nf90_def_var ( nc_id, "wdoti1_rz", NF90_REAL, &
       (/ nR_id, nz_id /), wdoti1_rz_id ) )
            call check ( &
       nf90_def_var ( nc_id, "wdoti2_rz", NF90_REAL, &
       (/ nR_id, nz_id /), wdoti2_rz_id ) )
            call check ( &
       nf90_def_var ( nc_id, "jeDotE", NF90_REAL, &
       (/ nR_id, nz_id /), jeDotE_id ) )
            call check ( &
       nf90_def_var ( nc_id, "eAlpha", NF90_REAL, &
       (/ nR_id, nz_id /), eAlpha_id ) )
            call check ( &
       nf90_def_var ( nc_id, "eBeta", NF90_REAL, &
       (/ nR_id, nz_id /), eBeta_id ) )
            call check ( &
       nf90_def_var ( nc_id, "eParallel", NF90_REAL, &
       (/ nR_id, nz_id /), eParallel_id ) )
            call check ( &
       nf90_def_var ( nc_id, "ex", NF90_REAL, &
       (/ nR_id, nz_id /), ex_id ) )
            call check ( &
       nf90_def_var ( nc_id, "ey", NF90_REAL, &
       (/ nR_id, nz_id /), ey_id ) )
            call check ( &
       nf90_def_var ( nc_id, "ez", NF90_REAL, &
       (/ nR_id, nz_id /), ez_id ) )




       call check ( nf90_def_var ( nc_id, "ePlus_real", NF90_REAL, &
       (/ nR_id, nz_id /), ePlus_real_id ) )
       call check ( nf90_def_var ( nc_id, "ePlus_imag", NF90_REAL, &
       (/ nR_id, nz_id /), ePlus_imag_id ) )
       call check ( nf90_def_var ( nc_id, "eMinu_real", NF90_REAL, &
       (/ nR_id, nz_id /), eMinu_real_id ) )
       call check ( nf90_def_var ( nc_id, "eMinu_imag", NF90_REAL, &
       (/ nR_id, nz_id /), eMinu_imag_id ) )
 
       call check ( nf90_def_var ( nc_id, "bx_wave_real", NF90_REAL, &
       (/ nR_id, nz_id /), bx_wave_real_id ) )
       call check ( nf90_def_var ( nc_id, "bx_wave_imag", NF90_REAL, &
       (/ nR_id, nz_id /), bx_wave_imag_id ) )
       call check ( nf90_def_var ( nc_id, "bz_wave_real", NF90_REAL, &
       (/ nR_id, nz_id /), bz_wave_real_id ) )
       call check ( nf90_def_var ( nc_id, "bz_wave_imag", NF90_REAL, &
       (/ nR_id, nz_id /), bz_wave_imag_id ) )
 
       call check ( nf90_def_var ( nc_id, "density", NF90_REAL, &
       (/ nR_id, nz_id /), density_id ) )
       call check ( nf90_def_var ( nc_id, "janty", NF90_REAL, &
       (/ nR_id, nz_id /), janty_id ) )
       call check ( nf90_def_var ( nc_id, "jantx", NF90_REAL, &
       (/ nR_id, nz_id /), jantx_id ) )
 
       call check ( nf90_def_var ( nc_id, "mask", NF90_INT, &
       (/ nR_id, nz_id /), mask_id ) )
 
       call check ( nf90_def_var ( nc_id, "xkperp_real", NF90_REAL, &
       (/ nR_id, nz_id /), xkperp_real_id ) )
       call check ( nf90_def_var ( nc_id, "xkperp_imag", NF90_REAL, &
       (/ nR_id, nz_id /), xkperp_imag_id ) )
 
       call check ( nf90_def_var ( nc_id, "bmod", NF90_REAL, &
       (/ nR_id, nz_id /), bmod_id ) )
 
            call check ( &
       nf90_def_var ( nc_id, "pscale", NF90_REAL, &
       scalar_id, pscale_id ) )

            call check ( nf90_enddef ( nc_id ) )

            call check ( &
       nf90_put_var ( nc_id, capR_id, capR(1:nnodex) ) )
            call check ( &
       nf90_put_var ( nc_id, y_id, y(1:nnodey) ) )
            call check ( &
       nf90_put_var ( nc_id, rho_id, &
       rhon_half(1:nnoderho_half) ) )
            call check ( &
       nf90_put_var ( nc_id, wdote_id, &
       wdoteavg(1:nnoderho_half) ) )
            call check ( &
       nf90_put_var ( nc_id, wdoti1_id, &
       wdoti1avg(1:nnoderho_half) ) )
            call check ( &
       nf90_put_var ( nc_id, wdoti2_id, &
       wdoti2avg(1:nnoderho_half) ) )
            call check ( &
       nf90_put_var ( nc_id, wdote_rz_id, &
       wdote(1:nnodex,1:nnodey) ) )
            call check ( &
       nf90_put_var ( nc_id, wdoti1_rz_id, &
       wdoti1(1:nnodex,1:nnodey) ) )
            call check ( &
       nf90_put_var ( nc_id, wdoti2_rz_id, &
       wdoti2(1:nnodex,1:nnodey) ) )
            call check ( &
       nf90_put_var ( nc_id, pscale_id, &
       pscale ) )

        call check ( nf90_put_var ( nc_id, mask_id, mask(1:nnodex,1:nnodey) ) )




        !
        !--------------------------------------
!
!--plot plasma profiles
!
      call PGSCH (1.5)


      title = 'bmod_mid surfaces'
      call ezconc(capr, y, bmod_mid, ff, nnodex, nnodey, 50, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)


      title = 'bratio surfaces'
      call ezconc(capr, y, bratio, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      title = 'ipsi'
      call ezconc(capr, y, ipsi, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
          nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      title= 'I(psi) '
      titll= 'I(psi) m-T'
      titlr='       '
      call ezplot0(title, titll, titlr, rhon, ipsi_avg, &
          nnoderho, nrhomax)


      title = 'capr_bpol_mid2 surfaces'
      call ezconc(capr, y, capr_bpol_mid2, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)


      title= 'Flux average bmod_mid'
      titll= 'bmod_mid (T)'
      titlr='       '

      call ezplot0(title, titll, titlr, rhon_half, bmod_midavg, &
          nnoderho_half, nrhomax)

      title= 'parallel gradient of B'
      titll= 'gradprlb_avg (T/m)'
      titlr='       '

      call ezplot0(title, titll, titlr, rhon_half, gradprlb2_avg, &
          nnoderho_half, nrhomax)


      title= 'Flux average capr_bpol'
      titll= 'capr_bpol (mT)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, capr_bpol_midavg, &
          nnoderho_half, nrhomax)

      call a1mnmx(capr, nxmx, nnodex, caprmin, caprmax)
      caprminp = caprmin * 1.01
      caprmaxp = caprmax * .99
      call a2mnmx(xnmid, nxmx, nnodex, &
         capr, caprminp, caprmaxp, xnmin, xnmax)

      xnmax = 2.0 * xnmax

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (caprmin, caprmax, xnmin, xnmax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BNST', 0.0, 0)
      call pgmtxt('b', 3.2, 0.5, 0.5, 'R (m)')
      call pgmtxt('t', 2.0, 0.5, 0.5, 'Plasma profiles')

      CALL PGSCI(nblue)
      call pgmtxt('l', 2.0, 0.5, 0.5, ' n (m-3)')
      call pgline(nnodex, capr, xnmid)

      call a1mnmx(xktemid, nxmx, nnodex, temin, temax)
      call a1mnmx(xktimid, nxmx, nnodex, timin, timax)
      tmax = max(timax, temax)
      tmin = 0.0

      CALL PGSWIN (caprmin, caprmax, tmin, tmax)
      CALL PGSCI(nblack)
      CALL PGBOX  (' ', 0.0, 0, 'CMST', 0.0, 0)

      CALL PGSCI(nred)
      call pgsls(2)
      call pgline(nnodex, capr, xktemid)
      call PGMTXT ('r', 2.0, 0.5, 0.5, 'T (eV)')

      CALL PGSCI(ngreen)
      call pgsls(3)
      call pgline(nnodex, capr, xktimid)

      CALL PGSCI(nblack)
      call pgsls(1)


      title= 'Flux surface average density'
      titll= 'density (m-3)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xnavg, &
         nnoderho_half, nrhomax)

      title= 'Differential volume element'
      titll= 'dVol (m3)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, dvol, &
                                                nnoderho_half, nrhomax)


      title= 'Integrated volume'
      titll= 'volume (m3)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon, volume, nnoderho, nrhomax)



      do n = 1, nnoderho
         xkteavg(n) = xkteavg(n) / q
         xktiavg(n) = xktiavg(n) / q
         xkti2avg(n) = xkti2avg(n) / q
         xkti3avg(n) = xkti3avg(n) / q
      end do

      title= 'Flux average electron temperature'
      titll= 'kTe (eV)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xkteavg, &
          nnoderho_half, nrhomax)

      title= 'Flux average ion temperature'
      titll= 'kTi (eV)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xktiavg, &
          nnoderho_half, nrhomax)


      title= 'Flux average beam density'
      titll= 'density (m-3)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xn3avg, &
         nnoderho_half, nrhomax)


      title= 'Flux average beam temperature'
      titll= 'kTi_fast (eV)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xkti3avg, &
          nnoderho_half, nrhomax)


      title= 'Flux average minority density'
      titll= 'density (m-3)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xn2avg, &
         nnoderho_half, nrhomax)

      title= 'Flux average minority temperature'
      titll= 'kT2 (eV)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xkti2avg, &
          nnoderho_half, nrhomax)

      title = 'Contour plot of density'
      call ezconc(capr, y, xn, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      call check ( nf90_put_var ( nc_id, density_id, &
      xn(1:nnodex,1:nnodey) ) )
  
      title= 'Flux surface average redotj'
      titll= 'redotj (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon_half, &
         redotj1avg, redotj2avg, redotj3avg, redotj4avg, redotj5avg, &
         redotj6avg, redotjeavg,  nnoderho_half, nrhomax)

      title= 'Integrated redotj'
      titll= 'P (Watts)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon, &
         redotj1avg_int, redotj2avg_int, redotj3avg_int, redotj4avg_int, &
         redotj5avg_int, redotj6avg_int, redotjeavg_int, &
         nnoderho, nrhomax)

      write(15, *)
      write(15, *) 'Flux surface driven current'
      write(15, *)
      write(15, *) '        n      rho       J (A/m2)     I (A)'
      write(15, *)

!      do n = 1, nnoderho_half
!       write (15, 1312) n, rhon_half(n), xjprlavg(n), xjprl_int(n)
!            write (6, 1312) n, rhon_half(n), xjprlavg(n), xjprl_int(n)
!      end do

      title= 'Flux surface driven current'
      titll= 'xjprl (Amps/m2)'
      titlr= 'I (Amps)'
      titlb= 'rho'

      call ezplot2(title, titll, titlr, titlb, rhon_half, xjprlavg, &
          xjprl_int, nnoderho_half, nrhomax)

      if(xjprl_int(nnoderho) .lt. 0.0) then


         titll= '-xjprl (Amps/m2)'
         titlr= '-I (Amps)'

         xjprlavg = - xjprlavg
         xjprl_int = - xjprl_int

         call ezplot2(title, titll, titlr, titlb, rhon_half, xjprlavg, &
            xjprl_int, nnoderho_half, nrhomax)

         titll= 'xjprl (MA/m2/MW)'
       xjprlavg = xjprlavg  / prfinMOD

         call ezplot1(title, titll, titlr, rhon_half, xjprlavg, &
            nnoderho_half, nrhomax)

         titll= 'I (kA)'
       xjprl_int  = xjprl_int / 1.0e+03
         call ezplot1(title, titll, titlr, rhon_half, xjprl_int, &
            nnoderho_half, nrhomax)

      end if

      if(xjprl_int(nnoderho) .gt. 0.0) then

         titll= 'xjprl (Amps/m2)'
         titlr= 'I (Amps)'

         xjprlavg = xjprlavg
         xjprl_int = xjprl_int

         titll= 'xjprl (MAmps/m2/MW)'
       xjprlavg = xjprlavg  / prfinMOD

         call ezplot1(title, titll, titlr, rhon_half, xjprlavg, &
            nnoderho_half, nrhomax)

         titll= 'I (kA)'
       xjprl_int  = xjprl_int / 1.0e+03
         call ezplot1(title, titll, titlr, rhon_half, xjprl_int, &
            nnoderho_half, nrhomax)

      end if

      title= 'Flux surface average wdot'
      titll= 'wdot (Watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon_half, &
         wdoti1avg, wdoti2avg, wdoti3avg, wdoti4avg, wdoti5avg, &
         wdoti6avg, wdoteavg, nnoderho_half, nrhomax)

      title= 'Integrated wdot'
      titll= 'P (Watts)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon, &
         wdoti1avg_int, wdoti2avg_int, wdoti3avg_int, wdoti4avg_int, &
         wdoti5avg_int, wdoti6avg_int, wdoteavg_int, nnoderho, nrhomax)

      title= 'Flux average Toroidal force'
      titll= 'force (Nt/m3)'
      titlr='       '

      call ezplot7(title, titll, titlr, titlb, rhon_half, &
          fz0i1avg, fz0i2avg, fz0i3avg, fz0i4avg, fz0i5avg, &
          fz0i6avg, fz0eavg, nnoderho_half, nrhomax)

      title= 'Flux average toroidal force'
      titll= 'force (Nt/m3)'
      titlr= 'integrated force (Nt/m)'
      titlb= 'rho'

      call ezplot2(title, titll, titlr, titlb, rhon_half, fz0avg, &
          fz0_int, nnoderho_half, nrhomax)

      title= 'Flux surface average Vy'
      titll= 'Vy (m/s)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, vyavg, &
          nnoderho_half, nrhomax)

!     -------------------
!     plot G(psi) in 2-D
!     -------------------

      title = 'Toroidal angular speed G(psi)'
      call ezconc(capr, y, gpsi, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

!     -----------------------
!     plot G(psi) in midplane
!     -----------------------

      do i = 1, nnodex
         fmidre(i) = gpsi(i, jmid)
      end do

      title = 'Toroidal angular speed G(psi)'
      titll = 'Gpsi (sec-1)'
      titlr = '      '
      titlb = 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmidre, &
         nnodex, nxmx)

!     -------------------
!     plot K(psi) in 2-D
!     -------------------

      title = 'Poloidal speed / B_pol = K(psi)'
      call ezconc(capr, y, kpsi, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

!     -----------------------
!     plot K(psi) in midplane
!     -----------------------

      do i = 1, nnodex
         fmidre(i) = kpsi(i, jmid)
      end do

      title = 'Poloidal speed / B_pol = K(psi)'
      titll = 'K(psi) (m/s/T)'
      titlr = '      '
      titlb = 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmidre, &
         nnodex, nxmx)


!     -----------------
!     plot uzeta in 2-D
!     -----------------

      title = 'Toroidal velocity (uzeta)'
      call ezconc(capr, y, uzeta, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

!     -----------------------
!     plot uzeta in midplane
!     -----------------------

      do i = 1, nnodex
         fmidre(i) = uzeta(i, jmid)
      end do

      title = 'Toroidal velocity (uzeta)'
      titll = 'uzeta (m/sec)'
      titlr = '  '
      titlb = 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmidre, &
         nnodex, nxmx)




!     -----------------
!     plot utheta in 2-D
!     -----------------

      title = 'Poloidal velocity (utheta)'
      call ezconc(capr, y, utheta, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

!     -----------------------
!     plot utheta in midplane
!     -----------------------
      do i = 1, nnodex
         fmidre(i) = utheta(i, jmid)
      end do

      title = 'Poloidal velocity (utheta)'
      titll = 'utheta (m/sec)'
      titlr = '  '
      titlb = 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmidre, &
         nnodex, nxmx)









!     -------------------
!     plot omgexb in 2-D
!     -------------------

      title = 'ExB shearing rate (omgexb)'
      call ezconc(capr, y, omgexb, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

!     -----------------------
!     plot omgexb in midplane
!     -----------------------

      do i = 1, nnodex
         fmidre(i) = omgexb(i, jmid)
      end do

      title = 'ExB shearing rate (omgexb)'
      titll = 'omgexb (sec-1)'
      titlr = '      '
      titlb = 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmidre, &
         nnodex, nxmx)



      title= 'xjhat(rho)'
      titll= 'xjhat (kg m-2 s-1)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xjhat, &
          nnoderho_half, nrhomax)


      title= 'xkhat(rho)'
      titll= 'xkhat (kg m-2 s-1)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, xkhat, &
          nnoderho_half, nrhomax)

      title= 'poloidal speed <K(rho)>'
      titll= 'K(rho) (m s-1 T-1)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, kpsi_avg, &
          nnoderho_half, nrhomax)


      title= 'toroidal angular speed <G(rho)>'
      titll= 'G(rho) (s-1)'
      titlr='       '
      call ezplot1(title, titll, titlr, rhon_half, gpsi_avg, &
          nnoderho_half, nrhomax)

      title= 'normalized viscosity <mu_hat>'
      titll= 'mu_hat (kg/m**3/s-1)'
      titlr='       '
      call ezplot0(title, titll, titlr, rhon_half, muhat_avg, &
          nnoderho_half, nrhomax)


      title= 'collisionality (nu_star)'
      titll= 'nu_star'
      titlr='       '
      call ezplot0(title, titll, titlr, rhon_half, nu_star_avg, &
          nnoderho_half, nrhomax)




!     -------------------
!     plot nu_star in 2-D
!     -------------------

      title = 'collisionality (nu_star)'
      call ezconc(capr, y, nu_star, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)



!
!--plot edotj in midplane
!
      do i = 1, nnodex
         fmidre(i) = redotje(i,jmid)
         fmidim(i) = redotj2(i,jmid)
         fmid1(i)  = redotj1(i,jmid)
         fmid2(i)  = redotj2(i,jmid)
         fmid3(i)  = redotj3(i,jmid)
         fmid4(i)  = redotj4(i,jmid)
       fmid5(i)  = redotj5(i,jmid)
         fmidt(i)  = redotjt(i,jmid)
      end do

      title= 'J dot E in midplane'
      titll= 'Re edotj (W/m3)'
      titlr= ' '
      titlb= 'R (m)'

      call ezplot7(title, titll, titlr, titlb, capr, &
          fmid1, fmid2, fmid3, fmid4, fmid5, fmid6, fmidre, &
          nnodex, nxmx)

!
!--plot wdot at jcut
!

      dy = y(2) - y(1)
      jcut = int((ycut - y(1) - dy / 2.0) / dy) + 1

      if (ycut .eq. 0.0)jcut = jmid

      do i = 1, nnodex
         fmidre(i) = wdote(i, jcut)
         fmidim(i) = wdoti2(i, jcut)
         fmid1(i)  = wdoti1(i, jcut)
         fmid2(i)  = wdoti2(i, jcut)
         fmid3(i)  = wdoti3(i, jcut)
       fmid4(i)  = wdoti4(i, jcut)
       fmid5(i)  = wdoti5(i, jcut)
         fmidt(i)  = wdott(i, jcut)
      end do

      title= 'Midplane Wdot'
      titll= 'Wdot (W/m3)'
      titlr= ' '
      titlb= 'R (m)'

      call ezplot70(title, titll, titlr, titlb, capr, &
          fmid1, fmid2, fmid3, fmid4, fmid5, fmid6, fmidre, &
          nnodex, nxmx)


      do i = 1, nnodex
         fmid1(i) = bmod(i, jcut)
      end do

      title= 'Midplane Mod B'
      titll= 'Mod B (T)'
      titlr= ' '
      titlb= 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmid1, &
         nnodex, nxmx)



!
!--plot F_toroidal at y_cut
!
      dy = y(2) - y(1)
      jcut = int((ycut - y(1) - dy / 2.0) / dy) + 1

      if (ycut .eq. 0.0)jcut = jmid

      do i = 1, nnodex
         fmidre(i) = fz0e(i, jcut)
         fmidim(i) = fz0i2(i, jcut)
         fmid1(i)  = fz0i1(i, jcut)
         fmid2(i)  = fz0i2(i, jcut)
         fmid3(i)  = fz0i3(i, jcut)
       fmid4(i)  = fz0i4(i, jcut)
       fmid5(i)  = fz0i5(i, jcut)
         fmidt(i)  = fz0(i, jcut)
      end do

      title= 'Toroidal force'
      titll= 'F_toroidal (Nt/m3)'
      titlr= ' '
      titlb= 'R (m)'

      call ezplot7(title, titll, titlr, titlb, capr, &
          fmid1, fmid2, fmid3, fmid4, fmid5, fmid6, fmidre, &
          nnodex, nxmx)



!--Contour plots using 16 contour levels:

      if (ibackground.eq.1)then
!        black background
         ncollin = ncyan
         ncolln2 = nblueviolet
      end if

      if (ibackground.eq.0)then
!        white background
         ncollin = nblue
         ncolln2 = nyellow
      end if

      titx = 'R (m)'
      tity = 'Z (m)'

!
!--plot q profile in midplane
!
      call a2mnmx(qmid, nxmx, nnodex, &
         capr, caprminp, caprmaxp, qmin, qmax)
      qmin = 0.0
      qmax = 5.0

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (caprmin, caprmax, qmin, qmax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BNST', 0.0, 0)
      call pgmtxt('b', 3.2, 0.5, 0.5, 'R (m)')
      call pgmtxt('t', 2.0, 0.5, 0.5, 'q profile in midplane')

      CALL PGSCI(nred)
      call pgmtxt('l', 2.0, 0.5, 0.5, 'q midplane')
      call pgline(nnodex, capr, qmid)

      CALL PGSCI(nblue)
      call pgsls(2)
      call pgline(nnodex, capr, xiotamid)

      call a1mnmx(bpmid, nxmx, nnodex, bpmin, bpmax)
      bpmin = 0.0

      CALL PGSWIN (caprmin, caprmax, bpmin, bpmax)
      CALL PGSCI(nblack)
      CALL PGBOX  (' ', 0.0, 0, 'CMST', 0.0, 0)

      CALL PGSCI(ngreen)
      call pgsls(3)
      call pgline(nnodex, capr, bpmid)
      call PGMTXT ('r', 2.0, 0.5, 0.5, 'btau')

      CALL PGSCI(nblack)
      call pgsls(1)


      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = rho(i,j)**2
          psi(i,j) = rho(i,j)**2
         end do
      end do

      title = 'Contour plot of psi'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
          nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

!     --------------------
!     plot psi in midplane
!     --------------------

      do i = 1, nnodex
         fmidre(i) = fmod(i, jmid)
      end do

      title = 'psi = rho**2 in midplane'
      titll = 'psi'
      titlr = '      '
      titlb = 'R (m)'

      call ezplot1q(title, titll, titlr, titlb, capr, fmidre, &
         nnodex, nxmx)

!     ---------------------
!     plot Mod 1/B surfaces
!     ---------------------

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = 1.0 / bmod(i,j)
         end do
      end do

      title = 'Mod 1/B surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

     call check ( nf90_put_var ( nc_id, bmod_id, &
        bmod(1:nnodex,1:nnodey) ) )
 



      title = 'Contour plot of theta'
      call ezconc(capr, y, theta, ff, nnodex, nnodey, 28, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)


!--   Plot xkperp_cold:
      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(xkperp_cold(i,j))
            fimag(i,j) = imag(xkperp_cold(i,j))
         end do
      end do


      title = 'Re xkperp_cold'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
          nnodey, 1, &
         nxmx, nymx, nlevmax, title, titx, tity)

      title = 'Im xkperp_cold'
      call ezconc(capr, y, fimag, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
          nnodey, 1, &
         nxmx, nymx, nlevmax, title, titx, tity)

        call check ( nf90_put_var ( nc_id, xkperp_real_id, &
            freal(1:nnodex,1:nnodey) ) )
        call check ( nf90_put_var ( nc_id, xkperp_imag_id, &
            fimag(1:nnodex,1:nnodey) ) )



!     --------------------------------------
!     plot xkperp_flux_plot(i,j) in midplane
!     --------------------------------------
      do i = 1, nnodex
         fmidre(i) = real(xkperp_cold(i,jmid))
         fmidim(i) = imag(xkperp_cold(i,jmid))
      end do

      title = 'Re xkperp_cold'
      titll = 'Re xkperp (1/m)'
      titlr = 'Im xkperp (1/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim, &
         nnodex, nxmx)




!--   Plot Acold:
      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(acold(i,j))
            fimag(i,j) = imag(acold(i,j))
         end do
      end do


      title = 'Re Acold'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, 1, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
          nnodey, 1, &
         nxmx, nymx, nlevmax, title, titx, tity)

!      call ezconz(capr, y, freal, ff, nnodex, nnodey, 1,
!     1   nxmx, nymx, nlevmax, title, titx, tity)


!     -----------------------
!     plot resonance surfaces
!     -----------------------

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = reomg1a(i,j)
         end do
      end do

      title = 'Majority resonant surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, 10, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
        nxmx, nymx, nlevmax, title, titx, tity)


      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = reomg2a(i,j)
         end do
      end do

      title = 'Minority resonant surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, 10, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)


      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = reomg3a(i,j)
         end do
      end do

      title = 'Species 3 resonant surfaces'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, 10, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)




      title = 'Contour plot of Janty'
      call ezconc(capr, y, xjy, ff, nnodex, nnodey, 2, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
          nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      title = 'Contour plot of Jantx'
      call ezconc(capr, y, xjx, ff, nnodex, nnodey, 2, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
       if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
          nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

    !DLG: write antenna current

    call check ( nf90_put_var ( nc_id, janty_id, &
        xjy(1:nnodex,1:nnodey) ) )
    call check ( nf90_put_var ( nc_id, jantx_id, &
        xjx(1:nnodex,1:nnodey) ) )

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(ealpha(i,j))
            fimag(i,j) = imag(ealpha(i,j))
         end do
      end do


      title = 'Real E alpha'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

    call check ( nf90_put_var ( nc_id, eAlpha_id, &
        freal(1:nnodex,1:nnodey) ) )
    call check ( nf90_put_var ( nc_id, eParallel_id, &
        real (eb(1:nnodex,1:nnodey)) ) )
 
!      call ezconz(capr, y, freal, ff, nnodex, nnodey, numb,
!     .   nxmx, nymx, nlevmax, title, titx, tity)


!     --------------------------------------
!     Write 2D.vtk file "structured points"
!     --------------------------------------
      number_points = nnodex * nnodey
      dx = capr(2) - capr(1)
      dy = y(2) - y(1)

      write(140, 2840)
 2840 format('# vtk DataFile Version 2.0')

      write(140, 2845)
 2845 format('Real E_alpha')

      write(140, 2846)
 2846 format('ASCII')

      write(140, 2847)
 2847 format('DATASET STRUCTURED_POINTS')

      write(140, 2841) nnodex,  nnodey, 1
 2841 format('DIMENSIONS', 3i8)

      write(140, 2842) capr(1), y(1), 0
 2842 format('ORIGIN', 2f9.3, 1i8)

      write(140, 2843) dx, dy, 1
 2843 format('SPACING', 2f9.3, 1i8)

      write(140, 2844) number_points
 2844 format('POINT_DATA', i10)

      write(140, 2848)
 2848 format('SCALARS Real_E_alpha float 1')

      write(140, 2849)
 2849 format('LOOKUP_TABLE default')

      write (140, 3411) ((freal(i,j), i = 1, nnodex),  j = 1, nnodey)

 3410 format(1p4e10.2)
 3411 format(6f16.4)

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(acold(i,j))
         end do
      end do

      write(140, 3853)
 3853 format('SCALARS Acold float 1')
      write(140, 2849)
      write (140, 3411) ((freal(i,j), i = 1, nnodex),  j = 1, nnodey)


      title = 'Imag E alpha'
      call ezconc(capr, y, fimag, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(ebeta(i,j))
            fimag(i,j) = imag(ebeta(i,j))
         end do
      end do
      title = 'Real Ebeta'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)


      title = 'Imag Ebeta'
      call ezconc(capr, y, fimag, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      call check ( nf90_put_var ( nc_id, eBeta_id, &
      freal(1:nnodex,1:nnodey) ) )




      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(eb(i,j))
            fimag(i,j) = imag(eb(i,j))
         end do
      end do
      title = 'Real Eb'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)



      write(140, 3966)
 3966 format('SCALARS Real_eb float 1')
      write(140, 2849)
      write(140, 3411) ((freal(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3967)
 3967 format('SCALARS Imag_eb float 1')
      write(140, 2849)
      write(140, 3411) ((fimag(i,j), i = 1, nnodex),  j = 1, nnodey)


      title = 'Imag Eb'
      call ezconc(capr, y, fimag, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)




!
!--plot E alpha in midplane
!
      do i = 1, nnodex
         fmidre(i) = real(ealpha(i,jmid))
         fmidim(i) = imag(ealpha(i,jmid))
      end do

      title = 'Ealpha'
      titll = 'Re Ealpha (V/m)'
      titlr = 'Im Ealpha (V/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim, &
         nnodex, nxmx)

!
!--plot Ebeta in midplane
!
      do i = 1, nnodex
         fmidre(i) = real(ebeta(i,jmid))
         fmidim(i) = imag(ebeta(i,jmid))
      end do

      title = 'Ebeta'
      titll = 'Re Ebeta (V/m)'
      titlr = 'Im Ebeta (V/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim, &
         nnodex, nxmx)


!     ----------------------
!     plot Eplus in midplane
!     ----------------------
      ISQ2 = SQRT(0.5)
      do i = 1, nnodex
         do j = 1, nnodey
            eplus(i,j)  = isq2 * (ealpha(i,j) + zi * ebeta(i,j))

            fmod(i,j)   = conjg(eplus(i,j)) * eplus(i,j)
            mod_Eplus(i,j) = sqrt(fmod(i,j))
         end do
      end do


      do i = 1, nnodex
         fmidre(i) = real(eplus(i,jmid))
         fmidim(i) = imag(eplus(i,jmid))
         mod_Eplus_mid(i) = mod_Eplus(i,jmid)
      end do

      title = 'Eplus'
      titll = 'Re Eplus (V/m)'
      titlr = 'Im Eplus (V/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim, &
         nnodex, nxmx)



!     -----------------------
!     plot Eminus in midplane
!     -----------------------

      do i = 1, nnodex
         do j = 1, nnodey
            eminus(i,j) = isq2 * (ealpha(i,j) - zi * ebeta(i,j))

            fmod(i,j)   = conjg(eminus(i,j)) * eminus(i,j)
            mod_Eminus(i,j) = sqrt(fmod(i,j))
         end do
      end do

      do i = 1, nnodex
         fmidre(i) = real(eminus(i,jmid))
         fmidim(i) = imag(eminus(i,jmid))
         mod_Eminus_mid(i) = mod_Eminus(i,jmid)
      end do

      title = 'Eminus'
      titll = 'Re Eminus (V/m)'
      titlr = 'Im Eminus (V/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim, &
         nnodex, nxmx)


!     ----------------------
!     plot Eb in midplane
!     ----------------------

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = conjg(eb(i,j)) * eb(i,j)
            mod_Eb(i,j) = sqrt(fmod(i,j))
!            if(i .eq. 70 .and. j .eq. 64)then
!               write(6, *) 'R(70) = ', capr(i)
!               write(6, *) 'Z(64) = ', y(j)
!               write(6, *) 'mod_Eb(70, 64) = ', mod_Eb(i,j)
!            end if
         end do
      end do

      do i = 1, nnodex
         fmidre(i) = real(eb(i,jmid))
         fmidim(i) = imag(eb(i,jmid))
         mod_Eb_mid(i) = mod_Eb(i,jmid)
      end do

      title = 'E parallel'
      titll = 'Re E parallel (V/m)'
      titlr = 'Im E parallel (V/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim, &
         nnodex, nxmx)


!     -----------------------
!     write bharvey (66) file
!     ------------------------
      write(66,309) nnodex, nnodey
      write(66,310) (capr(i), i = 1, nnodex)
      write(66,310) (y(j), j = 1, nnodey)

      write(66,310) ((eplus(i, j), i = 1, nnodex), j = 1, nnodey)
      write(66,310) ((eminus(i, j), i = 1, nnodex), j = 1, nnodey)
      write(66,310) ((eb(i, j), i = 1, nnodex), j = 1, nnodey)

      write(66,310) ((ex(i, j), i = 1, nnodex), j = 1, nnodey)
      write(66,310) ((ey(i, j), i = 1, nnodex), j = 1, nnodey)
      write(66,310) ((ez(i, j), i = 1, nnodex), j = 1, nnodey)

!     -----------------------
!     write zowens (68) file
!     ------------------------
      write(68,309) nnodex, nnodey
      write(68,310) (capr(i), i = 1, nnodex)
      write(68,310) (y(j), j = 1, nnodey)

      write(68,310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
      write(68,310) ((bratio(i, j), i = 1, nnodex), j = 1, nnodey)



      ncollin = nblue
      ncolln2 = nblueviolet




      write(140, 3949)
 3949 format('SCALARS redotj2 float 1')
      write(140, 2849)
      write(140, 3411) ((redotj2(i,j), i = 1, nnodex),  j = 1, nnodey)


!
!***********************************************************
!
      ncollin = nred
      ncolln2 = nmagenta

      title = '1/2 Re JedotE'
      call ezconc(capr, y, redotje, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      call check ( nf90_put_var ( nc_id, jeDotE_id, &
      redotje(1:nnodex,1:nnodey) ) )
 


!
!***********************************************************
!
      ncollin = ncyan
      ncolln2 = nblueviolet

      title = '1/2 Re J1dotE'
      call ezconc(capr, y, redotj1, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)



      title = '1/2 Re J2dotE'

      call ezconc(capr, y, redotj2, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)

      if (iflag .eq. 0) call boundary  (capr, y, rho, ff, &
         nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

!
!**********************************************************
!
      ncollin = npink
      ncolln2 = ngrey

      title = '1/2 Re J3dotE'
      call ezconc(capr, y, redotj3, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)


!
!***********************************************************
!
      ncollin = ngreen
      ncolln2 = nyellow

      title = '1/2 Re J4dotE'
      call ezconc(capr, y, redotj4, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      title = '1/2 Re J5dotE'
      call ezconc(capr, y, redotj5, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      title = '1/2 Re J6dotE'
      call ezconc(capr, y, redotj6, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

!
!***********************************************************
!

      ncollin = ngreen
      ncolln2 = nyellow

      title = '1/2 Re JsdotE'
      call ezconc(capr, y, redotjs, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

!
!***********************************************************
!
      if (ibackground.eq.1)then !  black background
         ncollin = nblue
         ncolln2 = nblueviolet
      end if

      if (ibackground.eq.0)then ! white background
         ncollin = nblue
         ncolln2 = nyellow
      end if

      title = '1/2 Re JdotE'
      call ezconc(capr, y, redotjt, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)




      ncollin=nred
      ncolln2=nmagenta

      title = 'wdote'
      call ezconc(capr, y,wdote, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)



      write(140, 3852)
 3852 format('SCALARS wdote float 1')
      write(140, 2849)
      write (140, 3411) ((wdote(i,j), i = 1, nnodex),  j = 1, nnodey)



      if (ibackground.eq.1)then
!        black background
         ncollin = ncyan
         ncolln2 = nblueviolet
      end if

      if (ibackground.eq.0)then
!        white background
         ncollin = nblue
         ncolln2 = nyellow
      end if

      title = 'wdoti1'
      call ezconc(capr, y,wdoti1, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)


      title = 'wdoti2'
      call ezconc(capr, y,wdoti2, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)




      write(140, 3849)
 3849 format('SCALARS wdoti2 float 1')
      write(140, 2849)
      write (140, 3411) ((wdoti2(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3850)
 3850 format('SCALARS wdoti1 float 1')
      write(140, 2849)
      write (140, 3411) ((wdoti1(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3851)
 3851 format('SCALARS rho float 1')
      write(140, 2849)
      write (140, 3411) ((rho(i,j), i = 1, nnodex),  j = 1, nnodey)





      title = 'wdoti3'
      call ezconc(capr, y,wdoti3, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      title = 'wdoti4'
      call ezconc(capr, y,wdoti4, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      title = 'wdoti5'
      call ezconc(capr, y,wdoti5, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)


      title = 'wdot_tot'
      call ezconc(capr, y,wdott, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)



      title = 'Toroidal force'
      call ezconc(capr, y, fz0, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)



      write(140, 3857)
 3857 format('SCALARS fz0 float 1')
      write(140, 2849)
      write (140, 3411) ((fz0(i,j), i = 1, nnodex),  j = 1, nnodey)



      titx = 'kx (m-1)'
      tity = 'ky (m-1)'
      title = 'Mod E alpha(kx, ky)'
      call ezconc(xkxsav, xkysav, exkmod, ff, nkxplt, nkyplt, numb, &
         nkpltdim, mkpltdim, nlevmax, title, titx, tity, iflag)

      title='3-D plot of Mod E alpha(kx,ky)'
      titz='Mod E alpha(kx,ky)'
      call ezcon3d(xkxsav, xkysav, exkmod, ff, nkxplt, nkyplt, numb, &
         nkpltdim, mkpltdim, nlevmax, title, titx, tity, titz)


      nmid = nkxplt / 2
      do m = 1, nkyplt
         fmodm(m) = exklog(nmid, m)
      end do
!
!--plot exklog vs ky
!
      title= 'Mod E alpha(ky) Spectrum'
      titll= 'Mod E alpha (V/m)'
      titlr= ''

      call ezlog1(title, titll, titlr, tity, xkysav, fmodm, &
          nkyplt, mkpltdim, exkmin, exkmax)







      mmid = nkyplt / 2
      do n = 1, nkxplt
         fmodm(n) = exklog(n, mmid)
      end do

!
!--plot exklog vs kx
!
      title= 'Mod E alpha(kx) Spectrum'
      titll= 'Mod E alpha (V/m)'
      titlr= ''

      call ezlog1(title, titll, titlr, titx, xkxsav, fmodm, &
          nkxplt, nkpltdim, exkmin, exkmax)




!--   plot log of E alpha

      titx = 'kx (m-1)'
      tity = 'ky (m-1)'
      title = 'log E alpha(kx, ky)'
      call ezconc(xkxsav, xkysav, exklog, ff, nkxplt, nkyplt, numb, &
         nkpltdim, mkpltdim, nlevmax, title, titx, tity, iflag)

      title='3-D plot of log E alpha(kx,ky)'
      titz='log E alpha(kx,ky)'
      call ezcon3d(xkxsav, xkysav, exklog, ff, nkxplt, nkyplt, numb, &
         nkpltdim, mkpltdim, nlevmax, title, titx, tity, titz)



      title = 'Mod E beta(kx, ky)'
      call ezconc(xkxsav, xkysav, eykmod, ff, nkxplt, nkyplt, numb, &
         nkpltdim, mkpltdim, nlevmax, title, titx, tity, iflag)

      title='3-D plot of Mod E beta(kx,ky)'
      titz='Mod E beta(kx,ky)'
      call ezcon3d(xkxsav, xkysav, eykmod, ff, nkxplt, nkyplt, numb, &
         nkpltdim, mkpltdim, nlevmax, title, titx, tity, titz)


      do m = 1, nkyplt
         fmodm(m) = eyklog(nmid, m)
      end do
!
!--plot eyklog vs ky
!
      title= 'Mod E beta(ky) Spectrum'
      titll= 'Mod E beta (V/m)'
      titlr= ''

      call ezlog1(title, titll, titlr, tity, xkysav, fmodm, &
          nkyplt, mkpltdim, exkmin, exkmax)



      mmid = nkyplt / 2
      do n = 1, nkxplt
         fmodm(n) = eyklog(n, mmid)
      end do

!
!--plot eyklog vs kx
!
      title= 'Mod E beta(kx) Spectrum'
      titll= 'Mod E beta (V/m)'
      titlr= ''

      call ezlog1(title, titll, titlr, titx, xkxsav, fmodm, &
          nkxplt, nkpltdim, exkmin, exkmax)



!--   plot log of E beta

      title = 'log E beta(kx, ky)'
      call ezconc(xkxsav, xkysav, eyklog, ff, nkxplt, nkyplt, numb, &
         nkpltdim, mkpltdim, nlevmax, title, titx, tity, iflag)

      title='3-D plot of log E beta(kx,ky)'
      titz='log E beta(kx,ky)'
      call ezcon3d(xkxsav, xkysav, eyklog, ff, nkxplt, nkyplt, numb, &
         nkpltdim, mkpltdim, nlevmax, title, titx, tity, titz)



      title = 'Mod E b(kx, ky)'
      call ezconc(xkxsav, xkysav, ezkmod, ff, nkxplt, nkyplt, numb, &
         nkpltdim, mkpltdim, nlevmax, title, titx, tity, iflag)

      title='3-D plot of Mod E b(kx,ky)'
      titz='Mod E b(kx,ky)'
      call ezcon3d(xkxsav, xkysav, ezkmod, ff, nkxplt, nkyplt, numb, &
         nkpltdim, mkpltdim, nlevmax, title, titx, tity, titz)


      do m = 1, nkyplt
         fmodm(m) = ezklog(nmid, m)
      end do
!
!--plot ezklog vs ky
!
      title= 'Mod Eb(ky) Spectrum'
      titll= 'Mod Eb (V/m)'
      titlr= ''

      call ezlog1(title, titll, titlr, tity, xkysav, fmodm, &
          nkyplt, mkpltdim, exkmin, exkmax)



      mmid = nkyplt / 2
      do n = 1, nkxplt
         fmodm(n) = ezklog(n, mmid)
      end do

!
!--plot ezklog vs kx
!
      title= 'Mod Eb(kx) Spectrum'
      titll= 'Mod Eb (V/m)'
      titlr= ''

      call ezlog1(title, titll, titlr, titx, xkxsav, fmodm, &
          nkxplt, nkpltdim, exkmin, exkmax)


!--   plot log of Eb

      title = 'log E b(kx, ky)'
      call ezconc(xkxsav, xkysav, ezklog, ff, nkxplt, nkyplt, numb, &
         nkpltdim, mkpltdim, nlevmax, title, titx, tity, iflag)

      title='3-D plot of log E b(kx,ky)'
      titz='Mod E b(kx,ky)'
      call ezcon3d(xkxsav, xkysav, ezklog, ff, nkxplt, nkyplt, numb, &
         nkpltdim, mkpltdim, nlevmax, title, titx, tity, titz)





      titx = 'R (m)'
      tity = 'Z (m)'
      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = conjg(ealpha(i,j)) * ealpha(i,j)
         end do
      end do

      title = 'Mod E alpha'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)


      title='3-D plot of Mod E alpha'
      titz='Mod E alpha'
      call ezcon3d(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, titz)

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j) = conjg(ebeta(i,j)) * ebeta(i,j)
         end do
      end do

      title = 'Mod E beta'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      title='3-D plot of Mod E beta'
      titz='Mod Ebeta'
      call ezcon3d(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, titz)

!     ------------
!     plot mod(Eb)
!     ------------


      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = conjg(eb(i,j)) * eb(i,j)
         end do
      end do

      title = ' Mod Eb'
      call ezconc(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      title='3-D plot of Mod Eb'
      titz='Mod Eb'
      call ezcon3d(capr, y, fmod, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, titz)

      write(140, 3956)
 3956 format('SCALARS mod_eb float 1')
      write(140, 2849)
      write(140, 3411) ((fmod(i,j), i = 1, nnodex),  j = 1, nnodey)



!     ----------------
!     plot mod(Eplus)
!     ----------------

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = conjg(eplus(i,j)) * eplus(i,j)
            mod_Eplus(i,j) = sqrt(fmod(i,j))
!            if(i .eq. 70 .and. j .eq. 64)then
!               write(6, *) 'R(70) = ', capr(i)
!               write(6, *) 'Z(64) = ', y(j)
!               write(6, *) 'mod_Eplus(70, 64) = ', mod_Eplus(i,j)
!            end if
         end do
      end do

      title = ' Mod(E_plus)'
      call ezconc(capr, y, mod_Eplus, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)


      title='3-D plot of Mod Eplus'
      titz='Mod Eplus'
      call ezcon3d(capr, y, mod_Eplus, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, titz)


      write(140, 3954)
 3954 format('SCALARS mod_eplus float 1')
      write(140, 2849)
      write(140, 3411) ((fmod(i,j), i = 1, nnodex),  j = 1, nnodey)

!     ----------------
!     plot mod(Eminus)
!     ----------------

      do i = 1, nnodex
         do j = 1, nnodey
            fmod(i,j)   = conjg(eminus(i,j)) * eminus(i,j)
            mod_Eminus(i,j) = sqrt(fmod(i,j))
!            if(i .eq. 70 .and. j .eq. 64)then
!               write(6, *) 'R(70) = ', capr(i)
!               write(6, *) 'Z(64) = ', y(j)
!               write(6, *) 'mod_Eminus(70, 64) = ', mod_Eminus(i,j)
!            end if
         end do
      end do

      title = ' Mod(E_minus)'
      call ezconc(capr, y, mod_Eminus, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)


      title='3-D plot of Mod Eminus'
      titz='Mod Eminus'
      call ezcon3d(capr, y, mod_Eminus, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, titz)


      write(140, 3955)
 3955 format('SCALARS mod_eminus float 1')
      write(140, 2849)
      write(140, 3411) ((fmod(i,j), i = 1, nnodex),  j = 1, nnodey)


      write(64,309) nnodex, nnodey
      write(64,310) (capr(i), i = 1, nnodex)
      write(64,310) (y(j), j = 1, nnodey)
      write(64,310) ((ex(i, j), i = 1, nnodex), j = 1, nnodey)
      write(64,310) ((ey(i, j), i = 1, nnodex), j = 1, nnodey)
      write(64,310) ((ez(i, j), i = 1, nnodex), j = 1, nnodey)


!     ---------------
!     plot real Eplus
!     ---------------

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j)   = real(eplus(i,j))
          fimag(i,j)   = imag(eplus(i,j))
         end do
      end do

!DLG
      call check ( nf90_put_var ( nc_id, ePlus_real_id, &
      freal(1:nnodex,1:nnodey) ) )
      call check ( nf90_put_var ( nc_id, ePlus_imag_id, &
      fimag(1:nnodex,1:nnodey) ) )

    call check ( nf90_put_var ( nc_id, ex_id, &
      real(ex(1:nnodex,1:nnodey)) ) )
    call check ( nf90_put_var ( nc_id, ey_id, &
      real(ey(1:nnodex,1:nnodey)) ) )
    call check ( nf90_put_var ( nc_id, ez_id, &
      real(ez(1:nnodex,1:nnodey)) ) )
  



      title = 'real(E_plus)'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      write(140, 3950)
 3950 format('SCALARS re_eplus float 1')
      write(140, 2849)
      write(140, 3411) ((freal(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3951)
 3951 format('SCALARS im_eplus float 1')
      write(140, 2849)
      write(140, 3411) ((fimag(i,j), i = 1, nnodex),  j = 1, nnodey)



!     ----------------
!     plot real Eminus
!     ----------------

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j)   = real(eminus(i,j))
          fimag(i,j)   = imag(eminus(i,j))
         end do
      end do

      call check ( nf90_put_var ( nc_id, eMinu_real_id, &
      freal(1:nnodex,1:nnodey) ) )
      call check ( nf90_put_var ( nc_id, eMinu_imag_id, &
      fimag(1:nnodex,1:nnodey) ) )
 
      title = 'Real(E_minus)'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      write(140, 3952)
 3952 format('SCALARS re_eminus float 1')
      write(140, 2849)
      write(140, 3411) ((freal(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3953)
 3953 format('SCALARS im_eminus float 1')
      write(140, 2849)
      write(140, 3411) ((fimag(i,j), i = 1, nnodex),  j = 1, nnodey)


      write(140, 3854)
 3854 format('SCALARS omgexb float 1')
      write(140, 2849)
      write (140, 3411) ((omgexb(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3855)
 3855 format('SCALARS uzeta float 1')
      write(140, 2849)
      write (140, 3411) ((uzeta(i,j), i = 1, nnodex),  j = 1, nnodey)


      title = 'B toroidal'
      call ezconc(capr, y, bzeta, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

!      title='3-D plot of B toroidal'
!      titz='Btor'
!      call ezcon3d(capr, y, bzeta, ff, nnodex, nnodey, numb,
!     1   nxmx, nymx, nlevmax, title, titx, tity, titz)

!      do i = 1, nnodex
!         write(6, 2163)i, capr(i), x(i), btau(i, 16),
!     1      bzeta(i,16)
!      end do
 2163 format(i5, 1p11e12.3)



      title = 'B poloidal'
      call ezconc(capr, y, btau, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      title='3-D plot of B poloidal'
      titz='Bpol'
      call ezcon3d(capr, y, btau, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, titz)



!      title='3-D plot of safety factor'
!      titz='q'
!      call ezcon3d(capr, y, qsafety, ff, nnodex, nnodey, numb,
!     1   nxmx, nymx, nlevmax, title, titx, tity, titz)





!      title = 'Spx'
!      call ezconc(capr, y, spx, ff, nnodex, nnodey, numb,
!     1   nxmx, nymx, nlevmax, title, titx, tity, iflag)
!      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
!     .   nnodey, numb,
!     1   nxmx, nymx, nlevmax, title, titx, tity)

!      title = 'Spy'
!      call ezconc(capr, y, spy, ff, nnodex, nnodey, numb,
!     1   nxmx, nymx, nlevmax, title, titx, tity, iflag)
!      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
!     .   nnodey, numb,
!     1   nxmx, nymx, nlevmax, title, titx, tity)

!      title = 'Spz'
!      call ezconc(capr, y, spz, ff, nnodex, nnodey, numb,
!     1   nxmx, nymx, nlevmax, title, titx, tity, iflag)
!      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex,
!     .   nnodey, numb,
!     1   nxmx, nymx, nlevmax, title, titx, tity)


!     --------------
!     read from rf2x
!     --------------

      read(38, 309) nnodex, nnodey
      read(38, 310) rhoplasm
      read(38, 310) (x(i), i = 1, nnodex)
      read(38, 310) (y(j), j = 1, nnodey)
      read(38, 310) (capr(i), i = 1, nnodex)

      read(38, 310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((theta(i, j), i = 1, nnodex), j = 1, nnodey)


      read(38, 309) nnoderho2
      read(38, 310) (rhon(n), n = 1, nnoderho2)
!      read(38, 310) (gradprlb2_avg(n), n = 1, nnoderho2)

      read(38, 310) (dvol(n), n = 1, nnoderho2)
      read(38, 310) (volume(n), n = 1, nnoderho2)

      read(38, 310) (redotjeavg(n), n = 1, nnoderho2)
      read(38, 310) (redotj1avg(n), n = 1, nnoderho2)
      read(38, 310) (redotj2avg(n), n = 1, nnoderho2)
      read(38, 310) (redotj3avg(n), n = 1, nnoderho2)
      read(38, 310) (redotj4avg(n), n = 1, nnoderho2)
      read(38, 310) (redotj5avg(n), n = 1, nnoderho2)
      read(38, 310) (redotj6avg(n), n = 1, nnoderho2)
      read(38, 310) (redotjiavg(n), n = 1, nnoderho2)

      read(38, 310) (xjprlavg(n), n = 1, nnoderho2)
      read(38, 310) (xjprl_int(n), n = 1, nnoderho2)

      read(38, 310) (wdoteavg(n), n = 1, nnoderho2)
      read(38, 310) (wdoti1avg(n), n = 1, nnoderho2)
      read(38, 310) (wdoti2avg(n), n = 1, nnoderho2)
      read(38, 310) (wdoti3avg(n), n = 1, nnoderho2)
      read(38, 310) (wdoti4avg(n), n = 1, nnoderho2)
      read(38, 310) (wdoti5avg(n), n = 1, nnoderho2)
      read(38, 310) (wdoti6avg(n), n = 1, nnoderho2)



      read(38, 310) (redotje_int(n), n = 1, nnoderho2)
      read(38, 310) (redotj1_int(n), n = 1, nnoderho2)
      read(38, 310) (redotj2_int(n), n = 1, nnoderho2)
      read(38, 310) (redotj3_int(n), n = 1, nnoderho2)
      read(38, 310) (redotj4_int(n), n = 1, nnoderho2)
      read(38, 310) (redotj5_int(n), n = 1, nnoderho2)
      read(38, 310) (redotj6_int(n), n = 1, nnoderho2)

      read(38, 310) (wdote_int(n), n = 1, nnoderho2)
      read(38, 310) (wdot1_int(n), n = 1, nnoderho2)
      read(38, 310) (wdot2_int(n), n = 1, nnoderho2)
      read(38, 310) (wdot3_int(n), n = 1, nnoderho2)
      read(38, 310) (wdot4_int(n), n = 1, nnoderho2)
      read(38, 310) (wdot5_int(n), n = 1, nnoderho2)
      read(38, 310) (wdot6_int(n), n = 1, nnoderho2)

      read(38, 310) (redotje_dvol(n), n = 1, nnoderho2)
      read(38, 310) (redotj1_dvol(n), n = 1, nnoderho2)
      read(38, 310) (redotj2_dvol(n), n = 1, nnoderho2)
      read(38, 310) (redotj3_dvol(n), n = 1, nnoderho2)
      read(38, 310) (redotj4_dvol(n), n = 1, nnoderho2)
      read(38, 310) (redotj5_dvol(n), n = 1, nnoderho2)
      read(38, 310) (redotj6_dvol(n), n = 1, nnoderho2)


      read(38, 310) (wdote_dvol(n),  n = 1, nnoderho2)
      read(38, 310) (wdoti1_dvol(n), n = 1, nnoderho2)
      read(38, 310) (wdoti2_dvol(n), n = 1, nnoderho2)
      read(38, 310) (wdoti3_dvol(n), n = 1, nnoderho2)
      read(38, 310) (wdoti4_dvol(n), n = 1, nnoderho2)
      read(38, 310) (wdoti5_dvol(n), n = 1, nnoderho2)
      read(38, 310) (wdoti6_dvol(n), n = 1, nnoderho2)



!     --------------
!     read from main
!     --------------

      read(38, 310) (wdote_ql(n),  n = 1, nnoderho)
      read(38, 310) (wdoti1_ql(n), n = 1, nnoderho)
      read(38, 310) (wdoti2_ql(n), n = 1, nnoderho)
      read(38, 310) (wdoti3_ql(n), n = 1, nnoderho)
      read(38, 310) (wdoti4_ql(n), n = 1, nnoderho)
      read(38, 310) (wdoti5_ql(n), n = 1, nnoderho)
      read(38, 310) (wdoti6_ql(n), n = 1, nnoderho)

      read(38, 310) (dldbavg(n), n = 1, nnoderho)

      read(38, 310) ((bxwave(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((bywave(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((bzwave(i, j), i = 1, nnodex), j = 1, nnodey)

      read(38, 310) ((fpsi0(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((ftheta0(i, j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((pressi(i,j), i = 1, nnodex), j = 1, nnodey)
      read(38, 310) ((dldb_tot12(i, j), i = 1, nnodex), &
                                        j = 1, nnodey)

      read (38, 309) ndisti1, ndisti2, ndisti3

      read (38, 310) ((eplus_flux_plot(i, j), i = 1, nnodex), &
                                               j = 1, nnodey)

      read (38, 310) ((eminus_flux_plot(i, j), i = 1, nnodex), &
                                               j = 1, nnodey)

      read (38, 310) ((xkperp_flux_plot(i, j), i = 1, nnodex), &
                                               j = 1, nnodey)

      read (38, 310) ((eplus_flux(n, m), n = 1, nnoderho), &
                                         m = 1, mnodetheta)
      read (38, 310) ((eminus_flux(n, m), n = 1, nnoderho), &
                                          m = 1, mnodetheta)

      read (38, 310) ((xkperp_flux(n, m), n = 1, nnoderho), &
                                          m = 1, mnodetheta)


      read (38, 310) ((capr_flux(n,m),   n = 1, nnoderho), &
                                         m = 1, mnodetheta)
      read (38, 310) ((capz_flux(n,m),   n = 1, nnoderho), &
                                         m = 1, mnodetheta)

      read (38, 310) ((bmod_flux(n,m),   n = 1, nnoderho), &
                                         m = 1, mnodetheta)




      write(66,310) ((bxwave(i, j), i = 1, nnodex), j = 1, nnodey)
      write(66,310) ((bywave(i, j), i = 1, nnodex), j = 1, nnodey)
      write(66,310) ((bzwave(i, j), i = 1, nnodex), j = 1, nnodey)

      write(66,310) rhoplasm
      write(66,310) ((rho(i, j), i = 1, nnodex), j = 1, nnodey)


      title= 'Integral dl/B'
      titll= 'Integral dl/B (m)'
      titlr='       '

      call ezplot1_0(title, titll, titlr, rhon_half_save, dldbavg, &
          nnoderho_half, nrhomax)



      title= 'Quasilinear wdot'
      titll= 'wdot (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon_half_save, &
          wdoti1_ql, wdoti2_ql, wdoti3_ql, wdoti4_ql, wdoti5_ql, &
          wdoti6_ql, wdote_ql, nnoderho_half, nrhomax)


      title= 'Integrated quasilinear wdot'
      titll= 'P (Watts)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon_half_save, &
         wdoti1_ql_int, wdoti2_ql_int, wdoti3_ql_int, wdoti4_ql_int, &
         wdoti5_ql_int, wdoti6_ql_int, wdote_ql_int, nnoderho_half, &
         nrhomax)

      title= 'Quasilinear wdoti1'
      titll= 'wdoti1 (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot1(title, titll, titlr, rhon_half_save, wdoti1_ql, &
          nnoderho_half, nrhomax)

      title= 'Quasilinear wdoti2'
      titll= 'wdoti2 (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot1(title, titll, titlr, rhon_half_save, wdoti2_ql, &
          nnoderho_half, nrhomax)


      title= 'Quasilinear wdoti3'
      titll= 'wdoti2 (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot1(title, titll, titlr, rhon_half_save, wdoti3_ql, &
          nnoderho_half, nrhomax)



!     ----------------------------------
!     Plot on fine mesh from rf2x_setup2
!     ----------------------------------

      nnoderho2_half = nnoderho2 - 1

      do n = 1, nnoderho2_half
       rhon_half(n) = (rhon(n) + rhon(n+1)) / 2.0
!         write(6, 1312)n, rhon(n), rhon_half(n)
      end do


      title= 'Flux surface average redotj'
      titll= 'redotj (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon_half, &
          redotj1avg, redotj2avg, redotj3avg, redotj4avg, redotj5avg, &
          redotj6avg, redotjeavg, nnoderho2_half, nrhomax)


      title= 'Integrated Real(EdotJ)'
      titll= 'P(watts)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon, &
          redotj1_int, redotj2_int, redotj3_int, redotj4_int, &
          redotj5_int, redotj6_int, redotje_int, nnoderho2, nrhomax)


      title= 'Real(EdotJ) * dvol'
      titll= 'P(watts)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon, &
          redotj1_dvol, redotj2_dvol, redotj3_dvol, redotj4_dvol, &
          redotj5_dvol, redotj6_dvol, redotje_dvol, nnoderho2, nrhomax)

      title= 'Flux surface average wdot'
      titll= 'wdot (watts/m3)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon_half, &
         wdoti1avg, wdoti2avg, wdoti3avg, wdoti4avg, wdoti5avg, &
         wdoti6avg, wdoteavg, nnoderho2_half, nrhomax)


      title= 'Integrated Wdot'
      titll= 'P (watts)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon, &
          wdot1_int, wdot2_int, wdot3_int, wdot4_int, wdot5_int, &
          wdot6_int, wdote_int, nnoderho2, nrhomax)

      title= 'Wdot * dvol'
      titll= 'P (watts)'
      titlr= '       '
      titlb= 'rho'

      call ezplot7(title, titll, titlr, titlb, rhon, &
          wdoti1_dvol, wdoti2_dvol, wdoti3_dvol, wdoti4_dvol, &
          wdoti5_dvol, wdoti6_dvol, wdote_dvol, nnoderho2, nrhomax)

      title= 'Flux surface average redotj'
      titll= 'redotj (watts/m3)'
      titlr='       '
      titlb= 'rho'

      call ezplot2p(title, titll, titlr, titlb, rhon_half, redotjeavg, &
         redotjiavg, nnoderho2_half, nrhomax)

      write(15, *)
      write(15, *) 'Flux surface driven current'
      write(15, *)
      write(15, *) '        n      rho       J (A/m2)     I (A)'
      write(15, *)

      do n = 1, nnoderho2_half
       write (15, 1312) n, rhon_half(n), xjprlavg(n), xjprl_int(n)
            write (6, 1312) n, rhon_half(n), xjprlavg(n), xjprl_int(n)
      end do

      title= 'Flux surface driven current'
      titll= 'xjprl (Amps/m2)'
      titlr= 'I (Amps)'
      titlb= 'rho'

      call ezplot2(title, titll, titlr, titlb, rhon_half, xjprlavg, &
          xjprl_int, nnoderho2_half, nrhomax)


      if(xjprl_int(nnoderho2) .lt. 0.0) then

         titll= '-xjprl (Amps/m2)'
         titlr= '-I (Amps)'

         xjprlavg = - xjprlavg
         xjprl_int = - xjprl_int

         call ezplot2(title, titll, titlr, titlb, rhon_half, xjprlavg, &
            xjprl_int, nnoderho2_half, nrhomax)

         titll= 'xjprl (MA/m2/MW)'
       xjprlavg = xjprlavg  / prfinMOD

         call ezplot1(title, titll, titlr, rhon_half, xjprlavg, &
            nnoderho2_half, nrhomax)

         titll= 'I (kA)'
       xjprl_int  = xjprl_int / 1.0e+03
         call ezplot1(title, titll, titlr, rhon_half, xjprl_int, &
            nnoderho2_half, nrhomax)

      end if

      if(xjprl_int(nnoderho2) .gt. 0.0) then

         titll= 'xjprl (Amps/m2)'
         titlr= 'I (Amps)'

         xjprlavg = xjprlavg
         xjprl_int = xjprl_int

         titll= 'xjprl (MAmps/m2/MW)'
       xjprlavg = xjprlavg  / prfinMOD

         call ezplot1(title, titll, titlr, rhon_half, xjprlavg, &
            nnoderho2_half, nrhomax)

         titll= 'I (kA)'
       xjprl_int  = xjprl_int / 1.0e+03
         call ezplot1(title, titll, titlr, rhon_half, xjprl_int, &
            nnoderho2_half, nrhomax)

      end if



      title= 'Volume element dvol'
      titll= 'dvol (m3)'
      titlr= '       '
      call ezplot1(title, titll, titlr, rhon_half, dvol, &
           nnoderho2_half, nrhomax)

      title= 'Integrated volume'
      titll= 'volume (m3)'
      titlr= '       '
      call ezplot1(title, titll, titlr, rhon, volume, &
           nnoderho2, nrhomax)

!     --------------------------------------
!     Reset rhon back to its original value
!     --------------------------------------
      do n = 1, nnoderho
       rhon(n) = rhon_save(n)
       rhon_half(n) = rhon_half_save(n)
      end do

!     ----------------------------------------------------------------
!     Read quasilinear diffusion coefficients from file out_cql3d.coef
!     ----------------------------------------------------------------


      if(ndisti2 .eq. 1)then

         open(unit=42,file='out_cql3d.coef2', status='unknown', &
                                                   form='formatted')

         read (42, 309) nuper
         read (42, 309) nupar
         read (42, 309) nnoderho

         allocate( UPERP(nuper) )
         allocate( UPARA(nupar) )
       allocate( f_cql_cart(nuper, nupar, nnoderho) )

       allocate( wperp_cql(nnoderho) )
       allocate( wpar_cql(nnoderho) )

         allocate( bqlavg_i1(nuper, nupar, nnoderho) )
         allocate( cqlavg_i1(nuper, nupar, nnoderho) )
         allocate( eqlavg_i1(nuper, nupar, nnoderho) )
         allocate( fqlavg_i1(nuper, nupar, nnoderho) )

         allocate( bqlavg_i1_2d(nupar, nuper) )
         allocate( E_kick_2d (nupar, nnoderho) )
       allocate( f_cql_cart_2d (nupar, nuper) )


         read (42, 3310) vc_cgs
         read (42, 3310) UminPara, UmaxPara

         read (42, 3310) (rhon(n), n = 1, nnoderho)
         read (42, 3310) (uperp(i_uperp), i_uperp = 1, nuper)
         read (42, 3310) (upara(i_upara), i_upara = 1, nupar)

         read (42, 3310) (((bqlavg_i1(i_uperp, i_upara, n), &
           i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         read (42, 3310) (((cqlavg_i1(i_uperp, i_upara, n), &
           i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         read (42, 3310) (((eqlavg_i1(i_uperp, i_upara, n), &
           i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         read (42, 3310) (((fqlavg_i1(i_uperp, i_upara, n), &
           i_uperp = 1, nuper), i_upara = 1, nupar), n = 1, nnoderho)

         read (42, 3310) xmi

         close (42)

!        --------------------------------------------
!       read data for plotting f(u_perp, u_parallel)
!        --------------------------------------------

         open(unit=238,file='out238',status='old',form='formatted')

       read  (238, 309) n_u
       read  (238, 309) n_psi
       read  (238, 310) vc

         read  (238, 310) (u(i_u), i_u = 1, n_u)
       read  (238, 309) (n_theta_(i_psi), i_psi = 1, n_psi)
         read  (238, 310) ((theta(i_theta, i_psi), &
                i_theta = 1, n_theta_(i_psi)), i_psi = 1, n_psi)

         read  (238, 310) (((f_cql(i_theta, i_u, i_psi), &
           i_theta = 1, n_theta_(i_psi)), i_u = 1, n_u), &
           i_psi = 1, n_psi)


       read (238, 309) nuper
       read (238, 309) nupar
       read (238, 309) nnoderho

!       write (6, *) "nuper =", nuper
!        write (6, *) "nupar =", nupar
!       write (6, *) "nnoderho =", nnoderho

         read (238, 310) (uperp(i_uperp), i_uperp = 1, nuper)
         read (238, 310)  (upara(i_upara), i_upara = 1, nupar)

         read (238, 310) (((f_cql_cart(i_uperp, i_upara, i_psi), &
           i_uperp = 1, nuper), i_upara = 1, nupar), &
           i_psi = 1, nnoderho)

         read (238, 310) (wperp_cql(i_psi), i_psi = 1, nnoderho)
         read (238, 310) (wpar_cql(i_psi),  i_psi = 1, nnoderho)

         close (238)

 4319    format(i10, 1p1e16.8)

         title= 'Wperp (keV)'
         titll= 'Wperp (keV)'
         titlr='       '
         call ezplot1(title, titll, titlr, rhon, wperp_cql, &
            nnoderho, nrhomax)

         title= 'Wpar (keV)'
         titll= 'Wpar (keV)'
         titlr='       '
         call ezplot1(title, titll, titlr, rhon, wpar_cql, &
            nnoderho, nrhomax)

!       ---------------------------------
!      1D plots of f(u, theta = const)
!       ---------------------------------
      titll = 'log f(u)'
      titlr = '    '
      titx = 'u'

      do i_psi_index = 1, 6

         if(i_psi_index .eq. 1)title = 'log f_psi1(u, theta = const)'
         if(i_psi_index .eq. 2)title = 'log f_psi2(u, theta = const)'
         if(i_psi_index .eq. 3)title = 'log f_psi3(u, theta = const)'
         if(i_psi_index .eq. 4)title = 'log f_psi4(u, theta = const)'
         if(i_psi_index .eq. 5)title = 'log f_psi5(u, theta = const)'
         if(i_psi_index .eq. 6)title = 'log f_psi6(u, theta = const)'

           i_psi = i_psi_array(i_psi_index)

         i_theta1 = 1
         i_theta2 = int(1./6. * n_theta_(i_psi))
         i_theta3 = int(2./6. * n_theta_(i_psi))
         i_theta4 = int(3./6. * n_theta_(i_psi))
         i_theta5 = int(4./6. * n_theta_(i_psi))
         i_theta6 = int(5./6. * n_theta_(i_psi))
         i_theta7 = n_theta_(i_psi)

           do i_u = 1, n_u
              f_cql_1d_1(i_u) = -40.0
            f_cql_1d_2(i_u) = -40.0
            f_cql_1d_3(i_u) = -40.0
            f_cql_1d_4(i_u) = -40.0
            f_cql_1d_5(i_u) = -40.0
            f_cql_1d_6(i_u) = -40.0
            f_cql_1d_7(i_u) = -40.0

              if (f_cql(i_theta1, i_u, i_psi) .gt. 0.0) &
              f_cql_1d_1(i_u) = alog10(f_cql(i_theta1, i_u, i_psi))
            if (f_cql(i_theta2, i_u, i_psi) .gt. 0.0) &
               f_cql_1d_2(i_u) = alog10(f_cql(i_theta2, i_u, i_psi))
            if (f_cql(i_theta3, i_u, i_psi) .gt. 0.0) &
               f_cql_1d_3(i_u) = alog10(f_cql(i_theta3, i_u, i_psi))
            if (f_cql(i_theta4, i_u, i_psi) .gt. 0.0) &
               f_cql_1d_4(i_u) = alog10(f_cql(i_theta4, i_u, i_psi))
            if (f_cql(i_theta5, i_u, i_psi) .gt. 0.0) &
               f_cql_1d_5(i_u) = alog10(f_cql(i_theta5, i_u, i_psi))
            if (f_cql(i_theta6, i_u, i_psi) .gt. 0.0) &
               f_cql_1d_6(i_u) = alog10(f_cql(i_theta6, i_u, i_psi))
            if (f_cql(i_theta7, i_u, i_psi) .gt. 0.0) &
               f_cql_1d_7(i_u) = alog10(f_cql(i_theta7, i_u, i_psi))
           end do


         call ezlog1_f(title, titll, titlr, titx, u, f_cql_1d_1, &
              f_cql_1d_2, f_cql_1d_3, f_cql_1d_4, f_cql_1d_5, &
              f_cql_1d_6, f_cql_1d_7, &
              n_u, n_u_max)

        end do




!        ------------------------------------------
!        2D plots of f_cql_cart(u_perp, u_parallel)
!        ------------------------------------------
         titx = 'u_parallel'
         tity = 'u_perp'

         write(6, *)"i_psi1 = ", i_psi1, "   rho1 = ", rhon(i_psi1)
         write(6, *)"i_psi2 = ", i_psi2, "   rho2 = ", rhon(i_psi2)
         write(6, *)"i_psi3 = ", i_psi3, "   rho3 = ", rhon(i_psi3)
         write(6, *)"i_psi4 = ", i_psi4, "   rho4 = ", rhon(i_psi4)
         write(6, *)"i_psi5 = ", i_psi5, "   rho5 = ", rhon(i_psi5)
         write(6, *)"i_psi6 = ", i_psi6, "   rho6 = ", rhon(i_psi6)


         write(15, *)"i_psi1 = ", i_psi1, "   rho1 = ", rhon(i_psi1)
         write(15, *)"i_psi2 = ", i_psi2, "   rho2 = ", rhon(i_psi2)
         write(15, *)"i_psi3 = ", i_psi3, "   rho3 = ", rhon(i_psi3)
         write(15, *)"i_psi4 = ", i_psi4, "   rho4 = ", rhon(i_psi4)
         write(15, *)"i_psi5 = ", i_psi5, "   rho5 = ", rhon(i_psi5)
         write(15, *)"i_psi6 = ", i_psi6, "   rho6 = ", rhon(i_psi6)

         numb = 20

       i_psi = i_psi1
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               f_cql_cart_2d(i_upara, i_uperp) = &
                           f_cql_cart(i_uperp, i_upara, i_psi)
            end do
         end do

         title = 'f_cql_psi1(u_perp, u_parallel)'
         call ezconcx(upara, uperp, f_cql_cart_2d, ff, nupar, nuper, 21, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)

       i_psi = i_psi2
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               f_cql_cart_2d(i_upara, i_uperp) = &
                           f_cql_cart(i_uperp, i_upara, i_psi)
            end do
         end do

         title = 'f_cql_psi2(u_perp, u_parallel)'
         call ezconcx(upara, uperp, f_cql_cart_2d, ff, nupar, nuper, 21, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)

        i_psi = i_psi3
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               f_cql_cart_2d(i_upara, i_uperp) = &
                           f_cql_cart(i_uperp, i_upara, i_psi)
            end do
         end do

         title = 'f_cql_psi3(u_perp, u_parallel)'
         call ezconcx(upara, uperp, f_cql_cart_2d, ff, nupar, nuper, 21, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)

        i_psi = i_psi4
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               f_cql_cart_2d(i_upara, i_uperp) = &
                           f_cql_cart(i_uperp, i_upara, i_psi)
            end do
         end do

         title = 'f_cql_psi4(u_perp, u_parallel)'
         call ezconcx(upara, uperp, f_cql_cart_2d, ff, nupar, nuper, 21, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)

         i_psi = i_psi5
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               f_cql_cart_2d(i_upara, i_uperp) = &
                           f_cql_cart(i_uperp, i_upara, i_psi)
            end do
         end do

         title = 'f_cql_psi5(u_perp, u_parallel)'
         call ezconcx(upara, uperp, f_cql_cart_2d, ff, nupar, nuper, 21, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)


        i_psi = i_psi6
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               f_cql_cart_2d(i_upara, i_uperp) = &
                           f_cql_cart(i_uperp, i_upara, i_psi)
            end do
         end do

         title = 'f_cql_psi6(u_perp, u_parallel)'
         call ezconcx(upara, uperp, f_cql_cart_2d, ff, nupar, nuper, 21, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)



!        --------------------------------------
!        2D plots of bqlavg(u_perp, u_parallel)
!        --------------------------------------

         i_psi = i_psi1
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlavg_i1_2d(i_upara, i_uperp) = &
                           bqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do

         title = 'bqlavg_psi1(u_perp, u_parallel)'
         call ezconc(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper, numb, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)


!        ---------------------------------------------
!        Write Bql_avg_2D.vtk file "structured points"
!        ---------------------------------------------
         number_points = nupar * nuper
         dx = upara(2) - upara(1)
         dy = uperp(2) - uperp(1)
         write(142, 2840)
         write(142, 4850)
 4850    format('Quasilinear diffusion coefficients')
         write(142, 2846)
         write(142, 2847)
         write(142, 2841) nupar,  nuper, 1
         write(142, 2842) upara(1), uperp(1), 0
         write(142, 2843) dx, dy, 1
         write(142, 2844) number_points

         write(142, 4851)
 4851    format('SCALARS bqlavg_psi1 float 1')
         write(142, 2849)
         write(142, 3411) ((bqlavg_i1_2d(i,j), i = 1, nupar), &
                                               j = 1, nuper)



         i_psi = i_psi2
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlavg_i1_2d(i_upara, i_uperp) = &
                           bqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do

         title = 'bqlavg_psi2(u_perp, u_parallel)'
         call ezconc(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper, numb, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)

         write(142, 4852)
 4852    format('SCALARS bqlavg_psi2 float 1')
         write(142, 2849)
         write(142, 3411) ((bqlavg_i1_2d(i,j), i = 1, nupar), &
                                               j = 1, nuper)






         i_psi = i_psi3
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlavg_i1_2d(i_upara, i_uperp) = &
                           bqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do

         title = 'bqlavg_psi3(u_perp, u_parallel)'
         call ezconc(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper, numb, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)

         write(142, 4853)
 4853    format('SCALARS bqlavg_psi3 float 1')
         write(142, 2849)
         write(142, 3411) ((bqlavg_i1_2d(i,j), i = 1, nupar), &
                                               j = 1, nuper)





         i_psi = i_psi4
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlavg_i1_2d(i_upara, i_uperp) = &
                           bqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do

         title = 'bqlavg_psi4(u_perp, u_parallel)'
         call ezconc(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper, numb, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)

         write(142, 4854)
 4854    format('SCALARS bqlavg_psi4 float 1')
         write(142, 2849)
         write(142, 3411) ((bqlavg_i1_2d(i,j), i = 1, nupar), &
                                               j = 1, nuper)






         i_psi = i_psi5
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlavg_i1_2d(i_upara, i_uperp) = &
                           bqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do

         title = 'bqlavg_psi5(u_perp, u_parallel)'
         call ezconc(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper, numb, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)

         write(142, 4855)
 4855    format('SCALARS bqlavg_psi5 float 1')
         write(142, 2849)
         write(142, 3411) ((bqlavg_i1_2d(i,j), i = 1, nupar), &
                                               j = 1, nuper)




         i_psi = i_psi6
         do  i_upara = 1, nupar
            do i_uperp = 1, nuper
               bqlavg_i1_2d(i_upara, i_uperp) = &
                           bqlavg_i1(i_uperp, i_upara, i_psi) / 1.0e+36
            end do
         end do

         title = 'bqlavg_psi6(u_perp, u_parallel)'
         call ezconc(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper, numb, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)

         write(142, 4856)
 4856    format('SCALARS bqlavg_psi6 float 1')
         write(142, 2849)
         write(142, 3411) ((bqlavg_i1_2d(i,j), i = 1, nupar), &
                                               j = 1, nuper)




!        --------------------------------------
!        2D plots energy kicks (eV)
!        --------------------------------------

         write(6,  *) 'energy kicks (eV)'
         write(15, *) 'energy kicks (eV)'

         title = 'Energy kick (eV)'
         titz='Energy kick (eV)'

         i_psi = i_psi1
         do  i_upara = 1, nupar
            vpara_cgs = abs(upara(i_upara) + .001) * vc_cgs
          vpara_mks = vpara_cgs * 1.0e-02
            do i_uperp = 1, nuper
             vperp_mks = uperp(i_uperp) * vc_cgs * 1.0e-02
             E_eV = 0.5 * xmi * (vperp_mks**2 + vpara_mks**2) / q
               bqlavg_i1_2d(i_upara, i_uperp) = (xmi / q) * &
                  sqrt(2.* bqlavg_i1(i_uperp, i_upara, i_psi)/vpara_cgs) &
                  * 1.0e-04
            end do
         end do

         call ezconcx(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper,numb, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)



         i_psi = i_psi2
         do  i_upara = 1, nupar
            vpara_cgs = abs(upara(i_upara) + .001) * vc_cgs
          vpara_mks = vpara_cgs * 1.0e-02
            do i_uperp = 1, nuper
             vperp_mks = uperp(i_uperp) * vc_cgs * 1.0e-02
             E_eV = 0.5 * xmi * (vperp_mks**2 + vpara_mks**2) / q
               bqlavg_i1_2d(i_upara, i_uperp) = (xmi / q) * &
                  sqrt(2.* bqlavg_i1(i_uperp, i_upara, i_psi)/vpara_cgs) &
                  * 1.0e-04
            end do
         end do

         call ezconcx(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper,numb, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)




         i_psi = i_psi3
         do  i_upara = 1, nupar
            vpara_cgs = abs(upara(i_upara) + .001) * vc_cgs
          vpara_mks = vpara_cgs * 1.0e-02
            do i_uperp = 1, nuper
             vperp_mks = uperp(i_uperp) * vc_cgs * 1.0e-02
             E_eV = 0.5 * xmi * (vperp_mks**2 + vpara_mks**2) / q
               bqlavg_i1_2d(i_upara, i_uperp) = (xmi / q) * &
                  sqrt(2.* bqlavg_i1(i_uperp, i_upara, i_psi)/vpara_cgs) &
                  * 1.0e-04
            end do
         end do

         call ezconcx(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper,numb, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)




         i_psi = i_psi4
         do  i_upara = 1, nupar
            vpara_cgs = abs(upara(i_upara) + .001) * vc_cgs
          vpara_mks = vpara_cgs * 1.0e-02
            do i_uperp = 1, nuper
             vperp_mks = uperp(i_uperp) * vc_cgs * 1.0e-02
             E_eV = 0.5 * xmi * (vperp_mks**2 + vpara_mks**2) / q
               bqlavg_i1_2d(i_upara, i_uperp) = (xmi / q) * &
                  sqrt(2.* bqlavg_i1(i_uperp, i_upara, i_psi)/vpara_cgs) &
                  * 1.0e-04
            end do
         end do

         call ezconcx(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper,numb, &
            NUPAR, NUPER, nlevmax, title, titx, tity, iflag)



         i_psi = i_psi5
         do  i_upara = 1, nupar
            vpara_cgs = abs(upara(i_upara) + .001) * vc_cgs
          vpara_mks = vpara_cgs * 1.0e-02
            do i_uperp = 1, nuper
             vperp_mks = uperp(i_uperp) * vc_cgs * 1.0e-02
             E_eV = 0.5 * xmi * (vperp_mks**2 + vpara_mks**2) / q
               bqlavg_i1_2d(i_upara, i_uperp) = (xmi / q) * &
                  sqrt(2.* bqlavg_i1(i_uperp, i_upara, i_psi)/vpara_cgs) &
                  * 1.0e-04
            end do
         end do

         call ezconcx(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper, &
            numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)




         i_psi = i_psi6
         do  i_upara = 1, nupar
            vpara_cgs = abs(upara(i_upara) + .001) * vc_cgs
          vpara_mks = vpara_cgs * 1.0e-02
            do i_uperp = 1, nuper
             vperp_mks = uperp(i_uperp) * vc_cgs * 1.0e-02
             E_eV = 0.5 * xmi * (vperp_mks**2 + vpara_mks**2) / q
               bqlavg_i1_2d(i_upara, i_uperp) = (xmi / q) * &
                  sqrt(2.* bqlavg_i1(i_uperp, i_upara, i_psi)/vpara_cgs) &
                  * 1.0e-04
            end do
         end do

         call ezconcx(upara, uperp, bqlavg_i1_2d, ff, nupar, nuper, &
            numb, NUPAR, NUPER, nlevmax, title, titx, tity, iflag)



!        --------------------------------------
!        2D plots energy kicks for 1keV ion (eV)
!        --------------------------------------
         write(6, *)
         write(6, *) 'energy kicks for 1 keV ion (eV)'

         write(15, *)
         write(15, *) 'energy kicks for 1 keV ion (eV)'

         title = 'Energy kick for 1 keV ion (eV)'
         titx = 'u_parallel'
         tity = 'rho'
         titz='Energy kick (eV)'
         numb = 15


         E_eV = 1000.
         vperp_mks = sqrt(2.0 * E_eV * q / xmi)
         vperp_cgs = vperp_mks * 100.
         uperp_1kev = vperp_cgs / vc_cgs
         duperp = uperp(3) - uperp(2)
         i_uperp = int(uperp_1kev / duperp) + 2

         uperp_1kev = (i_uperp - 1) * duperp
         vperp_cgs = uperp_1kev * vc_cgs
         vperp_mks = vperp_cgs / 100.

         E_eV = 0.5 * xmi * (vperp_mks**2) / q

         write (6, *) 'i_uperp =', i_uperp
         write (6, *) 'E = ', E_eV, ' eV'

         write (15, *) 'i_uperp =', i_uperp
         write (15, *) 'E = ', E_eV, ' eV'

         do i_upara = 1, nupar
            vpara_cgs = abs(upara(i_upara) + .001) * vc_cgs
          vpara_mks = vpara_cgs * 1.0e-02
            do i_psi = 1, nnoderho
             vperp_mks = uperp(i_uperp) * vc_cgs * 1.0e-02
             E_eV = 0.5 * xmi * (vperp_mks**2 + vpara_mks**2) / q
               E_kick_2d(i_upara, i_psi) = (xmi / q) * &
                  sqrt(2. * bqlavg_i1(i_uperp, i_upara, i_psi) &
                  / vpara_cgs) * 1.0e-04
            end do
         end do

         call ezconcx(upara, rhon, E_kick_2d, ff, nupar, nnoderho, numb, &
            NUPAR, nnoderho, nlevmax, title, titx, tity, iflag)
         call ezcon3dv(upara, rhon, E_kick_2d, ff, nupar, nnoderho, &
            numb, NUPAR, nnoderho, nlevmax, title, titx, tity, titz)


!        ---------------------------------------------
!        Write E_kicks_2D.vtk file "structured points"
!        ---------------------------------------------
         number_points = nupar * nnoderho
         dx = upara(2) - upara(1)
         dy = rhon(2) - rhon(1)

         write(141, 2840)
         write(141, 2852)
 2852    format('Energy kicks for 1 keV')
         write(141, 2846)
         write(141, 2847)
         write(141, 2841) nupar,  nnoderho, 1
         write(141, 2842) upara(1), rhon(1), 0
         write(141, 2843) dx, dy, 1
         write(141, 2844) number_points

         write(141, 3848)
 3848    format('SCALARS E_kick_2d float 1')
         write(141, 2849)
         write(141, 3411) ((E_kick_2d(i,j), i = 1, nupar), &
                                            j = 1, nnoderho)


      end if

      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(bxwave(i,j))
            fimag(i,j) = imag(bxwave(i,j))
         end do
      end do

      call check ( nf90_put_var ( nc_id, bx_wave_real_id, &
      freal(1:nnodex,1:nnodey) ) )
      call check ( nf90_put_var ( nc_id, bx_wave_imag_id, &
      fimag(1:nnodex,1:nnodey) ) )
 
      titx = 'R (m)'
      tity = 'Z (m)'

      title = 'Real Bx_wave'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)



      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(bzwave(i,j))
            fimag(i,j) = imag(bzwave(i,j))
         end do
      end do

      call check ( nf90_put_var ( nc_id, bz_wave_real_id, &
      freal(1:nnodex,1:nnodey) ) )
      call check ( nf90_put_var ( nc_id, bz_wave_imag_id, &
      fimag(1:nnodex,1:nnodey) ) )
 
      title = 'Real Bz_wave'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)



!
!--plot Bx_wave in midplane
!
      do i = 1, nnodex
         fmidre(i) = real(bxwave(i,jmid))
         fmidim(i) = imag(bxwave(i,jmid))
      end do

      title = 'Bx_wave'
      titll = 'Re Bx_wave (T)'
      titlr = 'Im Bx_wave (T)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim, &
         nnodex, nxmx)


!
!--plot Bz_wave in midplane
!
      do i = 1, nnodex
         fmidre(i) = real(bzwave(i,jmid))
         fmidim(i) = imag(bzwave(i,jmid))
      end do

      title = 'Bz_wave'
      titll = 'Re Bz_wave (T)'
      titlr = 'Im Bz_wave (T)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim, &
         nnodex, nxmx)


!     -------------------
!     plot fpsi in 2-D
!     -------------------

      title = 'Radial force (fpsi0)'
      call ezconc(capr, y, fpsi0, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      title = 'Poloidal force (ftheta0)'
      call ezconc(capr, y, ftheta0, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      title = 'pressure surfaces'
      call ezconc(capr, y, pressi, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)

      numb = 20

      title = 'dl/B surfaces'
      call ezconc(capr, y, dldb_tot12, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary(capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)


      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(eplus_flux_plot(i,j))
            fimag(i,j) = imag(eplus_flux_plot(i,j))
         end do
      end do


      title = 'E_plus_flux_plot'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)




      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(eminus_flux_plot(i,j))
            fimag(i,j) = imag(eminus_flux_plot(i,j))
         end do
      end do


      title = 'E_minus_flux_plot'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)





      do i = 1, nnodex
         do j = 1, nnodey
            freal(i,j) = real(xkperp_flux_plot(i,j))
            fimag(i,j) = imag(xkperp_flux_plot(i,j))
         end do
      end do


      title = 'xkperp_flux_plot'
      call ezconc(capr, y, freal, ff, nnodex, nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity, iflag)
      if (iflag .eq. 0) call boundary (capr, y, rho, ff, nnodex, &
         nnodey, numb, &
         nxmx, nymx, nlevmax, title, titx, tity)



!     --------------------------------
!     plot eplus_flux_plot(i,j) in midplane
!     --------------------------------
      do i = 1, nnodex
         fmidre(i) = real(eplus_flux_plot(i,jmid))
         fmidim(i) = imag(eplus_flux_plot(i,jmid))
      end do

      title = 'eplus_flux_plot'
      titll = 'Re Eplus (V/m)'
      titlr = 'Im Eplus (V/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim, &
         nnodex, nxmx)


!     --------------------------------
!     plot eminus_flux_plot(i,j) in midplane
!     --------------------------------
      do i = 1, nnodex
         fmidre(i) = real(eminus_flux_plot(i,jmid))
         fmidim(i) = imag(eminus_flux_plot(i,jmid))
      end do

      title = 'eminus_flux_plot'
      titll = 'Re Eminu s(V/m)'
      titlr = 'Im Eminus (V/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim, &
         nnodex, nxmx)



!     --------------------------------------
!     plot xkperp_flux_plot(i,j) in midplane
!     --------------------------------------
      do i = 1, nnodex
         fmidre(i) = real(xkperp_flux_plot(i,jmid))
         fmidim(i) = imag(xkperp_flux_plot(i,jmid))
      end do

      title = 'xkperp_flux_plot'
      titll = 'Re xkperp s(1/m)'
      titlr = 'Im xkperp (1/m)'
      titlb = 'R (m)'

      call ezplot2q(title, titll, titlr, titlb, capr, fmidre, fmidim, &
         nnodex, nxmx)



      do n = 1, nnoderho
!         bmod_plot(n) = bmod_flux(n, 1) / bmod_flux(1, 1)
          bmod_plot(n) = bmod_flux(n, 1)
      end do

      title= 'Mod B (theta = 0)'
      titll= 'bmod_flux (T)'
      titlr='       '
      call ezplot0(title, titll, titlr, rhon, bmod_plot, &
          nnoderho, nrhomax)





!DLG:   Change the species for mchoi files here (orbit-rf)

!     ---------------------
!     write mchoi (65) file
!     ---------------------
      write(65, 309) nnoderho, mnodetheta
      do n = 1, nnoderho
         do m = 1, mnodetheta
          write (65, 1313) n, m, rhon(n), thetam(m), &
               capr_flux(n, m), capz_flux(n, m), eplus_flux(n,m), &
                                                eminus_flux(n,m), &
                                                xkperp_flux(n,m)

          freal(n,m) = real(eplus_flux(n,m))
            fimag(n,m) = imag(eplus_flux(n,m))
         end do
      end do

!     ---------------------
!     write mchoi2 (67) file
!     ---------------------
      write(67, 309) nnoderho_half
      do n = 1, nnoderho_half
         write (67, 1312) n, rhon_half_save(n), wdoti2_ql(n)
      end do


      title = 'E_plus_flux'
      titx   = 'rho'
      tity = 'theta'
      call ezconc(rhon, thetam, freal, ff, nnoderho, mnodetheta, numb, &
         nrhomax, nthetamax, nlevmax, title, titx, tity, iflag)


!     ---------------------
!     write mchoi3 (69) file
!     ---------------------

      write(69, 309) nnodex, nnodey
      do i = 1, nnodex
         do j = 1, nnodey
          write (69, 1313) i, j, capr(i), y(j), wdoti2(i,j)
         end do
      end do



      do n = 1, nnoderho
         do m = 1, mnodetheta
          freal(n,m) = real(eminus_flux(n,m))
            fimag(n,m) = imag(eminus_flux(n,m))
         end do
      end do

      title = 'E_minus_flux'
      titx   = 'rho'
      tity = 'theta'
      call ezconc(rhon, thetam, freal, ff, nnoderho, mnodetheta, numb, &
         nrhomax, nthetamax, nlevmax, title, titx, tity, iflag)



      write(140, 3856)
 3856 format('SCALARS fpsi0 float 1')
      write(140, 2849)
      write (140, 3411) ((fpsi0(i,j), i = 1, nnodex),  j = 1, nnodey)

      write(140, 3858)
 3858 format('SCALARS ftheta0 float 1')
      write(140, 2849)
      write (140, 3411) ((ftheta0(i,j), i = 1, nnodex),  j = 1, nnodey)

 5000 continue

!DLG
    call check ( nf90_close ( nc_id ) )
! Close the graphics device.

      call pgclos

      close (140)


      close (38)

      close (141)
      close (142)

      close (51)
!      close (61)
      close (66)
      close (68)
      close (57)
      close (62)
      close (65)
      close (67)
      close (64)

      if (ndisti2 .eq. 1) then
         deallocate( UPERP )
         deallocate( UPARA )

         deallocate( bqlavg_i1 )
         deallocate( cqlavg_i1 )
         deallocate( eqlavg_i1 )
         deallocate( fqlavg_i1 )

         deallocate( bqlavg_i1_2d )

         deallocate( E_kick_2d  )
       deallocate( f_cql_cart)
       deallocate( f_cql_cart_2d)

       deallocate (wperp_cql)
       deallocate (Wpar_cql)
      end if


  310 format(1p6e12.4)
  309 format(10i10)
  311 format(1p10e12.4)
 3310 format(1p6e18.10)
 1312 format(i10,1p8e12.4)
 1313 format(2i10,1p10e12.4)
 9310 format(1p7e12.4)
   10 format(i10,1p4e10.3,i10,1pe10.3)

      return
      end subroutine fieldws
!
!********************************************************************
!

! Subroutine a1mnmx
!----------------------------------------------------------------------------!
! Minimum and the maximum elements of a 1-d array.
!----------------------------------------------------------------------------!

      subroutine a1mnmx(f, nxmax, nx, fmin, fmax)

      implicit none

      integer nx, i, nxmax
      real f(nxmax), fmin, fmax

      fmax=f(1)
      fmin=fmax
         do 23002 i=1, nx
            fmax=amax1(fmax, f(i))
            fmin=amin1(fmin, f(i))
23002    continue

      return
 2201 format(2i5,1p8e12.4)
      end subroutine a1mnmx

!
!********************************************************************
!
! Subroutine a2mnmx
!----------------------------------------------------------------------------!
! Minimum and the maximum elements of a 1-d array between x1 and x2
!----------------------------------------------------------------------------!

      subroutine a2mnmx(f, nxmax, nx, x, x1, x2, fmin, fmax)

      implicit none

      integer nx, i, nxmax
      real f(nxmax), fmin, fmax, x1, x2
      real x(nxmax)

      fmax=0.0
      fmin=fmax
      do i = 1, nx
         if(x(i).gt.x1.and.x(i).lt.x2)then
            fmax=amax1(fmax, f(i))
            fmin=amin1(fmin, f(i))
         endif
      end do

      return

      end subroutine a2mnmx

!
!********************************************************************
!
!----------------------------------------------------------------------------!
! Subroutine a2dmnmx_r4
!----------------------------------------------------------------------------!
! Minimum and the maximum elements of a 2-d array.
!----------------------------------------------------------------------------!

      subroutine a2dmnmx_r4(f, nxmax, nymax, nx, ny, fmin, fmax)

      implicit none

      integer nx, ny, i, j, nxmax, nymax
      real f(nxmax, nymax), fmin, fmax

      fmax=f(2, 1)
      fmin=fmax

      do 23000 j=1, ny
!         do 23002 i=1, nx
         do 23002 i=2, nx-1
            fmax=amax1(fmax, f(i, j))
            fmin=amin1(fmin, f(i, j))
23002    continue
23000 continue

      return
 2201 format(2i5,1p8e12.4)
      end subroutine a2dmnmx_r4
!
!********************************************************************
!
      subroutine ezconc_new(r, theta, f, flevel, nr, nth, nlevel, &
         nrmax, nthmax, nlevmax, title, titx, tity, iflag)

      implicit none

      integer ncolln2, ncollin, ncolbox, ncollab
      integer nxmx, imark, n0, i, nnode, ibackground, iant, nndm1

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax, &
         ypmin, theta, r, flevel, f, xmin, ymin, xmax
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr, &
         nrmax, nlevel, nthmax, nlevmax, &
         iflag

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character*8 xopt,yopt
      character*32 title
      character*32 titx
      character*32 tity

      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)

!--set up contour levels

      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)

      iflag = 0
       if(fmax .eq. 0.0 .and. fmin .eq. 0.0)then
         iflag = 1
         return
      end if



      df = abs(fmax - fmin) / float(nlevel)
      do i = 1, nlevel
!        flevel(i) = fmin + (i - 0.5) * df / 1000.
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if(nlevel.eq.1)flevel(1)=1.0e-03

      if(nlevel.eq.4)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
      endif

! Split contours into two parts, f > 0, and f < 0.
! Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
! Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0
      do i = 1, nlevel
         if(flevel(i) .gt. 0.) go to 100
         nlevlt = nlevlt + 1
      end do
  100 continue
      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt


!--Advance graphics frame and get ready to plot
!      call pgpage
!     call pladv(0)
!     call plcol(ncolbox)

!     call plenv(xmin, xmax, ymin, ymax, 1, 0)

!--Scale window to user coordinates
!--Make a bit larger so the boundary doesn't get clipped
      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps
!     call pllsty(1)
!     call plvpas(0.15,0.85,0.15,0.85,0.0)
!     call plwind(xpmin,xpmax,ypmin,ypmax)
      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0
!     call plbox(xopt,xtick,nxsub,yopt,ytick,nysub)

! Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then

!        call plwid(0.1)
!     call pllsty(4)
!     call plcol(ncolln2)
!        call plcol(nyellow)
!      call pgsci(nyellow)

!     call plcon1(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt),
!     &      nlevlt, r, theta)
      endif

      if(nlevgt .gt. 0) then

!        call plwid(10)
!     call pllsty(1)
!     call plcol(ncollin)
!        call plcol(nblue)
!      call pgsci(nblue)

!     call plcon1(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt),
!     &      nlevgt, r, theta)
      endif

!     call pllsty(1)
!     call plcol(ncollab)
!     call pllab(titx,tity,title)
      return
      end subroutine ezconc_new

!
!********************************************************************
!
      subroutine ezconc(r, theta, f, flevel, nr, nth, nlevel, &
         nrmax, nthmax, nlevmax, title, titx, tity, iflag)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua, &
         npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen, &
         nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5, &
         ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2, &
         ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant, &
         nndm1, norange

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax, &
         ypmin, theta, r, flevel, f, xmin, ymin, xmax
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr, &
         nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin, &
          ncolelec, iflag

      real tr(6), dx, dy

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character*8 xopt,yopt
      character*32 title
      character*32 titx
      character*32 tity

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      dx = r(2) - r(1)
      dy = theta(2) - theta(1)

      tr(1) = r(1) - dx
      tr(2) = dx
      tr(3) = 0.0
      tr(4) = theta(1) - dy
      tr(5) = 0.0
      tr(6) = dy


      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)

!--set up contour levels

      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)

!      write(6, *)"fmax = ", fmax, "   fmin = ", fmin
!      write(15,*)"fmax = ", fmax, "   fmin = ", fmin

      iflag = 0
       if(fmax .eq. 0.0 .and. fmin .eq. 0.0)then
!         write(6, *)"fmax = ", fmax
!        write(6, *)"fmin = ", fmin
         iflag = 1
         return
      end if

      df = abs(fmax - fmin) / float(nlevel)
      do i = 1, nlevel
!        flevel(i) = fmin + (i - 0.5) * df / 1000.
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if(nlevel.eq.1)flevel(1)=1.0e-03

      if(nlevel.eq.4)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
      endif


      if(nlevel.eq.10)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
       flevel(5)=5.0
         flevel(6)=6.0
         flevel(7)=7.0
         flevel(8)=8.0
       flevel(9)=9.0
         flevel(10)=10.0
      endif

! Split contours into two parts, f > 0, and f < 0.
! Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
! Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0
      do i = 1, nlevel
         if(flevel(i) .gt. 0.) go to 100
         nlevlt = nlevlt + 1
      end do
  100 continue
      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt


!--Scale window to user coordinates
!--Make a bit larger so the boundary doesn't get clipped

      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps

      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0

!--Advance graphics frame and get ready to plot

      call pgsci(nblack)
      call pgenv(xmin, xmax, ymin, ymax, 1, 0)

! Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then

       call pgsci(nyellow)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt), &
             nlevlt, tr)
      endif

      if(nlevgt .gt. 0) then

       call pgsci(nblue)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt), &
             nlevgt, tr)
      endif

      call pgsci(nblack)
      call pglab(titx, tity, title)

  310 format(1p6e12.4)
  312 format(i10, 1p6e12.4)


      return
      end subroutine ezconc
!
!*********************************************************************
!
      subroutine ezconz(r, theta, f, flevel, nr, nth, nlevel, &
         nrmax, nthmax, nlevmax, title, titx, tity)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua, &
         npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen, &
         nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5, &
         ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2, &
         ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant, &
         nndm1

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax, &
         ypmin, theta, r, flevel, f, xmin, ymin, xmax
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr, &
         nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin, &
         ncolelec, j
      real xmaxz, xminz, ymaxz, yminz

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character*8 xopt,yopt
      character*32 title
      character*32 titx
      character*32 tity

      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
          nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
          nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
          ncolion,ncolelec,ncollin,ncollab,ncolln2
      common/zoom/ xmaxz, xminz, ymaxz, yminz


      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)

!--set up contour levels

!      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)


      fmax=f(2, 1)
      fmin=fmax
      do j = 1, nth
         do i = 2, nr - 1
            if (r(i) .gt. xminz .and. r(i) .le. xmaxz .and. &
               theta(j) .gt. yminz .and. theta(j) .le. ymaxz) then
               fmax = amax1(fmax, f(i, j))
               fmin = amin1(fmin, f(i, j))
            end if
         end do
      end do

      if(fmax .eq. 0.0 .and. fmin .eq. 0.0)return


      df = abs(fmax - fmin) / float(nlevel)
      do i = 1, nlevel
!        flevel(i) = fmin + (i - 0.5) * df / 1000.
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if(nlevel .eq. 1)flevel(1) = 1.0e-03

      if(nlevel.eq.4)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
      endif

! Split contours into two parts, f > 0, and f < 0.
! Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
! Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0
      do i = 1, nlevel
         if(flevel(i) .gt. 0.) go to 100
         nlevlt = nlevlt + 1
      end do
  100 continue
      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt


!--Advance graphics frame and get ready to plot
!     call pladv(0)
!     call plcol(ncolbox)

      xmin = xminz
      xmax = xmaxz
      ymax = ymaxz
      ymin = yminz

!     call plenv(xmin, xmax, ymin, ymax, 1, 0)
!--Scale window to user coordinates
!--Make a bit larger so the boundary doesn't get clipped
      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps
!     call pllsty(1)
!     call plvpas(0.15,0.85,0.15,0.85,0.0)
!     call plwind(xpmin,xpmax,ypmin,ypmax)
      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0
!     call plbox(xopt,xtick,nxsub,yopt,ytick,nysub)

! Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then

!        call plwid(0.1)
!     call pllsty(4)
!     call plcol(ncolln2)
!        call plcol(nyellow)

!     call plcon1(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt),
!     &      nlevlt, r, theta)
      endif

      if(nlevgt .gt. 0) then

!        call plwid(10)
!     call pllsty(1)
!     call plcol(ncollin)
!        call plcol(nblue)

!     call plcon1(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt),
!     &      nlevgt, r, theta)
      endif

!     call plcol(ncollab)
!     call pllab(titx,tity,title)
      return
      end subroutine ezconz
!
!*********************************************************************
!
      subroutine ezcon3d(r,theta,f,flevel,nr,nth,nlevel, &
         nrmax,nthmax,nlevmax,title,titx,tity,titz)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua, &
         npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen, &
         nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5, &
         ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2, &
         ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant, &
         nndm1

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax, &
         ypmin, xmin, ymin, xmax
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr, &
         nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin, &
          ncolelec

      real r(nrmax), theta(nthmax)
      real f(nrmax, nthmax)
      real flevel(nlevmax)
      character*8 xopt,yopt
      character*32 title
      character*32 titx
      character*32 tity
      character*32 titz

      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
          nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
          nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
          ncolion,ncolelec,ncollin,ncollab,ncolln2

      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)

!--set up contour levels
      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)

      if(fmax .eq. 0.0 .and. fmin .eq. 0.0)return
      df = abs(fmax - fmin) / float(nlevel)

!      write(6, 100)fmin, fmax
  100 format(1p8e12.4)


!--Advance graphics frame and get ready to plot
!     call plcol(ncolbox)
!     call plenv(-2.5, 2.5, -2.5, 4.0, 0, -2)
!     call plenv(-2.3, 1.8, -1.5, 3.5, 0, -2)

!     call plw3d(2.,4.,3.,xmin,xmax,ymin,ymax,fmin,fmax,30.,30.)
!     call plschr(0.,.50)
!     call plbox3('binstu',titx,0.0,0,'binstu',
!    *     tity,0.0,0,'binstu',titz,0.0,0)
!     call plcol(ngreen)
!     call plmesh(r,theta,f,nr,nth,3,nrmax)
!     call plschr(0.,1.0)
!     call plcol(ncollab)
!     call plmtex('t',1.0,0.5,0.5,title)
      return
      end subroutine ezcon3d

!
!********************************************************************
!
      subroutine boundary(r, theta, f, flevel, nr, nth, nlevel, &
         nrmax, nthmax, nlevmax, title, titx, tity)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua, &
         npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen, &
         nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5, &
         ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2, &
         ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant, &
         nndm1, norange

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax, &
         ypmin, theta, r, flevel, f, xmin, ymin, xmax, rhoplasm
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr, &
         nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin, &
          ncolelec, nlevelb

      real tr(6), dx, dy

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character*8 xopt,yopt
      character*32 title
      character*32 titx
      character*32 tity

      common/boundcom/rhoplasm

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      dx = r(2) - r(1)
      dy = theta(2) - theta(1)

      tr(1) = r(1) - dx
      tr(2) = dx
      tr(3) = 0.0
      tr(4) = theta(1) - dy
      tr(5) = 0.0
      tr(6) = dy


      nlevelb = 1

      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)

!--set up contour levels

      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)

      if(nlevelb.eq.1)then
         flevel(1)= rhoplasm
      endif

! Split contours into two parts, f > 0, and f < 0.
! Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
! Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0
      do i = 1, nlevelb
         if(flevel(i) .gt. 0.) go to 100
         nlevlt = nlevlt + 1
      end do
  100 continue
      ilevgt = ilevlt + nlevlt
      nlevgt = nlevelb - nlevlt

!--Scale window to user coordinates
!--Make a bit larger so the boundary doesn't get clipped

      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps


      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0

!--Advance graphics frame and get ready to plot

      call pgsci(nblack)
      call pgsls(1)

      CALL PGVSTD
      CALL PGWNAD (xmin, xmax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)

! Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt), &
             nlevlt, tr)
      endif

      if(nlevgt .gt. 0) then
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt), &
             nlevgt, tr)
      endif

      return
      end subroutine boundary
!
!*********************************************************************
!
      subroutine ezcon3dv(r,theta,f,flevel,nr,nth,nlevel, &
         nrmax,nthmax,nlevmax,title,titx,tity,titz)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua, &
         npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen, &
         nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5, &
         ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2, &
         ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant, &
         nndm1

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax, &
         ypmin, xmin, ymin, xmax
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr, &
         nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin, &
          ncolelec

      real r(nrmax), theta(nthmax)
      real f(nrmax, nthmax)
      real flevel(nlevmax)
      character*8 xopt,yopt
      character*32 title
      character*32 titx
      character*32 tity
      character*32 titz

      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
          nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
          nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
          ncolion,ncolelec,ncollin,ncollab,ncolln2

      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)

!--set up contour levels
      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)

      if(fmax .eq. 0.0 .and. fmin .eq. 0.0)return
      df = abs(fmax - fmin) / float(nlevel)

!      write(6, 100)fmin, fmax
  100 format(1p8e12.4)


!--Advance graphics frame and get ready to plot
!     call plcol(ncolbox)
!     call plenv(-2.5, 2.5, -2.5, 4.0, 0, -2)
!        call plenv(-2.3, 1.8, -1.5, 3.5, 0, -2)
!     call plw3d(3.,3.,3.,xmin,xmax,ymin,ymax,fmin,fmax,30.,30.)

!     call plschr(0.,.50)
!     call plbox3('binstu',titx,0.0,0,'binstu',
!     *     tity,0.0,0,'binstu',titz,0.0,0)
!     call plcol(ngreen)
!     call plmesh(r,theta,f,nr,nth,3,nrmax)
!     call plschr(0.,1.0)
!     call plcol(ncollab)
!     call plmtex('t',1.0,0.5,0.5,title)
      return
      end subroutine ezcon3dv

!
!********************************************************************
!

      subroutine ezplot3(title, titll, titlr, x1, y1, y2, y3, y4, &
          nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax), y2(nrmax), y3(nrmax), y4(nrmax)
      real y1max, y2max, y3max, y4max
      real y1min, y2min, y3min, y4min
      real ymin, ymax

      character*32 title
      character*32 titll
      character*32 titlr
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
       nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
       nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
       ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen


      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      call a1mnmx(y3, nrmax, nr, y3min, y3max)
      call a1mnmx(y4, nrmax, nr, y4min, y4max)

      ymax = amax1(y1max, y2max, y3max, y4max)
!      ymax = 1.0e+06

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax=ymax*1.1
!      ymax = .1
      ymin=0.0



!     call pladv(0)
!     call pllsty(1)
!     call plcol(ncolbox)
!     call plvpor(0.15,0.85,0.15,0.85)
!     call plwind(rhomin,rhomax,ymin,ymax)
!     call plbox('bcnst',0.0,0,'bnstv',0.0,0)

!     call plcol(nblue)
!     call plline(nr,x1,y1)
!     call plptex(0.65*rhomax,.95*ymax,1.0,0.0,0.0,'minority ions')

!     call plcol(nred)
!     call plline(nr,x1,y2)
!     call plptex(0.65*rhomax,.89*ymax,1.0,0.0,0.0,'electrons')

!     call plcol(ncyan)
!     call plline(nr,x1,y3)
!     call plptex(0.65*rhomax,.83*ymax,1.0,0.0,0.0,'majority ions')

!     call plcol(ngreen)
!     call plline(nr, x1, y4)
!     call plptex(0.65*rhomax,.77*ymax,1.0,0.0,0.0,'ion3')


!     call plcol(ncollin)
!     call plmtex('l',5.0,0.5,0.5,titll)
!     call plmtex('b',3.2,0.5,0.5,'rho')
!     call plmtex('t',2.0,0.5,0.5,title)

!     call pllsty(1)
!     call plcol(ncyan)
!     call plpoin(n12,x2,y3,4)
!     call pllsty(1)
!     call plcol(ncolbox)
!     call plwind(rhomin,rhomax,ymin,ymax)
!     call plbox(' ',0.0,0,'cstv',0.0,0)
!     call plcol(3)
!     call pllsty(1)
!     call plline(nr,x1,y2)
!     call plpoin(n12,x2,y4,4)

!     call pllsty(1)
!     call plmtex('r',5.0,0.5,0.5,titlr)
  300 format (1p9e11.3)
      return
      end subroutine ezplot3

!
!***************************************************************************
!

      subroutine ezplot5p(title, titll, titlr, x1, y1, &
         y2, y3, y4, y5, &
         nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax), y2(nrmax), y3(nrmax), y4(nrmax), &
           y5(nrmax)
      real y1max, y2max, y3max, y4max, y5max
      real y1min, y2min, y3min, y4min, y5min
      real ymin, ymax

      character*32 title
      character*32 titll
      character*32 titlr
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
       nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
       nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
       ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen


      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      call a1mnmx(y3, nrmax, nr, y3min, y3max)
      call a1mnmx(y4, nrmax, nr, y4min, y4max)
      call a1mnmx(y5, nrmax, nr, y5min, y5max)

      ymax = amax1(y1max, y2max, y3max, y4max, y5max)

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax=ymax*1.1
      ymin=0.0



!     call pladv(0)
!     call pllsty(1)
!     call plcol(ncolbox)
!     call plvpor(0.15,0.85,0.15,0.85)
!     call plwind(rhomin,rhomax,ymin,ymax)
!     call plbox('bcnst',0.0,0,'bnstv',0.0,0)

!     call plcol(nblue)
!     call plline(nr,x1,y1)
!     call plptex(0.65*rhomax,.95*ymax,1.0,0.0,0.0,'minority ions')

!     call plcol(nred)
!     call plline(nr,x1,y2)
!     call plptex(0.65*rhomax,.89*ymax,1.0,0.0,0.0,'electrons')

!     call plcol(ncyan)
!     call plline(nr,x1,y3)
!     call plptex(0.65*rhomax,.83*ymax,1.0,0.0,0.0,'majority ions')

!     call plcol(ngreen)
!     call plline(nr, x1, y4)
!     call plptex(0.65*rhomax,.77*ymax,1.0,0.0,0.0,'ion3')

!     call plcol(ngreen)
!     call plline(nr, x1, y5)
!     call plptex(0.65*rhomax,.71*ymax,1.0,0.0,0.0,'fast ions')


!     call plcol(ncollin)
!     call plmtex('l',5.0,0.5,0.5,titll)
!     call plmtex('b',3.2,0.5,0.5,'rho')
!     call plmtex('t',2.0,0.5,0.5,title)

!     call pllsty(1)
!     call plcol(ncyan)
!     call plpoin(n12,x2,y3,4)
!     call pllsty(1)
!     call plcol(ncolbox)
!     call plwind(rhomin,rhomax,ymin,ymax)
!     call plbox(' ',0.0,0,'cstv',0.0,0)
!     call plcol(3)
!     call pllsty(1)
!     call plline(nr,x1,y2)
!     call plpoin(n12,x2,y4,4)

!     call pllsty(1)
!     call plmtex('r',5.0,0.5,0.5,titlr)
  300 format (1p9e11.3)
      return
      end subroutine ezplot5p

!
!***************************************************************************
!


      subroutine ezplot7(title, titll, titlr, titlb, x1, &
         y1, y2, y3, y4, y5, y6, ye, &
         nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax), y2(nrmax), y3(nrmax), y4(nrmax), &
           y5(nrmax), y6(nrmax), ye(nrmax)
      real y1max, y2max, y3max, y4max, y5max, y6max, yemax
      real y1min, y2min, y3min, y4min, y5min, y6min, yemin
      real ymin, ymax

      character*32 title
      character*32 titll
      character*32 titlr
      character*32 titlb

      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd

      integer nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      call a1mnmx(y3, nrmax, nr, y3min, y3max)
      call a1mnmx(y4, nrmax, nr, y4min, y4max)
      call a1mnmx(y5, nrmax, nr, y5min, y5max)
      call a1mnmx(y6, nrmax, nr, y6min, y6max)
      call a1mnmx(ye, nrmax, nr, yemin, yemax)

      ymax = amax1(y1max, y2max, y3max, y4max, y5max, y6max, yemax)
      ymin = amin1(y1min, y2min, y3min, y4min, y5min, y6min, yemin)
!      ymin=0.0

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax=ymax*1.1

!      ymax = 2.8e+06

! Advance plotter to a new page, define coordinate range of graph and draw axes
!      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)

! Label the axes (note use of \u and \d for raising exponent).

      call pglab(titlb, titll, title)

! Plot the line graph.

      call pgsci(nred)
      call pgline(nr, x1, ye)

      call pgsci(ncyan)
      call pgline(nr, x1, y1)

      call pgsci(nblue)
      call pgline(nr, x1, y2)

      call pgsci(ngreen)
      call pgline(nr, x1, y3)

      call pgsci(nyellow)
      call pgline(nr, x1, y4)

      call pgsci(nmagenta)
      call pgline(nr, x1, y5)

      call pgsci(norange)
      call pgline(nr, x1, y6)

      call pgsci(nblack)

  300 format (1p9e11.3)

      return
      end subroutine ezplot7

!
!***************************************************************************
!

      subroutine ezplot70(title, titll, titlr, titlb, x1, &
         y1, y2, y3, y4, y5, y6, ye, &
         nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax), y2(nrmax), y3(nrmax), y4(nrmax), &
           y5(nrmax), y6(nrmax), ye(nrmax)
      real y1max, y2max, y3max, y4max, y5max, y6max, yemax
      real y1min, y2min, y3min, y4min, y5min, y6min, yemin
      real ymin, ymax

      character*32 title
      character*32 titll
      character*32 titlr
      character*32 titlb

      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd

      integer nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      call a1mnmx(y3, nrmax, nr, y3min, y3max)
      call a1mnmx(y4, nrmax, nr, y4min, y4max)
      call a1mnmx(y5, nrmax, nr, y5min, y5max)
      call a1mnmx(y6, nrmax, nr, y6min, y6max)
      call a1mnmx(ye, nrmax, nr, yemin, yemax)

      ymax = amax1(y1max, y2max, y3max, y4max, y5max, y6max, yemax)
      ymin = amin1(y1min, y2min, y3min, y4min, y5min, y6min, yemin)
      ymin=0.0

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax=ymax*1.1

!      ymax = 2.8e+06

! Advance plotter to a new page, define coordinate range of graph and draw axes
!      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)

! Label the axes (note use of \u and \d for raising exponent).

      call pglab(titlb, titll, title)

! Plot the line graph.

      call pgsci(nred)
      call pgline(nr, x1, ye)

      call pgsci(ncyan)
      call pgline(nr, x1, y1)

      call pgsci(nblue)
      call pgline(nr, x1, y2)

      call pgsci(ngreen)
      call pgline(nr, x1, y3)

      call pgsci(nyellow)
      call pgline(nr, x1, y4)

      call pgsci(nmagenta)
      call pgline(nr, x1, y5)

      call pgsci(norange)
      call pgline(nr, x1, y6)

      call pgsci(nblack)

  300 format (1p9e11.3)

      return
      end subroutine ezplot70

!
!***************************************************************************
!
      subroutine ezplot2(title, titll, titlr, titlb, &
                                             x1, y1, y2, nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax),y1(nrmax), y2(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real ymin,ymax

      character*32 title
      character*32 titll
      character*32 titlr
      character*32 titlb

      integer nr, nrmax, n

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      call a1mnmx(y1, nrmax, nr, y1min, y1max)

!      write(6, 300) y2min, y2max
      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return

      ymax = y1max
      ymin = y1min

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax = ymax * 1.1
      ymin = ymin

      if (ymin .le. 0.0) ymin = ymin * 1.1

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BNST', 0.0, 0)
      call pgmtxt('b', 3.2, 0.5, 0.5, titlb)
      call pgmtxt('t', 2.0, 0.5, 0.5, title)

      CALL PGSCI(nblue)
      call pgmtxt('l', 2.0, 0.5, 0.5, titll)


      call pgline(nr, x1, y1)

      call a1mnmx(y2, nrmax, nr, y2min, y2max)
      y2min = y2min * 1.1
      y2max = y2max * 1.1


      CALL PGSWIN (rhomin, rhomax, y2min, y2max)
      CALL PGSCI(nblack)
      CALL PGBOX  (' ', 0.0, 0, 'CMST', 0.0, 0)

      CALL PGSCI(ngreen)
      call pgline(nr, x1, y2)
      call PGMTXT ('r', 2.0, 0.5, 0.5, titlr)

      CALL PGSCI(nblack)


  300 format (1p9e11.3)
      return
      end subroutine ezplot2
!
!***************************************************************************
!

      subroutine ezplot2p(title, titll, titlr, titlb, &
                                             x1, y1, y2, nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax),y1(nrmax), y2(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real ymin,ymax

      character*32 title
      character*32 titll
      character*32 titlr
      character*32 titlb

      integer nr, nrmax, n

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      call a1mnmx(y2, nrmax, nr, y2min, y2max)

      ymax = amax1(y1max, y2max)
      ymin = amin1(y1min, y2min)

      if(ymax .eq. 0.0 .and. ymin .eq. 0.0)return

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax = ymax * 1.1
      ymin = ymin

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)

      call pgmtxt('b', 3.2, 0.5, 0.5, titlb)
      call pgmtxt('t', 2.0, 0.5, 0.5, title)
      call pgmtxt('l', 2.0, 0.5, 0.5, titll)

      CALL PGSCI(nred)
      call pgline(nr, x1, y1)

      CALL PGSCI(nblue)
      call pgline(nr, x1, y2)

      CALL PGSCI(nblack)

  300 format (1p9e11.3)
      return
      end subroutine ezplot2p
!
!***************************************************************************
!


      subroutine ezplot2q(title, titll, titlr, titlb, &
                                             x1, y1, y2, nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax),y1(nrmax), y2(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real ymin,ymax

      character*32 title
      character*32 titll
      character*32 titlr
      character*32 titlb

      integer nr, nrmax, n

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      call a1mnmx(y2, nrmax, nr, y2min, y2max)

      ymax = amax1(y1max, y2max)
      ymin = amin1(y1min, y2min)

      if(ymax .eq. 0.0 .and. ymin .eq. 0.0)return

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax = ymax * 1.1
      ymin = ymin

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)

      call pgmtxt('b', 3.2, 0.5, 0.5, titlb)
      call pgmtxt('t', 2.0, 0.5, 0.5, title)
      call pgmtxt('l', 2.0, 0.5, 0.5, titll)

      CALL PGSCI(nblue)
      call pgline(nr, x1, y1)

      CALL PGSCI(ngreen)
      call pgline(nr, x1, y2)

      CALL PGSCI(nblack)

  300 format (1p9e11.3)
      return
      end subroutine ezplot2q
!
!***************************************************************************
!

      subroutine ezplot4(title, titll, titlr, x1, y1, y2, y3, y4, &
         nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax),y1(nrmax),y2(nrmax),y3(nrmax),y4(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real y4max, y4min
      real ymin,ymax

      character*32 title
      character*32 titll
      character*32 titlr
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
       nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
       nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
       ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen


      call a1mnmx(y1,nrmax,nr,y1min,y1max)
      call a1mnmx(y2,nrmax,nr,y2min,y2max)
      call a1mnmx(y3,nrmax,nr,y3min,y3max)
      call a1mnmx(y4,nrmax,nr,y4min,y4max)

      ymax = amax1(y1max,y2max,y3max,y4max)
      ymin = amin1(y1min,y2min,y3min,y4min)

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax = ymax * 1.1
      ymin = ymin * 1.1



!     call pladv(0)
!     call pllsty(1)
!     call plcol(ncolbox)
!     call plvpor(0.15,0.85,0.15,0.85)
!     call plwind(rhomin,rhomax,ymin,ymax)
!     call plbox('bcnst',0.0,0,'bnstv',0.0,0)

!     call plcol(nblue)
!     call plline(nr,x1,y1)
!     call plptex(0.65*rhomax,.95*ymax,1.0,0.0,0.0,'minority ions')

!     call plcol(nred)
!     call plline(nr,x1,y2)
!     call plptex(0.65*rhomax,.89*ymax,1.0,0.0,0.0,'electrons')

!     call plcol(ncyan)
!     call plline(nr,x1,y3)
!     call plptex(0.65*rhomax,.83*ymax,1.0,0.0,0.0,'majority ions')

!     call plcol(ngreen)
!     call plline(nr,x1,y4)
!     call plptex(0.65*rhomax,.77*ymax,1.0,0.0,0.0,'fast ions')

!     call plcol(ncollin)
!     call plmtex('l',5.0,0.5,0.5,titll)
!     call plmtex('b',3.2,0.5,0.5,'rho')
!     call plmtex('t',2.0,0.5,0.5,title)

!     call pllsty(1)
!     call plcol(ncyan)
!     call plpoin(n12,x2,y3,4)
!     call pllsty(1)
!     call plcol(ncolbox)
!     call plwind(rhomin,rhomax,ymin,ymax)
!     call plbox(' ',0.0,0,'cstv',0.0,0)
!     call plcol(3)
!     call pllsty(1)
!     call plline(nr,x1,y2)
!     call plpoin(n12,x2,y4,4)

!     call pllsty(1)
!     call plmtex('r',5.0,0.5,0.5,titlr)
  300 format (1p9e11.3)
      return
      end subroutine ezplot4
!
!***************************************************************************
!
      subroutine ezplot5(title,titll,titlr,x1, &
         y1,y2,y3,y4,y5,y6,nr,nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax),y1(nrmax),y2(nrmax),y3(nrmax)
      real y4(nrmax),y5(nrmax),y6(nrmax)
      real y1max,y2max,y3max,y4max,y5max,y6max
      real y1min,y2min,y3min,y4min,y5min,y6min
      real ymin,ymax

      character*32 title
      character*32 titll
      character*32 titlr
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd
      integer ncolln4,ncolln5,ncolln6
      integer nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
       nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
       nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
       ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen
      ncolln4=nred
      ncolln5=ngrey
      ncolln6=nbrown


      call a1mnmx(y1,nrmax,nr,y1min,y1max)
      call a1mnmx(y2,nrmax,nr,y2min,y2max)
      call a1mnmx(y3,nrmax,nr,y3min,y3max)
      call a1mnmx(y4,nrmax,nr,y4min,y4max)
      call a1mnmx(y5,nrmax,nr,y5min,y5max)
      call a1mnmx(y6,nrmax,nr,y6min,y6max)
      ymax = amax1(y1max,y2max,y3max,y4max,y5max,y6max)

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax=ymax*1.1
      ymin=0.0



!     call pladv(0)
!     call pllsty(1)
!     call plcol(ncolbox)
!     call plvpor(0.15,0.85,0.15,0.85)
!     call plwind(rhomin,rhomax,ymin,ymax)
!     call plbox('bcnst',0.0,0,'bnstv',0.0,0)

!     call plcol(ncollin)
!     call pllsty(3)
!     call plline(nr,x1,y1)
!     call plptex(0.6*rhomax,.95*ymax,1.0,0.0,0.0,'Mode -2')

!     call plcol(ncolln2)
!     call pllsty(2)
!     call plline(nr,x1,y2)
!     call plptex(0.6*rhomax,.90*ymax,1.0,0.0,0.0,'Mode -1')

!     call plcol(ncolln3)
!     call pllsty(1)
!     call plline(nr,x1,y3)
!     call plptex(0.6*rhomax,.85*ymax,1.0,0.0,0.0,'Mode 0')

!     call plcol(ncolln4)
!     call pllsty(2)
!     call plline(nr,x1,y4)
!     call plptex(0.6*rhomax,.80*ymax,1.0,0.0,0.0,'Mode 1')

!     call plcol(ncolln5)
!     call pllsty(3)
!     call plline(nr,x1,y5)
!     call plptex(0.6*rhomax,.75*ymax,1.0,0.0,0.0,'Mode 2')



!     call plcol(ncollin)
!     call plmtex('l',5.0,0.5,0.5,titll)
!     call plmtex('b',3.2,0.5,0.5,'rho')
!     call plmtex('t',2.0,0.5,0.5,title)

!     call pllsty(1)
!     call plcol(ncyan)
!     call plpoin(n12,x2,y3,4)
!     call pllsty(1)
!     call plcol(ncolbox)
!     call plwind(rhomin,rhomax,ymin,ymax)
!     call plbox(' ',0.0,0,'cstv',0.0,0)
!     call plcol(3)
!     call pllsty(1)
!     call plline(nr,x1,y2)
!     call plpoin(n12,x2,y4,4)

!     call pllsty(1)
!     call plmtex('r',5.0,0.5,0.5,titlr)
  300 format (1p9e11.3)
      return
      end subroutine ezplot5

!
!***************************************************************************
!
      subroutine ezplot1(title, titll, titlr, x1, y1, nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real ymin,ymax

      character*32 title
      character*32 titll
      character*32 titlr
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
       nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
       nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
       ncolion,ncolelec,ncollin,ncolln2,ncollab

      ncolln3=ngreen

      call a1mnmx(y1,nrmax,nr,y1min,y1max)
      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return

      ymax = y1max
      ymin = y1min

      rhomax = x1(nr)
      rhomin = x1(1)
      ymax = ymax * 1.1
      ymin = ymin

      if (ymin .le. 0.0) ymin = ymin * 1.1

! Advance plotter to a new page, define coordinate range of graph and draw axes

!      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)


! Label the axes (note use of \u and \d for raising exponent).

      call pglab('rho', titll, title)

! Plot the line graph.

      call pgline(nr, x1, y1)

  300 format (1p9e11.3)

      return
      end subroutine ezplot1
!
!***************************************************************************
!

      subroutine ezplot0(title, titll, titlr, x1, y1, nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real ymin,ymax

      character*32 title
      character*32 titll
      character*32 titlr
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
       nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
       nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
       ncolion,ncolelec,ncollin,ncolln2,ncollab

      ncolln3=ngreen

      call a1mnmx(y1,nrmax,nr,y1min,y1max)
      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return

      ymax = y1max
      ymin = y1min

      rhomax = x1(nr)
      rhomin = x1(1)
!      ymax = ymax * 1.1
      ymin = ymin

      if (ymin .le. 0.0) ymin = ymin * 1.1

!      ymin = 0.0

! Advance plotter to a new page, define coordinate range of graph and draw axes

!      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)


! Label the axes (note use of \u and \d for raising exponent).

      call pglab('rho', titll, title)

! Plot the line graph.

      call pgline(nr, x1, y1)

  300 format (1p9e11.3)

      return
      end subroutine ezplot0
!
!***************************************************************************
!


      subroutine ezplot1q(title, titll, titlr, titlb, x1, y1, nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real ymin,ymax

      character*32 title
      character*32 titll
      character*32 titlr
      character*32 titlb

      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
       nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
       nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
       ncolion,ncolelec,ncollin,ncolln2,ncollab

      ncolln3=ngreen

      call a1mnmx(y1,nrmax,nr,y1min,y1max)
      if(y1max .eq. 0.0 .and. y1min .eq. 0.0)return

      ymax = y1max
      ymin = y1min

      rhomax = x1(nr)
      rhomin = x1(1)
      ymax = ymax * 1.1
      ymin = ymin

      if (ymin .le. 0.0) ymin = ymin * 1.1

! Advance plotter to a new page, define coordinate range of graph and draw axes

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)
      CALL PGBOX  ('BCNST', 0.0, 0, 'BCNST', 0.0, 0)


! Label the axes (note use of \u and \d for raising exponent).

      call pglab(titlb, titll, title)

! Plot the line graph.

      call pgline(nr, x1, y1)

  300 format (1p9e11.3)

      return
      end subroutine ezplot1q
!
!***************************************************************************
!


      subroutine ezplot1p(title, titll, titlr, x1, y1, nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax),y1(nrmax)
      real y1max,y2max,y3max,y1min,y2min,y3min
      real ymin,ymax

      character*32 title
      character*32 titll
      character*32 titlr
      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd
      integer nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3
      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
       nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
       nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
       ncolion,ncolelec,ncollin,ncolln2,ncollab


      ncolln3=ngreen


      call a1mnmx(y1, nrmax, nr, y1min, y1max)
      if(y1max .eq. 0.0 .and. y2min .eq. 0.0)return

      ymax = y1max
      ymin = y1min

      rhomax=x1(nr)
      rhomin=x1(1)
      ymax = ymax * 1.1
      ymin = ymin

      if (ymin .le. 0.0) ymin = ymin * 1.1

!     call pladv(0)
!     call pllsty(1)
!     call plcol(ncolbox)
!     call plvpor(0.15,0.85,0.15,0.85)
!     call plwind(rhomin,rhomax,ymin,ymax)
!     call plbox('bcnst',0.0,0,'bnstv',0.0,0)

!     call plcol(ncollin)
!     call plline(nr,x1,y1)


!     call plcol(ncollin)
!     call plmtex('l',5.0,0.5,0.5,titll)
!     call plmtex('b',3.2,0.5,0.5,'rho')
!     call plmtex('t',2.0,0.5,0.5,title)

!     call pllsty(1)
!     call plcol(ncyan)

!     call pllsty(1)
!     call plcol(ncolbox)
!     call plwind(rhomin,rhomax,ymin,ymax)
!     call plbox(' ',0.0,0,'cstv',0.0,0)


  300 format (1p9e11.3)
      return
      end subroutine ezplot1p
!
!***************************************************************************
!
      subroutine ezconcx(r, theta, f, flevel, nr, nth, nlevel0, &
         nrmax, nthmax, nlevmax, title, titx, tity, iflag)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua, &
         npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen, &
         nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5, &
         ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2, &
         ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant, &
         nndm1, norange

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax, &
         ypmin, theta, r, flevel, f, xmin, ymin, xmax, fact
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr, &
         nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin, &
          ncolelec, iflag, nlevel0

      real tr(6), dx, dy

      dimension r(nrmax),theta(nthmax)
      dimension f(nrmax, nthmax)
      dimension flevel(nlevmax)
      character*8 xopt,yopt
      character*32 title
      character*32 titx
      character*32 tity

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      dx = r(2) - r(1)
      dy = theta(2) - theta(1)

      tr(1) = r(1) - dx
      tr(2) = dx
      tr(3) = 0.0
      tr(4) = theta(1) - dy
      tr(5) = 0.0
      tr(6) = dy


      call a1mnmx(r, nrmax, nr, xmin, xmax)
      call a1mnmx(theta, nthmax, nth, ymin, ymax)

!--set up contour levels

      call a2dmnmx_r4(f, nrmax, nthmax, nr, nth, fmin, fmax)

      write(6, *)"fmax = ", fmax, "   fmin = ", fmin
      write(15,*)"fmax = ", fmax, "   fmin = ", fmin

      iflag = 0
       if(fmax .eq. 0.0 .and. fmin .eq. 0.0)then
!         write(6, *)"fmax = ", fmax
!        write(6, *)"fmin = ", fmin
         iflag = 1
         return
      end if

      df = abs(fmax - fmin) / float(nlevel0)
      do i = 1, nlevel0
!        flevel(i) = fmin + (i - 0.5) * df / 1000.
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if(nlevel0.eq.1)flevel(1)=1.0e-03

      if(nlevel0.eq.4)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
      endif





      if (nlevel0 .eq. 21) then

!         nlevel = nlevel0 * 2.0

!         fact = 2.0
!         i = 1
!         flevel(i) = fmin

!         do i = 2, nlevel / 2
!            fact = fact * 1.2
!            flevel(i) = flevel(i-1) / fact
!         end do

!         fact = 2.0
!         i = nlevel / 2 + 1
!         flevel(i) = fmax

!         do i = nlevel / 2 + 2, nlevel
!            fact = fact * 1.2
!            flevel(i) = flevel(i-1) / fact
!         end do

         nlevel = nlevel0
         fact = 2.0
         i = 1
         flevel(i) = fmax

         do i = 2, nlevel
            fact = fact * 1.2
            flevel(i) = flevel(i-1) / fact
         end do

      else

         nlevel = nlevel0

      end if








! Split contours into two parts, f > 0, and f < 0.
! Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
! Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0
      do i = 1, nlevel
         if(flevel(i) .gt. 0.) go to 100
         nlevlt = nlevlt + 1
      end do
  100 continue
      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt


!--Scale window to user coordinates
!--Make a bit larger so the boundary doesn't get clipped

      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps

      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0

!--Advance graphics frame and get ready to plot

      call pgsci(nblack)
      call pgenv(xmin, xmax, ymin, ymax, 1, 0)

! Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then

       call pgsci(nyellow)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevlt), &
             nlevlt, tr)
      endif

      if(nlevgt .gt. 0) then

       call pgsci(nblue)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth, flevel(ilevgt), &
             nlevgt, tr)
      endif

      call pgsci(nblack)
      call pglab(titx, tity, title)

  310 format(1p6e12.4)
  312 format(i10, 1p6e12.4)


      return
      end subroutine ezconcx
!
!*********************************************************************
!

      subroutine ezlog1(title, titll, titlr, titlb, x1, y1, &
         nr, nrmax, ymin, ymax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax)
      real y1max, y2max, y3max, y4max, y5max, y6max
      real y1min, y2min, y3min, y4min, y5min, y6min
      real ymin, ymax, tmin, tmax

      character*32 title
      character*32 titll
      character*32 titlr
      character*32 titlb

      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd

      integer nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      rhomax=x1(nr)
      rhomin=x1(1)

      call a1mnmx(y1, nrmax, nr, ymin, ymax)

! Advance plotter to a new page, define coordinate range of graph and draw axes
!      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)

      CALL PGBOX  ('BCNST', 0.0, 0, 'LBCNST', 0.0, 0)

      call pgmtxt('b', 3.2, 0.5, 0.5, titlb)
      call pgmtxt('t', 2.0, 0.5, 0.5, title)
      call pgmtxt('l', 2.0, 0.5, 0.5, titll)

      call pgsci(nblue)
      call pgline(nr, x1, y1)



      CALL PGSCI(nblack)
      call pgsls(1)

  300 format (1p9e11.3)
      return
      end subroutine ezlog1

!
!***************************************************************************
!

      subroutine ezlog1_f(title, titll, titlr, titlb, x1, y1, &
         y2, y3, y4, y5, y6, y7, &
         nr, nrmax)

      implicit none

      real xzmax,xzmin,xnmin,xnmax,rhomin,rhomax
      real x1(nrmax), y1(nrmax)
      real y2(nrmax), y3(nrmax), y4(nrmax)
      real y5(nrmax), y6(nrmax), y7(nrmax)
      real y1max, y2max, y3max, y4max, y5max, y6max
      real y1min, y2min, y3min, y4min, y5min, y6min
      real ymin, ymax, tmin, tmax

      character*32 title
      character*32 titll
      character*32 titlr
      character*32 titlb

      integer nr,nrmax

      integer nplot1,ncollab, ncolion,ncolbox, ncyan, &
          ncolelec, ncolln2, ncollin, ncolbrd

      integer nblack,nred,nyellow, ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan1, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolln3, norange

      nwhite = 0
      nblack = 1
      nred = 2
      ngreen = 3
      nblue = 4
      ncyan = 5
      nmagenta = 6
      nyellow = 7
      norange = 8

      rhomax=x1(nr)
      rhomin=x1(1)

      call a1mnmx(y1, nrmax, nr, ymin, ymax)

      ymin = -10.0

! Advance plotter to a new page, define coordinate range of graph and draw axes
!      call pgenv(rhomin, rhomax, ymin, ymax, 0, 0)

      CALL PGPAGE
      CALL PGSVP (0.15,0.85,0.15,0.85)
      CALL PGSWIN (rhomin, rhomax, ymin, ymax)

      CALL PGBOX  ('BCNST', 0.0, 0, 'LBCNST', 0.0, 0)

      call pgmtxt('b', 3.2, 0.5, 0.5, titlb)
      call pgmtxt('t', 2.0, 0.5, 0.5, title)
      call pgmtxt('l', 2.0, 0.5, 0.5, titll)

      call pgsci(nblue)
      call pgline(nr, x1, y1)
      call pgline(nr, x1, y2)
      call pgline(nr, x1, y3)
      call pgline(nr, x1, y4)
      call pgline(nr, x1, y5)
      call pgline(nr, x1, y6)
      call pgline(nr, x1, y7)


      CALL PGSCI(nblack)
      call pgsls(1)

  300 format (1p9e11.3)
      return
      end subroutine ezlog1_f

!
!***************************************************************************
!



      subroutine ezconpx(r, theta, f, flevel, nr, nth, nlevel0, &
         nrmax, nthmax, nlevmax, xg, yg, title, r0)

      implicit none

      integer nxmx, ncolln10, imark, ncolln9, nwheat, ngrey, naqua, &
         npink, nblueviolet, ncyan, nbrown, nblue, nyellow, ngreen, &
         nblack, nred, nturquoise, ncolln6, ncolln7, ncolln4, ncolln5, &
         ncolln8, nmodes, nwhite, ncolbox, nmagenta, nsalmon, ncolln2, &
         ncolln3, ncolbrd, ncolln1, n0, i, nnode, ibackground, iant, &
         nndm1

      real fmin, ymax, df, fmax, eps,xtick, ytick, xpmax, xpmin, ypmax, &
         ypmin, xmin, ymin, xmax
      integer nlevlt, ilevlt, nlevgt, ilevgt, nxsub, nysub, nth, nr, &
         nrmax, nlevel, nthmax, ncollab, ncolion, nlevmax, ncollin, &
         ncolelec, j, n, m, nth1, nlevel0
      real xmaxz, xminz, ymaxz, r0, fact

      real r(nrmax),theta(nthmax)
      real f(nrmax,nthmax)
      real xg(nrmax,nthmax)
      real yg(nrmax,nthmax)
      real flevel(nlevmax)

      character*8 xopt,yopt
      character*32 title


      common/colcom/nblack,nred,nyellow,ngreen,naqua,npink, &
         nwheat,ngrey,nbrown,nblue,nblueviolet,ncyan, &
         nturquoise,nmagenta,nsalmon,nwhite,ncolbox,ncolbrd, &
         ncolion,ncolelec,ncollin,ncollab,ncolln2



!--set up transformation
      do n = 1, nr
         do m = 1, nth
            xg(n,m) = r(n) * cos(theta(m)) + r0
            yg(n,m) = r(n) * sin(theta(m))
         end do
         xg(n, nth+1) = xg(n,1)
         yg(n, nth+1) = yg(n,1)
      end do
      nth1 = nth + 1
      call a2dmnmx(xg, nrmax, nthmax, nr, nth, xmin, xmax)
      call a2dmnmx(yg, nrmax, nthmax, nr, nth, ymin, ymax)



!--set up contour levels
      call a2dmnmx(f, nrmax, nthmax, nr, nth, fmin, fmax)
      df=abs(fmax-fmin)/float(nlevel0)

!        write (6, *) "df = ", df

      if(df .eq. 0.0)return

      do i = 1, nlevel0
         flevel(i) = fmin + (i - 0.5) * df
      end do

      if(nlevel0.eq.1)flevel(1)=1.0e-03
      if(nlevel0.eq.4)then
         flevel(1)=1.0
         flevel(2)=2.0
         flevel(3)=3.0
         flevel(4)=4.0
      endif

      if (nlevel0 .eq. 20) then

         nlevel = nlevel0 * 2.0

         fact = 2.0
         i = 1
         flevel(i) = fmin

         do i = 2, nlevel / 2
            fact = fact * 1.2
            flevel(i) = flevel(i-1) / fact
         end do


         fact = 2.0
         i = nlevel / 2 + 1
         flevel(i) = fmax

         do i = nlevel / 2 + 2, nlevel
            fact = fact * 1.2
            flevel(i) = flevel(i-1) / fact
         end do

      else
         nlevel = nlevel0
      end if


! Split contours into two parts, f > 0, and f < 0.
! Dashed contours will be at levels 'ilevlt' through 'ilevlt+nlevlt'.
! Solid  contours will be at levels 'ilevgt' through 'ilevgt+nlevgt'.

      ilevlt = 1
      nlevlt = 0

      do i = 1, nlevel
         if(flevel(i) .gt. 0.) go to 100
         nlevlt = nlevlt + 1
      end do

  100 continue

      ilevgt = ilevlt + nlevlt
      nlevgt = nlevel - nlevlt



!--Advance graphics frame and get ready to plot

!      call plcol(ncolbox)
!      call plenv(xmin, xmax, ymin, ymax, 0, 0)

      call pgsci(nblack)
      call pgenv(xmin, xmax, ymin, ymax, 1, 0)

!--Scale window to user coordinates
!--Make a bit larger so the boundary doesn't get clipped
      eps=0.05
      xpmin = xmin - abs(xmin) * eps
      xpmax = xmax + abs(xmax) * eps
      ypmin = ymin - abs(ymin) * eps
      ypmax = ymax + abs(ymax) * eps


!      call pllsty(1)
      xopt='bstn'
      yopt='bstn'
      xtick=0.
      nxsub=0
      ytick=0.
      nysub=0

! Call plotter once for f < 0 (dashed), once for f > 0 (solid lines).

      if(nlevlt .gt. 0) then
!         call plwid(1)
!         call pllsty(2)
!         call plcol(ncolln2)
         call pgsci(nyellow)
         call pgcont(f, nrmax, nthmax, 1, nr, 1, nth1, flevel(ilevlt), &
             nlevlt, xg, yg)
      endif

      if(nlevgt .gt. 0) then
!         call plwid(1)
!         call pllsty(1)
!         call plcol(ncollin)
!         call plcol(nblue)
!         call plcon2(f,nrmax,nthmax,1,nr,1,nth1,flevel(ilevgt),
!     &   nlevgt, xg, yg)
      endif
!      call plwid(1)


!      call plcol(ncollab)
!      call plcol(nblue)
!      call pllab('u_parallel','u_perp', title)


  310 format(1p6e12.4)
  311 format(10i10)

      return
      end subroutine ezconpx
!
!*********************************************************************
!

end module plot_aorsa2dps
