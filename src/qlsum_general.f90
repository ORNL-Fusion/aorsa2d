module qlsum_general

contains

    subroutine qlSum ( k_uper, b_sum, c_sum, e_sum, f_sum, &
        sum_wdot, sum_fx0, sum_fy0, W, ZSPEC, ASPEC, BMAG, &
        lmax, ENORM, UPARMIN, UPARMAX, &
        NUPAR, NUPER, UPER, UPAR, DFDUPER, DFDUPAR,  &
        ealphak, ebetak, ebk, nkdim1, nkdim2, mkdim1, mkdim2,   &
        nkx1, nkx2, nky1, nky2, &
        uxx, uxy, uxz, &
        uyx, uyy, uyz, &
        uzx, uzy, uzz, &
        nxdim, nydim, xkxsav, xkysav, xkphi, xx, yy, i_global, j_global, &
        lmaxdim, ndist, nzeta, &
        gradprlb, bmod, omgc, alpha, xm, upshift, xk_cutoff, &
        maxwellian )

    use qlsum_deps
    use constants

    implicit none

    logical, optional, intent(IN) :: maxwellian
    logical :: doMaxwellian

    integer, intent(IN):: NUPAR, NUPER, lmax
    real, intent(IN):: W, ZSPEC, ASPEC, BMAG
    real, intent(IN):: ENORM, UPARMIN, UPARMAX
    real, dimension(NUPER), intent(IN):: UPER
    real, dimension(NUPAR), intent(IN):: UPAR
    real, dimension(NUPER,NUPAR), intent(IN):: DFDUPER, DFDUPAR
 
    integer i_global, j_global, ier, nxdim, nydim, lmaxdim, ndist
    integer nkx1, nkx2, nky1, nky2, k_uper, l
    integer nkdim1, nkdim2, mkdim1, mkdim2
    integer:: NHARM, IHARM, M, N, i, nzeta
    integer i_uprl, upshift
    integer ires, iresmax

    complex, dimension(:,:), allocatable :: zbeta, zBeta_iHarm
    
    real  y, y0, alpha, xm, akprl, sgn_kprl, omgc, bmod, xkprl_eff,  &
      descrim, fgam, gradprlb, gammab, xk_cutoff

    real  uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz
    real  xkphi, sinth, facte, sinth_inv
    real  xkxsav(nkdim1 : nkdim2), xkysav(mkdim1 : mkdim2)
    real  xkperpn, xkperpni, xkrhon, xketan, xkprln
    real:: W2, WCW, RRP, RRM, WC, WCI
    real:: MUT0, SQMUT0, PISQMUT0, SQMUT0I
    real:: ISQ2, SQ2, NWCW, DFACTPAR, DFACTPER, U0
    real:: UPAR0, dfdupar0, dfduper0, du, dui, p 
    real:: WI, uperpk, uperpk2
    real:: dzeta, dzetai, zetamax, zetamin, zeta0
    real:: A1, A2, A3, u
    real:: factor
    
    real, dimension(:),     allocatable :: zetai
    real, dimension(:,:),   allocatable :: Jni
    real, dimension(:,:,:), allocatable :: Jn
    real, dimension(:,:),   allocatable :: NPARA_sav
    
    integer, dimension(:),  allocatable :: nres, mres
        
    complex, dimension(:),  allocatable :: &
        sumb_11, sumb_31, sumc_11, sumc_31, &
        sume_11, sume_31, sumf_11, sumf_31
    
    complex sumf_11_nm, sumf_31_nm, sume_11_nm, &
        sume_31_nm, sumc_11_nm, sumc_31_nm, &
        sumb_11_nm, sumb_31_nm
 
    logical, dimension(:,:), allocatable :: is_resonance_nm      

    complex epsx, epsy, epsz
    complex ealphak(nkdim1 : nkdim2, mkdim1 : mkdim2), &
            ebetak(nkdim1 : nkdim2, mkdim1 : mkdim2), &
               ebk(nkdim1 : nkdim2, mkdim1 : mkdim2)
    complex zeta

    complex xx(nkdim1 : nkdim2, 1 : nxdim),   &
           yy(mkdim1 : mkdim2, 1 : nydim)

    complex cexp1, cexp2, cexp0
    complex sumwdot_11_nm, sumwdot_31_nm
    
    complex sumwdot_11, sumwdot_31, sumwdotkx_11, &
        sumwdotkx_31, sumwdotky_11, sumwdotky_31, & 
        sum2_1, sum2_2, sum2_3, sumkx2_1, sumkx2_2, sumkx2_3, & 
        sumky2_1, sumky2_2, sumky2_3        
    
    complex, allocatable :: b_sum(:), c_sum(:), e_sum(:), f_sum(:)
    complex sum_wdot, sum_fx0, sum_fy0
    complex, allocatable :: bbb(:)

    real :: cosbeta_n_m, sinbeta_n_m
    
 
    if ( present ( maxwellian ) ) then
        doMaxwellian = maxwellian
    else
        doMaxwellian = .false.
    endif 

    !   Allocate arrays

    allocate ( b_sum ( nuPar ), c_sum ( nuPar ), &
        e_sum ( nuPar ), f_sum ( nuPar) )

    allocate ( bbb ( lMax + 3) )
    allocate( zbeta(nkx1:nkx2,nky1:nky2) )
    allocate( zbeta_iharm(nkx1:nkx2,nky1:nky2) )

    allocate(zetai(nzeta + 1) )
    allocate(Jni(-lmaxdim : lmaxdim, nzeta + 1) )
    allocate( Jn(-lmaxdim : lmaxdim, nkdim1 : nkdim2, mkdim1 : mkdim2))    
    allocate(NPARA_sav(nkdim1 : nkdim2, mkdim1 : mkdim2) ) 
    
    allocate ( nres(nxdim * nydim), mres(nxdim * nydim) )
    
    allocate(sumb_11(nupar), sumb_31(nupar), sumc_11(nupar), &
        sumc_31(nupar), sume_11(nupar), sume_31(nupar), &
        sumf_11(nupar), sumf_31(nupar) )

    allocate(is_resonance_nm(nkx1:nkx2,nky1:nky2))

!   initialize allocatable arrays to zero

    zbeta = 0.0
    zbeta_iharm = 0.0

    sumf_11_nm = 0.0
    sumf_31_nm = 0.0

    sume_11_nm = 0.0
    sume_31_nm = 0.0
    
    sumc_11_nm = 0.0
    sumc_31_nm = 0.0
    
    sumb_11_nm = 0.0
    sumb_31_nm = 0.0

    is_resonance_nm = .false.

    zetai = 0.0
    Jni = 0.0
    NPARA_sav = 0.0

    sumb_11 = 0.0
    sumb_31 = 0.0
      
    sumc_11 = 0.0
    sumc_31 = 0.0

    sume_11 = 0.0
    sume_31 = 0.0

    sumf_11 = 0.0
    sumf_31 = 0.0    
    
    uperpk = uper(k_uper)
    uperpk2 = uperpk**2

    W2 = W * W
    WI = 1.0 / W

    WCW = BMAG * ZSPEC * EOVERAMU / ASPEC / W
    WC = WCW * W
    WCI = 1.0 /WC

    MUT0 = 0.5 * MPC2 * ASPEC / ENORM
    SQMUT0 = SQRT(MUT0)
    SQMUT0I = 1.0 / SQMUT0
    PISQMUT0 = SQMUT0 * pi

    ISQ2 = SQRT(0.5)
    SQ2 = SQRT(2.0)
    NHARM = lmax
    
    du = (upar(nupar) - upar(1)) / (nupar - 1)
    dui = 1.0 / du
    
    if ( nzeta .eq. 1 ) then
    ! -------------------------------------------------------- !
    ! ---Don't interpolate: precalculate all Bessel functions- !
    ! -------------------------------------------------------- !

        kxLoop1: &
        do n = nkx1, nkx2
            kyLoop1: &
            do m = nky1, nky2
    
                xkrhon  = uxx * xkxsav(n) + uxy * xkysav(m) + uxz * xkphi
                xketan  = uyx * xkxsav(n) + uyy * xkysav(m) + uyz * xkphi
                xkprln  = uzx * xkxsav(n) + uzy * xkysav(m) + uzz * xkphi
                xkperpn = sqrt(xkrhon**2 + xketan**2)
         
                ! ------------------------------------
                ! Optional: leave out upshift in xkprl
                ! --------------------------------- --          
                if ( upshift .eq. 0 ) xkprln = uzz * xkphi
          
                if ( upshift .eq. -1 .and. xkperpn .gt. xk_cutoff ) then      
                    xkprln = uzz * xkphi
                endif
         
                sgn_kprl = sign ( 1.0, xkprln )
                akprl = abs ( xkprln )

                y0  = 1.5
                y   = y0
                l   = 1

                if ( xkprln .eq. 0 ) xkprln = 1.0e-06

                gammab = abs(l * omgc / (2.0 * alpha * xkprln**2)  &
                                         * gradprlb / bmod)

                if ( xm .eq. xme ) gammab = 0.0
                if ( abs ( gammab ) .lt. .01 ) gammab = .01

                !   Calculate kPar effective

                if ( sgn_kprl .ge. 0.0 ) then

                    fgam = 1.0

                    if ( gammab .gt. 1.0e-05 ) then
                        y = y0
                        fgam = (sqrt(1. +  4. * gammab * y) - 1.)   &
                            / (2. * gammab * y)
                    endif

                    xkprl_eff = xkprln / fgam 

                else

                    fgam = 1.0

                    if ( gammab .gt. 1.0e-05 ) then
                        descrim = 1. - 4. * gammab * y0
                        if (descrim .ge. 0.0) y =   y0
                        if (descrim .lt. 0.0) y = - y0
                        fgam = (1. - sqrt(1. -  4. * gammab * y) )  &
                            / (2. * gammab * y)
                    endif

                    xkprl_eff = xkprln / fgam 

                endif

                xkprln = xkprl_eff      
    
                NPARA_sav(n, m) = xkprln * C / W

                xkperpn = sqrt(xkrhon**2 + xketan**2)
                if ( xkperpn .eq. 0.0 ) xkperpn = 1.0e-08
    
                cosbeta_n_m  = xkrhon / xkperpn
                sinbeta_n_m  = xketan / xkperpn
                zbeta(n,m) = cmplx ( cosbeta_n_m , sinbeta_n_m  )

                zeta = xkperpn * uper(k_uper) * c * sqmut0i / wc

                call besjc ( zeta, nharm + 2, bbb, ier )
                if ( ier .ne. 0 ) write (6, *) "qlsum_general.f90:268 - besjc ier = ", ier

                do IHARM = 0, NHARM + 1
                    Jn(iharm,  n, m) = bbb(iharm + 1)
                    Jn(-iharm, n, m) = (-1.0)**iharm * Jn(iharm, n, m)
                enddo

          enddo kyLoop1
       enddo kxLoop1
    
    else
    ! -------------------------------------- !
    ! ---Interpolate; calculate zeta mesh -- !
    ! -------------------------------------- !
    
        zetamax = 0.0
        zetamin = 0.0
   
        kxLoop2: & 
        do n = nkx1, nkx2
            kyLoop2: &
            do m = nky1, nky2

                xkrhon = uxx * xkxsav(n) + uxy * xkysav(m) + uxz * xkphi
                xketan = uyx * xkxsav(n) + uyy * xkysav(m) + uyz * xkphi
                xkprln   = uzx * xkxsav(n) + uzy * xkysav(m) + uzz * xkphi

                xkperpn = sqrt(xkrhon**2 + xketan**2) + 1.0e-08
    
                zeta0 = xkperpn * uperpk * c * sqmut0i * wci
      
                if (zeta0 .gt. zetamax) zetamax = zeta0
                if (zeta0 .lt. zetamin) zetamin = zeta0

            enddo kyLoop2
        enddo kxLoop2

        if ( zetamax .eq. zetamin ) then
            zetamax =  1.0e-06
            zetamin = -1.0e-06
        endif

        dzeta = (zetamax - zetamin) / (nzeta - 1)
        dzetai = 1.0 / dzeta
    
        ! ------------------------------------------------- !
        ! ---Pre-calculate Bessel functions on zeta mesh -- !
        ! ------------------------------------------------- !
      
        do i = 1, nzeta + 1
            zetai(i) = zetamin + (i - 1) * dzeta
            zeta = cmplx(zetai(i), 0.0)
      
            call besjc(zeta, nharm + 2, bbb, ier)
      
            do iharm = 0, NHARM + 1
               Jni(iharm,  i) = bbb(iharm + 1)
               Jni(-iharm, i) = (-1.0)**iharm * bbb(iharm + 1)
            enddo
        enddo
     
       ! --------------------------------- !
       ! ---Interpolate Bessel functions-- !
       ! --------------------------------- !

        kxLoop3: &
        do n = nkx1, nkx2
            kyLoop3: &
            do m = nky1, nky2

                xkrhon = uxx * xkxsav(n) + uxy * xkysav(m) + uxz * xkphi
                xketan = uyx * xkxsav(n) + uyy * xkysav(m) + uyz * xkphi
                xkprln   = uzx * xkxsav(n) + uzy * xkysav(m) + uzz * xkphi
                xkperpn = sqrt(xkrhon**2 + xketan**2)
         
                ! ------------------------------------
                ! Optional: leave out upshift in xkprl
                ! --------------------------------- --          
                if(upshift .eq. 0)xkprln = uzz * xkphi
         
                if (upshift .eq. -1) then      
                    if (xkperpn  .gt. xk_cutoff) xkprln = uzz * xkphi
                endif
         
                sgn_kprl = sign(1.0, xkprln)
                akprl = abs(xkprln)

                y0 = 1.5
                y = y0
         
                l = 1
                if(xkprln .eq. 0)xkprln = 1.0e-06
                gammab = abs(l * omgc / (2.0 * alpha * xkprln**2)  &
                            * gradprlb / bmod)

                if(xm .eq. xme)gammab = 0.0
                if(abs(gammab) .lt. .01)gammab = .01

                if(sgn_kprl .ge. 0.0)then
                    fgam = 1.0

                    if ( gammab .gt. 1.0e-05 ) then
                        y = y0
                        fgam = (sqrt(1. +  4. * gammab * y) - 1.)   &
                            / (2. * gammab * y)
                    endif

                    xkprl_eff = xkprln / fgam 
                endif


                if(sgn_kprl .lt. 0.0)then
                    fgam = 1.0

                    if(gammab .gt. 1.0e-05)then
                        descrim = 1. - 4. * gammab * y0
                        if (descrim .ge. 0.0) y =   y0
                        if (descrim .lt. 0.0) y = - y0
                        fgam = (1. - sqrt(1. -  4. * gammab * y) )  &
                          / (2. * gammab * y)
                    endif

                    xkprl_eff = xkprln / fgam 

                endif

                xkprln = xkprl_eff
    
                NPARA_sav(n, m) = xkprln * C * WI

                xkperpn = sqrt(xkrhon**2 + xketan**2) + 1.0e-08
                xkperpni = 1.0 / xkperpn
    
                cosbeta_n_m = xkrhon * xkperpni
                sinbeta_n_m  = xketan * xkperpni
                zbeta(n,m) = cmplx( cosbeta_n_m , sinbeta_n_m )
      
                zeta0 = xkperpn * uperpk * c * sqmut0i * wci
      
                i = int((zeta0 - zetamin) * dzetai) + 1
                p = (zeta0 - zetai(i)) * dzetai
                A1 = 0.5 * P * (P - 1.)
                A2 = 1. - P * P
                A3 = 0.5 * P * (P + 1.)
     
                iHarmLoop2: & 
                do iharm = -NHARM - 1, NHARM + 1
      
                    Jn(iharm, n, m) = Jni(iharm, i)    &
                        + p * (Jni(iharm, i + 1) - Jni(iharm, i))

                    if ( i .ne. 1 ) then
                        Jn(iharm, n, m) = A1 * Jni(iharm, i - 1)     &
                                     + A2 * Jni(iharm, i)         &
                                     + A3 * Jni(iharm, i + 1)
                    endif
          
                enddo iHarmLoop2
      
          enddo kyLoop3
       enddo kxLoop3
    
    endif
                    

    ! ------------------------ !
    ! ---Sum over harmonics--- !
    ! ------------------------ !
    
    sum_wdot = 0.0
    sum_fx0  = 0.0
    sum_fy0  = 0.0    
    
    b_sum = 0.0
    c_sum = 0.0
    e_sum = 0.0
    f_sum = 0.0
      
    iHarmLoop: & 
    do IHARM = -NHARM, NHARM

       NWCW = real(IHARM) * WCW
       
       sumb_11 = 0.0
       sumb_31 = 0.0
      
       sumc_11 = 0.0
       sumc_31 = 0.0    
      
       sume_11 = 0.0
       sume_31 = 0.0
      
       sumf_11 = 0.0
       sumf_31 = 0.0      
              
       sumwdot_11 = 0.0
       sumwdot_31 = 0.0       
       
       sumwdotkx_11 = 0.0
       sumwdotkx_31 = 0.0
       
       sumwdotky_11 = 0.0
       sumwdotky_31 = 0.0       

       sum2_1 = 0.0
       sum2_2 = 0.0
       sum2_3 = 0.0
       
       sumkx2_1 = 0.0
       sumkx2_2 = 0.0
       sumkx2_3 = 0.0
       
       sumky2_1 = 0.0
       sumky2_2 = 0.0
       sumky2_3 = 0.0       
       
       call zpow((nkx2-nkx1+1)*(nky2-nky1+1), zbeta, iharm, zbeta_iharm) 
       
       ! ----------------------- !
       ! ---Find resonant modes--!
       ! ----------------------- !
       ires = 0
       kxLoop4: &
       do n = nkx1, nkx2
            kyLoop4: &
            do m = nky1, nky2     
                ! ------------------------ !
                ! -- Resonance relation -- !
                ! ------------------------ !
                RRP = 1.0 - NWCW - NPARA_sav(n, m) * (UPARMAX) * SQMUT0i
                RRM = 1.0 - NWCW - NPARA_sav(n, m) * (UPARMIN) * SQMUT0i

                is_resonance_nm(n,m) = (RRP * RRM .le. 0.0)
                if (is_resonance_nm(n,m)) then
                    ires = ires + 1
                    nres(ires) = n
                    mres(ires) = m
                endif
            enddo kyLoop4
        enddo kxLoop4

        iresmax = ires     

        ! ---------------------------- !
        ! ---Sum over resonant modes-- !
        ! ---------------------------- !

        !   wDot sums over resonant uPrl only (in contrast to dQL below)

        resonantModes: &
        do ires = 1, iresmax

            n = nres(ires)
            m = mres(ires)
                                                       
            cexp1 = xx(n, i_global) * yy(m, j_global) * zbeta_iharm(n,m) * zbeta(n,m)
            cexp2 = xx(n, i_global) * yy(m, j_global) * zbeta_iharm(n,m) / zbeta(n,m)
            cexp0 = xx(n, i_global) * yy(m, j_global) * zbeta_iharm(n,m)
        
            epsx = isq2 * (ealphak(n, m) - zi * ebetak(n, m)) * cexp1
            epsy = isq2 * (ealphak(n, m) + zi * ebetak(n, m)) * cexp2
            epsz = ebk(n, m) * cexp0
         
            sum2_1 = sum2_1 + conjg(epsx) * Jn(IHARM + 1, n, m)
            sum2_2 = sum2_2 + conjg(epsy) * Jn(IHARM - 1, n, m)
            sum2_3 = sum2_3 + conjg(epsz) * Jn(IHARM, n, m)
         
            sumkx2_1 = sumkx2_1 + xkxsav(n) * conjg(epsx) * Jn(IHARM + 1, n, m)
            sumkx2_2 = sumkx2_2 + xkxsav(n) * conjg(epsy) * Jn(IHARM - 1, n, m)
            sumkx2_3 = sumkx2_3 + xkxsav(n) * conjg(epsz) * Jn(IHARM, n, m)
         
            sumky2_1 = sumky2_1 + xkysav(m) * conjg(epsx) * Jn(IHARM + 1, n, m)
            sumky2_2 = sumky2_2 + xkysav(m) * conjg(epsy) * Jn(IHARM - 1, n, m)
            sumky2_3 = sumky2_3 + xkysav(m) * conjg(epsz) * Jn(IHARM, n, m)                       

            UPAR0   = SQMUT0 / NPARA_sav(n, m) * (1. - NWCW)
    
            u       = sqrt ( upar0**2 + uperpk2 ) + 1.0e-08
            sinth   = uperpk / u + 1.0e-08
            sinth_inv = 1.0 / sinth
        
            facte = (nwcw - sinth**2) /  upar0
    
            i       = int((UPAR0 - UPAR(1)) * dui) + 1
            i_uprl  = i
            p       = (UPAR0 - UPAR(i)) * dui

            dfduper0 = dfduper(k_uper, NUPAR)

            if (i .ne. NUPAR) then
                dfduper0 = dfduper(k_uper, i) + (dfduper(k_uper, i+1) - dfduper(k_uper, i)) * p
            endif
     
            maxwellianOrNot: & 
            if ( doMaxwellian ) then

                   U0   = DFDUPER0

            else 

                DFACTPAR = NPARA_sav(n, m) * UPAR0 * SQMUT0I
                DFACTPER = NPARA_sav(n, m) * UPER(k_uper) * SQMUT0I
      
                dfdupar0 = dfdupar(k_uper, NUPAR)
                
                if ( i .ne. NUPAR ) then
                   dfdupar0 = dfdupar(k_uper, i) + &
                                (dfdupar(k_uper, i+1) - dfdupar(k_uper, i)) * p
                endif      
    
                U0 = (1. - DFACTPAR) * DFDUPER0 + DFACTPER * DFDUPAR0   

            endif maxwellianOrNot
            
            factor = PISQMUT0 / abs(NPARA_sav(n, m)) 
            
            sumb_11_nm = UPER(k_uper) * UPER(k_uper)  * Jn(IHARM + 1, n, m) * epsx     &
                    + UPER(k_uper) * UPER(k_uper)  * Jn(IHARM - 1, n, m) * epsy     &
                    + SQ2 * UPER(k_uper) * UPAR0  * Jn(IHARM, n, m)     * epsz

            sumb_31_nm = SQ2 * UPER(k_uper) * UPAR0  * Jn(IHARM + 1, n, m) * epsx     &
                    + SQ2 * UPER(k_uper) * UPAR0  * Jn(IHARM - 1, n, m) * epsy     &
                    + 2.0 * UPAR0 * UPAR0 * Jn(IHARM, n, m)     * epsz
     
            sumb_11_nm = sumb_11_nm * factor 
            sumb_31_nm = sumb_31_nm * factor 
                  
            sume_11_nm = sumb_11_nm * facte
            sume_31_nm = sumb_31_nm * facte
            
            sumc_11_nm = sume_11_nm * sinth_inv
            sumc_31_nm = sume_31_nm * sinth_inv     
     
            sumf_11_nm = sumc_11_nm * facte
            sumf_31_nm = sumc_31_nm * facte
                  
            sumwdot_11_nm = sumb_11_nm * u0
            sumwdot_31_nm = sumb_31_nm * u0   
                                                   
            sumf_11(i_uprl) = sumf_11(i_uprl) + sumf_11_nm
            sumf_31(i_uprl) = sumf_31(i_uprl) + sumf_31_nm

            sume_11(i_uprl) = sume_11(i_uprl) + sume_11_nm
            sume_31(i_uprl) = sume_31(i_uprl) + sume_31_nm

            sumc_11(i_uprl) = sumc_11(i_uprl) + sumc_11_nm
            sumc_31(i_uprl) = sumc_31(i_uprl) + sumc_31_nm

            sumb_11(i_uprl) = sumb_11(i_uprl) + sumb_11_nm
            sumb_31(i_uprl) = sumb_31(i_uprl) + sumb_31_nm
        
            sumwdot_11 = sumwdot_11 + sumwdot_11_nm
            sumwdot_31 = sumwdot_31 + sumwdot_31_nm
              
            sumwdotkx_11 = sumwdotkx_11 + xkxsav(n) * sumwdot_11_nm
            sumwdotkx_31 = sumwdotkx_31 + xkxsav(n) * sumwdot_31_nm   
              
            sumwdotky_11 = sumwdotky_11 + xkysav(m) * sumwdot_11_nm
            sumwdotky_31 = sumwdotky_31 + xkysav(m) * sumwdot_31_nm                           
              
        enddo resonantModes
       
        sum_wdot = sum_wdot + sum2_1 * sumwdot_11 &
                   + sum2_2 * sumwdot_11 &
                   + sum2_3 * sumwdot_31 
     
        sum_fx0 = sum_fx0  + sumkx2_1 * sumwdot_11 &
                  + sumkx2_2 * sumwdot_11 &
                  + sumkx2_3 * sumwdot_31 &
                  + sum2_1 * sumwdotkx_11 &
                  + sum2_2 * sumwdotkx_11 &
                  + sum2_3 * sumwdotkx_31
     
        sum_fy0 = sum_fy0  + sumky2_1 * sumwdot_11 &
                 + sumky2_2 * sumwdot_11 &
                 + sumky2_3 * sumwdot_31 &
                 + sum2_1 * sumwdotky_11 &
                 + sum2_2 * sumwdotky_11 &
                 + sum2_3 * sumwdotky_31    
      
        !   QL diffusion coefficients, need to sum over all uPrl

        uPrlLoop: & 
        do i_uprl = 1, nupar

            b_sum(i_uprl) = b_sum(i_uprl) + sum2_1 * sumb_11(i_uprl) &
                                          + sum2_2 * sumb_11(i_uprl) &
                                          + sum2_3 * sumb_31(i_uprl)
            c_sum(i_uprl) = c_sum(i_uprl) + sum2_1 * sumc_11(i_uprl) &
                                          + sum2_2 * sumc_11(i_uprl) &
                                          + sum2_3 * sumc_31(i_uprl)
            e_sum(i_uprl) = e_sum(i_uprl) + sum2_1 * sume_11(i_uprl) &
                                          + sum2_2 * sume_11(i_uprl) &
                                          + sum2_3 * sume_31(i_uprl)
            f_sum(i_uprl) = f_sum(i_uprl) + sum2_1 * sumf_11(i_uprl) &
                                          + sum2_2 * sumf_11(i_uprl) &
                                          + sum2_3 * sumf_31(i_uprl)     
        enddo uPrlLoop
       
    enddo iHarmLoop
    
    deallocate( zbeta, zbeta_iharm, nres, mres, &
        zetai, Jni, Jn, NPARA_sav, sumb_11, sumb_31, &
        sumc_11, sumc_31, sume_11, sume_31, sumf_11, &
        sumf_31, is_resonance_nm )
    
  end subroutine qlSum

end module qlSum_general
