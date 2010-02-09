module orbit_mod

contains

      subroutine f(x, y, dy)
      use parameters
      use fitpack

      implicit none

      integer fcount, nnodex, nnodey

      real x, phi, y(:), dy(:), br, bz, bphi, modb, x_extint, sgn_vprl
      real bratio_phi
      real xprime, yprime, dxdphi, dydphi, dr, dz, rt, capr, xwleft
      real bxn(nmodesmax, mmodesmax), byn(nmodesmax, mmodesmax), &
         bzn(nmodesmax, mmodesmax), bmod(nmodesmax, mmodesmax), &
         bratio(nmodesmax, mmodesmax)

      real sigma, surf2
      real zbxn (nmodesmax, mmodesmax, 3)
      real zbyn (nmodesmax, mmodesmax, 3)
      real zbzn (nmodesmax, mmodesmax, 3)
      real zbmod(nmodesmax, mmodesmax, 3)
      real zbratio(nmodesmax, mmodesmax, 3)
      real xprimea(nmodesmax), yprimea(mmodesmax)

      common/fcom/fcount, bxn, byn, bzn, bmod, bratio, &
         dr, dz, &
         nnodex, nnodey, rt, xwleft, sgn_vprl, modb, bratio_phi, &
         dxdphi, dydphi, capr

      common/spline_com/sigma, zbxn, zbyn, zbzn, zbmod, zbratio, &
         xprimea, yprimea

      phi = x
      xprime = y(1)
      yprime = y(2)

      br = surf2 (xprime, yprime, nnodex, nnodey, xprimea, yprimea, &
         bxn, nxmx, zbxn, sigma)

      bz = surf2 (xprime, yprime, nnodex, nnodey, xprimea, yprimea, &
         byn, nxmx, zbyn, sigma)

      bphi =surf2(xprime, yprime, nnodex, nnodey, xprimea, yprimea, &
         bzn, nxmx, zbzn, sigma)

      modb =surf2(xprime, yprime, nnodex, nnodey, xprimea, yprimea, &
         bmod, nxmx, zbmod, sigma)

      bratio_phi=surf2(xprime, yprime, nnodex, nnodey, xprimea, yprimea, &
         bratio, nxmx, zbratio, sigma)

      x_extint = xprime + xwleft
      capr = x_extint + rt

      dxdphi = sgn_vprl * capr * br / bphi
      dydphi = sgn_vprl * capr * bz / bphi


      dy(1) = dxdphi
      dy(2) = dydphi

      fcount = fcount + 1

 1312 format(1p8e12.4)

      return
      end subroutine f

!
!***************************************************************************
!


      subroutine fdtau(dtau, nxdim, nydim, len_x, modb_x, &
                     n_theta_max, n_psi_max, norb_dim, sinth2_init, &
                     modb_init, n_theta_, &
                     i_psi, i, j, nphi_enter, nphi_exit, sgn_vprl)

      implicit none

      integer n_theta_max, n_psi_max, norb_dim, n_theta, i_psi
      integer nxdim, nydim

      real dtau(n_theta_max)
      real len_x(norb_dim), modb_x(norb_dim)

      real sinth2_init(n_theta_max, n_psi_max)

      integer n_theta_(n_psi_max), i, j, nphi_enter, nphi_exit

      real mri, b_ip1, b_i, dl, l_ip1, l_i, sgn_vprl, argi, modb_init
      real costh_i, tanth2_i, sinth2_i, costh2_i

!     -------------------------------------------------
!     This subroutine calculates the time (dtau)
!     that a particle stays in the (i, j) spatial cell
!     -------------------------------------------------

      l_i =   len_x(nphi_enter)
      l_ip1 = len_x(nphi_exit)

      b_i =   modb_x(nphi_enter)
      b_ip1 = modb_x(nphi_exit)

      mri = b_ip1 / b_i
      dl = l_ip1 - l_i
      dl = dl * sgn_vprl

      do n_theta = 1, n_theta_(i_psi)

!         thetai = theta_(n_theta, i_psi)
!         costhi = cos(thetai)
!         sinthi = sin(thetai)
!         tanthi = tan(thetai)

!         write(6, 100)modb_init
!  100    format(1p8e12.4)
!         call exit

         dtau(n_theta) = 0.0

         sinth2_i = b_i / modb_init * sinth2_init(n_theta, i_psi)
       costh2_i = 1.0 - sinth2_i
       tanth2_i = sinth2_i / costh2_i
       costh_i = sqrt(costh2_i)

       if(sinth2_i .gt. 1.0) go to 300

         argi = 1.0 - (mri - 1.0) * tanth2_i


!        -------
!        Passing
!        -------
         if(argi .ge. 0.0) then          ! passing

            dtau(n_theta) = dl / sgn_vprl &
                  * 2.0 / abs(costh_i) / (1.0 + sqrt(argi))
         end if



!        -------
!        Trapped
!        -------
         if(argi .lt. 0.0) then          ! trapped

            dtau(n_theta) = dl / sgn_vprl &
                 * 2.0 * abs(costh_i)  / (mri - 1.0) / sinth2_i
         end if

  300    continue

      end do




      return
      end subroutine fDtau


!
!***************************************************************************
!



      SUBROUTINE EXTINT (NMAX, X, Y, F, H0, MMAX, ERROR)
!     IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      INTEGER NMAX, MMAX
!     REAL*8 X, H0, Y(NMAX)
      REAL   X, H0, Y(NMAX)
!
!          A DEFERRED-LIMIT INTEGRATOR (J.P. BORIS AND N.K. WINSOR)
!
!         THIS SUBROUTINE INTEGRATES UP TO 100 SIMULTANEOUS FIRST ORDER
!     ORDINARY DIFFERENTIAL EQUATIONS FROM X TO X + H0 BY REPEATED EX-
!     TRAPOLATIONS ON A MIDPOINT RULE. UP TO 10 EXTRAPOLATIONS MAY BE
!     REQUESTED BEFORE REDUCTION OF THE INITIAL STEPSIZE IS CARRIED OUT.
!     (REF. R. BULIRSCH AND J. STOER -  NUMERISCHE MATHEMATIK 8,1 (1966)
!
!     NMAX         THE TOTAL NUMBER OF DEPENDENT VARIABLES BEING INTE-
!                  GRATED. THERE WILL BE ONE FIRST ORDER EQUATION FOR
!                  EACH DEPENDENT VARIABLE.
!
!     X            THE INDEPENDENT VARIABLE, X IS TREATED AS REAL AND
!                  MONOTONIC DURING THE INTEGRATION STEP. THE VALUE OF
!                  X ON RETURN FROM ''EXT INT'' CONTAINS THE VALUE WHICH
!                  IS APPROPRIATE TO THE LENGTH OF THE INTEGRATION ACTU-
!                  ALLY PERFORMED. IF CONVERGENCE HAS BEEN OBSERVED IN
!                  M MAX OF FEWER EXTRAPOLATIONS, X (AT EXIT) = X (AT
!                  ENTRY) + H0. ON ENTRY X MUST BE SET TO THE INITIAL
!                  VALUE OF THE INDEPENDENT VARIABLE FOR THE INTEGRATION
!                  STEP BEING CONTEMPLATED.
!
!     Y            THE DEPENDENT VARIABLES. EACH DEPENDENT VARIABLE Y(N)
!                  (FOR N = 1, 2, .., NMAX) IS INTEGRATED FROM X TO X+H0
!                  IF ADEQUATE CONVERGENCE, AS DEFINED BY THE ''ERROR''
!                  SUBROUTINE, OCCURS IN M MAX OR FEWER EXTRAPOLATIONS.
!
!     F            THE DERIVATIVE SUBROUTINE SUPPLIED BY THE USER FOR
!                  HIS PARTICULAR PROBLEM. IT MUST BE OF THE FORM
!                  F (X, Y, DY) WHERE THE ARRAY DY(N) (FOR N = 1,2, ...,
!                  NMAX) IS RETURNED CONTAINING THE NMAX DERIVATIVES
!                  (DY(N)/DX)(X,Y).
!
!     H0           IS THE BASIC STEPSIZE OF THE INTEGRATION. IF CONVER-
!                  GENCE OCCURS WITHIN MMAX EXTRAPOLATIONS, X RETURNS
!                  FROM EXT INT WITH THE VALUE X + H0. THE VALUES IN THE
!                  ARRAY Y ARE THE VALUES OF THE DEPENDENT VARIABLES AT
!                  THIS VALUE OF X. IF CONVERGENCE DOES NOT OCCUR, H0 IS
!                  HALVED AND THE ENTIRE EXTRAPOLATION PROCEDURE IS
!                  REPEATED AND REPEATED AGAIN UNTIL CONVERGENCE OCCURS.
!                  AN ATTEMPT HAS BEEN MADE TO UTILIZE AS MUCH PREVIOUS-
!                  LY COMPUTED INFORMATION AS POSSIBLE WHEN H0 MUST BE
!                  HALVED FOR CONVERGENCE.
!
!     M MAX        CONTAINS THE NUMBER OF TIMES EXTRAPOLATION IS ATTEM-
!                  PTED BEFORE H0 IS HALVED. THIS VALUE WILL VARY WITH
!                  COMPUTER ROUND-OFF ERROR AND WITH THE TYPE AND NUMBER
!                  OF EQUATIONS BEING INTEGRATED. IN ALL CASES, HOWEVER,
!                  ONE SHOULD SPECIFY AN M MAX 'GQ' 2. THE STEPSIZE FOR
!                  FUTURE ITERATIONS IS SELECTED AND RETURNED IN H0.
!                  THIS NEW VALUE IS CHOSEN SO THAT CONVERGENCE WILL BE
!                  OBSERVED AT ABOUT THE 0.66*MMAX-TH EXTRAPOLATION.
!                  AS WRITTEN, MMAX.LE.10.
!
!     ERROR        IS A SUBROUTINE WHICH IS USED TO DETERMINE THE SATIS-
!                  FACTORY CONVERGENCE OF EACH INDIVIDUAL EQUATION BEING
!                  INTEGRATED. AN EXAMPLE IS GIVEN WHICH CORRESPONDS TO
!                  THE ERROR CRITERION USED BY BULIRSCH AND STOER.
!                      THE CALLING SEQUENCE IS
!                  ERROR (M, DY, CONV, FINISH) WHERE M IS THE ORDER OF
!                  EXTRAPOLATION, DY IS THE VECTOR OF INCREMENTS TO THE
!                  DEPENDENT VARIABLES Y, CONV IS A VECTOR OF LOGICALS
!                  WHICH ARE .TRUE. COMPONENTWISE WHEN THE CORRESPONDING
!                  DEPENDENT VARIABLE HAS CONVERGED, AND FINISH IS
!                  RETURNED .TRUE. WHEN THE ENTIRE SYSTEM OF EQUATIONS
!                  HAS SATISFIED THE CONVERGENCE CRITERION.
!
      LOGICAL LATERL, CONV(100), PREVIN, FINISH
!     REAL*8 STEPFC, X0, U, SUM, YM, BETA, H, DEN, SQRT2, YP
      REAL   STEPFC, X0, U, SUM, YM, BETA, H, DEN, SQRT2, YP
      INTEGER J, K, L, M, N, LMAX, KASIDE, PTS, MM, MMAXP
      INTEGER KMIN
!     REAL*8 HM(11), S(11), P(11), YBAR(100,11), Y0(100)
      REAL   HM(11), S(11), P(11), YBAR(100,11), Y0(100)
!     REAL*8 Y NEW(100), Y OLD(100), DY(100), DY0(100), Y HOLD(7,10,100)
      REAL   YNEW(100), YOLD(100), DY(100), DY0(100), YHOLD(7,10,100)
!
!
! INITIALIZE..100
      MMAXP = MMAX + 1
      SQRT2 =  SQRT(2.0)
      FINISH = .FALSE.
      LATERL = .FALSE.
      X0 = X
      DO 100 N = 1,NMAX
  100   Y0(N) = Y(N)
      LMAX = (MMAX + 1)/2 + 1
      CALL F (X0, Y, DY0)
!
! START A NEW LEVEL..200
  204 X = X0 + H0
      KMIN = 1
      STEPFC = 2.0**(MMAX/3.0+0.5)
      DO 205 N = 1, NMAX
  205   CONV(N) = .FALSE.
      IF (.NOT.LATERL .OR. (MMAX.LT.1)) GO TO 203
!     ELSE BEGIN SHIFTING OLD INFORMATION INTO POSITION..
      DO 202 N = 1, NMAX
        Y(N) = YHOLD(2,1,N)
        YBAR(N,1) = Y(N)
        DO 202 L = 1, LMAX
          MM2 = MMAX - 2
          MMIN = IABS(2*L-3)
          DO 202 M = MMIN, MM2
  202       YHOLD(L,M,N) = YHOLD (L+1, M+2, N)
!
! COMPUTING BETA AND ASSOCIATED QUANTITIES..300
  203 H = H0/2
      HM(1) = H0/2.0
      HM(2) = H0/4.0
      HM(3) = H0/6.0
      BETA = 0.25/(H0*H0)
!
! EXTRAPOLATIONS OF HIGHER ORDER EXPANSIONS..400
      DO 400 MM = 1, MMAXP
!     BEGINNING THE LOOP OVER MMAX EXTRAPOLATIONS IN THIS LEVEL..
        M = MM - 1
        STEPFC = STEPFC/SQRT2
        PREVIN = .FALSE.
        IF (LATERL.AND.(M.LT.MMAX - 1)) PREVIN = .TRUE.
        KASIDE = 2
        IF (2*(M/2).EQ.M) KASIDE = 3
        L = (M + 1)/2 + 1
        IF (M.GT.2) HM(MM) = HM(MM-2)/2
        H = HM(MM)
!       S(MM) = 1 - DEXP(-BETA*H*H)
        S(MM) = 1 -  EXP(-BETA*H*H)
        IF (PREVIN) GO TO 503
!       ELSE GENERATE THE M-TH MIDPOINT INTEGRAL..
          DO 404 N = 1, NMAX
            YOLD(N) = Y0(N)
  404       YNEW(N) = Y0(N) + H*DY0(N)
          CALL F (X0+H, YNEW, DY)
          PTS = (H0*1.000001)/H
          DO 405 K = 2, PTS
!         BEGIN THE MIDPOINT INTEGRATION..
            DO 406 N = 1, NMAX
              U = YOLD(N) + 2*H*DY(N)
              YOLD(N) = YNEW (N)
  406         YNEW (N) = U
            CALL F (X0 + K*H, YNEW, DY)
            IF ((K.NE.KASIDE).OR.(L.LT.2)) GO TO 405
!           ELSE BEGIN PUTTING INFORMATION ASIDE..
              DO 408 N = 1, NMAX
  408           YHOLD (L, M, N) = (YNEW(N) + YOLD(N) + H*DY(N))/2.0
              L = L - 1
              KASIDE = 2*KASIDE
!           END OF PUTTING INFORMATION ASIDE
  405     CONTINUE
!         END OF THE MIDPOINT INTEGRATION
!
! NOW ADVANCE THE DEPENDENT VARIABLES..500
  503   IF (M.GT.0) GO TO 504
        DO 505 N = 1, NMAX
          YBAR(N,1) = (YNEW(N) + YOLD(N) + H*DY(N))/2
  505     Y(N) = YBAR(N,1)
        IF (MMAX.EQ.0) GO TO 700
        GO TO 400
!       NOW DETERMINE THE INTERPOLATIONAL POLYNOMIALS..
  504   IF (STEPFC.LT.1.1) KMIN = KMIN + 1
        DEN = 1
        DO 401 K = KMIN, M
          P(K) = ((H/HM(K))**2)
          DO 410 J = KMIN, M
            IF (J.NE.K) P(K) = P(K)*(S(J) - S(MM))/(S(J) - S(K))
  410     CONTINUE
  401     DEN = DEN - P(K)
!       END DETERMINATION OF THE INTERPOLATIONAL POLYNOMIALS
        DO 500 N=1, NMAX
          IF (CONV(N)) GO TO 500
          YP = Y(N)
          YM = YHOLD (1,M,N)
          IF (.NOT.PREVIN) YM = (YNEW(N) + YOLD(N) + H*DY(N))/2.0
          SUM = 0.0
          IF (M.LT.2) GO TO 501
          DO 502 J = KMIN, M
  502       SUM = SUM + (YBAR(N,J) - YP)*P(J)
  501     DY(N) = 0.0
          IF (DEN.NE.0.0) DY(N) = ((YM - YP) - SUM)/DEN
          Y(N) = YP + DY(N)
          YBAR(N,MM) = YM
  500   CONTINUE
        CALL ERROR (M, DY, CONV, FINISH)
        IF (FINISH) GO TO 700
  400 CONTINUE
!
! PREPARE FOR THE NEXT LEVEL IF NECESSARY..600
      LATERL = .TRUE.
      H0 = H0/2
      GO TO 204
!
! RETURN..700
  700 H0 = H0 * STEPFC
      RETURN
      END subroutine extInt

!
!***************************************************************************
!

      SUBROUTINE ERROR (M, DY, CONV, FINISH)
!
!         THIS VERSION OF THE ERROR ROUTINE CORRESPONDS CLOSELY TO THE
!     VERSION EMPLOYED BY BULIRSCHE AND STOER. THE DIFFERENCES, HOWEVER,
!     ARE NOTEWORTHY. THIS VERSION IS CALLED ONLY AFTER ALL UNCONVERGED
!     DEPENDENT VARIABLES HAVE BEEN EXTRAPOLATED AT ORDER M. THUS AN
!     EXCESSIVE NUMBER OF SUBROUTINE REFERENCES ARE MADE UNNECESSARY.
!         NEW VALUES FOR ANY PARTICULAR DEPENDENT VARIABLE ARE NOT CAL-
!     CULATED AFTER CONVERGENCE, AS DEFINED BY THIS SUBROUTINE, HAS BEEN
!     OBSERVED FOR THAT VARIABLE. THUS POSSIBLE INSTABILITIES ASSOCIATED
!     WITH ATTEMPTS AT OVERCONVERGENCE CAN BE ELIMINATED. ONE OF THESE
!     INSTABILITIES, ARISING FROM RESETTING THE MAGNITUDE VECTOR S AT
!     EVERY EXTRAPOLATION AS DOES THE B - S PROGRAM, IS ELIMINATED BY
!     RESETTING S ONLY AFTER THE CORRESPONDING DEPENDENT VARIABLE HAS
!     CONVERGED ACCORDING TO THE CRITERIA CHOSEN.
!     IMPLICIT REAL*8 (A-H,O-Z), INTEGER (I-N)
      INTEGER M
!     REAL*8 DY(M)
      REAL   DY(M)
      LOGICAL CONV(M), FINISH
      COMMON /ERRCOM/ EPS, S(100), Y(100), NMAX
!     REAL*8 EPS, S, Y
      REAL   EPS, S, Y
      INTEGER NMAX, NTIMES (100)
!
      IF (M.NE.1) GO TO 1
        DO 3 N = 1, NMAX
    3     NTIMES(N) = 0
        NCONV = 0
    1 DO 2 N = 1, NMAX
!       IF (.NOT.(DABS(DY(N))/S(N).LT.EPS).OR. CONV(N))  GO TO 2
        IF (.NOT.( ABS(DY(N))/S(N).LT.EPS).OR. CONV(N))  GO TO 2
        NTIMES(N) = NTIMES(N) + 1
        IF (NTIMES(N).EQ. 1) NCONV = NCONV + 1
        IF (NTIMES(N).EQ. 2) CONV(N) = .TRUE.
!       IF (DABS(Y(N)).GT. S(N)) S(N) = DABS(Y(N))
        IF ( ABS(Y(N)).GT. S(N)) S(N) =  ABS(Y(N))
    2 CONTINUE
      IF (NCONV.EQ.NMAX) FINISH = .TRUE.
      RETURN
      END subroutine error

!
!***************************************************************************
!

end module orbit_mod
