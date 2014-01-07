
	ShrinkFac = 0.25

	freq = 53.0e6
	nphi = -27
	b0 = 4.3
	r0 = 6.0 * ShrinkFac

	x0 = r0
	y0 = 0.0 * ShrinkFac

	bField_eqdsk = 0
    bField_gaussian = 1
    bField_flat = 0

	eqdskFileName = 'Scen4_bn2.57_129x129.dlgMod'

    br_flat = 1d-5
	bt_flat = 0.0
	bz_flat = 0.0

	bpmax_b0 = 0.4

	flux_profiles = 0
	flat_profiles = 1
	gaussian_profiles = 0

	atomicZ	= [-1,2]
	amu = [me/mi,4]

	; Flux function is : a[0] + ( a[1] - a[0] ) * ( 1.0 - x^a[3] )^a[2]
	;	a[0] = limiter, a[1] = center, a[2] = alpha, a[3] = beta
   	;		limiter,	center,		alp,	bet

	nn = [	[0.0,		0.0,		9.0,	10.0],$ ; electrons are spec 0
			[1.0d19,	2.0d19,		9.0,	10.0] ]

	tt = [	[00.5d-1,	00.5d3,		6.0,	4.0],$
			[00.5d-1,	00.5d3,		6.0,	4.0] ]

	psi_xsig = 2 * ShrinkFac
	psi_ysig = 4 * ShrinkFac

	Density_xsig=1.1 * ShrinkFac
	Density_ysig=2.2 * ShrinkFac

	Temp_xsig=1.0 * ShrinkFac
	Temp_ysig=2.0 * ShrinkFac

	DensityMin = 1.0d16
	TempMin = 1e-3

	; Grid

	nR = 512
	nZ = 512
	rMin = 4.0*ShrinkFac
	rMax = 8.4*ShrinkFac
	rMax = 12*ShrinkFac

	zMin = -4.7*ShrinkFac
	zMax = 4.7*ShrinkFac

    rDomainBox = [rMin,rMax,rMax,rMin,rMin]
    zDomainBox = [zMin,zMin,zMax,zMax,zMin]

    @ar2_coupling_common

    rLim = LeftSide_rLim
    zLim = LeftSide_zLim
