
	ShrinkFac = 0.25

	freq = 53.0e6
	nphi = -27
	b0 = 4.3
	r0 = 6.0 * ShrinkFac

	x0 = r0
	y0 = 0.0 * ShrinkFac

	eqdsk = 0
	eqdskFileName = 'Scen4_bn2.57_129x129.dlgMod'
	flux_profiles = 0
	flat_profiles = 0
	gaussian_profiles = 1
	br_flat = 1d-5
	bt_flat = 0.0
	bz_flat = 0.0

	bpmax_b0 = 0.4

	atomicZ	= [-1,2]
	amu = [me/mi,4]

	; Flux function is : a[0] + ( a[1] - a[0] ) * ( 1.0 - x^a[3] )^a[2]
	;	a[0] = limiter, a[1] = center, a[2] = alpha, a[3] = beta
   	;		limiter,	center,		alp,	bet

	nn = [	[0.0,		0.0,		9.0,	10.0],$ ; electrons are spec 0
			[1.0d17,	2.0d19,		9.0,	10.0] ]

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

	nR = 64
	nZ = 128
	rMin = 4.0*ShrinkFac
	rMax = 8.4*ShrinkFac
	rMax = 12*ShrinkFac

	zMin = -4.7*ShrinkFac
	zMax = 4.7*ShrinkFac

	rlim = [ $
			4.066, $
			4.325, $
			5.000, $
			5.779, $
			6.951, $
			7.578, $
			7.993, $
			8.3, $
			8.3, $
			8.3, $
			7.924, $
			7.318, $
			6.315, $
			5.796, $
			4.516, $
			4.083, $
			4.083, $
			4.066 ]*ShrinkFac

	zlim = ([ $
			3.581, $
			4.256, $
			4.671, $
			4.498, $
			3.599, $	
			2.976, $
			2.215, $
			1.626, $
			0.588, $
			-0.519, $
			-1.401, $
			-2.318, $
			-3.253, $
			-3.374, $
			-3.235, $
			-2.561, $
			0.000, $
			3.581 ]-0.5)*ShrinkFac


	VorpalBox_rOffSet = 2.0
	VorpalBox_r = [0,0.8,0.8,0.0,0.0]+VorpalBox_rOffset
	VorpalBox_zOffSet = -0.36/2.0
	VorpalBox_z = [0,0,0.36,0.36,0.0]+VorpalBox_zOffset
