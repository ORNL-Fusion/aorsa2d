
	freq = 30.0d6
	nphi = 12

	eqdsk = 1
	eqdskFileName = 'g130608.00355.EFIT02.mds.corrected.qscale_1.00000.dlgMod'
	bTorFactor = 1.0

	flux_profiles = 1
	flat_profiles = 0

	atomicZ	= [-1,1]
	amu = [me/mi,2]

	; Flux function is : a[0] + ( a[1] - a[0] ) * ( 1.0 - x^a[3] )^a[2]
	;	a[0] = limiter, a[1] = center, a[2] = alpha, a[3] = beta
   	;		limiter,	center,		alp,	bet

	nn = [	[0.0,		0.0,		0.2,	2.0],$ ; electrons are spec 0
			[0.3d19,	3.0d19,		0.2,	2.0] ]

	tt = [	[00.1d3,	1.0d3,		0.5,	2.0],$
			[00.2d3,	8.0d3,		2.0,	2.0] ]

	DensityMin = 0.2e18
	TempMin = 1e-3

	; Grid
	nR = 600
	nZ = 600
	rMin = 0.05
	rMax = 1.8
	zMin = -1.8
	zMax = +1.8

	
