
	freq = 56.0d6
	nphi = -27
	eqdsk = 1
	eqdskFileName = 'Scen4_bn2.57_129x129.dlgMod'
	flux_profiles = 1
	atomicZ	= [-1,1,1]
	amu = [me/mi,2,3]

	; Flux function is : a[0] + ( a[1] - a[0] ) * ( 1.0 - x^a[3] )^a[2]
	;	a[0] = limiter, a[1] = center, a[2] = alpha, a[3] = beta
   	;		limiter,	center,		alp,	bet

	nn = [	[0.0,		0.0,		9.0,	10.0],$ ; electrons are spec 0
			[0.5d17,	2.0d19,		9.0,	10.0],$
			[0.5d17,	2.0d19,		9.0,	10.0] ]

	tt = [	[00.2d3,	20.0d3,		6.0,	4.0],$
			[00.3d3,	20.0d3,		6.0,	4.0],$
			[00.3d3,	20.0d3,		6.0,	4.0] ]

	; Grid
	nR = 600
	nZ = 600
	rMin = 3.5
	rMax = 9.0
	zMin = -5.0
	zMax = +5.0


