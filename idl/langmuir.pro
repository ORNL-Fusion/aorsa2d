
	freq = 1.15d8
	nphi = 0 
	eqdsk = 0
	b0 = 0.0001
	br_frac = 100.0
	bz_frac = 0.0
	r0 = 100.0
	flux_profiles = 0
	atomicZ	= [-1,1]
	amu = [me/mi,1]

	; Flux function is : a[0] + ( a[1] - a[0] ) * ( 1.0 - x^a[3] )^a[2]
	;	a[0] = limiter, a[1] = center, a[2] = alpha, a[3] = beta
   	;		limiter,	center,		alp,	bet

	nn = [	[0.0,		0.0,		9.0,	10.0],$ ; electrons are spec 0
			[1.1d14,	1.1d14,		9.0,	10.0] ]

	tt = [	[00.5d3,	00.5d-1,		6.0,	4.0],$
			[00.5d3,	00.5d-1,		6.0,	4.0] ]

	; Grid

	nR = 600
	nZ = 10
	rMin = 95.0 
	rMax = 105.0
	zMin = -0.1
	zMax = +0.1


