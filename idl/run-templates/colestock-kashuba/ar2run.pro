
	freq = 42e6
	nphi = 0 

	bfield_eqdsk = 0
    b0 = 2.9
    r0 = 1.32
    br_frac = 0
    bz_frac = 0

	flux_profiles = 0
    gaussian_profiles = 0
    parabolic_profiles = 1
    parabolic_half_length = 0.35

    ; [e,D,H]
	atomicZ	= [-1,1,1]
	amu = [me/mi,2,1]

	; Flux function is : a[0] + ( a[1] - a[0] ) * ( 1.0 - x^a[3] )^a[2]
	;	a[0] = limiter, a[1] = center, a[2] = alpha, a[3] = beta
   	;		limiter,	center,		alp,	bet

    ne0 = 4d19
    nD = 0.95 * ne0
    nH = 0.05 * ne0 

	nn = [	[0.0,	0.0,	0.0,	0.0],$ ; electrons are spec 0
			[0.0,	nD,	   	0.0,	0.0],$
            [0.0,	nH,		0.0,	0.0]]

	tt = [	[0,	2d3,		0.0,	0.0],$
			[0,	2d3,		0.0,	0.0],$
            [0,	2d3,		0.0,	0.0] ]

	DensityMin = 5d17
	TempMin = 2d3 

	; Grid

	nR = 201
	nZ = 201
    rMinor = 0.4
	rMin = r0-rMinor 
	rMax = r0+rMinor
	zMin = -0.3
	zMax = +0.3

	
