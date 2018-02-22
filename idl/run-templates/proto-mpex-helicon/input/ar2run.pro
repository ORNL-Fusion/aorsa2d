
	freq = 13.56d6
    dispersionAxialWaveLength = 0.2
    r0 = 100
    radius = 0.06

	nphi = fix(2*!pi*r0/dispersionAxialWaveLength)

    bField_flat = 1
	bField_eqdsk = 0

    br_flat = 0.0
	bt_flat = 0.03
	bz_flat = 0.0

	gaussian_profiles = 0
	flux_profiles = 0
   	numeric_flux_profiles = 0
	fred_namelist_input = 0
    numeric_profiles = 1

    fileName = 'input/exp.data.1'
    nLines = file_lines(fileName)
    data = replicate({x:0.0,y:0.0}, nLines)
    openr, lun, fileName, /get_lun
    readf, lun, data
    free_lun, lun 
    data.x = data.x*1e-2*2.4; convert to m and add expand
    dx = data[1].x-data[0].x
    iiUse = where(data.x ge 0,nUse) ; choose one or the other radial profile
    numeric_range = 0.1
    numeric_decay_length = 0.03
    numeric_n = fix(numeric_range / dx)
    numeric_r = fIndGen(numeric_n)*dx
    data_r = data[iiUse].x
    data_ne = data[iiUse].y
    iiSort = sort(abs(data_r))
    data_r = data_r[iiSort]
    data_ne = data_ne[iiSort]
    numeric_max_r = data_r[-1]
    numeric_max_r_value = data_ne[-1] 
    numeric_norm = exp(-(numeric_max_r)/numeric_decay_length)
    numeric_ne = exp(-(numeric_r)/numeric_decay_length)/numeric_norm*numeric_max_r_value
    numeric_ne[0:nUse-1] = data_ne
    numeric_r = numeric_r 

    x0 = r0
    y0 = 0
    rCenter = r0
    zCenter = y0

    density_xsig = 0.03
    density_ysig = 0.03

    temp_xsig = 0.01
    temp_ysig = 0.01

	atomicZ	= [-1,1]
	amu = [me/mi,2]

	nn = [	[0.0,		0.0,		9.0,	10.0],$ ; electrons are spec 0
			[0.0,	    5.0d18,		9.0,	10.0] ]

	tt = [	[00.5d-1,	00.5d3,		6.0,	4.0],$
			[00.5d-1,	00.5d3,		6.0,	4.0] ]

	DensityMin = [1e12, 1e12]
	TempMin = 1e-3

	; Grid
	nR = 91
	nZ = 91
	rMin = r0-radius
	rMax = r0+radius
	zMin = -radius
	zMax = +radius

    rDomainBox = [rMin,rMax,rMax,rMin,rMin]
    zDomainBox = [zMin,zMin,zMax,zMax,zMin]

    rLim = rDomainBox
    zLim = zDomainBox

    angle = fIndGen(100)/99*2*!pi
    rlcfs = (0.9*radius*cos(angle))+r0
    zlcfs = (0.9*radius*sin(angle)) 


