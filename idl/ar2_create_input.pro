;
; Create ar2 input file such that the following requirements are met:
; 
; 1. The problem is periodic in all variables.
; 2. An absorbing SOL is created.
; 3. The transition to metal is made soft.
;



function MakePeriodic, Arr, MaskIn, look = look

	StartSmoothWidth = n_elements(Arr[*,0])/5.0;20.0
	FinalSmoothWidth = n_elements(Arr[*,0])/50.0;2.0
	nSmooth = 50.0

	; Create new expanded and smoothed mask

	nR = n_elements(Arr[*,0])
	nZ = n_elements(Arr[0,*])

	zSlice = nZ/2
	rSlice = nR/2
	if keyword_set(look) then pr=plot(Arr[*,zSlice],thick=3)
	if keyword_set(look) then pz=plot(Arr[rSlice,*],thick=3)
	ArrTmp = Arr*MaskIn
	ArrTmp = smooth ( ArrTmp, StartSmoothWidth )
	iiMaskInSmooth_lim = where (abs(ArrTmp) gt 0)

	if keyword_set(look) then p=plot(ArrTmp[*,zSlice],over=pr,thick=2,color='blue')
	if keyword_set(look) then p=plot(ArrTmp[rSlice,*],over=pz,thick=2,color='blue')

	for i=0,nSmooth-1 do begin

		smoothWidth = StartSmoothWidth + (FinalSmoothWidth-StartSmoothWidth)/nSmooth*i

		ArrTmp[iiMaskInSmooth_lim] = Arr[iiMaskInSmooth_lim]
		tmp = [[ArrTmp,ArrTmp,ArrTmp],[ArrTmp,ArrTmp,ArrTmp],[ArrTmp,ArrTmp,ArrTmp]]
		tmp = smooth (tmp, smoothWidth )
		ArrTmp = tmp[nR:nR+nR-1,nZ:nZ+nZ-1]
		if keyword_set(look) then p=plot(ArrTmp[*,zSlice],over=pr,thick=1,color='red')
		if keyword_set(look) then p=plot(ArrTmp[rSlice,*],over=pz,thick=1,color='red')

	endfor

	return, ArrTmp

end

pro ar2_create_input

	@constants

	freq = 56.0d6
	nphi = -27

	eqdskFileName = 'Scen4_bn2.57_129x129.dlgMod'

	z	= [-1,1,1]
	amu = [me/mi,2,3]

	; Flux function is : a[0] + ( a[1] - a[0] ) * ( 1.0 - x^a[3] )^a[2]
	;	a[0] = limiter, a[1] = center, a[2] = alpha, a[3] = beta
   	;		limiter,	center,		alp,	bet

	nn = [	[0.0,		0.0,		1.0,	32.0],$ ; electrons are spec 0
			[2.5d19,	3.7d19,		1.0,	32.0],$
			[0.5d19,	3.6d19,		1.0,	32.0] ]

	tt = [	[00.2d3,	23.9d3,		6.0,	4.0],$
			[00.3d3,	25.2d3,		6.0,	4.0],$
			[00.3d3,	25.2d3,		6.0,	4.0] ]

	nSpec = n_elements ( amu )
	w	= freq * 2d0 * !dpi

	for n=1,nSpec-1 do begin
		nn[0,0]	+= z[n]*nn[0,n] 
		nn[1,0]	+= z[n]*nn[1,n] 
	endfor

	; Grid

	nR = 100
	nZ = 200

	rMin = 3.5
	rMax = 9.0

	zMin = -5.0
	zMax = +5.0

	r = fIndGen(nR)/(nR-1)*(rMax-rMin)+rMin
	z = fIndGen(nZ)/(nZ-1)*(zMax-zMin)+zMin

	r2d = rebin ( r, nR, nZ )
	z2d = transpose ( rebin ( z, nZ, nR ) )

	; Load eqdsk file

	g = readgeqdsk ( eqdskFileName, /noToroidalFlux)

	; Get b field values on this new grid

	br = interpolate ( g.br, ( r2d - min(g.r) ) / (max(g.r)-min(g.r)) * (n_elements(g.r)-1), $
		( z2d - min ( g.z ) ) / (max(g.z)-min(g.z)) * (n_elements(g.z)-1) )
	bt = interpolate ( g.bPhi, ( r2d - min(g.r) ) / (max(g.r)-min(g.r)) * (n_elements(g.r)-1), $
		( z2d - min ( g.z ) ) / (max(g.z)-min(g.z)) * (n_elements(g.z)-1) )
	bz = interpolate ( g.bz, ( r2d - min(g.r) ) / (max(g.r)-min(g.r)) * (n_elements(g.r)-1), $
		( z2d - min ( g.z ) ) / (max(g.z)-min(g.z)) * (n_elements(g.z)-1) )

	bMag = sqrt ( br^2 + bt^2 + bz^2 )


	; Get psi (poloidal flux) for profile generation

	psi = interpolate ( g.psizr, ( r2d - min(g.r) ) / (max(g.r)-min(g.r)) * (n_elements(g.r)-1), $
		( z2d - min ( g.z ) ) / (max(g.z)-min(g.z)) * (n_elements(g.z)-1) )
	psiNorm = (psi-g.simag) / (g.siBry - g.siMag)


	; Create masks for outside the LCFS and Limiting structure

	print, 'Creating masks (may take a while) ...'

	oversample_boundary, g.rbbbs, g.zbbbs, rbbbs_os, zbbbs_os 
	oversample_boundary, g.rlim, g.zlim, rlim_os, zlim_os 

	mypoly=obj_new('IDLanROI',rbbbs_os,zbbbs_os,type=2)
	mask_bbb = mypoly->ContainsPoints(r2d[*],z2d[*])
	mask_bbb = reform(mask_bbb,nR,nZ)
	iiMaskIn_bbb = where ( mask_bbb eq 1 )
	iiMaskOut_bbb = where ( mask_bbb eq 0 )

	mypoly=obj_new('IDLanROI',rlim_os,zlim_os,type=2)
	mask_lim = mypoly->ContainsPoints(r2d[*],z2d[*])
	mask_lim = reform(mask_lim,nR,nZ)
	iiMaskIn_lim = where ( mask_lim eq 1 )
	iiMaskOut_lim = where ( mask_lim eq 0 )


	print, 'DONE'


	; Create a distance from these surfaces

	print, 'Creating distances from surfaces ... '

	d_bbb = fltarr ( nR, nZ )
	d_lim = fltarr ( nR, nZ )

	for i=0, nR-1 do begin
		for j=0, nZ-1 do begin

			distbbb	= sqrt ( (z2d[i,j] - zbbbs_os)^2 + (r2d[i,j] - rbbbs_os)^2 )
			tmp	= min ( distbbb, iiMindistbbb )
			distbbb = sqrt ( (zbbbs_os[iiMindistbbb] - g.zmaxis)^2 $
				+ (rbbbs_os[iiMindistbbb] - g.rmaxis)^2 )
			d_bbb[i,j]	= $
					sqrt ( (z2d[i,j] - g.zmaxis)^2 $
						+ (r2d[i,j] - g.rmaxis)^2 ) - distbbb

			distlim	= sqrt ( (z2d[i,j] - zlim_os)^2 + (r2d[i,j] - rlim_os)^2 )
			tmp	= min ( distlim, iiMindistlim )
			distlim = sqrt ( (zlim_os[iiMindistlim] - g.zmaxis)^2 $
				+ (rlim_os[iiMindistlim] - g.rmaxis)^2 )
			d_lim[i,j]	= $
					sqrt ( (z2d[i,j] - g.zmaxis)^2 $
						+ (r2d[i,j] - g.rmaxis)^2 ) - distlim

		endfor
	endfor


	br = MakePeriodic ( br, mask_lim);, /look )
	bt = MakePeriodic ( bt, mask_lim);, /look )
	bz = MakePeriodic ( bz, mask_lim);, /look )

	bMag = sqrt ( br^2 + bt^2 + bz^2 )


	; Create profiles

	density_m3	= fltArr ( nR, nZ, nSpec )
	temp_eV	= fltArr ( nR, nZ, nSpec )

	ar2_create_profiles, nn, tt, nR, nZ, PsiNorm, Mask_bbb, d_bbb, $
			Density_m3, Temp_eV

	; Write netCdf file

	outFileName	= 'ar2_bFieldInput.nc'
	nc_id	= nCdf_create ( outFileName, /clobber )
	nCdf_control, nc_id, /fill
	
	nR_id	= nCdf_dimDef ( nc_id, 'nR', nR )
	nz_id	= nCdf_dimDef ( nc_id, 'nZ', nZ )
	nSpec_id	= nCdf_dimDef ( nc_id, 'nSpec', nSpec )
	nbbbs_id	= nCdf_dimDef ( nc_id, 'nbbbs', n_elements(g.rbbbs) )
	nlim_id	= nCdf_dimDef ( nc_id, 'nlim', n_elements(g.rlim) )

	scalar_id	= nCdf_dimDef ( nc_id, 'scalar', 1 )
	
	r_id = nCdf_varDef ( nc_id, 'r', [ nR_id ], /float )
	z_id = nCdf_varDef ( nc_id, 'z', [ nz_id ], /float )
	br_id = nCdf_varDef ( nc_id, 'br', [nR_id, nz_id], /float )
	bt_id = nCdf_varDef ( nc_id, 'bt', [nR_id, nz_id], /float )
	bz_id = nCdf_varDef ( nc_id, 'bz', [nR_id, nz_id], /float )

	Density_id = nCdf_varDef ( nc_id, 'density_m3', [nR_id, nz_id, nSpec_id], /float )
	Temp_id = nCdf_varDef ( nc_id, 'temp_eV', [nR_id, nz_id, nSpec_id], /float )

	rbbbs_id = nCdf_varDef ( nc_id, 'rbbbs', [nbbbs_id], /float )
	zbbbs_id = nCdf_varDef ( nc_id, 'zbbbs', [nbbbs_id], /float )

	rlim_id = nCdf_varDef ( nc_id, 'rlim', [nlim_id], /float )
	zlim_id = nCdf_varDef ( nc_id, 'zlim', [nlim_id], /float )

	nCdf_control, nc_id, /enDef
	
	nCdf_varPut, nc_id, r_id, r 
	nCdf_varPut, nc_id, z_id, z
	nCdf_varPut, nc_id, br_id, br
	nCdf_varPut, nc_id, bt_id, bt 
	nCdf_varPut, nc_id, bz_id, bz

	nCdf_varPut, nc_id, Density_id, Density_m3
	nCdf_varPut, nc_id, Temp_id, Temp_eV

	nCdf_varPut, nc_id, rbbbs_id, g.rbbbs[*]
	nCdf_varPut, nc_id, zbbbs_id, g.zbbbs[*]

	nCdf_varPut, nc_id, rlim_id, g.rlim[*]
	nCdf_varPut, nc_id, zlim_id, g.zlim[*]

	nCdf_close, nc_id

	stop
end
