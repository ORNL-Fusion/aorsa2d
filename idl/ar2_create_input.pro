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

	bField_eqdsk = 1
	bField_gaussian = 0
    bField_flat = 0

	br_flat = 0.0
	bt_flat = 1.0
	bz_flat = 0.0

	@constants

	;@gorden_bell
	;@gorden_bell_b
	;@ar2_run_langmuir
	;@ar2_run_nstxslow
	;@ar2_run_ar_vo_bench
    @ar2_run_coupling_right_simple
    ;@ar2_run_coupling_left_simple

	nSpec = n_elements ( amu )
	wrf	= freq * 2d0 * !dpi

	for n=1,nSpec-1 do begin
		nn[0,0]	+= atomicZ[n]*nn[0,n] 
		nn[1,0]	+= atomicZ[n]*nn[1,n] 
	endfor

	r = fIndGen(nR)/(nR-1)*(rMax-rMin)+rMin
	z = fIndGen(nZ)/(nZ-1)*(zMax-zMin)+zMin

	r2d = rebin ( r, nR, nZ )
	z2d = transpose ( rebin ( z, nZ, nR ) )

	x2d = r2d
	y2d = z2d


	if bField_eqdsk eq 1 then begin

		; Load eqdsk file

		g = readgeqdsk ( eqdskFileName, /noToroidalFlux, bTorFactor = bTorFactor)

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

	endif else if bField_flat eq 1 then begin

		mask_bbb = FltArr(nR,nZ)+1	
		mask_lim = FltArr(nR,nZ)+1	

		br = fltArr(nR,nZ)+br_flat
		bt = fltArr(nR,nZ)+bt_flat
		bz = fltArr(nR,nZ)+bz_flat

	endif else if bField_gaussian eq 1 then begin

		oversample_boundary, rlim, zlim, rlim_os, zlim_os 

		mypoly=obj_new('IDLanROI',rlim_os,zlim_os,type=2)
		mask_lim = mypoly->ContainsPoints(r2d[*],z2d[*])
		mask_lim = reform(mask_lim,nR,nZ)
		iiMaskIn_lim = where ( mask_lim eq 1 )
		iiMaskOut_lim = where ( mask_lim eq 0 )

		bt = r0 * b0 / r2d

		psi = 1 - exp(-( (x2d-x0)^2/psi_xsig^2+(y2d-y0)^2/psi_ysig^2 ))

		bx2d = fltArr(nR,nZ)
		by2d = fltArr(nR,nZ)

		for j=0,nZ-1 do begin
			by2d[*,j] = deriv(psi[*,j])
		endfor
		for i=0,nR-1 do begin
			bx2d[i,*] = deriv(psi[i,*])
		endfor

		bp_mag = sqrt(bx2d^2+by2d^2)

		bx2d = bx2d/max(bp_mag)*b0*bpmax_b0
		by2d = by2d/max(bp_mag)*b0*bpmax_b0

		br = bx2d
		bz = by2d

		p=plot(r,bt[*,nZ/2])
		p=plot(r,br[*,nZ/2],/over,color='b')
		p=plot(r,bz[*,nZ/2],/over,color='r')

	endif else begin

		mask_bbb = FltArr(nR,nZ)+1	
		mask_lim = FltArr(nR,nZ)+1	

		bt = r0 * b0 / r2d
		br = bt * br_frac
		bz = bt * bz_frac

	endelse

	bMag = sqrt ( br^2 + bt^2 + bz^2 )

	; Create profiles

	Density_m3	= fltArr ( nR, nZ, nSpec )
	Temp_eV	= fltArr ( nR, nZ, nSpec )
    nuOmg = fltArr ( nR, nZ, nSpec )

	if flux_profiles eq 1 then begin

		ar2_create_flux_profiles, nSpec, nn, tt, nR, nZ, PsiNorm, Mask_bbb, d_bbb, $
			Density_m3, Temp_eV, DensityMin = DensityMin, TempMin = TempMin

	endif else if gaussian_profiles eq 1 then begin

		Density_m3[*] = 0

		SmoothWidth=min([nR,nZ])*0.2
		for s=1,nSpec-1 do begin

			Density_m3[*,*,s] = nn[1,s]*exp(-((x2d-x0)^2/Density_xsig^2+(y2d-y0)^2/Density_ysig^2 ))
			Density_m3[*,*,s] = Density_m3[*,*,s]>DensityMin
			Density_m3[*,*,s] = Smooth(Density_m3[*,*,s],SmoothWidth,/edge_mirror)

			Density_m3[*,*,0] = Density_m3+Density_m3[*,*,s]*atomicZ[s]


		endfor

		for s=0,nSpec-1 do begin

			Temp_eV[*,*,s] = tt[1,s]*exp(-((x2d-x0)^2/Temp_xsig^2+(y2d-y0)^2/Temp_ysig^2 ))
			Temp_eV[*,*,s] = Temp_eV[*,*,s]>TempMin
			Temp_eV[*,*,s] = Smooth(Temp_eV[*,*,s],SmoothWidth,/edge_mirror)

		endfor


	endif else begin

		for s=0,nSpec-1 do begin
			Density_m3[*,*,s] = nn[0,s]
			Temp_ev[*,*,s] = tt[0,s]
		endfor

	endelse
    
	p=plot(r,Density_m3[*,nZ/2,0],title='Density [1/m3]',/ylog)
	p=plot(r,Temp_eV[*,nZ/2,0],title='Temp [eV]')


    ; Create nuOmg profiles

    @ar2_run_coupling_right_simple_nuomg
    ;@ar2_run_coupling_left_simple_nuomg

    p=plot(r,nuOmg[*,nZ/2,0],title='nuOmg [electrons]')

   ;nSmooth = 5 
    ;for s=0,nSpec-1 do begin
    ;    for n=0,nSmooth-1 do begin
    ;        nuOmg[*,*,s] = smooth(nuOmg[*,*,s]>MinNuOmg,min([nR,nZ])*0.05,/edge_truncate)
    ;    endfor
    ;endfor

	; Look at the dispersion relation for these data

	n_nPhi = 101 
	nPhiMin = -n_nPhi/2
	nPhiArray = IndGen(n_nPhi)+nPhiMin

	kPerp_F = complexArr(nR,n_nPhi)
	kPerp_S = complexArr(nR,n_nPhi)

	kPerp_F2D_avg = complexArr(nR,nZ)
	kPerp_S2D_avg = complexArr(nR,nZ)

	StixP_nPhi = FltArr(nR,n_nPhi)
	StixS_nPhi = FltArr(nR,n_nPhi)
	StixR_nPhi = FltArr(nR,n_nPhi)
	StixL_nPhi = FltArr(nR,n_nPhi)

	for nphi_i=0,n_nPhi-1 do begin

	ar2_input_dispersion, wrf, amu, atomicZ, nn, nPhiArray[nphi_i], nSpec, nR, nZ, $
			Density_m3, bMag, r2D, resonances = resonances, $
			IonIonHybrid_res_freq=IonIonHybrid_res_freq, Spec1=1.0,Spec2=2.0, $
			kPerSq_F=kPerpSq_F,kPerSq_S=kPerpSq_S, $
			StixP=StixP,StixL=StixL,StixR=StixR,StixS=StixS
		
		slice = nZ/4
		kPerp_F[*,nphi_i] = sqrt(kPerpSq_F[*,slice])
		kPerp_S[*,nphi_i] = sqrt(kPerpSq_S[*,slice])

		kPerp_F2D_avg = kPerp_F2D_avg + sqrt(kPerpSq_F)
		kPerp_S2D_avg = kPerp_S2D_avg + sqrt(kPerpSq_S)

	endfor

	kPerp_F2D_avg = kPerp_F2D_avg/n_nPhi
	kPerp_S2D_avg = kPerp_S2D_avg/n_nPhi

	nLevs=21
	range=20
	levels = (findGen(nLevs)+1)/nLevs*range
	colors = 256-(bytScl(levels,top=254)+1)
	c = contour(kPerp_F,r,nPhiArray,c_value=levels,rgb_indices=colors,$
            rgb_table=1,/fill, title='kPerp_F')

	nLevs=31
	range=500
	levels = (findGen(nLevs)+1)/nLevs*range
	levels = [1,2,3,10,20,30,100,200,300,1e3,2e3,4e3]
	levels_map = findgen(n_elements(levels))
	colors = 256-(bytScl(levels_map,top=254)+1)
	c = contour(kPerp_S,r,nPhiArray,c_value=levels,rgb_indices=colors,$
            rgb_table=3,/fill, title='kPerp_S')

	c = contour(stixp,r,z,layout=[2,2,1],n_levels=21,title='StixP')
	c = contour(stixs,r,z,layout=[2,2,2],n_levels=21,/current,title='StixS')
	c = contour(stixr,r,z,layout=[2,2,3],n_levels=21,/current,title='StixR')
	c = contour(stixl,r,z,layout=[2,2,4],n_levels=21,/current,title='StixL')

	nLevs=31
	range=20
	levels = (findGen(nLevs)+1)/nLevs*range
	levels = [1,2,3,10,20,30,100,200,300,1e3,2e3,4e3]
	levels_map = findgen(n_elements(levels))
	colors = 256-(bytScl(levels_map,top=254)+1)
	c = contour(kPerp_F2D_avg,r,z,layout=[2,1,1],$
			c_value=levels,rgb_indices=colors,rgb_table=1,$
            /fill,aspect_ratio=1.0, title='kPerp_F2D_avg')
	p=plot(rlim,zlim,/over,thick=2)

	p=plot(RightSide_rlim,RightSide_zlim,/over,thick=2)
	p=plot(LeftSide_rlim,LeftSide_zlim,/over,thick=2)

	;p=plot(VorpalBox_r,VorpalBox_z,/over,thick=2,color='b')
	
	nLevs=31
	range=500
	levels = (findGen(nLevs)+1)/nLevs*range
	levels = [1,2,3,10,20,30,100,200,300,1e3,2e3,4e3]
	levels_map = findgen(n_elements(levels))
	colors = 256-(bytScl(levels_map,top=254)+1)
	c = contour(kPerp_S2D_avg,r,z,layout=[2,1,2],/current,$
			c_value=levels,rgb_indices=colors,rgb_table=3,$
            /fill,aspect_ratio=1.0, title='kPerp_S2D_avg')
	p=plot(rlim,zlim,/over,thick=2)
	;p=plot(VorpalBox_r,VorpalBox_z,/over,thick=1,color='b',transparency=50)
    p=plot(rDomainBox,zDomainBox,/over,thick=1,linestyle='dash')
	
	; Plot up resonance locations too.

	if bField_eqdsk eq 1 and nZ gt 1 then begin
		p = plot(g.rbbbs,g.zbbbs,thick=2,aspect_ratio=1.0)
		p = plot(g.rlim,g.zlim,thick=2,/over)

		for s=1,nSpec-1 do begin
			c=contour(resonances[*,*,s],r,z,c_value=fIndGen(5)/4.0*0.01,/over)
		endfor
		if nSpec gt 2 then $
		c=contour(1/(abs(IonIonHybrid_res_freq mod wrf)/wrf),r,z,c_value=fIndGen(25)*10,/over)
	endif else if bField_eqdsk eq 0 and nZ gt 1 then begin

		c=contour(psi,r,z,aspect_ratio=1.0, title='psi')	
		p=plot(rlim,zlim,/over,thick=2)
		;p=plot(VorpalBox_r,VorpalBox_z,/over,thick=2,color='b')
		for s=1,nSpec-1 do begin
			c=contour(resonances[*,*,s],r,z,c_value=fIndGen(5)/4.0*0.01,/over)
		endfor
		if nSpec gt 2 then $
		c=contour(1/(abs(IonIonHybrid_res_freq mod wrf)/wrf),r,z,c_value=fIndGen(25)*10,/over)

	endif else if bField_eqdsk eq 0 and nZ eq 0 then begin
		for s=1,nSpec-1 do begin
			p = plot(r,resonances[*,0,s],/over)
		endfor
	endif

	; Write netCdf file

	save, freq, nphi, bField_eqdsk, eqdskFileName, flux_profiles, atomicZ, amu, $
			nn, tt, nR, nZ, rMin, rMax, zMin, zMax, fileName = 'ar2RunCreationParameters.sav'

	outFileName	= 'ar2Input.nc'
	nc_id	= nCdf_create ( outFileName, /clobber )
	nCdf_control, nc_id, /fill
	
	nR_id	= nCdf_dimDef ( nc_id, 'nR', nR )
	nz_id	= nCdf_dimDef ( nc_id, 'nZ', nZ )
	nSpec_id	= nCdf_dimDef ( nc_id, 'nSpec', nSpec )
	nlim_id	= nCdf_dimDef ( nc_id, 'nlim', n_elements(rlim) )

	scalar_id	= nCdf_dimDef ( nc_id, 'scalar', 1 )

	AtomicZ_id = nCdf_varDef ( nc_id, 'AtomicZ', [nSpec_id], /short )
	amu_id = nCdf_varDef ( nc_id, 'amu', [nSpec_id], /short )

	rMin_id = nCdf_varDef ( nc_id, 'rMin', [scalar_id], /float )
	rMax_id = nCdf_varDef ( nc_id, 'rMax', [scalar_id], /float )
	zMin_id = nCdf_varDef ( nc_id, 'zMin', [scalar_id], /float )
	zMax_id = nCdf_varDef ( nc_id, 'zMax', [scalar_id], /float )

	r_id = nCdf_varDef ( nc_id, 'r', [ nR_id ], /float )
	z_id = nCdf_varDef ( nc_id, 'z', [ nz_id ], /float )
	br_id = nCdf_varDef ( nc_id, 'br', [nR_id, nz_id], /float )
	bt_id = nCdf_varDef ( nc_id, 'bt', [nR_id, nz_id], /float )
	bz_id = nCdf_varDef ( nc_id, 'bz', [nR_id, nz_id], /float )

	Density_id = nCdf_varDef ( nc_id, 'Density_m3', [nR_id, nz_id, nSpec_id], /float )
	Temp_id = nCdf_varDef ( nc_id, 'Temp_eV', [nR_id, nz_id, nSpec_id], /float )
	nuOmg_id = nCdf_varDef ( nc_id, 'nuOmg', [nR_id, nz_id, nSpec_id], /float )

	LimMask_id = nCdf_varDef ( nc_id, 'LimMask', [nR_id,nz_id], /short )
	;BbbMask_id = nCdf_varDef ( nc_id, 'BbbMask', [nR_id,nz_id], /short )

	Lim_r_id = nCdf_varDef ( nc_id, 'Lim_r', [nlim_id], /float )
	Lim_z_id = nCdf_varDef ( nc_id, 'Lim_z', [nlim_id], /float )

	nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, rMin_id, rMin
	nCdf_varPut, nc_id, rMax_id, rMax
	nCdf_varPut, nc_id, zMin_id, zMin
	nCdf_varPut, nc_id, zMax_id, zMax

	nCdf_varPut, nc_id, r_id, r 
	nCdf_varPut, nc_id, z_id, z
	nCdf_varPut, nc_id, br_id, br
	nCdf_varPut, nc_id, bt_id, bt 
	nCdf_varPut, nc_id, bz_id, bz

	nCdf_varPut, nc_id, AtomicZ_id, AtomicZ 
	nCdf_varPut, nc_id, amu_id, amu 

	nCdf_varPut, nc_id, Density_id, Density_m3
	nCdf_varPut, nc_id, Temp_id, Temp_eV
	nCdf_varPut, nc_id, nuOmg_id, nuOmg

	nCdf_varPut, nc_id, LimMask_id, mask_lim 
	;nCdf_varPut, nc_id, BbbMask_id, mask_bbb

	nCdf_varPut, nc_id, Lim_r_id, rlim
	nCdf_varPut, nc_id, Lim_z_id, zlim

	nCdf_close, nc_id

	; Write VORPAL-AORSA Coupling text file

	for s=0,nSpec-1 do begin

		VorpalFileName = 'VorpalProfiles_'+string(s,format='(i1)')+'.txt'

		Vorpal_nX = 64
		Vorpal_nY = 32
		Vorpal_nZ = 16

		Vorpal_xDim = 0.8
		Vorpal_yDim = 0.4
		Vorpal_zDim = 0.36

		Vorpal_xMin = min(VorpalBox_r)
		Vorpal_xMax = max(VorpalBox_r)
		Vorpal_x_grid = fIndGen(Vorpal_nX)*(Vorpal_xMax-Vorpal_xMin)/(Vorpal_nX-1)+Vorpal_xMin

		Vorpal_yMin = -Vorpal_yDim/2.0
		Vorpal_yMax = +Vorpal_yDim/2.0
		Vorpal_y_grid = fIndGen(Vorpal_nY)*(Vorpal_yMax-Vorpal_yMin)/(Vorpal_nY-1)+Vorpal_yMin

		Vorpal_zMin = min(VorpalBox_z)
		Vorpal_zMax = max(VorpalBox_z)
		Vorpal_z_grid = fIndGen(Vorpal_nZ)*(Vorpal_zMax-Vorpal_zMin)/(Vorpal_nz-1)+Vorpal_zMin

		Vx3D = rebin(Vorpal_x_grid, Vorpal_nX, Vorpal_nY, Vorpal_nZ)
		Vy3D = transpose(rebin(Vorpal_y_grid, Vorpal_nY, Vorpal_nZ, Vorpal_nX),[2,0,1])
		Vz3D = transpose(rebin(Vorpal_z_grid, Vorpal_nZ, Vorpal_nX, Vorpal_nY),[1,2,0])

		Vr3D = sqrt(Vx3D^2+Vy3D^2)
		Vt3D = atan(Vy3D,Vx3D)

		openw, lun, VorpalFileName, /get_lun

		printf, lun, 'nX: '+string(Vorpal_nX, format='(i4.4)')
		printf, lun, 'nY: '+string(Vorpal_nY, format='(i4.4)')
		printf, lun, 'nZ: '+string(Vorpal_nZ, format='(i4.4)')
		printf, lun, 'amu: '+string(amu[s],format='(f12.10)')
		printf, lun, 'AtomicZ: ',+string(AtomicZ[s],format='(f+6.2)')

		i3D = (Vr3D - min(r))/(max(r)-min(r))*(nR-1)
		j3D = (Vz3D - min(z))/(max(z)-min(z))*(nZ-1)

		V_Br = reform(interpolate(br,i3D[*],j3D[*],cubic=-0.5),Vorpal_nX,Vorpal_nY,Vorpal_nZ)
		V_Bt = reform(interpolate(bt,i3D[*],j3D[*],cubic=-0.5),Vorpal_nX,Vorpal_nY,Vorpal_nZ)
		V_Bz = reform(interpolate(bz,i3D[*],j3D[*],cubic=-0.5),Vorpal_nX,Vorpal_nY,Vorpal_nZ)

		V_Bx = cos(Vt3D)*V_Br-sin(Vt3D)*V_Bt
		V_By = sin(Vt3D)*V_Br+cos(Vt3D)*V_Bt

		V_T = reform(interpolate(Temp_eV[*,*,s],i3D[*],j3D[*],cubic=-0.5),Vorpal_nX,Vorpal_nY,Vorpal_nZ)
		V_n = reform(interpolate(Density_m3[*,*,s],i3D[*],j3D[*],cubic=-0.5),Vorpal_nX,Vorpal_nY,Vorpal_nZ)

		printf, lun, 'X, Y, Z, Bx[T], By[T], Bz[T], T[eV], n[m^-3]'

		for i=0,Vorpal_nX-1 do begin
			for j=0,Vorpal_nY-1 do begin
				for k=0,Vorpal_nZ-1 do begin

					printf, lun, Vx3D[i,j,k], Vy3D[i,j,k], Vz3D[i,j,k], $
						V_Bx[i,j,k], V_By[i,j,k], V_Bz[i,j,k], $
						V_T[i,j,k], V_n[i,j,k], format='(7(f10.3,1x),e12.3)'

				endfor
			endfor
		endfor

		close, lun

	endfor

	stop
end
