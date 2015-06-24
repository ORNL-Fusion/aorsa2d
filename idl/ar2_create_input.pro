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

    gaussian_profiles = 0
    parabolic_profiles = 0
    flux_profiles = 0
	fred_namelist_input = 0
	flat_profiles = 1

	br_flat = 0.0
	bt_flat = 1.0
	bz_flat = 0.0

    generate_vorpal_input = 0

	@constants
    @ar2run

	wrf	= freq * 2d0 * !dpi

	r = fIndGen(nR)/(nR-1)*(rMax-rMin)+rMin
	z = fIndGen(nZ)/(nZ-1)*(zMax-zMin)+zMin

	r2d = rebin ( r, nR, nZ )
	z2d = transpose ( rebin ( z, nZ, nR ) )

	x2d = r2d
	y2d = z2d
    if gaussian_profiles or parabolic_profiles or bField_gaussian then x0 = r0

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

		; Get psi (poloidal flux) for profile generation

		psi = interpolate ( g.psizr, ( r2d - min(g.r) ) / (max(g.r)-min(g.r)) * (n_elements(g.r)-1), $
			( z2d - min ( g.z ) ) / (max(g.z)-min(g.z)) * (n_elements(g.z)-1) )
		psiNorm = (psi-g.simag) / (g.siBry - g.siMag)


        rlim = g.rlim[*]
        zlim = g.zlim[*]

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

		;br = MakePeriodic ( br, mask_lim);, /look )
		;bt = MakePeriodic ( bt, mask_lim);, /look )
		;bz = MakePeriodic ( bz, mask_lim);, /look )

	endif else if bField_flat eq 1 then begin

		br = fltArr(nR,nZ)+br_flat
		bt = fltArr(nR,nZ)+bt_flat
		bz = fltArr(nR,nZ)+bz_flat

	endif else if bField_gaussian eq 1 then begin

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


	endif else begin

		bt = r0 * b0 / r2d
		br = bt * br_frac
		bz = bt * bz_frac

	endelse

	; Create masks for outside the LCFS and Limiting structure

	print, 'Creating masks (may take a while) ...'

	oversample_boundary, rlim, zlim, rlim_os, zlim_os 

	mypoly=obj_new('IDLanROI',rlim_os,zlim_os,type=2)
	mask_lim = mypoly->ContainsPoints(r2d[*],z2d[*])
	mask_lim = reform(mask_lim,nR,nZ)<1
	iiMaskIn_lim = where ( mask_lim eq 1 )
	iiMaskOut_lim = where ( mask_lim eq 0 )

	print, 'DONE'

    layout = [5,3]
    !x.margin = !x.margin / 3
    !y.margin = !y.margin / 2
	ScreenSize = get_screen_size()
	dimensions = ScreenSize*0.8
    plotpos = 1
	p=plot(r,bt[*,nZ/2],layout=[layout,plotpos], title='bField',dimensions=dimensions)
	p=plot(r,br[*,nZ/2],/over,color='b')
	p=plot(r,bz[*,nZ/2],/over,color='r')
    ++plotpos
	plotFile = 'inputCreationPlots.pdf'

	bMag = sqrt ( br^2 + bt^2 + bz^2 )
	nSpec = n_elements ( amu )

	; Create profiles

	if fred_namelist_input then begin

		@aorsa2d.in.ejf.idl
		
		_nRho = fred.s_nRho_n
		_nRhoFile = n_elements(fred.s_rho_n_grid)
		_nS = n_elements(fred.s_s_name)
		_nSFile = n_elements(fred.s_m_s)
		_rho = fred.s_rho_n_grid[0:_nRho-1]


		_NumericData_t_eV = reform(fred.s_T_s,_nRhoFile,_nSFile)
		_NumericData_t_eV = _NumericData_t_eV[0:_nRho-1,0:_nS-1]	
		_NumericData_t_eV = [[_rho],[_NumericData_t_ev*1e3]]

		_NumericData_n_m3 = reform(fred.s_n_s,_nRhoFile,_nSFile)
		_NumericData_n_m3 = _NumericData_n_m3[0:_nRho-1,0:_nS-1]	
		_NumericData_n_m3 = [[_rho],[_NumericData_n_m3]]

		_atomicZ = fred.s_q_s[0:_nS-1]/_e
		_amu = fred.s_m_s[0:_nS-1]/_mi

		amu = _amu
		atomicZ = _atomicZ

        _nS = nS_numeric
        amu = amu[0:_nS-1]
        atomicZ = atomicZ[0:_nS-1]

		nSpec = n_elements ( amu )
		NumericData_n_m3 = _NumericData_n_m3[*,0:_nS-1+1]
		NumericData_T_eV = _NumericData_T_eV[*,0:_nS-1+1]
		nS = _nS
		nRho = _nRho
	
		print, 'Read from Fred namelist ...'
		print, 'nS: ', nS
		print, 'amu: ', amu
		print, 'Z: ', atomicZ

	endif else begin
    
		for n=1,nSpec-1 do begin
			nn[0,0]	+= atomicZ[n]*nn[0,n] 
			nn[1,0]	+= atomicZ[n]*nn[1,n] 
		endfor

	endelse


	Density_m3	= fltArr ( nR, nZ, nSpec )
	Temp_eV	= fltArr ( nR, nZ, nSpec )
    nuOmg = fltArr ( nR, nZ, nSpec )

	if flux_profiles eq 1 then begin

			ar2_create_flux_profiles, nSpec, nn, tt, nR, nZ, PsiNorm, Mask_bbb, d_bbb, $
			Density_m3, Temp_eV, DensityMin = DensityMin, TempMin = TempMin, $
            NumericProfiles = Numeric_flux_profiles, NumericData_n_m3 = NumericData_n_m3, $
            NumericData_T_eV = NumericData_T_eV, r2d=r2d, z2d=z2d

	endif else if gaussian_profiles eq 1 then begin

		Density_m3[*] = 0

		SmoothWidth=min([nR,nZ])*0.2
		for s=1,nSpec-1 do begin

			Density_m3[*,*,s] = nn[1,s]*exp(-((x2d-x0)^2/Density_xsig^2+(y2d-y0)^2/Density_ysig^2 ))
			Density_m3[*,*,s] = Density_m3[*,*,s]>DensityMin
			Density_m3[*,*,s] = Smooth(Density_m3[*,*,s],SmoothWidth,/edge_mirror)

			Density_m3[*,*,0] = Density_m3[*,*,0]+Density_m3[*,*,s]*atomicZ[s]


		endfor

		for s=0,nSpec-1 do begin

			Temp_eV[*,*,s] = tt[1,s]*exp(-((x2d-x0)^2/Temp_xsig^2+(y2d-y0)^2/Temp_ysig^2 ))
			Temp_eV[*,*,s] = Temp_eV[*,*,s]>TempMin
			Temp_eV[*,*,s] = Smooth(Temp_eV[*,*,s],SmoothWidth,/edge_mirror)

		endfor


    endif else if parabolic_profiles eq 1 then begin

        Density_m3[*] = 0
		SmoothWidth=min([nR,nZ])*0.2

        ; Density ^2
        x = x2d-x0 
        l = parabolic_half_length
        for s=1,nSpec-1 do begin
            Density_m3[*,*,s] = (nn[1,s]*(1-x^2/l^2))>DensityMin
			Density_m3[*,*,s] = Smooth(Density_m3[*,*,s],SmoothWidth,/edge_mirror)

			Density_m3[*,*,0] = Density_m3[*,*,0]+Density_m3[*,*,s]*atomicZ[s]
        endfor

        ; Temp ^4
        for s=0,nSpec-1 do begin
            Temp_eV[*,*,s] = (tt[1,s]*(1-(x^2/l^2)^2))>TempMin
			Temp_eV[*,*,s] = Smooth(Temp_eV[*,*,s],SmoothWidth,/edge_mirror)
        endfor

	endif else if flat_profiles then begin

		for s=0,nSpec-1 do begin
			Density_m3[*,*,s] = nn[0,s]
			Temp_ev[*,*,s] = tt[0,s]
		endfor

	endif

    yRange = [zMin,zMax]
    xRange = [rMin,rMax]

	if not flat_profiles then begin 
    s_ne = contour(density_m3[*,*,0],r,z, layout=[layout,plotpos],$
            /current,title='density',aspect_ratio=1.0, xRange=xRange, yRange=yRange )
    if bField_eqdsk then p=plot(rLim,zLim,/over)
    ++plotpos

   	s_te = contour(temp_eV[*,*,0],r,z, layout=[layout,plotpos],$
            /current,title='temp',aspect_ratio=1.0, xRange=xRange, yRange=yRange)
    if bField_eqdsk then p=plot(rLim,zLim,/over)
    ++plotpos
	endif

    _c = ['b','g','r','c','m','y','k']
	densityRange=[min(Density_m3),max(Density_m3)]
	p=plot(r,Density_m3[*,nZ/2,0],$
			title='Density [1/m3]',thick=2,$
			layout=[layout,plotpos],/current,yRange=densityRange,/yLog)
    _p = [p]
    for s=1,nSpec-1 do begin
		p=plot(r,Density_m3[*,nZ/2,s],/over,color=_c[s-1])
        _p = [_p,p]
    endfor
    plotpos++	

    _c = ['b','g','r','c','m','y','k']
	densityRange=[min(Density_m3),max(Density_m3)]
	p=plot(r,Density_m3[*,nZ/2,0],$
			title='Density [1/m3]',thick=2,$
			layout=[layout,plotpos],/current,yRange=densityRange)
    _p = [p]
    for s=1,nSpec-1 do begin
		p=plot(r,Density_m3[*,nZ/2,s],/over,color=_c[s-1])
        _p = [_p,p]
    endfor
    plotpos++	

	p=plot(r,Temp_eV[*,nZ/2,0],title='Temp [eV]', thick=2,layout=[layout,plotpos],/current)
    for s=1,nSpec-1 do begin
	    p=plot(r,Temp_eV[*,nZ/2,s],/over)
    endfor
    plotpos++	

    ; Create nuOmg profiles

    @ar2nuomg

    p=plot(r,nuOmg[*,nZ/2,0],title='nuOmg [electrons]',layout=[layout,plotpos],/current)
    ++plotpos

   	@ar2jant

	; Set the jAnt 2-D function

	if fancy_antenna then begin

		if bField_eqdsk then begin	
			rCenter = g.rcentr
			zCenter = 0.0	
		endif

		;get angular points on LCFS with respect to center core

		theta = atan(zlcfs-zCenter, rlcfs-rCenter)*!radeg	

		; Get antenna location on LCFS, the shift away from LCFS
		ii = where(theta gt theta_ant1 and theta le theta_ant2, iiCnt)
		if iiCnt lt 2 then begin
				print, 'ERROR: start and end of the antenna positions are too close together'
				stop
		endif
		antr = rlcfs(ii)
		antz = zlcfs(ii)  
		antr = antr + rshift
		antz = antz + zshift
		
		n_ant = n_elements(antr)
		; Get antenna current: A Guassian peaked at the antenna
		_d = FltArr(nR,nZ,n_ant)
		for k=0,n_ant-1 do begin
		  _d[*,*,k] = sqrt ( ( r2d - antr(k) )^2 + ( z2d - antz(k) )^2 )
		endfor
		d = min(_d,dim=3)
		jAnt =  exp(-(d^2)/(sigma_ant^2))
		
		n = floor(n_ant/2)
		Ju_r = (antr(n) - antr(n+1))/norm([antr(n+1), antz(n+1)] - [antr(n), antz(n)])
		Ju_z = (antz(n) - antz(n+1))/norm([antr(n+1), antz(n+1)] - [antr(n), antz(n)])
		
		jAnt_r = jAnt*Ju_r
		jAnt_z = jAnt*Ju_z
		jAnt_t = jAnt*0.0

	endif else begin

		jAnt =  exp(-(  (r2d-rAnt)^2/antSig_r^2 + (z2d-zAnt)^2/antSig_z^2 ) )

		jAnt_r = jAnt*jRmag
		jAnt_t = jAnt*jTMag
		jAnt_z = jAnt*jZMag

	endelse

   nlevs=10
   scale = 1.1
   levels = (fIndGen(nlevs)+1)/nlevs*scale
   colors=bytscl(levels)
   c1 = contour(jant, r, z,title='Antenna Current', $
   aspect_ratio = 1.0,layout=[layout,plotpos],/fill,/current,$
   c_value=levels,rgb_table=51,rgb_indices=colors, xRange=xRange, yRange=yRange)
   
   p1 = plot(rlim, zlim,/over)
   if fancy_antenna then p1 = plot(antr, antz,/over)
   ++plotpos

	; Get dispersion solution for some nPhi 
	ar2_input_dispersion, wrf, amu, atomicZ, nn, nPhi, nSpec, nR, nZ, $
			Density_m3, bMag, r2D, resonances = resonances, $
			IonIonHybrid_res_freq=IonIonHybrid_res_freq, Spec1=1.0,Spec2=2.0, $
			kPerSq_F=kPerpSq_F,kPerSq_S=kPerpSq_S, $
			StixP=StixP,StixL=StixL,StixR=StixR,StixS=StixS
	
	kPer_F_2D = sqrt(kPerpSq_F)
	kPer_S_2D = sqrt(kPerpSq_S)
	kPer_F = kPer_F_2D[*,nZ/2]
	kPer_S = kPer_S_2D[*,nZ/2]

    ; plot up the ion cyclortron resonances

    for s=1,nSpec-1 do begin
        if s eq 1 then p=plot(r,resonances[*,nZ/2,s], $
                title='Ion cyclotron resonsances',layout=[layout,PlotPos],/current) $
                else p=plot(r,resonances[*,nZ/2,s],/over)
    endfor
    ++PlotPos

	p = plot(r,kPer_F, title='kPer_F',layout=[layout,PlotPos],/current)
    ++PlotPos

	p = plot(r,kPer_S, title='kPer_S',layout=[layout,PlotPos],/current)
    ++PlotPos

	nLevs=31
	range=20
	levels = (findGen(nLevs)+1)/nLevs*range
	levels = [1,2,3,10,20,30,100,200,300,1e3,2e3,4e3]
	levels_map = findgen(n_elements(levels))
	colors = 256-(bytScl(levels_map,top=254)+1)
    if nZ gt 1 then begin
	    c = contour(kPer_F_2D,r,z,layout=[layout,plotpos],$
			c_value=levels,rgb_indices=colors,rgb_table=1,$
            /fill,aspect_ratio=1.0, title='kPer Fast Branch',/current, xRange=xRange, yRange=yRange )
            ++plotpos
            c_zero_set = contour(kPerpSq_F,r,z,/over,c_value=0.001,color='r',C_LABEL_SHOW=0,c_thick=2)
    endif 

	if bfield_eqdsk then p=plot(g.rlim,g.zlim,/over,thick=2)

	nLevs=31
	range=10
	levels = (findGen(nLevs)+1)/nLevs*range
	levels = [1,2,3,10,20,30,100,200,300,1e3,2e3,4e3]
	levels_map = findgen(n_elements(levels))
	colors = 256-(bytScl(levels_map,top=254)+1)
    if nZ gt 1 then begin
	    c = contour(kPer_S_2D,r,z,layout=[layout,plotpos],/current,$
			c_value=levels,rgb_indices=colors,rgb_table=3,$
            /fill,aspect_ratio=1.0, title='kPer Slow Branch', xRange=xRange, yRange=yRange)
        ++PlotPos
		;c.save, plotFile, /append, /close
    endif 

    if bfield_eqdsk eq 0 and size(rLim,/type) eq 0 then begin
            rLim = [rMin, rMax, rMax, rMin, rMin]
            zLim = [zMin, zMin, zMax, zMax, zMin]
    endif

	if bField_eqdsk then p=plot(g.rlim,g.zlim,/over,thick=2)

	; Plot up resonance locations too.

	if bField_eqdsk eq 1 and nZ gt 1 then begin
		p = plot(g.rbbbs,g.zbbbs,thick=2,aspect_ratio=1.0,layout=[layout,PlotPos],$
                /current,title='resonances', xRange=xRange, yRange=yRange)
		p = plot(g.rlim,g.zlim,thick=2,/over)
        ++PlotPos

		for s=1,nSpec-1 do begin
			c=contour(resonances[*,*,s],r,z,c_value=fIndGen(5)/4.0*0.01,/over)
		endfor
		if nSpec gt 2 then $
		c=contour(1/(abs(IonIonHybrid_res_freq mod wrf)/wrf),r,z,c_value=fIndGen(25)*10,/over)
	endif else if flux_profiles eq 1 and nZ gt 1 then begin

		c=contour(psinorm,r,z,aspect_ratio=1.0, title='psi', xRange=xRange, yRange=yRange)	
		p=plot(rlim,zlim,/over,thick=2)
		;p=plot(VorpalBox_r,VorpalBox_z,/over,thick=2,color='b')
		for s=1,nSpec-1 do begin
			c=contour(resonances[*,*,s],r,z,c_value=fIndGen(5)/4.0*0.01,/over)
		endfor
		if nSpec gt 2 then $
		c=contour(1/(abs(IonIonHybrid_res_freq mod wrf)/wrf),r,z,c_value=fIndGen(25)*10,/over)

	endif else if nZ gt 1 then begin
		;for s=1,nSpec-1 do begin
		;	p = plot(r,resonances[*,0,s],/over)
		;endfor
        for s=1,nSpec-1 do begin
			c=contour(resonances[*,*,s],r,z,c_value=fIndGen(5)/4.0*0.01,/over)
		endfor
		if nSpec gt 2 then $
		c=contour(1/(abs(IonIonHybrid_res_freq mod wrf)/wrf),r,z,c_value=fIndGen(25)*10,/over)


    endif

	; Write netCdf file

	;save, freq, nphi, bField_eqdsk, eqdskFileName, flux_profiles, atomicZ, amu, $
	;		nn, tt, nR, nZ, rMin, rMax, zMin, zMax, fileName = 'ar2RunCreationParameters.sav'

	outFileName	= 'ar2Input.nc'
	nc_id	= nCdf_create ( outFileName, /clobber )
	nCdf_control, nc_id, /fill
	
	nR_id	= nCdf_dimDef ( nc_id, 'nR', nR )
	nz_id	= nCdf_dimDef ( nc_id, 'nZ', nZ )
	nSpec_id = nCdf_dimDef ( nc_id, 'nSpec', nSpec )
	nlim_id	= nCdf_dimDef ( nc_id, 'nlim', n_elements(rlim) )

	scalar_id	= nCdf_dimDef ( nc_id, 'scalar', 1 )

	AtomicZ_id = nCdf_varDef ( nc_id, 'AtomicZ', [nSpec_id], /float )
	amu_id = nCdf_varDef ( nc_id, 'amu', [nSpec_id], /float )

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
	jant_id = nCdf_varDef ( nc_id, 'jAnt', [nR_id, nz_id], /float )

	jant_r_id = nCdf_varDef ( nc_id, 'jAnt_r', [nR_id, nz_id], /float )
	jant_z_id = nCdf_varDef ( nc_id, 'jAnt_z', [nR_id, nz_id], /float )
	jant_t_id = nCdf_varDef ( nc_id, 'jAnt_t', [nR_id, nz_id], /float )

	LimMask_id = nCdf_varDef ( nc_id, 'LimMask', [nR_id,nz_id], /short )

	Lim_r_id = nCdf_varDef ( nc_id, 'Lim_r', [nlim_id], /float )
	Lim_z_id = nCdf_varDef ( nc_id, 'Lim_z', [nlim_id], /float )

	kPerSq_F_id = nCdf_varDef ( nc_id, 'kPerSq_F', [nR_id, nz_id], /float )
	kPerSq_S_id = nCdf_varDef ( nc_id, 'kPerSq_S', [nR_id, nz_id], /float )

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
	
	nCdf_varPut, nc_id, jant_id, jAnt
    nCdf_varPut, nc_id, jant_r_id, jAnt_r
    nCdf_varPut, nc_id, jant_z_id, jAnt_z
    nCdf_varPut, nc_id, jant_t_id, jAnt_t

	nCdf_varPut, nc_id, LimMask_id, mask_lim 

    nCdf_varPut, nc_id, kPerSq_F_id, kPerpSq_F
    nCdf_varPut, nc_id, kPerSq_S_id, kPerpSq_S

	nCdf_varPut, nc_id, Lim_r_id, rlim
	nCdf_varPut, nc_id, Lim_z_id, zlim


	nCdf_close, nc_id

	; Write VORPAL-AORSA Coupling text file

    if generate_vorpal_input then begin

	for s=0,nSpec-1 do begin

		VorpalFileName = 'VorpalProfiles_'+string(s,format='(i1)')+'.txt'

		Vorpal_nX = 64
		Vorpal_nY = 36 
		Vorpal_nZ = 64

		Vorpal_xDim = rMax-rMin
		Vorpal_yDim = 0.1 * Vorpal_xDim
		Vorpal_zDim = zMax-zMin

		Vorpal_xMin = rMin
		Vorpal_xMax = rMax
		Vorpal_x_grid = fIndGen(Vorpal_nX)*(Vorpal_xMax-Vorpal_xMin)/(Vorpal_nX-1)+Vorpal_xMin

		Vorpal_yMin = -Vorpal_yDim/2.0
		Vorpal_yMax = +Vorpal_yDim/2.0
		Vorpal_y_grid = fIndGen(Vorpal_nY)*(Vorpal_yMax-Vorpal_yMin)/(Vorpal_nY-1)+Vorpal_yMin

		Vorpal_zMin = zMin
		Vorpal_zMax = zMax
		Vorpal_z_grid = fIndGen(Vorpal_nZ)*(Vorpal_zMax-Vorpal_zMin)/(Vorpal_nz-1)+Vorpal_zMin

		openw, lun, VorpalFileName, /get_lun

		printf, lun, 'nX: '+string(Vorpal_nX, format='(i4.4)')
		printf, lun, 'nY: '+string(Vorpal_nY, format='(i4.4)')
		printf, lun, 'nZ: '+string(Vorpal_nZ, format='(i4.4)')
		printf, lun, 'amu: '+string(amu[s],format='(f12.10)')
		printf, lun, 'AtomicZ: ',+string(AtomicZ[s],format='(f+6.2)')

		printf, lun, 'X, Y, Z, Bx[T], By[T], Bz[T], T[eV], n[m^-3], nuOmg'

		for i=0,Vorpal_nX-1 do begin
			for j=0,Vorpal_nY-1 do begin
				for k=0,Vorpal_nZ-1 do begin

                    thisX = Vorpal_x_grid[i]
                    thisY = Vorpal_y_grid[j]
                    thisZ = Vorpal_z_grid[k]

                    thisI = (thisX-rMin)/(rMax-rMin)*(nR-1)
                    thisJ = (thisZ-zMin)/(zMax-zMin)*(nZ-1)

		            this_Br = interpolate(br,thisI,thisJ,cubic=-0.5)
		            this_Bt = interpolate(bt,thisI,thisJ,cubic=-0.5)
		            this_Bz = interpolate(bz,thisI,thisJ,cubic=-0.5)

                    this_Bx = this_Br
                    this_By = this_Bt
                    this_Bz = this_Bz

		            this_nuOmg = interpolate(nuOmg[*,*,s],thisI,thisJ,cubic=-0.5)
		            this_T_eV = interpolate(Temp_eV[*,*,s],thisI,thisJ,cubic=-0.5)
		            this_n_m3 = interpolate(density_m3[*,*,s],thisI,thisJ,cubic=-0.5)

					printf, lun, thisX, thisY, thisZ, $
						this_Bx, this_By, this_Bz, this_T_eV, this_n_m3, this_nuOmg, $
                        format='(7(f10.3,1x),2(e12.3,1x))';, $
						;V_T[i,j,k], V_n[i,j,k], format='(7(f10.3,1x),e12.3)'

				endfor
			endfor
		endfor

		close, lun

	endfor

    endif
	p.save, plotFile
	stop
exit ; this is here for OMFit - don't remove
end
