pro compare_ar2_genray

	gFileName = 'param_along_ray_nphi71.txt'
	;gFileName = 'param_along_ray_nphi94.txt'

	g = read_genray(gFileName)
	__n = n_elements(g)
	downSampleFac = 10
	__i = IndGen(__n/downSampleFac)*downSampleFac
	g = g[__i]
   	traj_r = g.r
    traj_z = g.z

	nPhi = -71
	;nPhi = -94
	freq = 500e6

    ;nPhi = -21
    ;freq = 30e6

	eqdskFileName = 'g122976.03021'
    ;eqdskFileName = 'NSTX_130608A03_trxpl_0.410_plasma_state.geq'
	eqdsk = readgeqdsk(eqdskFileName,bInterpS=bS)

	;flipGenRayParDir = 0
	;if flipGenRayParDir then begin
	;	g.npar = -g.npar
	;endif

 	useAR2File = 0
	if useAR2File then begin
		a = ar2_read_solution('./',1)
		ar2_read_ar2input,'ar2Input.nc',ar2=d
		r = ar2_read_rundata('./',1)

		_r = a.r
		_z = a.z

		_er = a.e_r
		_et = a.e_t
		_ez = a.e_z

		nPhi = r.nPhi
		freq = r.freq
	endif

	useFredFile = 1	
	if useFredFile then begin
		eLabFrameFileName = 'E_lab_frame'
		;eLabFrameFileName = 'E_lab_frame-damp-500a'
		;vtkFileName = 'Efield_2D_PARKMURAKAMI.vtk'

		f = read_e_lab_frame(eLabFrameFileName)
		;f = read_vtk(vtkFileName)

		_r = f.r
		_z = f.z
		_er = f.er
		_et = f.et
		_ez = f.ez
	endif

	eMag = sqrt(_er^2+_et^2+_ez^2)

	useMaxFieldTrajectory = 0
	if useMaxFieldTrajectory then begin
        _ii = where(g.r gt 1.7)
        g = g[_ii]
        traj_r = g.r - (g.distance*0.0005)^2
        traj_z = g.z + (g.distance*0.003)^2

	endif

	ScreenSize = get_screen_size()
	nLevs = 11
	scale = max(eMag)/6
	levels = (fIndGen(nLevs)+1)/(nLevs)*scale
	colors = bytScl(levels,top=253)+1
	c=contour(eMag,_r,_z,aspect_ratio=1.0,layout=[2,1,1], $
			/current,xtitle='R [m]',yTitle='z [m]', $
			c_value=levels, rgb_table=55, rgb_indices=colors, /fill,$
			dimensions=ScreenSize*0.8,title='AORSA |E| '+' nPhi: '+string(nPhi))
	
	p=plot(eqdsk.rlim,eqdsk.zlim,/over)
	p=plot(eqdsk.rbbbs,eqdsk.zbbbs,/over)
	if useMaxFieldTrajectory then p=plot(traj_r[*],traj_z[*],/over,thick=2)
	p=plot(g.r[*],g.z[*],/over, color='b',thick=2)

	nGenRay = n_elements(traj_r)

	kParMin = -150
	kParMax = +150
	nkPar = 80
	kParGrid = fIndgen(nkpar)/(nkPar-1)*(kParMax-kParMin)+kParMin
	kParCnt = intarr(nkpar)
	pPar = fltArr(nkpar,nGenRay)

	kPerMin = 0
	kPerMax = 1400
	nkPer = 80
	kPerGrid = fIndgen(nkper)/(nkPer-1)*(kPerMax-kPerMin)+kPerMin
	kPerCnt = intarr(nkper)
	pPer = fltArr(nkper,nGenRay)

	_nf = 150
	__pPar=fltArr(nGenRay,_nF*2-1)
	__pPer=fltArr(nGenRay,_nF*2-1)

	for n=0,nGenRay-1 do begin

		this_r = traj_r[n]
		this_z = traj_z[n]
		i = (this_r-eqdsk.r[0])/(eqdsk.r[-1]-eqdsk.r[0])*(n_elements(eqdsk.r)-1)
		j = (this_z-eqdsk.z[0])/(eqdsk.z[-1]-eqdsk.z[0])*(n_elements(eqdsk.z)-1)

		this_br = interpolate(eqdsk.br,i,j,cubic=-0.5)
		this_bt = interpolate(eqdsk.bphi,i,j,cubic=-0.5)
		this_bz = interpolate(eqdsk.bz,i,j,cubic=-0.5)

		bMag = sqrt(this_br^2+this_bt^2+this_bz^2)
		_bru = this_br/bMag
		_btu = this_br/bMag
		_bzu = this_br/bMag

		rSpan = 80
		zSpan = 80

		ia = where(abs(_r-this_r) eq min(abs(_r-this_r)))	
		ja = where(abs(_z-this_z) eq min(abs(_z-this_z)))	
		ja = ja[0]

		this_r_span_r = ((ia+rSpan)<(n_elements(_r)-1))-ia
		this_r_span_l = -(((ia-rSpan)>0)-ia)
		this_r_span = min([this_r_span_r,this_r_span_l])

		this_z_span_r = ((ja+zSpan)<(n_elements(_z)-1))-ja
		this_z_span_l = -(((ja-zSpan)>0)-ja)
		this_z_span = min([this_z_span_r,this_z_span_l])

		e_r_r = _er[ia-this_r_span:ia+this_r_span,ja]
		e_t_r = _et[ia-this_r_span:ia+this_r_span,ja]
		e_z_r = _ez[ia-this_r_span:ia+this_r_span,ja]

		e_r_z = _er[ia,ja-this_z_span:ja+this_z_span]
		e_t_z = _et[ia,ja-this_z_span:ja+this_z_span]
		e_z_z = _ez[ia,ja-this_z_span:ja+this_z_span]

		_nR = n_elements(e_r_r)
		_nZ = n_elements(e_r_z)
		_dR = _r[1]-_r[0]
		_dZ = _z[1]-_z[0]

		han = hanning(_nR)
		this_E = han*e_r_r
		p_r = abs(fft(this_E,/center))^2
		k_r = fIndGen(_nR)*2*!pi/(_nR*_dR)-_nR*!pi/(_nR*_dR)

		han = hanning(_nZ)
		this_E = han*e_r_z
		p_z = abs(fft(this_E,/center))^2
		k_z = fIndGen(_nZ)*2*!pi/(_nZ*_dZ)-_nZ*!pi/(_nZ*_dZ)

		;if n eq 100 then stop

		p_t = 0
		k_t = nPhi/this_r

		kParCnt[*] = 0
		kPerCnt[*] = 0

		; try an explicit field line trace and interp to get kPar
		_pt = [this_r,0,this_z]
		_ds = 0.005
		
		_line1 = dlg_fieldlinetrace ( bS, _pt, direction = 1, ds = _ds, nS = _nF, $ 
			B_FieldLine_CYL = _line_b1, perp = 0 )

		_line2 = dlg_fieldlinetrace ( bS, _pt, direction = -1, ds = _ds, nS = _nF, $ 
			B_FieldLine_CYL = _line_b2, perp = 0 )

		fLine = [[reverse(_line1[*,0:-2],2)],[_line2[*,1:-2]]]
		bLine = [[reverse(_line_b1[*,0:-2],2)],[_line_b2[*,1:-2]]]

		_i = (fLine[0,*]-min(_r))/(max(_r)-min(_r))*(n_elements(_r)-1)
		_j = (fLine[2,*]-min(_z))/(max(_z)-min(_z))*(n_elements(_z)-1)
		_ii = complex(0,1)
		__er = interpolate(_er,_i,_j,cubic=-0.5)*exp(-_ii*nPhi*fLine[1,*])
		this_E = hanning(n_elements(__er))*__er
		;this_E = __er
		_nPar = n_elements(this_E)
		__pPar[n,*] = abs(fft(this_E,/center))^2
		__kPar = fIndGen(_nPar)*2*!pi/(_nPar*_dS)-_nPar*!pi/(_nPar*_dS)

		; now do the same for some perp direction

		_ds = 0.001

		_line1 = dlg_fieldlinetrace ( bS, _pt, direction = 1, ds = _ds, nS = _nF, $ 
			B_FieldLine_CYL = _line_b1, perp = 1 )

		_line2 = dlg_fieldlinetrace ( bS, _pt, direction = -1, ds = _ds, nS = _nF, $ 
			B_FieldLine_CYL = _line_b2, perp = 1 )

		fLine = [[reverse(_line1[*,0:-2],2)],[_line2[*,1:-2]]]
		bLine = [[reverse(_line_b1[*,0:-2],2)],[_line_b2[*,1:-2]]]

		_i = (fLine[0,*]-min(_r))/(max(_r)-min(_r))*(n_elements(_r)-1)
		_j = (fLine[2,*]-min(_z))/(max(_z)-min(_z))*(n_elements(_z)-1)
		_ii = complex(0,1)
		__er = interpolate(_er,_i,_j,cubic=-0.5)*exp(-_ii*nPhi*fLine[1,*])
		this_E = hanning(n_elements(__er))*__er
		;this_E = __er
		_nPer = n_elements(this_E)
		__pPer[n,*] = abs(fft(this_E,/center))^2
		__kPer = fIndGen(_nPer)*2*!pi/(_nPer*_dS)-_nPer*!pi/(_nPer*_dS)

		print, n	
		;;if n eq 50 then stop

		;for ii=0,_nR-1 do begin
		;	for jj=0,_nZ-1 do begin

		;		_kMag = sqrt(k_r[ii]^2+k_t^2+k_z[jj]^2)		
		;		_kPar = k_r[ii]*_bru+k_t*_btu*k_z[jj]*_bru		
		;		_kPer = sqrt(_kMag^2-_kPar^2)

		;		_pMag = p_r[ii]+p_z[jj]		
		;		_pPar = p_r[ii]*_bru+p_t*_btu*p_z[jj]*_bru		
		;		_pPer = sqrt(_pMag^2-_pPar^2)

		;		ikpar = (_kPar-kParMin)/(kParMax-kParMin)*(nkPar-1)
		;		if ikpar gt 0 and ikpar lt nkPar then begin
		;			pPar[ikpar,n] = pPar[ikpar,n]+_pMag
		;			++kParCnt[ikpar]
		;		endif

		;		ikper = (_kPer-kPerMin)/(kPerMax-kPerMin)*(nkPer-1)
		;		if ikper gt 0 and ikper lt nkPer then begin
		;			pPer[ikper,n] = pPer[ikper,n]+_pMag
		;			++kPerCnt[ikper]
		;		endif


		;	endfor
		;endfor			

		nz = where(kParCnt gt 0)
		pPar[nz,n] = pPar[nz,n]/kParCnt[nz]
		nz = where(kPerCnt gt 0)
		pPer[nz,n] = pPer[nz,n]/kPerCnt[nz]

	endfor

	w = 2*!pi*freq
	_c = 3e8
	nLevs = 31
	k2n = _c / w
	n2k = 1.0 / k2n

	__nPar = __kPar * k2n
	__nPer = __kPer * k2n

	;scale = max(ppar)/1.5
	;levels = (fIndGen(nLevs)+1)/(nLevs)*scale
	;colors = bytScl(levels,top=253)+1
	;c=contour(transpose(ppar),g.distance,kpargrid,layout=[2,2,2], $
	;		/current,xtitle='Genray Distance',yTitle='kPar [1/m]', $
	;		c_value=levels, rgb_table=55, rgb_indices=colors, /fill)

    margin=[0.05,0.05,0.05,0.05]*2

	range = 150*k2n
	plotThis = alog(__ppar)
	scale = max(plotThis)/1.0
	levels = (fIndGen(nLevs)+1)/(nLevs)*scale
	colors = bytScl(levels,top=253)+1
	c=contour(plotThis,g.distance*1e-2,__nPar,layout=[2,3,2], $
			/current,xtitle='Distance along GENRAY ray [m]',yTitle='nPar', $
			c_value=levels, rgb_table=55, rgb_indices=colors, /fill,$
			yRange=[-1,1]*range,xRange=[0,2],title=gFileName,margin=margin)
	
	p=plot(g.distance*1e-2,g.nPar,/over,color='b',thick=2,xRange=[0,2])

	ppower=plot(g.distance*1e-2,g.power,/current,layout=[2,3,4],$
            ytitle='Power',xTitle='Distance along GENRAY ray [m]',$
			margin=margin,xRange=[0,2])

	range = 1200*k2n
	plotThis = alog(__pper)
	scale = max(PlotThis)/1.
	levels = (fIndGen(nLevs)+1)/(nLevs)*scale
	colors = bytScl(levels,top=253)+1
	c=contour(PlotThis,g.distance*1e-2,__nPer,layout=[2,3,6], $
			/current,xtitle='Distance along GENRAY ray [m]',yTitle='nPer', $
			c_value=levels, rgb_table=55, rgb_indices=colors, /fill, $
			yRange=[-1,1]*range,margin=margin,xRange=[0,2])
	p=plot(g.distance*1e-2,g.nPer,/over,color='b',thick=2)

    p.save, 'aorsa-genray-comp.eps'
    p.save, 'aorsa-genray-comp.pdf', /bitmap, resolution=100

stop
end
