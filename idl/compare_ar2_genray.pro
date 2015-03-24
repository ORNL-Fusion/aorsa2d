pro compare_ar2_genray

	gFileName1 = 'param_along_ray_n-56.txt'
	gFileName2 = 'param_along_ray_n-71.txt'
	gFileName2 = 'param_along_ray_2.22_0.txt'
	eqdskFileName = 'g200201.00000.dlgMod'
	eqdsk = readgeqdsk(eqdskFileName)

	g1 = read_genray(gFileName1)
	g2 = read_genray(gFileName2)
	g2.z = g2.z-0.08
	__n = n_elements(g2)
	__i = IndGen(__n/2)*2
	g2 = g2[__i]

	flipGenRayParDir = 1
	if flipGenRayParDir then begin
		g1.npar = -g1.npar
		g2.npar = -g2.npar
	endif

	a = ar2_read_solution('./',1)
	ar2_read_ar2input,'ar2Input.nc',ar2=d
	r = ar2_read_rundata('./',1)

	_r = a.r
	_z = a.z

	_er = a.e_r
	_et = a.e_t
	_ez = a.e_z

	nPhi = r.nPhi

	useFredFile = 1	
	if useFredFile then begin
		eLabFrameFileName = 'E_lab_frame'
		eLabFrameFileName = 'E_lab_frame-damp-500a'
		vtkFileName = 'Efield_2D_PARKMURAKAMI.vtk'

		f = read_e_lab_frame(eLabFrameFileName)
		f = read_vtk(vtkFileName)

		_r = f.r
		_z = f.z
		_er = f.er
		_et = f.et
		_ez = f.ez
		nPhi = -71
	endif

	eMag = sqrt(_er^2+_et^2+_ez^2)

	g3 = g1
	for i=0,n_elements(g3)-1 do begin
		ir = (g3[i].r-_r[0])/(_r[-1]-_r[0])*(n_elements(_r)-1)
		tmpE = eMag[ir,*]

		iiPosZ = where(_z gt 0)
		tmpE[iiPosZ] = 0
		iiPosZ = where(_z lt -0.6)
		tmpE[iiPosZ] = 0

		iz = where(tmpE eq max(tmpE) )
		g3[i].z = _z[iz[0]]
	endfor

	ScreenSize = get_screen_size()
	nLevs = 11
	scale = max(eMag)/6
	levels = (fIndGen(nLevs)+1)/(nLevs)*scale
	colors = bytScl(levels,top=253)+1
	c=contour(eMag,_r,_z,aspect_ratio=1.0,layout=[2,1,1], $
			/current,xtitle='R [m]',yTitle='z [m]', $
			c_value=levels, rgb_table=55, rgb_indices=colors, /fill,$
			dimensions=ScreenSize*0.8)
	
	p=plot(d.lim_r,d.lim_z,/over)
	p=plot(eqdsk.rbbbs,eqdsk.zbbbs,/over)
	p=plot(g1.r[*],g1.z[*],/over)
	p=plot(g2.r[*],g2.z[*],/over,color='b')
	p=plot(g3.r[*],g3.z[*],/over,color='g')

	iiKeep_g3 = where(g3.r gt 1.6)
	g3 = g3[iiKeep_g3]
	g = g2
	nGenRay = n_elements(g)

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

	for n=0,nGenRay-1 do begin
		this_r = g[n].r
		this_z = g[n].z
		i = (this_r-d.rmin)/(d.rmax-d.rmin)*(n_elements(d.r)-1)
		j = (this_z-d.zmin)/(d.zmax-d.zmin)*(n_elements(d.z)-1)
		this_br = interpolate(d.br,i,j,cubic=-0.5)
		this_bt = interpolate(d.bt,i,j,cubic=-0.5)
		this_bz = interpolate(d.bz,i,j,cubic=-0.5)

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
		_g = readgeqdsk(eqdskFileName,fieldLineIn=_pt,fieldLine_CYL=_line1,$
			B_AlongFieldLine_CYL=_line_b1,fieldLineTraceDS=_ds,$
			fieldLineTraceNSteps=_nf,fieldLineTraceDir=1)
		_g = readgeqdsk(eqdskFileName,fieldLineIn=_pt,fieldLine_CYL=_line2,$
			B_AlongFieldLine_CYL=_line_b2,fieldLineTraceDS=_ds,$
			fieldLineTraceNSteps=_nf,fieldLineTraceDir=-1)

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

		print, n	
		;;if n eq 50 then stop

		for ii=0,_nR-1 do begin
			for jj=0,_nZ-1 do begin

				_kMag = sqrt(k_r[ii]^2+k_t^2+k_z[jj]^2)		
				_kPar = k_r[ii]*_bru+k_t*_btu*k_z[jj]*_bru		
				_kPer = sqrt(_kMag^2-_kPar^2)

				_pMag = p_r[ii]+p_z[jj]		
				_pPar = p_r[ii]*_bru+p_t*_btu*p_z[jj]*_bru		
				_pPer = sqrt(_pMag^2-_pPar^2)

				ikpar = (_kPar-kParMin)/(kParMax-kParMin)*(nkPar-1)
				if ikpar gt 0 and ikpar lt nkPar then begin
					pPar[ikpar,n] = pPar[ikpar,n]+_pMag
					++kParCnt[ikpar]
				endif

				ikper = (_kPer-kPerMin)/(kPerMax-kPerMin)*(nkPer-1)
				if ikper gt 0 and ikper lt nkPer then begin
					pPer[ikper,n] = pPer[ikper,n]+_pMag
					++kPerCnt[ikper]
				endif


			endfor
		endfor			

		nz = where(kParCnt gt 0)
		pPar[nz,n] = pPar[nz,n]/kParCnt[nz]
		nz = where(kPerCnt gt 0)
		pPer[nz,n] = pPer[nz,n]/kPerCnt[nz]

	endfor

	w = 2*!pi*r.freq
	_c = 3e8
	nLevs = 31

	;scale = max(ppar)/1.5
	;levels = (fIndGen(nLevs)+1)/(nLevs)*scale
	;colors = bytScl(levels,top=253)+1
	;c=contour(transpose(ppar),g.distance,kpargrid,layout=[2,2,2], $
	;		/current,xtitle='Genray Distance',yTitle='kPar [1/m]', $
	;		c_value=levels, rgb_table=55, rgb_indices=colors, /fill)
	plotThis = alog10(__ppar)
	scale = max(plotThis)/1.0
	levels = (fIndGen(nLevs)+1)/(nLevs)*scale
	colors = bytScl(levels,top=253)+1
	c=contour(plotThis,g.distance,__kpar,layout=[2,2,2], $
			/current,xtitle='Genray Distance',yTitle='kPar [1/m]', $
			c_value=levels, rgb_table=55, rgb_indices=colors, /fill,$
			yRange=[-150,150])
	
	p=plot(g.distance,g.nPar*w/_c,/over)

	ppower=plot(g.distance,g.power,/current,layout=[2,5,6])

	PlotThis = pper*transpose(rebin(g.distance,nGenRay,n_elements(pPer[*,0])))
	scale = max(PlotThis)/2.
	levels = (fIndGen(nLevs)+1)/(nLevs)*scale
	colors = bytScl(levels,top=253)+1
	c=contour(transpose(PlotThis),g.distance,kpergrid,layout=[2,2,4], $
			/current,xtitle='Genray Distance',yTitle='kPer [1/m]', $
			c_value=levels, rgb_table=55, rgb_indices=colors, /fill)
	p=plot(g.distance,g.nPer*w/_c,/over)

stop
end
