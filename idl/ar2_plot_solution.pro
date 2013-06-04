pro ar2_plot_solution, full = full, $
		scale1 = scale1, scale2 = scale2, scale3_ = scale3_, $
		sav = sav, NoRunData=NoRunData, $
		nPhi = _nPhi, RHS=_RHS, spec = _ThisSPEC

	if keyword_set(_RHS) then ThisRHS = _RHS else ThisRHS = 0
	if keyword_set(_ThisSPEC) then ThisSPEC = _ThisSPEC else ThisSPEC = 0

	ar2_read_namelist, ar2Input = ar2Input

	AR2InputFile = ar2Input['AR2InputFileName']

	ThisGridNo = 1
	GridNoStr = string(ThisGridNo,format='(i3.3)')

	if keyword_set(_nPhi) then ThisNPhi = _nPhi else ThisNPhi = ar2Input['nPhi']

	nPhiStr = string(ThisNPhi,format='(i+4.3)')

	SolutionFile = 'solution'+GridNoStr+nPhiStr+'.nc'
	RunDataFile = 'runData'+GridNoStr+nPhiStr+'.nc'

	print, 'AR2InputFile: ', AR2InputFile
	print, 'SolutionFile: ', SolutionFile
	print, 'RunDataFile: ', RunDataFile

	ar2_read_ar2input, AR2InputFile, $
			rLim=rLim,zLim=zLim,LimMask=LimMask,nPhi=ThisNPhi

	@constants

	xAll = 0.0
	yAll = 0.0
	nAll = 0
	eAlphaAll = complex (0,0)
	eBetaAll = complex (0,0)
	eBAll = complex (0,0)
	
	cdfId = ncdf_open ( RunDataFile, /noWrite ) 
		nCdf_varGet, cdfId, 'nPhi', nPhi 
		nCdf_varGet, cdfId, 'freq', freq 
		nCdf_varGet, cdfId, 'capR', x 
		nCdf_varGet, cdfId, 'y', y 
		nCdf_varGet, cdfId, 'jr_re', jr_re 
		nCdf_varGet, cdfId, 'jr_im', jr_im 
		nCdf_varGet, cdfId, 'jt_re', jt_re 
		nCdf_varGet, cdfId, 'jt_im', jt_im 
		nCdf_varGet, cdfId, 'jz_re', jz_re 
		nCdf_varGet, cdfId, 'jz_im', jz_im 
		nCdf_varGet, cdfId, 'bmod', bmod
		nCdf_varGet, cdfId, 'brU', brU
		nCdf_varGet, cdfId, 'btU', btU
		nCdf_varGet, cdfId, 'bzU', bzU
		nCdf_varGet, cdfId, 'densitySpec', densitySpec
		nCdf_varGet, cdfId, 'tempSpec', tempSpec
		nCdf_varGet, cdfId, 'xx_re', xx_re 
		nCdf_varGet, cdfId, 'xx_im', xx_im
		nCdf_varGet, cdfId, 'yy_re', yy_re 
		nCdf_varGet, cdfId, 'yy_im', yy_im
		nCdf_varGet, cdfId, 'drBfn_re', dRbFn_bFn_re
		nCdf_varGet, cdfId, 'drBfn_im', dRbFn_bFn_im
		nCdf_varGet, cdfId, 'dzBfn_re', dzbFn_bFn_re
		nCdf_varGet, cdfId, 'dzBfn_im', dzbFn_bFn_im
		nCdf_varGet, cdfId, 'nuOmg', nuOmg
	ncdf_close, cdfId

	xx	= complex ( xx_re, xx_im )
	yy	= complex ( yy_re, yy_im )

    nuOmg = nuOmg[*,*,ThisSpec]

	jA_r = complex ( jr_re[*,*,ThisRHS], jr_im[*,*,ThisRHS] )
	jA_t = complex ( jt_re[*,*,ThisRHS], jt_im[*,*,ThisRHS] )
	jA_z = complex ( jz_re[*,*,ThisRHS], jz_im[*,*,ThisRHS] )

	dRbFn_bFn	= complex ( dRbFn_bFn_re, dRbFn_bFn_im )
	dZbFn_bFn	= complex ( dZbFn_bFn_re, dZbFn_bFn_im )

	nN	= n_elements ( xx[*,0] )
	nM	= n_elements ( yy[*,0] )

	cdfId = ncdf_open ( SolutionFile, /noWrite ) 

		nCdf_varGet, cdfId, 'r', r
		nCdf_varGet, cdfId, 'z', z
		nCdf_varGet, cdfId, 'nPhi', nPhi 
		nCdf_varGet, cdfId, 'freqcy', freq 
	
		nCdf_varGet, cdfId, 'ealpha_re', ealpha_re 
		ealpha_re = temporary(ealpha_re[*,*,ThisRHS])
		nCdf_varGet, cdfId, 'ebeta_re', ebeta_re 
		ebeta_re = temporary(ebeta_re[*,*,ThisRHS])
		nCdf_varGet, cdfId, 'eB_re', eB_re 
		eb_re = temporary(eb_re[*,*,ThisRHS])

		nCdf_varGet, cdfId, 'ealpha_im', ealpha_im 
		ealpha_im = temporary(ealpha_im[*,*,ThisRHS])
		nCdf_varGet, cdfId, 'ebeta_im', ebeta_im 
		ebeta_im = temporary(ebeta_im[*,*,ThisRHS])
		nCdf_varGet, cdfId, 'eB_im', eB_im 
		eb_im = temporary(eb_im[*,*,ThisRHS])

		nCdf_varGet, cdfId, 'ealphak_re', ealphak_re 
		ealphak_re = temporary(ealphak_re[*,*,ThisRHS])
		nCdf_varGet, cdfId, 'ebetak_re', ebetak_re 
		ebetak_re = temporary(ebetak_re[*,*,ThisRHS])
		nCdf_varGet, cdfId, 'eBk_re', eBk_re 
		ebk_re = temporary(ebk_re[*,*,ThisRHS])

		nCdf_varGet, cdfId, 'ealphak_im', ealphak_im 
		ealphak_im = temporary(ealphak_im[*,*,ThisRHS])
		nCdf_varGet, cdfId, 'ebetak_im', ebetak_im 
		ebetak_im = temporary(ebetak_im[*,*,ThisRHS])
		nCdf_varGet, cdfId, 'eBk_im', eBk_im 
		ebk_im = temporary(ebk_im[*,*,ThisRHS])

		nCdf_varGet, cdfId, 'er_re', er_re 
		er_re = temporary(er_re[*,*,ThisRHS])
		nCdf_varGet, cdfId, 'et_re', et_re 
		et_re = temporary(et_re[*,*,ThisRHS])
		nCdf_varGet, cdfId, 'ez_re', ez_re 
		ez_re = temporary(ez_re[*,*,ThisRHS])

		nCdf_varGet, cdfId, 'er_im', er_im 
		er_im = temporary(er_im[*,*,ThisRHS])
		nCdf_varGet, cdfId, 'et_im', et_im 
		et_im = temporary(et_im[*,*,ThisRHS])
		nCdf_varGet, cdfId, 'ez_im', ez_im 
		ez_im = temporary(ez_im[*,*,ThisRHS])

		nCdf_varGet, cdfId, 'jalpha_re', jalpha_re 
		jalpha_re = temporary(jalpha_re[*,*,ThisSPEC,ThisRHS])
		nCdf_varGet, cdfId, 'jbeta_re', jbeta_re 
		jbeta_re = temporary(jbeta_re[*,*,ThisSPEC,ThisRHS])
		nCdf_varGet, cdfId, 'jB_re', jB_re 
		jb_re = temporary(jb_re[*,*,ThisSPEC,ThisRHS])

		nCdf_varGet, cdfId, 'jalpha_im', jalpha_im 
		jalpha_im = temporary(jalpha_im[*,*,ThisSPEC,ThisRHS])
		nCdf_varGet, cdfId, 'jbeta_im', jbeta_im 
		jbeta_im = temporary(jbeta_im[*,*,ThisSPEC,ThisRHS])
		nCdf_varGet, cdfId, 'jB_im', jB_im 
		jb_im = temporary(jb_im[*,*,ThisSPEC,ThisRHS])

		nCdf_varGet, cdfId, 'jP_r_re', jPr_re 
		jpr_re = temporary(jpr_re[*,*,ThisSpec,ThisRHS])
		nCdf_varGet, cdfId, 'jP_t_re', jPt_re 
		jpt_re = temporary(jpt_re[*,*,ThisSpec,ThisRHS])
		nCdf_varGet, cdfId, 'jP_z_re', jPz_re 
		jpz_re = temporary(jpz_re[*,*,ThisSpec,ThisRHS])

		nCdf_varGet, cdfId, 'jP_r_im', jPr_im 
		jpr_im = temporary(jpr_im[*,*,ThisSpec,ThisRHS])
		nCdf_varGet, cdfId, 'jP_t_im', jPt_im 
		jpt_im = temporary(jpt_im[*,*,ThisSPEC,ThisRHS])
		nCdf_varGet, cdfId, 'jP_z_im', jPz_im 
		jpz_im = temporary(jpz_im[*,*,ThisSPEC,ThisRHS])

		nCdf_varGet, cdfId, 'jouleHeating', jouleHeating 

		ealpha	= complex ( ealpha_re, ealpha_im )
		ebeta	= complex ( ebeta_re, ebeta_im )
		eb	= complex ( eb_re, eb_im )

		ealphak	= complex ( ealphak_re, ealphak_im )
		ebetak	= complex ( ebetak_re, ebetak_im )
		ebk	= complex ( ebk_re, ebk_im )

		e_r	= complex ( er_re, er_im )
		e_t	= complex ( et_re, et_im )
		e_z	= complex ( ez_re, ez_im )

		jP_r	= complex ( jPr_re, jPr_im )
		jP_t	= complex ( jPt_re, jPt_im )
		jP_z	= complex ( jPz_re, jPz_im )

		jPAlpha = complex ( jAlpha_re, jAlpha_im )
		jPBeta = complex ( jBeta_re, jBeta_im )
		jPB = complex ( jB_re, jB_im )

	ncdf_close, cdfId

	x = r
	y = z
	nx = n_elements(r)
	ny = n_elements(z)
	x2D	= rebin ( r, nX, nY )
	if nY gt 1 then begin
		y2D = transpose(rebin ( z, nY, nX ))
	endif else begin
		y2D = fltArr(1)+z
	endelse

	; 1D catch
	; --------

	@load_colors

	; Plot all grid solutions
	; -----------------------

	if size(y,/dim) eq 0 then begin

		iix = n_elements(x)/2
		wrf = freq * 2 * !pi

		eRange = max(abs([e_r,e_t,e_z]))
		p_r = plot ( x, e_r, layout=[1,3,1],$
				title='Er',yRange=[-eRange,eRange],ytitle='Er [V/m]',name='Re',window_title='aorsa',font_size=16)
		p_i = plot ( x, imaginary(e_r), color='red',/over,name='Im')
		;l = legend(target=[p_r,p_i],position=[0.98,0.95],/norm,font_size=10,horizontal_alignment='RIGHT')

		p_r = plot ( x, e_t, layout=[1,3,2],/current,$
				title='Et',yRange=[-eRange,eRange],ytitle='Et [V/m]',name='Re')
		p_i = plot ( x, imaginary(e_t), color='red',/over,name='Im')
		;l = legend(target=[p_r,p_i],position=[0.98,0.62],/norm,font_size=10,horizontal_alignment='RIGHT')

		p_r = plot ( x, e_z, layout=[1,3,3],/current,$
				title='Ez',yRange=[-eRange,eRange],ytitle='Ez [V/m]',name='Re')
		p_i = plot ( x, imaginary(e_z), color='red',/over,name='Im')
		;l = legend(target=[p_r,p_i],position=[0.98,0.28],/norm,font_size=10,horizontal_alignment='RIGHT')

		nSpec	= n_elements ( jp_r[0,0,*] )

		; jDotE
		; -----

		jDotE	= jp_r[*,*,0] * conj(e_r) $
					+ jp_t[*,*,0] * conj(e_t) $
					+ jp_z[*,*,0] * conj(e_z)

		p = plot (x,jDotE,color='b',thick=3,transparency=50,$
				title='J dot E',name='jDote_0',font_size=10,$
				layout=[1,2,1],window_title='aorsa')

		p_array = !NULL
		p_array = [p_array,p]

		for s=1,nSpec-1 do begin

			jDotE	= jp_r[*,*,s] * conj(e_r) $
					+ jp_t[*,*,s] * conj(e_t) $
					+ jp_z[*,*,s] * conj(e_z)

			p = plot ( x, jDotE, /over,name='jDotE_'+strTrim(string(s),2) )
			p_array = [p_array,p]

		endfor

		; jP
		; --

		jpRange = max(abs([abs(jp_r),abs(jp_t),abs(jp_z)]))

		s = 0
		;p_array = !NULL
		p = plot (x,jp_r[*,0,s],thick=2,transparency=50,$
				title='jPr',name='jPr_re_'+strTrim(string(s),2),font_size=10,$
				layout=[1,3,1],yRange=[-jpRange,jpRange],window_title='aorsa')
		;p_array = [p_array,p]
		p = plot (x,imaginary(jp_r[*,0,s]),/over,name='jPr_im_'+strTrim(string(s),2),color='r',thick=2)
		;p_array = [p_array,p]
		for s=1,nSpec-1 do begin
			p = plot ( x, jp_r[*,0,s],/over,name='jPr_re_'+strTrim(string(s),2),thick=2,lineStyle=s)
			;p_array = [p_array,p]
			p = plot ( x, imaginary(jp_r[*,0,s]),/over,name='jPr_im_'+strTrim(string(s),2),thick=2,lineStyle=s,color='r')
			;p_array = [p_array,p]
		endfor
	   	;l = legend(target=p_array,position=[0.99,0.9],/norm,font_size=10,horizontal_alignment='RIGHT')

		s = 0
		;p_array = !null
		p = plot (x,jp_t[*,0,s],thick=2,transparency=50,$
				title='jPt',name='jPt_re_'+strtrim(string(s),2),font_size=10,$
				layout=[1,3,2],yrange=[-jPrange,jPrange],/current)
		;p_array = [p_array,p]
		p = plot (x,imaginary(jp_t[*,0,s]),/over,name='jPt_im_'+strtrim(string(s),2),color='r',thick=2)
		;p_array = [p_array,p]
		for s=1,nspec-1 do begin
			p = plot ( x, jp_t[*,0,s],/over,name='jPt_re_'+strtrim(string(s),2),thick=2,linestyle=s)
			;p_array = [p_array,p]
			p = plot ( x, imaginary(jp_t[*,0,s]),/over,name='jPt_im_'+strtrim(string(s),2),thick=2,linestyle=s,color='r')
			;p_array = [p_array,p]
		endfor
	   	;l = legend(target=p_array,position=[0.99,0.6],/norm,font_size=10,horizontal_alignment='right')

		s = 0
		;p_array = !null
		p = plot (x,jP_z[*,0,s],thick=2,transparency=50,$
				title='jPz',name='jPz_re_'+strtrim(string(s),2),font_size=10,$
				layout=[1,3,3],yrange=[-jPrange,jPrange],/current)
		;p_array = [p_array,p]
		p = plot (x,imaginary(jP_z[*,0,s]),/over,name='jPz_im_'+strtrim(string(s),2),color='r',thick=2)
		;p_array = [p_array,p]
		for s=1,nspec-1 do begin
			p = plot ( x, jP_z[*,0,s],/over,name='jPz_re_'+strtrim(string(s),2),thick=2,linestyle=s)
			;p_array = [p_array,p]
			p = plot ( x, imaginary(jP_z[*,0,s]),/over,name='jPz_im_'+strtrim(string(s),2),thick=2,linestyle=s,color='r')
			;p_array = [p_array,p]
		endfor
	   	;l = legend(target=p_array,position=[0.99,0.25],/norm,font_size=10,horizontal_alignment='right')

	endif else begin ; Now 2D plotting

        dimensions = [500,600]
        scaleFac = 0.05
        nPhiString = ' (nPhi: '+string(ThisNPHI,format='(i+4.3)')+')'

		scale = max(abs([jP_r,jP_t,jP_z]))*scaleFac
        scale = 0.16
        print, 'jP scale: ', scale

		nLevs = 11
		thisField = jP_r[*,*]*LimMask
        title = 'jP_r'+ nPhiString
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		c = contour ( thisField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
                aspect_ratio=1.0, title=title, layout=[1,3,1], dimensions=dimensions )
		c = contour ( -thisField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
		p = plot ( rLim, zLim, /over )

		thisField = jP_t[*,*]*LimMask
        title = 'jP_t'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		c = contour ( thisField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
                aspect_ratio=1.0, title=title, layout=[1,3,2], /current )
		c = contour ( -thisField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
		p = plot ( rLim, zLim, /over )

		thisField = jP_z[*,*]*LimMask
        title = 'jP_z'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		c = contour ( thisField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
                aspect_ratio=1.0, title=title, layout=[1,3,3], /current )
		c = contour ( -thisField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
		p = plot ( rLim, zLim, /over )

	    dimensions = [800,600]
        scaleFac = 0.3
	    scale = max(abs([e_r,e_t,e_z]))*scaleFac 

		thisField = e_r[*,*]
        title = 'E_r'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
            aspect_ratio=1.0, layout=[2,3,1], dimensions=dimensions, title=title )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
		p = plot ( rlim, zlim, /over )

		thisField = abs(e_r[*,*])
        title = 'abs(E_r)'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
            aspect_ratio=1.0, layout=[2,3,2], /current, title=title )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
		p = plot ( rlim, zlim, /over )


		thisField = e_t[*,*]
        title = 'E_t'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
                aspect_ratio=1.0, layout=[2,3,3], /current, title=title )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
		p = plot ( rlim, zlim, /over )

		thisField = abs(e_t[*,*])
        title = 'abs(E_t)'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
                aspect_ratio=1.0, layout=[2,3,4], /current, title=title )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
		p = plot ( rlim, zlim, /over )


		thisField = e_z[*,*]
        title = 'E_z'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
                aspect_ratio=1.0, layout=[2,3,5], /current, title=title )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
		p = plot ( rlim, zlim, /over )

		thisField = abs(e_z[*,*])
        title = 'abs(E_z)'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
                aspect_ratio=1.0, layout=[2,3,6], /current, title=title )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
		p = plot ( rlim, zlim, /over )



		scale = max(abs([jA_r,jA_t,jA_z]))*scaleFac

		thisField = jA_r[*,*]
        title = 'jA_r'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
                aspect_ratio=1.0, title=title, layout=[1,3,1], dimensions=dimensions )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
		p = plot ( rlim, zlim, /over )

		thisField = jA_t[*,*]
        title = 'jA_t'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
                aspect_ratio=1.0, title=title, layout=[1,3,2], /current )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
		p = plot ( rlim, zlim, /over )

		thisField = jA_z[*,*]
        title = 'jA_z'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
                aspect_ratio=1.0, title=title, layout=[1,3,3], /current )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
		p = plot ( rlim, zlim, /over )

		scale = max(abs([nuOmg]))*scaleFac

		thisField = nuOmg[*,*]
        title = 'nuOmg'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
                aspect_ratio=1.0, title=title, layout=[1,3,1], dimensions=dimensions )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
		p = plot ( rlim, zlim, /over )


	endelse

	if keyword_set(sav) then begin

		xorig = x
		eAlphaOrig = eAlpha
		eBetaOrig = eBeta
		eBOrig = eB

		save, xorig, eAlphaOrig, eBetaOrig, eBOrig, $
				file = 'soln.sav'
	endif

	stop
end
