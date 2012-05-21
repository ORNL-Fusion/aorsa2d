pro ar2_plot_solution_71, NoRunData = NoRunData

	@constants

	; Load data from all grids
	; ------------------------

	fileList = file_search ( 'solution*.nc' )
	fileListData = file_search ( 'runData*.nc' )

	xAll = 0.0
	yAll = 0.0
	nAll = 0
	eAlphaAll = complex (0,0)
	eBetaAll = complex (0,0)
	eBAll = complex (0,0)

	if not keyword_set(NoRunData) then begin

		cdfId = ncdf_open ( fileListData[0], /noWrite ) 
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
		ncdf_close, cdfId

	endif

		cdfId = ncdf_open ( fileList[0], /noWrite ) 

			nCdf_varGet, cdfId, 'ealpha_re', ealpha_re 
			nCdf_varGet, cdfId, 'ebeta_re', ebeta_re 
			nCdf_varGet, cdfId, 'eB_re', eB_re 
			nCdf_varGet, cdfId, 'ealpha_im', ealpha_im 
			nCdf_varGet, cdfId, 'ebeta_im', ebeta_im 
			nCdf_varGet, cdfId, 'eB_im', eB_im 

			nCdf_varGet, cdfId, 'ealphak_re', ealphak_re 
			nCdf_varGet, cdfId, 'ebetak_re', ebetak_re 
			nCdf_varGet, cdfId, 'eBk_re', eBk_re 
			nCdf_varGet, cdfId, 'ealphak_im', ealphak_im 
			nCdf_varGet, cdfId, 'ebetak_im', ebetak_im 
			nCdf_varGet, cdfId, 'eBk_im', eBk_im 

			nCdf_varGet, cdfId, 'er_re', er_re 
			nCdf_varGet, cdfId, 'et_re', et_re 
			nCdf_varGet, cdfId, 'ez_re', ez_re 
			nCdf_varGet, cdfId, 'er_im', er_im 
			nCdf_varGet, cdfId, 'et_im', et_im 
			nCdf_varGet, cdfId, 'ez_im', ez_im 

			nCdf_varGet, cdfId, 'jalpha_re', jalpha_re 
			nCdf_varGet, cdfId, 'jbeta_re', jbeta_re 
			nCdf_varGet, cdfId, 'jB_re', jB_re 
			nCdf_varGet, cdfId, 'jalpha_im', jalpha_im 
			nCdf_varGet, cdfId, 'jbeta_im', jbeta_im 
			nCdf_varGet, cdfId, 'jB_im', jB_im 

			nCdf_varGet, cdfId, 'jP_r_re', jPr_re 
			nCdf_varGet, cdfId, 'jP_t_re', jPt_re 
			nCdf_varGet, cdfId, 'jP_z_re', jPz_re 
			nCdf_varGet, cdfId, 'jP_r_im', jPr_im 
			nCdf_varGet, cdfId, 'jP_t_im', jPt_im 
			nCdf_varGet, cdfId, 'jP_z_im', jPz_im 

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

	
	nR = n_elements ( e_r[*,0] )
	nZ = n_elements ( e_r[0,*] )

	rMin = 3.5
	rMax = 9.0
	zMin = -5.0
	zMax = 5.0

	r = fIndGen(nR)/(nR-1)*(rMax-rMin)+rMin
	z = fIndGen(nZ)/(nZ-1)*(zMax-zMin)+zMin

	g = readgeqdsk('Scen4_bn2.57_129x129.dlgMod',/noTor)

	nLevs = 11
	scale = 0.05e4 
	levels = fIndGen(nLevs)/(nLevs-1)*scale
	colors = reverse(bytScl(levels, top=253)+1)
	PlotField = (real_part(eAlpha)<max(levels))>min(-levels)
	device, decomposed = 0
	window, 0, xSize=400, ySize=700
	!p.background = 255
	loadct, 3
	contour, PlotField, r, z, levels=levels,c_colors=colors,/fill,/iso,color=0
	loadct, 1
	contour, -PlotField, r, z, levels=levels,c_colors=colors,/fill,/iso,/over,color=0
	loadct, 0
	oplot, g.rlim, g.zlim, color=0,thick=2
	oplot, g.rbbbs, g.zbbbs, color=0,thick=1

	outfname = 'ar2_solution.eps'
	old_dev = !D.name
	set_plot, 'ps'
	device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
		xsize=9, ysize=16,xoffset=.1, yoffset=.1, /encapsul
	loadct, 3
	contour, PlotField, r, z, levels=levels,c_colors=colors,/fill,/iso,color=0
	loadct, 1
	contour, -PlotField, r, z, levels=levels,c_colors=colors,/fill,/iso,/over,color=0
	loadct, 0
	oplot, g.rlim, g.zlim, color=0,thick=2
	oplot, g.rbbbs, g.zbbbs, color=0,thick=1
	device, /close
	set_plot, old_dev

	;c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, aspect_ratio=1.0 )
	;c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
	;p = plot ( g.rlim, g.zlim, /over )
	;p = plot ( g.rbbbs, g.zbbbs, /over )
stop
end 
