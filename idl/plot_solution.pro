pro plot_solution

	cdfId = ncdf_open ( 'runData.nc', /noWrite ) 

		nCdf_varGet, cdfId, 'x', x 
		nCdf_varGet, cdfId, 'y', y 
		nCdf_varGet, cdfId, 'jy_re', jy_re 
		nCdf_varGet, cdfId, 'jy_im', jy_im 
		nCdf_varGet, cdfId, 'bmod', bmod 
		nCdf_varGet, cdfId, 'kx', kx 
		nCdf_varGet, cdfId, 'ky', ky 

	ncdf_close, cdfId


	cdfId = ncdf_open ( 'solution.nc', /noWrite ) 

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

		ealpha	= complex ( ealpha_re, ealpha_im )
		ebeta	= complex ( ebeta_re, ebeta_im )
		eb	= complex ( eb_re, eb_im )

	ncdf_close, cdfId

	window, 0
	!p.multi = [0,3,2]
	scale = max ( abs ( [ealpha[*],ebeta[*],eb[*]] ) ) 
	nLevs	=101 
	device, decomposed = 0
	loadct, 13, file = 'davect.tbl'
	levels	= (fIndGen(nLevs)/(nLevs-1)-0.5) * 2.0 * scale * 1.1
	colors	= bytScl ( levels, top = 253 )+1
	contour, ealpha,x, y, levels = levels, c_colors=colors, /fill
	contour, ebeta,x, y, levels = levels, c_colors=colors, /fill
	contour, eb,x, y, levels = levels*1e-2, c_colors=colors, /fill
	contour, imaginary(ealpha),x, y, levels = levels, c_colors=colors, /fill
	contour, imaginary(ebeta),x, y, levels = levels, c_colors=colors, /fill
	contour, imaginary(eb),x, y, levels = levels*1e-2, c_colors=colors, /fill

	window, 1
	!p.multi = [0,3,2]
	scale = max ( abs ( [ealphak_re[*],ebetak_re[*],ebk_re[*]] ) ) 
	scalePar = max ( abs ( [ebk_re[*]] ) ) 
	nLevs	= 101
	device, decomposed = 0
	loadct, 13, file = 'davect.tbl'
	levels	= (fIndGen(nLevs)/(nLevs-1)-0.5) * 2.0 * scale * 1.1
	levelsPar	= (fIndGen(nLevs)/(nLevs-1)-0.5) * 2.0 * scalePar * 1.1
	colors	= bytScl ( levels, top = 253 )+1
	contour, ealphak_re,kx, ky, levels = levels, c_colors=colors, /fill
	contour, ebetak_re,kx, ky, levels = levels, c_colors=colors, /fill
	contour, ebk_re,kx, ky, levels = levelsPar, c_colors=colors, /fill
	contour, ealphak_im,kx, ky, levels = levels, c_colors=colors, /fill
	contour, ebetak_im,kx, ky, levels = levels, c_colors=colors, /fill
	contour, ebk_im,kx, ky, levels = levelsPar, c_colors=colors, /fill


stop
end
