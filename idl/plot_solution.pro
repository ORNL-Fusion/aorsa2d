pro plot_solution

	cdfId = ncdf_open ( 'runData.nc', /noWrite ) 
		nCdf_varGet, cdfId, 'x', x 
		nCdf_varGet, cdfId, 'y', y 
		nCdf_varGet, cdfId, 'jy_re', jy_re 
		nCdf_varGet, cdfId, 'jy_im', jy_im 
		nCdf_varGet, cdfId, 'bmod', bmod 
	ncdf_close, cdfId


	cdfId = ncdf_open ( 'solution.nc', /noWrite ) 
		nCdf_varGet, cdfId, 'ealpha_re', ealpha_re 
		nCdf_varGet, cdfId, 'ebeta_re', ebeta_re 
		nCdf_varGet, cdfId, 'eB_re', eB_re 
		nCdf_varGet, cdfId, 'ealpha_im', ealpha_im 
		nCdf_varGet, cdfId, 'ebeta_im', ebeta_im 
		nCdf_varGet, cdfId, 'eB_im', eB_im 
	ncdf_close, cdfId

	!p.multi = [0,3,2]
	scale = max ( abs ( [ealpha_re[*],ebeta_re[*],eb_re[*]] ) ) 
	nLevs	=21 
	device, decomposed = 0
	loadct, 13, file = 'davect.tbl'
	levels	= (fIndGen(nLevs)-nLevs/2)/nLevs/2 * scale
	colors	= bytScl ( levels, top = 253 )+1
	contour, ealpha_re,x, y, levels = levels, c_colors=colors
	contour, ebeta_re,x, y, levels = levels, c_colors=colors
	contour, eb_re,x, y, levels = levels*1e-2, c_colors=colors
	contour, ealpha_im,x, y, levels = levels, c_colors=colors
	contour, ebeta_im,x, y, levels = levels, c_colors=colors
	contour, eb_im,x, y, levels = levels*1e-2, c_colors=colors


stop
end
