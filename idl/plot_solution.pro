pro contour_field, field, x, y, nLevs, scale

	loadct, 1
	levels	= (fIndGen(nLevs)+1)/(nLevs-1) * scale * 1.1
	colors	= 256-(bytScl ( levels, top = 253 )+1)
	contour, field,x, y, levels = levels, c_colors=colors, color = 0, /fill
	loadct, 3
	levels	= (fIndGen(nLevs)+1)/(nLevs-1) * scale * 1.1
	colors	= 256-(bytScl ( levels, top = 253 )+1)
	contour, -field,x, y, levels = levels, c_colors=colors, /over, /fill

end

pro plot_solution

	cdfId = ncdf_open ( 'runData.nc', /noWrite ) 

		nCdf_varGet, cdfId, 'capR', x 
		nCdf_varGet, cdfId, 'y', y 
		nCdf_varGet, cdfId, 'xjy_re', jy_re 
		nCdf_varGet, cdfId, 'xjy_im', jy_im 
		nCdf_varGet, cdfId, 'bmod', bmod 
		nCdf_varGet, cdfId, 'xkxsav', kx 
		nCdf_varGet, cdfId, 'xkysav', ky 

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

		ealphak	= complex ( ealphak_re, ealphak_im )
		ebetak	= complex ( ebetak_re, ebetak_im )
		ebk	= complex ( ebk_re, eb_im )

	ncdf_close, cdfId
stop
	window, 0
	!p.multi = [0,3,2]
	!p.background = 255
	scale = max ( abs ( [ealpha[*],ebeta[*],eb[*]] ) ) * 0.4 
	nLevs	= 1001 
	device, decomposed = 0

	contour_field, ealpha, x, y, nLevs, scale
	contour_field, ebeta, x, y, nLevs, scale
	contour_field, eb, x, y, nLevs, scale

	contour_field, imaginary(ealpha), x, y, nLevs, scale
	contour_field, imaginary(ebeta), x, y, nLevs, scale
	contour_field, imaginary(eb), x, y, nLevs, scale

	window, 1
	!p.multi = [0,3,2]
	scale = max ( abs ( [ealphak_re[*],ebetak_re[*],ebk_re[*]] ) ) * 1.5
	scalePar = max ( abs ( [ebk_re[*]] ) ) 
	nLevs	= 1001
	contour_field, ealphak,	kx, ky, nLevs, scale
	contour_field, ebetak,	kx, ky, nLevs, scale 
	contour_field, ebk,		kx, ky, nLevs, scale 
	contour_field, imaginary(ealphak),	kx, ky, nLevs, scale 
	contour_field, imaginary(ebetak),	kx, ky, nLevs, scale 
	contour_field, imaginary(ebk),		kx, ky, nLevs, scale 


stop
end
