pro plot_rundata

	cdfId = ncdf_open ( 'runData.nc', /noWrite ) 
		nCdf_varGet, cdfId, 'capR', x 
		nCdf_varGet, cdfId, 'y', y 
		nCdf_varGet, cdfId, 'xjy_re', jy_re 
		nCdf_varGet, cdfId, 'xjy_im', jy_im 
		nCdf_varGet, cdfId, 'bmod', bmod
		nCdf_varGet, cdfId, 'bxn', bxn
		nCdf_varGet, cdfId, 'byn', byn
		nCdf_varGet, cdfId, 'bzn', bzn
		nCdf_varGet, cdfId, 'densitySpec', densitySpec
		nCdf_varGet, cdfId, 'tempSpec', tempSpec
		nCdf_varGet, cdfId, 'tempSpec', tempSpec
		nCdf_varGet, cdfId, 'xx_re', xx_re 
		nCdf_varGet, cdfId, 'xx_im', xx_im
		nCdf_varGet, cdfId, 'yy_re', yy_re 
		nCdf_varGet, cdfId, 'yy_im', yy_im
		nCdf_varGet, cdfId, 'kxsav',  kxsav
		nCdf_varGet, cdfId, 'kysav', kysav 
	ncdf_close, cdfId

	xx	= complex ( xx_re, xx_im )
	yy	= complex ( yy_re, yy_im )

	bx	= bxn * bmod
	by	= byn * bmod
	bz	= bzn * bmod

	window, 0, xSize = 800, ySize = 800
	!p.multi = [0,3,3]
	contour, jy_im, x, y
	contour, bmod, x, y, nLev = 30
	contour, bx, x, y
	contour, by, x, y
	contour, bz, x, y
	contour, densitySpec[*,*,0], x, y
	contour, tempSpec[*,*,0], x, y

end
