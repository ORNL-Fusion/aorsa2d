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

	nR	= n_elements ( capR )
	nz	= n_elements ( y )
	nSpec	= n_elements ( densitySpec[0,0,*] )

	iPlot, x, bMod[*,nR/2], $
			view_grid = [2,2], $
			title = 'z=0 equil B field'
	iPlot, x, bxn[*,nR/2]*bMod, /over
	iPlot, x, byn[*,nR/2]*bMod, /over
	iPlot, x, bzn[*,nR/2]*bMod, /over

	iPlot, x, densitySpec[*,nz/2,0], $
			/view_next, $
			title = 'z=0 Density'
	for i=1,nSpec-1 do $
		iPlot, x, densitySpec[*,nz/2,i], /over

	iPlot, x, tempSpec[*,nz/2,0], $
			/view_next, $
			title = 'z=0 Temp'
	for i=1,nSpec-1 do $
		iPlot, x, tempSpec[*,nz/2,i], /over


	bxn_	= conGrid ( bxn, 20, 20, /center )
	byn_	= conGrid ( byn, 20, 20, /center )
	x_		= conGrid ( x, 20, 20, /center )
	y_		= conGrid ( y, 20, 20, /center )

	colors	= 256 - ( bytScl ( sqrt ( bxn_^2+byn_^2 ), top = 253 ) + 1 )

	iVector, bxn_, byn_, x_, y_, $
			vector_colors = colors, $
			rgb_table = 1, $
			scale_isotropic = 1, $
			/view_next

stop

end
