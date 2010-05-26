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
		nCdf_varGet, cdfId, 'omgc', omgc 
		nCdf_varGet, cdfId, 'omgp2', omgp2
	
		nCdf_varGet, cdfId, 'xx_re', xx_re 
		nCdf_varGet, cdfId, 'xx_im', xx_im
		nCdf_varGet, cdfId, 'yy_re', yy_re 
		nCdf_varGet, cdfId, 'yy_im', yy_im

		nCdf_varGet, cdfId, 'drUrr', drUrr
		nCdf_varGet, cdfId, 'drUrt', drUrt
		nCdf_varGet, cdfId, 'drUrz', drUrz

		nCdf_varGet, cdfId, 'drUtr', drUtr
		nCdf_varGet, cdfId, 'drUtt', drUtt
		nCdf_varGet, cdfId, 'drUtz', drUtz

		nCdf_varGet, cdfId, 'drUzr', drUzr
		nCdf_varGet, cdfId, 'drUzt', drUzt
		nCdf_varGet, cdfId, 'drUzz', drUzz

		nCdf_varGet, cdfId, 'dzUrr', dzUrr
		nCdf_varGet, cdfId, 'dzUrt', dzUrt
		nCdf_varGet, cdfId, 'dzUrz', dzUrz

		nCdf_varGet, cdfId, 'dzUtr', dzUtr
		nCdf_varGet, cdfId, 'dzUtt', dzUtt
		nCdf_varGet, cdfId, 'dzUtz', dzUtz

		nCdf_varGet, cdfId, 'dzUzr', dzUzr
		nCdf_varGet, cdfId, 'dzUzt', dzUzt
		nCdf_varGet, cdfId, 'dzUzz', dzUzz

		nCdf_varGet, cdfId, 'Urr', Urr
		nCdf_varGet, cdfId, 'Urt', Urt
		nCdf_varGet, cdfId, 'Urz', Urz

		nCdf_varGet, cdfId, 'Utr', Utr
		nCdf_varGet, cdfId, 'Utt', Utt
		nCdf_varGet, cdfId, 'Utz', Utz

		nCdf_varGet, cdfId, 'Uzr', Uzr
		nCdf_varGet, cdfId, 'Uzt', Uzt
		nCdf_varGet, cdfId, 'Uzz', Uzz

	ncdf_close, cdfId

	xx	= complex ( xx_re, xx_im )
	yy	= complex ( yy_re, yy_im )

	bx	= bxn * bmod
	by	= byn * bmod
	bz	= bzn * bmod

	nR	= n_elements ( x )
	nz	= n_elements ( y )
	nSpec	= n_elements ( densitySpec[0,0,*] )

	iPlot, x, bMod[*,nz/2.0], $
			view_grid = [2,2], $
			title = 'z=0 equil B field'
	iPlot, x, bx[*,nz/2.0], /over
	iPlot, x, by[*,nz/2.0], /over
	iPlot, x, bz[*,nz/2.0], /over

	iPlot, x, densitySpec[*,nz/2.0,0], $
			/view_next, $
			title = 'z=0 Density'
	for i=1,nSpec-1 do $
		iPlot, x, densitySpec[*,nz/2,i], /over

	iPlot, x, tempSpec[*,nz/2,0], $
			/view_next, $
			title = 'z=0 Temp'
	for i=1,nSpec-1 do $
		iPlot, x, tempSpec[*,nz/2,i], /over

	if (nR gt 1) and (nZ gt 1) then begin

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

	endif
stop

end
