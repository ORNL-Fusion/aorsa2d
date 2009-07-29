pro plot_aorsa

	eqDskFileName	= 'g123435.00400'
	eqdsk	= readGEQDSK ( eqDskFileName )
	cdfId = ncdf_open ( 'output/plotData.nc', /noWrite ) 
	
	ncdf_varGet, cdfId, 'rho', rho 
	ncdf_varGet, cdfId, 'wdote', wdote
	ncdf_varget, cdfId, 'wdoti1', wdoti1 
	ncdf_varget, cdfId, 'wdoti2', wdoti2 
	ncdf_varget, cdfId, 'capR', capR 
	ncdf_varget, cdfId, 'zLoc', zLoc 
	ncdf_varGet, cdfId, 'wdote_rz', wdote_rz
	ncdf_varget, cdfId, 'wdoti1_rz', wdoti1_rz
	ncdf_varget, cdfId, 'wdoti2_rz', wdoti2_rz
	nCdf_varGet, cdfId, 'pscale', pscale
	nCdf_varGet, cdfId, 'ePlus_real', ePlus_real
	nCdf_varGet, cdfId, 'ePlus_imag', ePlus_imag
	nCdf_varGet, cdfId, 'eMinu_real', eMinu_real
	nCdf_varGet, cdfId, 'eMinu_imag', eMinu_imag
	nCdf_varGet, cdfId, 'bx_wave_real', bx_wave_real
	nCdf_varGet, cdfId, 'bx_wave_imag', bx_wave_imag
	nCdf_varGet, cdfId, 'bz_wave_real', bz_wave_real
	nCdf_varGet, cdfId, 'bz_wave_imag', bz_wave_imag

	ncdf_close, cdfId


	loadct, 13, file = 'davect.tbl', /silent
	device, decomposed = 0

	window, 0, xSize = 1200, ySize = 600
	!p.multi = [0,3,2]
	!p.charSize = 2.6

	nLevs	= 21

	range	= max ( ePlus_real ) / 2.0
	brange	= max ( bx_wave_real ) / 2.0

	levels	= ( fIndGen ( nLevs ) - nLevs / 2.0 ) / ( nLevs / 2.0 ) * range 
	colors	= bytScl ( levels, top = 253 ) + 1

	blevels	= ( fIndGen ( nLevs ) - nLevs / 2.0 ) / ( nLevs / 2.0 ) * brange 
	bcolors	= bytScl ( blevels, top = 253 ) + 1

	mod_nLevs	= 11
	mod_levels	= fIndGen ( mod_nLevs ) / mod_nLevs * range 
	mod_colors	= 255 - ( bytScl ( mod_levels, top = 253 ) + 1 )

	mod_blevels	= fIndGen ( mod_nLevs ) / mod_nLevs * brange 
	mod_bcolors	= 255 - ( bytScl ( mod_blevels, top = 253 ) + 1 )

	contour, (ePlus_real<range)>(-range), capR, zLoc, $
			color = 255, $
			levels = levels, $
			c_colors = colors, $
			title = 'ePlus_real', $
			/fill
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 255 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 255 

   	contour, (ePlus_imag<range)>(-range), capR, zLoc, $
			color = 255, $
			levels = levels, $
			c_colors = colors, $
			title = 'ePlus_imag', $
			/fill
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 255 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 255 

   	loadct, 3, /silent
   	contour, sqrt ( ePlus_imag^2 + ePlus_real^2 ), capR, zLoc, $
			color = 0, $
			levels = mod_levels, $
			c_colors = mod_colors, $
			title = 'mod ( ePlus )', $
			/fill
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 0 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 0 

   	loadct, 13, file = 'davect.tbl', /silent
    contour, (eMinu_real<range)>(-range), capR, zLoc, $
			color = 255, $
			levels = levels, $
			c_colors = colors, $
			title = 'ePlus_real', $
			/fill
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 255 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 255 

   	contour, (eMinu_imag<range)>(-range), capR, zLoc, $
			color = 255, $
			levels = levels, $
			c_colors = colors, $
			title = 'ePlus_imag', $
			/fill
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 255 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 255 

   	loadct, 3, /silent
   	contour, sqrt ( eMinu_imag^2 + eMinu_real^2 ), capR, zLoc, $
			color = 0, $
			levels = mod_levels, $
			c_colors = mod_colors, $
			title = 'mod ( eMinu )', $
			/fill
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 0 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 0 

   	loadct, 13, file = 'davect.tbl', /silent
   	window, 1, xSize = 1200, ySize = 600

	contour, (bx_wave_real<brange)>(-brange), capR, zLoc, $
			color = 255, $
			levels = blevels, $
			c_colors = bcolors, $
			title = 'bx_wave_real', $
			/fill
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 255 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 255 

   	contour, (bx_wave_imag<brange)>(-brange), capR, zLoc, $
			color = 255, $
			levels = blevels, $
			c_colors = bcolors, $
			title = 'bx_wave_imag', $
			/fill
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 255 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 255 

   	loadct, 3, /silent
   	contour, (sqrt(bx_wave_imag^2+bx_wave_real^2)<brange)>(-brange), capR, zLoc, $
			color = 0, $
			levels = mod_blevels, $
			c_colors = mod_bcolors, $
			title = 'mod ( bx_wave )', $
			/fill
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 0 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 0 

   	loadct, 13, file = 'davect.tbl', /silent
   	contour, (bz_wave_real<brange)>(-brange), capR, zLoc, $
			color = 255, $
			levels = blevels, $
			c_colors = bcolors, $
			title = 'bz_wave_real', $
			/fill
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 255 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 255 

   	contour, (bz_wave_imag<brange)>(-brange), capR, zLoc, $
			color = 255, $
			levels = blevels, $
			c_colors = bcolors, $
			title = 'bz_wave_imag', $
			/fill
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 255 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 255 

   	loadct, 3, /silent
   	contour, (sqrt(bz_wave_imag^2+bz_wave_real^2)<brange)>(-brange), capR, zLoc, $
			color = 0, $
			levels = mod_blevels, $
			c_colors = mod_bcolors, $
			title = 'mod ( bz_wave )', $
			/fill
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 0 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 0 

stop

end 
