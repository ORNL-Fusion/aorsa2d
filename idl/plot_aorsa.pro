pro plot_aorsa, $
		range = range, $
		brange = brange, $
		prange = prange, $
		flat = flat

	;eqDskFileName	= 'g123435.00400'
	;eqDskFileName	= 'eqdsk.122993'
	eqDskFileName = 'eqdsk'

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
	nCdf_varGet, cdfId, 'density', density 
	nCdf_varGet, cdfId, 'mask', mask
	nCdf_varGet, cdfId, 'janty', janty 
	nCdf_varGet, cdfId, 'jantx', jantx
	nCdf_varGet, cdfId, 'jeDotE', jeDotE
	nCdf_varGet, cdfId, 'eAlpha', eAlpha
	nCdf_varGet, cdfId, 'eBeta', eBeta
	nCdf_varGet, cdfId, 'eParallel', eParallel
	nCdf_varGet, cdfId, 'ex', ex
	nCdf_varGet, cdfId, 'ey', ey
	nCdf_varGet, cdfId, 'ez', ez
	ncdf_close, cdfId

	;	create an interpolated limiter boundary

	newR	= (eqdsk.rLim)[0]
	newZ	= (eqdsk.zLim)[0]

	for i=0,n_elements(eqdsk.rLim)-2 do begin

		;	get slope

		m	= ( (eqdsk.zLim)[i+1]-(eqdsk.zLim)[i] ) $
				/ ( (eqdsk.rLim)[i+1] - (eqdsk.rLim)[i] )
		b	= (eqdsk.zLim)[i] - m * (eqdsk.rLim)[i]

		;	distance

		d	= sqrt ( ( (eqdsk.rLim)[i+1] - (eqdsk.rLim)[i] )^2 $
				+ ( (eqdsk.zLim)[i+1] - (eqdsk.zLim)[i] )^2 )

		;dMin	= ( capR[0] - capR[1] ) / 2.0

		;if d gt abs(dMin) then begin

			nExtra	= 10;fix ( d / abs(dMin) )
			dStep	= ((eqdsk.rLim)[i+1] - (eqdsk.rLim)[i]) / nExtra

			for j = 0, nExtra - 1 do begin

				if dStep ne 0 then begin
					newR	= [ newR, (eqdsk.rLim)[i] + dStep*j ]
					newZ	= [ newZ, m * ((eqdsk.rLim)[i] + dStep*j) + b ]
				endif

			endfor

		;endif

	endfor

	;	calculate div jant to test for div . jant = 0

	R2D	= rebin ( capR, n_elements ( capR ), n_elements ( zLoc ) )
	RStep	= capR[1]-capR[0]
	zStep	= zLoc[1] - zLoc[0]
	divJ	= 1.0 / R2D * dlg_pderiv ( R2D * jantx, 1, RStep ) $
			+ dlg_pderiv ( janty, 1, zStep )

	loadct, 13, file = 'davect.tbl', /silent

	set_plot, 'ps'
	device, fileName = 'output/aorsa_dlg.ps', $
			/color, $
			bits_per_pixel = 8

	!p.multi = [0,3,2]
	!p.charSize = 1.0

	nLevs	= 21

	if(not keyword_set(range)) then range = max ( ePlus_real ) / 2.0
	if(not keyword_set(brange)) then brange	= max ( bx_wave_real ) / 2.0
	if(not keyword_set(prange)) then prange	= max ( jedote ) / 5.0

	levels	= ( fIndGen ( nLevs ) - nLevs / 2.0 ) / ( nLevs / 2.0 ) * range 
	colors	= bytScl ( levels, top = 253 ) + 1

	blevels	= ( fIndGen ( nLevs ) - nLevs / 2.0 ) / ( nLevs / 2.0 ) * brange 
	bcolors	= bytScl ( blevels, top = 253 ) + 1

	plevels	= ( fIndGen ( nLevs ) - nLevs / 2.0 ) / ( nLevs / 2.0 ) * prange 
	pcolors	= bytScl ( plevels, top = 253 ) + 1


	xRange	= [ min ( eqdsk.rlim ), max ( eqdsk.rlim ) ]
	yRange	= [ min ( eqdsk.zlim ), max ( eqdsk.zlim ) ]


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
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
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
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
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
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 0 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 0 
   	loadct, 12, /silent
   	contour, janty, capR, zLoc, $
		   c_color = 1*16-1, $
		   /overplot, $ 
   			levels = fIndGen(10)
	
   	loadct, 13, file = 'davect.tbl', /silent
    contour, (eMinu_real<range)>(-range), capR, zLoc, $
			color = 255, $
			levels = levels, $
			c_colors = colors, $
			title = 'eMinu_real', $
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
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
			title = 'eMinu_imag', $
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
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
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 0 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 0 

   	loadct, 13, file = 'davect.tbl', /silent
   	;window, 1, xSize = 1200, ySize = 600

	contour, (ex<range)>(-range), capR, zLoc, $
			color = 255, $
			levels = levels, $
			c_colors = colors, $
			title = 'ex', $
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 255 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 255 
	contour, (ey<range)>(-range), capR, zLoc, $
			color = 255, $
			levels = levels, $
			c_colors = colors, $
			title = 'ey', $
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 255 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 255 
	contour, (ez<range)>(-range), capR, zLoc, $
			color = 255, $
			levels = levels, $
			c_colors = colors, $
			title = 'ez', $
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 255 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 255 
	
   	contour, (eAlpha<range)>(-range), capR, zLoc, $
			color = 255, $
			levels = levels, $
			c_colors = colors, $
			title = 'eAlpha', $
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 255 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 255 

	contour, (eBeta<range)>(-range), capR, zLoc, $
			color = 255, $
			levels = levels, $
			c_colors = colors, $
			title = 'eBeta', $
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 255 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 255 

   	contour, ((eParallel*100)<range)>(-range), capR, zLoc, $
			color = 255, $
			levels = levels, $
			c_colors = colors, $
			title = 'eParallel x100', $
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 255 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 255 


	contour, (bx_wave_real<brange)>(-brange), capR, zLoc, $
			color = 255, $
			levels = blevels, $
			c_colors = bcolors, $
			title = 'bx_wave_real', $
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
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
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
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
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
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
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
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
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
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
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 0 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 0 

	!p.multi = 0   
	loadct, 0, /silent
   	contour, mask, capR, zLoc, $
			title = 'mask', $
			/fill
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 0 
   	loadct, 12, /silent
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   psym = -4, $
		   symSize = 2.0, $
		   color = 8*16-1
   	loadct, 0, /silent
	plots, rebin ( capR, n_elements ( capR ), n_elements ( zLoc ) ), $
			transpose ( rebin ( zLoc, n_elements ( zLoc ), n_elements ( capR ) ) ), $
			psym = 1, $
			symSize = 0.5
	loadct, 12, /silent
	plots, newR, newZ, psym = 5, color = 12*16-1

	if (not keyword_set(flat)) then begin
	!p.multi = [0,2,2]
	contour, density, capR, zLoc, $
		   nlev = 20, $
		   title = 'density + janty'
   	contour, janty, capR, zLoc, $
		   c_color = 8*16-1, $
		   /overplot, $ 
   			levels = fIndGen(10)*10
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 0 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   color = 8*16-1

   	!p.multi=[7,2,4]
	plot, capR, density[*,n_elements(density[0,*])/2], $
			xTitle = 'R[m]', $
			charSize = 2, $
			xCharSize = 1.0
	axis, xaxis=1, $
			xrange=[-2,2], $
		   	/save, $
			xTitle = 'z[m]', $
			color = 8*16-1, $
			charSize = 1.5
	oplot, zloc, density[n_elements(density[*,0])/2,*], $
			color = 8*16-1
   	!p.multi=[5,2,4]
	plot, capR, density[*,n_elements(density[0,*])/2], $
			xTitle = 'R[m]', $
			charSize = 2, /yLog, $
			yRange = [1e16,1e20], $
			yStyle = 1, $
			min_val = 1e12
	axis, xaxis=1, $
			xrange=[-2,2], $
		   	/save, $
			xTitle = 'z[m]', $
			color = 8*16-1, $
			charSize = 1.5 
	oplot, zloc, density[n_elements(density[*,0])/2,*], $
			color = 8*16-1
	
	!p.multi = [2,2,2]

	nxGrid	= 25
	nyGrid	= 50
	densityGrid	= congrid ( density, nxGrid, nyGrid, /interp )
	capRGrid	= congrid ( capR, nxGrid, /interp )
	zLocGrid	= congrid ( zLoc, nyGrid, /interp )

	loadct, 3, /silent
   	surface, densityGrid, capRGrid, zLocGrid, $
			title = 'density', $
			charSize = 1.0, $
			ax = 80, $
			az = -15, font=0, /save, $
			shades = (255 - ( bytScl ( densityGrid, top = 253, min = 1e17, max = 1e20) + 1 ))<220, $
			color = 0, $
			zRange = [1e16,5e19]
	loadct, 12, /silent
	plots, eqdsk.rbbbs, eqdsk.zbbbs, eqdsk.rbbbs*0+2e18, $
			color = 12*16-1, $
		   	/t3d, $
			thick = 5
	plots, eqdsk.rlim, eqdsk.zlim, eqdsk.rbbbs*0+2e18, $
			color = 8*16-1, $
		   	/t3d, $
			thick = 5
	
   	surface, density, capR, zLoc, $
			ax = 40, $
			az = -80, $
			title = 'density [log]', $
			font = 0, $
			charSize = 1.0, $
			zRange = [1.0e16,1000.0e17], $
			zStyle = 1, $
			min_val = 1.0e12, $
		   	/zlog, $
			/save
	plots, eqdsk.rbbbs, eqdsk.zbbbs, eqdsk.rbbbs*0+2e18, $
			color = 12*16-1, $
		   	/t3d, $
			thick = 3
	plots, eqdsk.rlim, eqdsk.zlim, eqdsk.rbbbs*0+2e18, $
			color = 8*16-1, $
		   	/t3d, $
			thick = 3
   endif	

   loadct, 13, file = 'davect.tbl'	
   !p.multi = [0,3,2]
	contour, (jedote<prange)>(-prange), capR, zLoc, $
			color = 255, $
			levels = plevels, $
			c_colors = pcolors, $
			title = 'jedote', $
			/fill, $
			yRange = yRange, $
			xRange = xRange, $
			xSty = 1, ySty = 1
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
		   thick = 2, $
		   color = 255 
	oPlot, eqdsk.rLim, eqdsk.zLim, $
		   thick = 2, $
		   color = 255 
	!p.multi = 0
  	
   device, /close_file
stop

end 
