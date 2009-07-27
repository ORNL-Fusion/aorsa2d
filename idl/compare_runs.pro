pro compare_runs

	runs	= ['128x128/60MHz/num_60keV_max_flat', $
				'128x128/60MHz/num_60keV_max_flat_ba' ]

	loadct, 12
	device, decomposed = 0
	window, 0, xSize = 600, ySize = 800
	!p.multi = [0,1,3]
	!p.background = 255
	!p.charSize = 2

	firstRun	= 1
	color	= 1
	for i=0,n_elements(runs)-1 do begin

		cdfId	= ncdf_open ( runs[i]+'/output/plotData.nc', /noWrite )
		glob	= ncdf_inquire ( cdfId )
		
		ncdf_varGet, cdfId, 'rho', rho 
		ncdf_varGet, cdfId, 'wdote', wdote
	
		ncdf_close, cdfId

		if firstRun eq 1 then begin

			plot, rho, wdote/max(wdote), $
				thick = 2.0, $
				xTitle = 'rho', $
				yTitle = 'norm power', color = 0

			firstRun = 0
			color = color + 3

		endif else begin

			oPlot, rho, wdote/max(wdote), thick = 2.0, $
				color = color * 16 - 1
			color = color + 3

		endelse

	endfor


	firstRun	= 1
	color	= 1
	for i=0,n_elements(runs)-1 do begin

		cdfId	= ncdf_open ( runs[i]+'/output/plotData.nc', /noWrite )
		glob	= ncdf_inquire ( cdfId )
		
		ncdf_varGet, cdfId, 'rho', rho 
		ncdf_varget, cdfId, 'wdoti1', wdoti1 
	
		ncdf_close, cdfId

		if firstRun eq 1 then begin

			plot, rho, wdoti1/max(wdoti1), $
				thick = 2.0, $
				xTitle = 'rho', $
				yTitle = 'norm power', color = 0

			firstRun = 0
			color = color + 3

		endif else begin

			oPlot, rho, wdoti1/max(wdoti1), thick = 2.0, $
				color = color * 16 - 1
			color = color + 3

		endelse

	endfor

	firstRun = 1
	color = 1
	for i=0,n_elements(runs)-1 do begin

		cdfId	= ncdf_open ( runs[i]+'/output/plotData.nc', /noWrite )
		glob	= ncdf_inquire ( cdfId )
		
		ncdf_varGet, cdfId, 'rho', rho 
		ncdf_varget, cdfId, 'wdoti2', wdoti2 
		;ncdf_varget, cdfId, 'capR', capR 
		;ncdf_varget, cdfId, 'zLoc', zLoc 
		;ncdf_varGet, cdfId, 'wdote_rz', wdote_rz
		;ncdf_varget, cdfId, 'wdoti1_rz', wdoti1_rz
		;ncdf_varget, cdfId, 'wdoti2_rz', wdoti2_rz
	
		ncdf_close, cdfId

		if firstRun eq 1 then begin

			plot, rho, wdoti2, $
				thick = 2.0, $
				xTitle = 'rho', $
				yTitle = 'norm power', color = 0

			firstRun = 0
			color = color + 3

		endif else begin

			oPlot, rho, wdoti2, thick = 2.0, $
				color = color * 16 - 1
			color = color + 3

		endelse

	endfor


stop
end
