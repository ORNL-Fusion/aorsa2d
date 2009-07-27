pro plot_wdot

	if strCmp ( getEnv ( 'MACHINE' ), 'jaguar' ) then scratchDir = getEnv ( 'SYSTEM_USERDIR' )+'/'
	if strCmp ( getEnv ( 'MACHINE' ), 'franklin' ) then scratchDir = getEnv ( 'SCRATCH' )+'/'
	if strCmp ( getEnv ( 'MACHINE' ), 'dlghp' ) then scratchDir = '/home/dg6/scratch/'

	runList	= [ 'AORSA2D/AORSA-15-DLG/CMOD_H_MINORITY_2.4MW/t_00ms_ana/', $
		'AORSA2D/AORSA-15-DLG/CMOD_H_MINORITY_2.4MW/t_00ms/', $
		'AORSA2D/AORSA-15-DLG/CMOD_H_MINORITY_2.4MW/t_01ms/',$ 
		'AORSA2D/AORSA-15-DLG/CMOD_H_MINORITY_2.4MW/t_02ms/' ]
	;eqDskFileName	= 'g122080.03100'
	eqDskFileName	= 'g129x129_1051206002.01120'
	eqdsk	= readGEQDSK ( eqDskFileName )
	
	;	Plot the flux surface averaged wdot arrays

	loadct, 12	
	device, decomposed = 0
	!p.charSize = 2.0
	!p.background = 255
	!p.thick = 2.0
	;window, 0, xSize = 1200, ySize = 400
	old_dev = !D.name
	set_plot, 'ps'
	outfname	= 'output/wdot_1D.eps'
	device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
		xsize=40, ysize=10,xoffset=.1, yoffset=.1, /encapsul
	

	for i = 0, n_elements ( runList ) - 1 do begin
	
	if i eq 0 then $
		cdfId = ncdf_open ( 'output/plotData.nc', /noWrite ) $
	else $
		cdfId = ncdf_open ( scratchDir + runList[i] + 'output/plotData.nc', /noWrite )

	glob	= ncdf_inquire ( cdfId )
	
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

	ncdf_close, cdfId


	;	color	= [0,0,0,$
	;				14*16-1,14*16-1,14*16-1,$
	;				12*16-1,12*16-1,12*16-1,$
	;				8*16-1,8*16-1,$
	;				4*16-1,4*16-1,4*16-1,4*16-1]
		color	= ([0,7,2,12, 11,11,12,13] * 16 - 1)>0
		style	= [0,0,0,0, 0,0,0,0]
		thick	= [2,1,1,1, 1,2,2,2]

	;if i gt 3 then loadct, 1
	
	print, color
	print, pscale
	!p.charSize = 2.0
	;pscale =  median ( wdoti1 )
	pscale = 1.0
	;print, pscale

	!p.multi = [3,3,1]
	if i eq 0 then begin
		yRange_e	= [0,  1.5*max (wdote/abs(pscale))]
		plot, rho, wdote/abs(pscale), $
			color = color[i], lineStyle = style[i], xTitle = 'rho', yTitle = 'wDote', $
			yRange = yRange_e, thick = thick[i]
		xSave_e	= rho
		ySave_e	= wdote/abs(pscale)
	endif else begin
		plot, xSave_e, ySave_e, /noData, $
			yRange = yRange_e, color = 0, /noErase
		oplot, rho, wdote/abs(pscale), color = color[i], lineStyle = style[i], thick = thick[i] 
	endelse

	!p.multi = [2,3,1]
	if i eq 0 then begin
		yRange_i1 = [0, 1.5*max (wdoti1/abs(pscale))]
		plot, rho, wdoti1/abs(pscale), $
			color = color[i], xTitle = 'rho', yTitle = 'wDoti1', $
			yRange = yRange_i1, lineStyle = style[i], thick = thick[i]
		xSave_i1	= rho
		ySave_i1	= wdoti1/abs(pscale)
		print, max(wdoti1/abs(pscale))
	endif else begin
		plot, xSave_i1, ySave_i1, /noData, $
			yRange = yRange_i1, color = 0, /noErase
		oplot, rho, wdoti1/abs(pscale), color = color[i], lineStyle = style[i], thick = thick[i] 
		print, max(wdoti1/abs(pscale))
	endelse

	!p.multi = [1,3,1]
	if i eq 0 then begin
		yRange_i2 = [0, 1.5*max (wdoti2/abs(pscale))*1d-6]
		plot, rho, wdoti2/abs(pscale)*1d-6, $
			color = color[i], xTitle = 'rho', yTitle = 'power [watts/m!U3!N] x10!U6!N', $
			yRange = yRange_i2, lineStyle = style[i], thick = thick[i]
		xSave_i2	= rho
		ySave_i2	= wdoti2/abs(pscale)
	endif else begin
		plot, xSave_i2, ySave_i2*1d-6, /noData, $
			yRange = yRange_i2, color = 0, /noErase 
		oplot, rho, wdoti2/abs(pscale)*1d-6, color = color[i], lineStyle = style[i], thick = thick[i]
	endelse

	xyOuts, 0.95, 0.85-i*0.05, file_baseName ( runList[i] ),$
		color = color[i], /norm, align = 1.0, charSize = 1.0

	endfor

	device, /close
	!p.multi = 0
	!p.position = 0	
	outfname	= 'output/wdot_2D.eps'
	device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
		xsize=8, ysize=10,xoffset=.1, yoffset=.1, /encapsul
	
	nLevs	= 31
	levels	= (fIndGen(nLevs)-15.0)/nLevs*1e6
	colors	= bytScl ( levels, top = 253 ) + 1.0
	loadct, 13, file = 'davect.tbl'  

	contour, wdoti2_rz, capR, zLoc, $
		levels = levels, c_colors = colors, /fill, color = 255
	contour, wdoti2_rz, capR, zLoc, $
		levels = levels, c_colors = (colors*1.2)<254, /over
	loadct, 0
	plots, eqdsk.rbbbs, eqdsk.zbbbs, color = 0
	plots, zLoc * 0 + eqdsk.rmaxis, zLoc, color = 0, lineStyle = 1

	nLevs	= 10	
	levels	= fIndGen(nLevs)*( eqdsk.sibry - eqdsk.simag )	/ nLevs + eqdsk.simag
	contour, eqdsk.psizr, eqdsk.r, eqdsk.z, /over, $
		levels = levels, color = 200

	;	Resonant surfaces

	amu	= 1.0
	e_	= 1.60217e-19
	m	= 1.6726e-27
	q	= amu * e_
	freq	= 60.0e6 * 2.0 * !pi
	harm	= [3.0,4.0,5.0]

	B	= freq * m / ( q * harm/2.0 )

	thisR	= fltArr(n_elements(eqdsk.z))
	for i=0,n_elements(harm)-1 do begin
		for z=0,n_elements(eqdsk.bMag[0,*])-1 do begin
			rIndex	= where(abs(eqdsk.bMag[*,z]-B[i]) eq min(abs(eqdsk.bMag[*,z]-B[i])))
			thisR[z]	= eqdsk.r[rIndex]	
		endfor	
		plots, thisR, eqdsk.z, color = 0	
		xyOuts, mean(thisR), 1.3, string(harm[i],for='(i1.1)'),$
			charSize = 1.5, color = 0, align = 0.5, charThick = 2.0
		print, harm[i], mean(thisR), max(thisR)
	endfor

	device, /close
	set_plot, old_dev
	stop
end
