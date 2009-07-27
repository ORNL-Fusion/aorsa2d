pro plot_cql3d

	fileName	= 'cmod_scidac_test.0.nc' 
	cdfId	= ncdf_open ( fileName, /noWrite )
	glob	= ncdf_inquire ( cdfId )

	ncdf_varGet, cdfId, 'f', f
	ncdf_varGet, cdfId, 'rya', rya
	ncdf_varGet, cdfId, 'enorm', enorm
	ncdf_varGet, cdfId, 'vnorm', vnorm
	ncdf_varGet, cdfId, 'y', y
	ncdf_varGet, cdfId, 'x', x
	ncdf_varGet, cdfId, 'dy', dy
	ncdf_varGet, cdfId, 'dx', dx
	ncdf_varGet, cdfId, 'iy_', iy_ 
	ncdf_varGet, cdfId, 'density', density_file
	ncdf_varGet, cdfId, 'wperp', wperp
	ncdf_varGet, cdfId, 'wperp', wpar 
	ncdf_varGet, cdfId, 'time', time
	ncdf_varGet, cdfId, 'powers', powers 


	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'rdim' ), name, rdim
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'xdim' ), name, xdim
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'ydim' ), name, ydim
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'tdim' ), name, tdim 

	ncdf_close, cdfId

	x2D	= transpose ( rebin ( x, xdim, ydim ) )
	y2D	= ( rebin ( y[*,0], ydim, xdim ) )

	uPar	= cos ( y2D ) * x2D
	uPer	= sqrt ( x2D^2 - uPar^2 )

	uPar_max	= max ( abs ( uPar ) )
	uPer_max	= max ( uPer )
	uPer_nBins	= 128 
	uPer_reg	= (fIndGen ( uPer_nBins )) / ( uPer_nBins-1 ) * uPer_max
	uPar_nBins	= 256 
	uPar_reg	= fIndGen ( uPar_nBins ) / ( uPar_nBins - 1 ) * (2.0 * uPar_max ) - uPar_max
	uPer_binSize	= abs ( uPer_reg[0] - uPer_reg[1] )
	uPar_binSize	= abs ( uPar_reg[0] - uPar_reg[1] )

	uParTmp	= uPar - min ( uPar )
	uParTmp	= uParTmp / max ( uParTmp )
      
	f_rvv	= fltArr ( rdim, uPer_nBins, uPar_nBins ) 

	triangulate, uPar[*], uPer[*], tri, b
	for i=0,rdim-1 do begin

		;f_rvv[i,*,*]	= tsc ( (f[*,*,i])[*],$
		;	uPer[*]*(uPer_nBins-1), uPer_nBins, $
		;	(uParTmp[*]*uPar_nBins)<(uPar_nBins-1), uPar_nBins, /iso )

		f_rvv[i,*,*]	= triGrid ( uPer[*], uPar[*], (f[*,*,i])[*], tri, $
			[uPer_binSize, uPar_binSize], $
			[0, -uPar_max, uPer_max, uPar_max] )

	endfor

	device, decomposed = 0
	!p.background = 255
	window, 4, xSize = 1600, ySize = 200
	!p.multi = [0,8,1]
	!p.charSize = 2.0

	for i=0,7 do begin

		if i eq 0 then yTitleSize = 2.0 $
			else yTitleSize = 0.01
		
		contour, transpose(f_rvv[i*rdim/8,*,*]), uPar_reg, uPer_reg, $
			levels = 10d0^fIndGen(20), $
			xTitle = 'uPar', $
			yTitle = 'uPer', $
			color = 0, $
			yCharSize = yTitleSize	

		;contour, (f[*,*,25])[*], uPar[*], uPer[*], $
		;	levels = 10d0^fIndGen(20), $
		;	/irreg

	endfor
	

;	Convert f to SI units from  dammit!
;	vnorm^3/(cm^3*(cm/sec)^3) -> 1/(m^3*(m/s)^3)

	f_rvv	= f_rvv / ( vnorm^3) * (1d2)^6 

;	Get density profile

	vPer_binSize	= uPer_binSize * vnorm * 1d-2
	vPar_binSize	= uPar_binSize * vnorm * 1d-2
	vPer_reg	= uPer_reg * vnorm * 1d-2
	vPar_reg	= uPar_reg * vnorm * 1d-2

	vPer2D	= transpose ( rebin ( vPer_reg, uPer_nBins, rdim ) ) 
	vPar3D	= transpose ( rebin ( vPar_reg, uPar_nBins, uPer_nBins, rdim ) ) 
	density	= total ( total ( f_rvv, 3 ) * vPer2D, 2 ) $
		* 2.0 * !pi * vPer_binSize * vPar_binSize

;	Fit with 2D gaussian to get temperature profile

	k   = 1.3806504d-23
	e_   = 1.60217646d-19
	amu	= 1d0
	mi  = amu * 1.67262158d-27
	q   = 1d0 * e_
	c   = 3.0d8

;	iiDiscard	= where ( vPer2D gt 2.8e6, iiDisCnt )
;	vPer2D[iiDiscard]=0
;	iiDiscard	= where ( abs ( vPar3D ) gt 1.4e6, iiDisCnt )
;	f_rvv[iiDiscard] = 0

	wPerp_dlg	= total ( 0.5 * mi * vPer2D^2 * total ( f_rvv, 3 ) * vPer2D, 2 ) $
		* 2.0 * !pi * vPer_binSize * vPar_binSize / density / e_ * 1d-3
	wPar_dlg	= total ( total ( 0.5 * mi * vPar3D^2 * f_rvv, 3 ) * vPer2D, 2 ) $
		* 2.0 * !pi * vPer_binSize * vPar_binSize / density / e_ * 1d-3



	temp	= fltArr ( rdim )
	for i = 0, rdim-1 do begin
		print, i
		
			tmp	= [reverse(reform(f_rvv[i,*,*]),1),reform(f_rvv[i,*,*])]
			if mean ( tmp ) gt 0 then begin
				fit	= gauss2dFit ( tmp, A )
				temp[i]	= mi * ( $
					( A[2] * uPer_binsize * vnorm * 1d-2 )^2 + ( A[3] * uPar_binSize * vnorm * 1d-2 )^2 $
						) / ( 2.0 * 1e3 * e_ )  
			endif

	endfor

	!p.charSize = 3.0
	!p.multi	= [0,1,3]
	window, 1, xSize = 400, ySize = 750

	plot, rya, density/1d19, $
		psym = -4, $
		color = 0, $
		xTitle = 'rho', $
		yTitle = 'density [m^-3] x1e19', $
		yRange = [0, 2], $
		yStyle = 1, $
		title = 'CQL3D'

	for i=0,tdim-1 do begin

		if i lt tdim-1 then thick=1 $
			else thick=3	

		if i eq 0 then $
		plot, rya, wperp[*,i], $
			color = 0, $
			xTitle = 'rho', $
			yTitle = 'wPerp [keV]', $
			yRange = [0, 8], $
			yStyle = 1, $
			thick = thick

		if i gt 0 then $
		oPlot, rya, wperp[*,i], $
			color = 0, $
			thick = thick
	
	endfor

	loadct, 12
	oPlot, rya, wPerp_dlg, color = 8*16-1
	oPlot, rya, wPar_dlg, color = 12*16-1
	loadct, 0
	
	plot, rya, temp, $
		psym = -4, $
		color = 0, $
		xTitle = 'rho', $
		yTitle = 'temp [keV]', $
		yRange = [0, 5], $
		yStyle = 1

	window, 3, xSize = 400, ySize = 750
	
	powerMax	= 1d7

	plot, rya, powers[*,4,0,0]*1d6, $
			color = 0, $
			yRange=[0,powerMax], $
			xTitle = 'rho', $
			title = 'spec1'
	plot, rya, powers[*,4,0,1]*1d6, $
			color = 0, $
			yRange=[0,powerMax], $
			title = 'spec2'
	plot, rya, powers[*,4,0,2]*1d6, $
			color = 0, $
			yRange=[0,powerMax], $
			title = 'spec3'

	wPerp_cql3d	= wPerp[*,n_elements(wperp[0,*])-1]
	wPar_cql3d	= wPar_dlg
	rho_cql3d	= rya
	save, /variables, fileName = 'cql3dData.sav'

stop
end
