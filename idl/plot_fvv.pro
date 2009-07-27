pro plot_fvv

	k   = 1.3806504d-23
	e_   = 1.60217646d-19
	amu	= 1d0
	mi  = amu * 1.67262158d-27
	q   = 1d0 * e_
	c   = 3.0d8

	if strCmp ( getEnv ( 'MACHINE' ), 'jaguar' ) then scratchDir = getEnv ( 'SYSTEM_USERDIR' )
	if strCmp ( getEnv ( 'MACHINE' ), 'franklin' ) then scratchDir = getEnv ( 'SCRATCH' )
	if strCmp ( getEnv ( 'MACHINE' ), 'dlghp' ) then scratchDir = '/home/dg6/scratch' 

	cdfId	= ncdf_open ( 'output/fvv.nc', /noWrite )
	glob	= ncdf_inquire ( cdfId )
	
	ncdf_varGet, cdfId, 'f_vv', f_vv 
	ncdf_varGet, cdfId, 'dfduper', dfduper 
	ncdf_varGet, cdfId, 'dfdupar', dfdupar 
	ncdf_varGet, cdfId, 'i', i 
	ncdf_varGet, cdfId, 'j', j
	ncdf_varGet, cdfId, 'R', R
	ncdf_varGet, cdfId, 'z', z
	ncdf_close, cdfId

;	cdfId	= ncdf_open ( 'output/ql.nc', /noWrite )
;	glob	= ncdf_inquire ( cdfId )
;	
;	ncdf_varGet, cdfId, 'uper', uper 
;	ncdf_varGet, cdfId, 'upar', upar 
;	nCdf_varGet, cdfId, 'ql_b', ql_b
;	nCdf_varGet, cdfId, 'ql_c', ql_c
;	nCdf_varGet, cdfId, 'ql_e', ql_e
;	nCdf_varGet, cdfId, 'ql_f', ql_f
;	ncdf_close, cdfId

	cdfId	= ncdf_open ( 'output/f_vvp.nc', /noWrite )

	glob	= ncdf_inquire ( cdfId )
	
	ncdf_varGet, cdfId, 'uper', uper_ 
	ncdf_varGet, cdfId, 'upar', upar_ 
	nCdf_varGet, cdfId, 'rho', rho_
 	nCdf_varGet, cdfId, 'f_vvp', f_vvp 
	nCdf_varGet, cdfId, 'vc_mks', vc_mks
	ncdf_close, cdfId

	!p.multi = [0,1,3]
	!p.charSize = 1.5
	plot, dfdupar[*,n_elements(dfduper[0,*])/2],psym=-4
	plot, dfduper[*,n_elements(dfduper[0,*])/2],psym=-4
	plot, f_vv[*,n_elements(dfduper[0,*])/2],psym=-4

stop
	cql_per_keV	= ( uper_ * vc_mks )^2 * 0.5 * mi / e_ / 1d3
	cql_par_keV	= ( upar_ * vc_mks )^2 * 0.5 * mi / e_ / 1d3

	;cdfId	= ncdf_open ( 'fdis_cmod_t1000_2kev.delta.NFO.40x40.nc', /noWrite )
	cdfId	= ncdf_open ( scratchDir+'/p2f/cmod_t1000_1e6_delta/data/'$
		+ 'fdis.dav.nc', /noWrite )

	;cdfId	= ncdf_open ( 'fdis_D3D.nc', /noWrite )
	;cdfId	= ncdf_open ( 'fdis_heidbrink.nc', /noWrite )

	glob	= ncdf_inquire ( cdfId )
	
	ncdf_varGet, cdfId, 'f_rzvv', f_rzvv
	ncdf_varGet, cdfId, 'R_binCenters', R_binCenters
	ncdf_varGet, cdfId, 'z_binCenters', z_binCenters
	ncdf_varGet, cdfId, 'vPerp_binCenters', vPerp_binCenters
	ncdf_varGet, cdfId, 'vPar_binCenters', vPar_binCenters
	ncdf_varGet, cdfId, 'R_binEdges', R_binEdges
	ncdf_varGet, cdfId, 'z_binEdges', z_binEdges
	ncdf_varGet, cdfId, 'vPerp_binEdges', vPerp_binEdges
	ncdf_varGet, cdfId, 'vPar_binEdges', vPar_binEdges
	ncdf_varGet, cdfId, 'nP', nP

	ncdf_close, cdfId

	par_per_keV	= ( vPerp_binCenters )^2 * 0.5 * mi / e_ / 1d3
	par_par_keV	= ( vPar_binCenters )^2 * 0.5 * mi / e_ / 1d3


	;cdfId	= ncdf_open ( 'fdis_cmod_t1000_2kev.smooth.NFO.40x40.nc', /noWrite )
	;cdfId	= ncdf_open ( 'fdis_cmod_t1000_2kev.delta.NFO.40x40.nc', /noWrite )
	;cdfId	= ncdf_open ( 'fdis_D3D_NFO.nc', /noWrite )
	;cdfId	= ncdf_open ( 'fdis_heidbrink_smooth.nc', /noWrite )
	cdfId	= ncdf_open ( scratchDir+'/p2f/cmod_t1000_1e6_noAve/data/'$
		+ 'fdis.dav.nc', /noWrite )

	glob	= ncdf_inquire ( cdfId )
	
	ncdf_varGet, cdfId, 'f_rzvv', f_rzvv_s
	ncdf_varGet, cdfId, 'R_binCenters', R_binCenters_s
	ncdf_varGet, cdfId, 'z_binCenters', z_binCenters_s
	ncdf_varGet, cdfId, 'vPerp_binCenters', vPerp_binCenters_s
	ncdf_varGet, cdfId, 'vPar_binCenters', vPar_binCenters_s
	ncdf_varGet, cdfId, 'R_binEdges', R_binEdges_s
	ncdf_varGet, cdfId, 'z_binEdges', z_binEdges_s
	ncdf_varGet, cdfId, 'vPerp_binEdges', vPerp_binEdges_s
	ncdf_varGet, cdfId, 'vPar_binEdges', vPar_binEdges_s
	ncdf_varGet, cdfId, 'nP', nP_s

	ncdf_close, cdfId



	;!p.background = 255
	;device, decomposed = 0
	loadct, 0
	;window, 0, xSize = 1200, ySize =500 
	
	old_dev = !D.name
	set_plot, 'ps'
	outfname	= 'output/f_vvp.eps'
	device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
		xsize=36, ysize=15,xoffset=.1, yoffset=.1, /encapsul

	eqDskFileName	= 'g129x129_1051206002.01120.cmod'
	;eqDskFileName	= 'g122080.03100'
	;eqDskFileName	= 'eqdsk.122993'
	eqdsk	= readGEQDSK ( eqDskFileName )


	;iiArrR	= [32,32,32,32,32,32,32,32]
	;iiArrZ	= [19,21,23,25,27,29,31,32]
	iiArrR	= [33,35,37,39,41,43,45,47]/2.0
	iiArrZ	= [19,19,19,19,19,19,19,19]
	;iiArrR	= [32,30,28,26,24,22,20,18]
	;iiArrZ	= [19,19,19,19,19,19,19,19]
	;iiArrR	= [20,22,24,26,28,30,32,34]
	;iiArrZ	= [19,19,19,19,19,19,19,19]
	
	xRange	= [ min(upar_), max(upar_) ] / 3.0
	yRange	= [ 0, max ( uper_ ) ]
	;xRange	= [ min(vPar_bincenters), max(vPar_binCenters) ] / vc_mks 
	;yRange	= [ 0, max ( vPerp_binCenters ) ] / vc_mks 
	
	!p.charsize =0.01
	cql_duPer	= abs(uper_[0]-uper_[1])
	cql_duPar	= abs(upar_[0]-upar_[1])
	cql_uPer2D	= rebin ( uper_, n_elements ( uper_ ), n_elements ( upar_ ) )
	cqlNorm	= total ( f_vvp[*,*,0] * cql_duPer * cql_duPar * cql_uper2D )
	
	par_duPer	= abs(vPerp_binCenters[0]-vPerp_binCenters[1]) / vc_mks
	par_duPar	= abs(vPar_binCenters[0]-vPar_binCenters[1]) / vc_mks
	par_uPer2D	= rebin ( vPerp_binCenters, n_elements ( vPerp_binCenters ), n_elements ( vPar_binCenters ) ) / vc_mks
	iiR	= where ( abs ( R_binCenters - eqdsk.rmaxis ) eq min ( abs ( R_binCenters - eqdsk.rmaxis ) ) )
	iiz	= where ( abs ( z_binCenters - eqdsk.zmaxis ) eq min ( abs ( z_binCenters - eqdsk.zmaxis ) ) )
	parNorm	= total ( f_rzvv[iiR, iiz,*,*] * par_duPer * par_duPar * par_uper2D )
	
	par_duPer_s	= abs(vPerp_binCenters_s[0]-vPerp_binCenters_s[1]) / vc_mks
	par_duPar_s	= abs(vPar_binCenters_s[0]-vPar_binCenters_s[1]) / vc_mks
	par_uPer2D_s	= rebin ( vPerp_binCenters_s, n_elements ( vPerp_binCenters_s ),$
		 n_elements ( vPar_binCenters_s ) ) / vc_mks
	iiR_s	= where ( abs ( R_binCenters_s - eqdsk.rmaxis ) eq min ( abs ( R_binCenters_s - eqdsk.rmaxis ) ) )
	iiz_s	= where ( abs ( z_binCenters_s - eqdsk.zmaxis ) eq min ( abs ( z_binCenters_s - eqdsk.zmaxis ) ) )
	parNorm_s	= total ( f_rzvv_s[iiR_s, iiz_s,*,*] * par_duPer_s * par_duPar_s * par_uper2D_s )
	
	
	for i=0, n_elements ( iiArrR ) -1  do begin
			
		psiRange	= abs ( eqdsk.siMag - eqdsk.siBry )
		;if i eq 0 then begin
		;	psi  = interpolate ( eqdsk.psizr, ( eqdsk.rmaxis - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nW-1.0), $
    	;		( eqdsk.zmaxis - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nH-1.0) )
		;	iiArrR[i]	= where ( abs ( R_binCenters - eqdsk.rmaxis ) eq min ( abs ( R_binCenters - eqdsk.rmaxis ) ) )
		;	iiArrz[i]	= where ( abs ( z_binCenters - eqdsk.zmaxis ) eq min ( abs ( z_binCenters - eqdsk.zmaxis ) ) )
		;endif else begin	
			psi  = interpolate ( eqdsk.psizr, ( R_binCenters[iiArrR[i]] - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nW-1.0), $
   		 		( z_binCenters[iiArrZ[i]] - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nH-1.0) )
		;endelse 

		psiNorm	= ( psi - eqdsk.siMag ) / psiRange			
		rho	= sqrt ( psiNorm )
	
		iiCqlRho	= where ( abs(rho_ - rho) eq min ( abs ( rho_ - rho ) ) )

		print, R_binCenters[iiArrR[i]], z_binCenters[iiArrZ[i]],iiCqlRho, reform(rho), psiNorm

		!p.multi = [n_elements(iiArrR)*3-i,n_elements(iiArrR),3]
		cqlf	= 	f_vvp[*,*,iiCqlRho]
		cqlf	= cqlf / cqlNorm 
		title	= 'i: '+string(iiArrR[i],for='(i2.2)')$
			+', j: '+string(iiArrz[i],for='(i2.2)')$
			+', R:'+string(R_binCenters[iiArrR[i]],for='(f5.2)')$
			+', z:'+string(z_binCenters[iiArrZ[i]],for='(f5.2)')$
			+', rho:'+string(rho,for='(f5.2)')$
			+', max_keV: '+string((max([xRange,yRange])*vc_mks)^2*0.5*mi/e_/1d3,for='(f4.0)')
		levels	= 10.0^fIndgen(10)*1e-6
		colors	= ((reverse ( (bytScl ( fIndGen(10), top = 253 ) + 1 ))-100)>0)*2<254
		contour, transpose(cqlf),upar_, uper_, $
			levels=levels, $
			xRange = xRange, yRange = yRange, title = title, $
			color = 0, c_colors = colors, xStyle = 4, yStyle = 4, c_thick = 3.0, c_label = colors*0+1, c_charSize = 0.5
		xyOuts, 0.35,0.65, 'rho: '+string(rho,for='(f4.2)'), align = 1.0, charSize = 1.5, color = 0, charThick = 3.0
		!p.multi = [n_elements(iiArrR)*2-i,n_elements(iiArrR),3]
		parf	= reform(f_rzvv[iiArrR[i],iiArrZ[i],*,*])
		parf	= parf / parNorm 
		contour, transpose(parf),$
			vPar_binCenters/vc_mks,vPerp_binCenters/vc_mks, $
			levels=levels, xRange = xRange, yRange = yRange, $
			c_colors = colors, color = 0, xStyle = 4, yStyle = 4, c_thick = 3.0, c_label = colors*0+1, c_charSize = 0.5
		!p.multi = [n_elements(iiArrR)-i,n_elements(iiArrR),3]
		parf_s	= reform(f_rzvv_s[iiArrR[i],iiArrZ[i],*,*])
		parf_s	= parf_s / parNorm_s
		contour, transpose(parf_s),$
			vPar_binCenters_s/vc_mks,vPerp_binCenters_s/vc_mks, $
			levels=levels, xRange = xRange, yRange = yRange, $
			c_colors = colors, color = 0, xStyle = 4, yStyle = 4, c_thick = 3.0, c_label = colors*0+1, c_charSize = 0.5

	endfor

	device, /close

;	window, 1, xSize = 1700, ySize = 100
;
;	!p.multi = [0,8,1]
;	for i = 0, 7 do begin
;		cqlf	= 	f_vvp[*,*,i*8]
;		cqlf	= cqlf / cqlNorm 
;		title	= 'rho:'+string(rho_[i*8],for='(f5.2)')$
;			+', max_keV: '+string((max([xRange,yRange])*vc_mks)^2*0.5*mi/e_/1d3,for='(f4.0)')
;		contour, transpose(cqlf),upar_, uper_, $
;			levels=levels, $
;			xRange = xRange, yRange = yRange, title = title, $
;			c_colors = colors, color = 0
;	endfor

stop
end
