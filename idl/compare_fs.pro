pro compare_fs , lists = lists
	; Constants

	k   = 1.3806504d-23
	e_   = 1.60217646d-19
	amu	= 2d0
	mi  = amu * 1.67262158d-27
	q   = 1d0 * e_
	c   = 3.0d8

	nP	= 100000L

	data	= { psi : 0.0, $
					theta : 0.0, $
					R : 0.0, $
					z : 0.0, $
					E : 0.0, $
					lambda : 0.0, $
					weight : 0.0 }
	
	if keyword_set ( lists ) then begin

		mhc_data	= replicate ( data, nP )
		dlg_data	= replicate ( data, nP )
	
		mhc_file	= '~/data/particleLists/Max_25keV_nospatial_1mil'
		dlg_file	= '~/data/particleLists/fdis_25keV_H_6_g122080.03100.dav'
	
		openr, unit, mhc_file, /get_lun
		readf, unit, mhc_data
		close, unit
	
		openr, unit, dlg_file, /get_lun
		readf, unit, dlg_data
		close, unit
	
		mhc_v	= sqrt ( mhc_data.e * 1e3 * e_ * 2.0 / mi ) 
		dlg_v	= sqrt ( mhc_data.e * 1e3 * e_ * 2.0 / mi ) 
	
		mhc_eDist	= histogram ( mhc_data.e, loca = mhc_e )
		dlg_eDist	= histogram ( dlg_data.e, loca = dlg_e )
	
		mhc_binSize	= (mhc_e[1]-mhc_e[0])*1e3*e_
		dlg_binSize	= (mhc_e[1]-mhc_e[0])*1e3*e_
	
		!p.multi = [0,1,2]
		plot, mhc_e, mhc_eDist/mhc_binSize, xRange = [0,100]
		oplot, dlg_e, dlg_eDist/dlg_binSize
	
		T_keV    = 25.0 ; [keV]  
		T	= T_keV * 1e3 * e_ / k
		Ecoord	= double(dlg_e) * 1d3 * e_
		analytical_e	= 2.0 * sqrt ( Ecoord / ( !pi * (k*T)^3 ) ) * exp ( -Ecoord / (k*T) )
		oPlot, Ecoord / e_ / 1e3,  analytical_e * nP, thick = 2.0
	
		stop

	endif
	; Particle maxwellian
	
	cdfId	= ncdf_open ( 'output/p_f_aorsa_grid.nc', /noWrite )
	glob	= ncdf_inquire ( cdfId )

	nCdf_varGet, cdfId, 'p_f_vv', p_f_vv
	nCdf_varGet, cdfId, 'p_dfduPer', p_dfduPer
	nCdf_varGet, cdfId, 'p_dfduPar', p_dfduPar
	nCdf_varGet, cdfId, 'uPerp', p_uPerp
	nCdf_varGet, cdfId, 'uPara', p_uPara
	nCdf_varGet, cdfId, 'p_enorm', p_enorm

	ncdf_close, cdfId

	u_	= mi * c^2 / ( 2.0 * e_ * p_enorm )	
	vc_	= c / sqrt ( u_ )

	if keyword_set ( lists ) then begin

		mhc_u	= mhc_v / vc_
		dlg_u	= dlg_v / vc_
		
		mhc_uPar	= mhc_data.lambda * mhc_u
		mhc_uPerp	= sqrt ( mhc_u^2 - mhc_uPar^2 )
	
		dlg_uPar	= dlg_data.lambda * dlg_u
		dlg_uPerp	= sqrt ( dlg_u^2 - dlg_uPar^2 )
	
		mhc_uHist	= histogram ( mhc_u, locat = mhc_ugrid, nBins=100 )
		dlg_uHist	= histogram ( dlg_u, locat = dlg_ugrid, nBins=100 )
	
		dlg_uHist_cyl	= hist_2d ( dlg_uPerp, dlg_uPar, bin1 = p_uPerp[1]-p_uPerp[0],bin2=abs(p_uPara[1]-p_uPara[0]), $
			min1=min(p_uPerp),max1=max(p_uPerp),min2=min(p_uPara),max2=max(p_uPara) )
	
	endif
	
	; AORSA maxwellian

	cdfId	= ncdf_open ( 'output/aorsa_num_max.nc', /noWrite )
	glob	= ncdf_inquire ( cdfId )

	ncdf_varGet, cdfId, 'f_cql', f_cql
	ncdf_varGet, cdfId, 'fnorm', fNorm
	ncdf_varGet, cdfId, 'vc_mks', vc_mks
	ncdf_varGet, cdfId, 'xkt', xkt
	ncdf_varGet, cdfId, 'Enorm', Enorm_
	ncdf_varGet, cdfId, 'u', u
	
	nCdf_varGet, cdfId, 'f_cql_cart', f_cql_cart
	nCdf_varGet, cdfId, 'dfduper', dfduper 
	nCdf_varGet, cdfId, 'dfdupar', dfdupar 
	nCdf_varGet, cdfId, 'uPerp', uPerp
	nCdf_varGet, cdfId, 'uPara', uPara
	ncdf_close, cdfId


	print, 'vc_mks: ', vc_mks, vc_
	!p.multi = [0,3,3]
	window, ySize = 900
	device, decomposed = 0
	loadct, 12	

	levels	= 3.0^fIndGen(16)*1e-3
	dlevels	= 3.0^fIndGen(8)*1e0
	dLevels	= [-reverse(dLevels[1:*]),dLevels]

	colors	= bytScl ( fIndGen ( 16 ) )

	;contour, transpose ( f_vv ), uPar_binCenters, uPerp_binCenters,$
	;	levels = levels, c_colors = colors, xrange = [-0.5,0.5], $
	;	yrange = [0, 0.5]
	contour, transpose ( f_cql_cart[*,*,0] ), uPara, uPerp,$
		levels = levels, c_colors = colors, xrange = [-0.5,0.5], $
		yrange = [0, 0.5]
	contour, transpose ( dfduPer[*,*,0] ), uPara, uPerp,$
		levels = dlevels, c_colors = colors, xrange = [-0.5,0.5], $
		yrange = [0, 0.5]
	contour, transpose ( dfduPar[*,*,0] ), uPara, uPerp,$
		levels = dlevels, c_colors = colors, xrange = [-0.5,0.5], $
		yrange = [0, 0.5]

	contour, transpose ( p_f_vv ), p_uPara, p_uPerp,$
		levels = levels, c_colors = colors, xrange = [-0.5,0.5], $
		yrange = [0, 0.5]
	contour, transpose ( p_dfduPer ), p_uPara, p_uPerp,$
		levels = dlevels, c_colors = colors, xrange = [-0.5,0.5], $
		yrange = [0, 0.5]
	contour, transpose ( p_dfduPar ), p_uPara, p_uPerp,$
		levels = dlevels, c_colors = colors, xrange = [-0.5,0.5], $
		yrange = [0, 0.5]

	plot, upara, f_cql_cart[4,*,0]
	oPlot, p_uPara, p_f_vv[4,*], color = 8*16-1, thick = 2

	plot, upara, dfduPer[4,*,0]
	oPlot, p_uPara, p_dfduPer[4,*], color = 8*16-1, thick = 2

	plot, upara, dfduPar[4,*,0]
	oPlot, p_uPara, p_dfduPar[4,*], color = 8*16-1, thick = 2


stop
end
