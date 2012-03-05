
pro dlg_deriv1, x, y, dydx = dydx, dy2dx2 = dy2dx2 

		h	= x[1]-x[0]
		dy2dx2	= y*0
		dy2dx2[0]	= ( 2*y[0]-5*y[1]+4*y[2]-1*y[3] ) / h^2
		dy2dx2[-1]	= ( 2*y[-1]-5*y[-2]+4*y[-3]-1*y[-4] ) / h^2
	
		for i=1,n_elements(x[*])-2 do begin
			dy2dx2[i]	= ( 1*y[i-1]-2*y[i]+1*y[i+1] ) / h^2
		endfor

		dydx	= y*0
		dydx[0]	= ( -3*y[0]+4*y[1]-1*y[2] ) / (2*h)
		dydx[-1]	= ( 3*y[-1]-4*y[-2]+1*y[-3] ) / (2*h)
	
		for i=1,n_elements(x[*])-2 do begin
			dydx[i]	= ( -1*y[i-1]+1*y[i+1] ) / (2*h)
		endfor

end



pro contour_field, field, x, y, nLevs, scale, $
		id = id, view = view, log = log, nolines = nolines

	levels	= (fIndGen(nLevs)+1)/(nLevs-1) * scale * 1.1
	colors	= 256-(bytScl ( levels, top = 253 )+1)

	if keyword_set ( log ) then begin

		levels	= 2.0^fIndGen(nLevs) * 1e-2
		colors	= 256-(bytScl ( findGen(nLevs), top = 253 )+1)

	endif

	iContour, field, x, y, $
		c_value = levels, $
		rgb_indices = colors, $
		/fill, $
		rgb_table = 1, $
		/zoom_on_resize, $
		xTickFont_size = 26.0, $
		yTickFont_size = 26.0, $
		overplot = id, $
		/iso, $
		view_number = view, /irreg

	if not keyword_set(noLines) then $
	iContour, field, x, y, $
		c_value = levels/2, $
		rgb_indices = colors, $
		/zoom_on_resize, $
		rgb_table = 1, $
		over = id, /irreg

	iContour, -field, x, y, $
		c_value = levels, $
		rgb_indices = colors, $
		/zoom_on_resize, $
		over = id, $
		/fill, $
		rgb_table = 3, /irreg

	if not keyword_set(noLines) then $
	iContour, -field, x, y, $
		c_value = levels/2, $
		rgb_indices = colors, $
		/zoom_on_resize, $
		rgb_table = 3, $
		over = id, /irreg

end

pro ar2_plot_solution, full = full, $
		scale1 = scale1, scale2 = scale2, scale3_ = scale3_, $
		sav = sav

	@constants

	; Load data from all grids
	; ------------------------

	fileList = file_search ( 'solution*.nc' )
	fileListData = file_search ( 'runData*.nc' )

	xAll = 0.0
	yAll = 0.0
	nAll = 0
	eAlphaAll = complex (0,0)
	eBetaAll = complex (0,0)
	eBAll = complex (0,0)

	for iii=0,n_elements(fileList)-1 do begin

		cdfId = ncdf_open ( fileListData[iii], /noWrite ) 
			nCdf_varGet, cdfId, 'nPhi', nPhi 
			nCdf_varGet, cdfId, 'freq', freq 
			nCdf_varGet, cdfId, 'capR', x 
			nCdf_varGet, cdfId, 'y', y 
			nCdf_varGet, cdfId, 'jr_re', jr_re 
			nCdf_varGet, cdfId, 'jr_im', jr_im 
			nCdf_varGet, cdfId, 'jt_re', jt_re 
			nCdf_varGet, cdfId, 'jt_im', jt_im 
			nCdf_varGet, cdfId, 'jz_re', jz_re 
			nCdf_varGet, cdfId, 'jz_im', jz_im 
			nCdf_varGet, cdfId, 'bmod', bmod
			nCdf_varGet, cdfId, 'brU', brU
			nCdf_varGet, cdfId, 'btU', btU
			nCdf_varGet, cdfId, 'bzU', bzU
			nCdf_varGet, cdfId, 'densitySpec', densitySpec
			nCdf_varGet, cdfId, 'tempSpec', tempSpec
			nCdf_varGet, cdfId, 'xx_re', xx_re 
			nCdf_varGet, cdfId, 'xx_im', xx_im
			nCdf_varGet, cdfId, 'yy_re', yy_re 
			nCdf_varGet, cdfId, 'yy_im', yy_im
			nCdf_varGet, cdfId, 'drBfn_re', dRbFn_bFn_re
			nCdf_varGet, cdfId, 'drBfn_im', dRbFn_bFn_im
			nCdf_varGet, cdfId, 'dzBfn_re', dzbFn_bFn_re
			nCdf_varGet, cdfId, 'dzBfn_im', dzbFn_bFn_im
		ncdf_close, cdfId

		xx	= complex ( xx_re, xx_im )
		yy	= complex ( yy_re, yy_im )

		jA_r = complex ( jr_re, jr_im )
		jA_t = complex ( jt_re, jt_im )
		jA_z = complex ( jz_re, jz_im )

		dRbFn_bFn	= complex ( dRbFn_bFn_re, dRbFn_bFn_im )
		dZbFn_bFn	= complex ( dZbFn_bFn_re, dZbFn_bFn_im )

		nX	= n_elements ( xx[0,*] )
		nY	= n_elements ( yy[0,*] )
		nN	= n_elements ( xx[*,0] )
		nM	= n_elements ( yy[*,0] )

		x2D	= rebin ( x, nX, nY )
		if nY gt 1 then begin
			y2D = transpose(rebin ( y, nY, nX ))
		endif else begin
			y2D = fltArr(1)+y
		endelse

	
		cdfId = ncdf_open ( fileList[iii], /noWrite ) 

			nCdf_varGet, cdfId, 'ealpha_re', ealpha_re 
			nCdf_varGet, cdfId, 'ebeta_re', ebeta_re 
			nCdf_varGet, cdfId, 'eB_re', eB_re 
			nCdf_varGet, cdfId, 'ealpha_im', ealpha_im 
			nCdf_varGet, cdfId, 'ebeta_im', ebeta_im 
			nCdf_varGet, cdfId, 'eB_im', eB_im 

			nCdf_varGet, cdfId, 'ealphak_re', ealphak_re 
			nCdf_varGet, cdfId, 'ebetak_re', ebetak_re 
			nCdf_varGet, cdfId, 'eBk_re', eBk_re 
			nCdf_varGet, cdfId, 'ealphak_im', ealphak_im 
			nCdf_varGet, cdfId, 'ebetak_im', ebetak_im 
			nCdf_varGet, cdfId, 'eBk_im', eBk_im 

			nCdf_varGet, cdfId, 'er_re', er_re 
			nCdf_varGet, cdfId, 'et_re', et_re 
			nCdf_varGet, cdfId, 'ez_re', ez_re 
			nCdf_varGet, cdfId, 'er_im', er_im 
			nCdf_varGet, cdfId, 'et_im', et_im 
			nCdf_varGet, cdfId, 'ez_im', ez_im 

			nCdf_varGet, cdfId, 'jalpha_re', jalpha_re 
			nCdf_varGet, cdfId, 'jbeta_re', jbeta_re 
			nCdf_varGet, cdfId, 'jB_re', jB_re 
			nCdf_varGet, cdfId, 'jalpha_im', jalpha_im 
			nCdf_varGet, cdfId, 'jbeta_im', jbeta_im 
			nCdf_varGet, cdfId, 'jB_im', jB_im 

			nCdf_varGet, cdfId, 'jP_r_re', jPr_re 
			nCdf_varGet, cdfId, 'jP_t_re', jPt_re 
			nCdf_varGet, cdfId, 'jP_z_re', jPz_re 
			nCdf_varGet, cdfId, 'jP_r_im', jPr_im 
			nCdf_varGet, cdfId, 'jP_t_im', jPt_im 
			nCdf_varGet, cdfId, 'jP_z_im', jPz_im 

			nCdf_varGet, cdfId, 'jouleHeating', jouleHeating 

			ealpha	= complex ( ealpha_re, ealpha_im )
			ebeta	= complex ( ebeta_re, ebeta_im )
			eb	= complex ( eb_re, eb_im )

			ealphak	= complex ( ealphak_re, ealphak_im )
			ebetak	= complex ( ebetak_re, ebetak_im )
			ebk	= complex ( ebk_re, ebk_im )

			e_r	= complex ( er_re, er_im )
			e_t	= complex ( et_re, et_im )
			e_z	= complex ( ez_re, ez_im )

			jP_r	= complex ( jPr_re, jPr_im )
			jP_t	= complex ( jPt_re, jPt_im )
			jP_z	= complex ( jPz_re, jPz_im )

			jPAlpha = complex ( jAlpha_re, jAlpha_im )
			jPBeta = complex ( jBeta_re, jBeta_im )
			jPB = complex ( jB_re, jB_im )

		ncdf_close, cdfId

		; 1D catch
		; --------

		@load_colors

		if size(y,/dim) eq 0 then begin

			if keyword_set(full) then begin

				ealpha_	= complexArr ( nX, nY )
				dRealpha	= complexArr ( nX, nY )
				d2Realpha	= complexArr ( nX, nY )
				d3Realpha	= complexArr ( nX, nY )
				d4Realpha	= complexArr ( nX, nY )

				ealpha_[*,*]	= 0 
				dRealpha[*,*]	= 0 
				d2Realpha[*,*]	= 0 
				d3Realpha[*,*]	= 0 
				d4Realpha[*,*]	= 0 

    			for i = 0, nX-1 do begin
    				for j = 0, nY-1 do begin
    			   		for n = 0, nN-1 do begin
    			        	for m = 0, nM-1 do begin

						          bfn = xx(n, i) * yy(m, j)
    			                  ealpha_(i,j) += ealphak(n,m) * bfn 

    			                  dRbFn = dRbFn_bFn(n,i) * xx(n, i) * yy(m, j)
    			                  dRealpha(i,j) += ealphak(n,m) * dRbFn 

    			        	endfor
    			   	 	endfor
    				endfor
  				endfor
			endif

			overlap = 2

			print, ''
			print, 'eAlpha: '
			print, '------- '
			print, 'E - ', eAlpha[1],eAlpha[nX-2]

			h1 = x[1]-x[0]
			h2 = x[2]-x[1]
			alpha = h2/h1
			dEdR_L = (-1*alpha^2*eAlpha[0]+1*eAlpha[1]*(alpha^2-1)+eAlpha[2])/(h1*(alpha+1)*alpha)

			h1 = x[nX-2]-x[nX-3]
			h2 = x[nX-1]-x[nX-2]
			alpha = h2/h1
			dEdR_R = (-1*alpha^2*eAlpha[nX-3]+1*eAlpha[nX-2]*(alpha^2-1)+eAlpha[nX-1])/(h1*(alpha+1)*alpha)

			print, 'dEdR - ', dEdR_L, dEdR_R;, dReAlpha[1], dReAlpha[nX-2]

			print, ''
			print, 'eBet: '
			print, '-----'
			print, eBeta[1],eBeta[nX-2]

			print, ''
			print, 'eB: ', eB[1],eB[nX-2]
			print, '---'



		endif else begin

			;; Reconstruct the fields using only a specific set of
			;; basis vectors
			;; ---------------------------------------------------

			;ealphak_	= ealphak * 0
			;ebetak_	= ebetak * 0
			;ebk_	= ebk * 0

			;for n = nN/6, nN-nN/6-1 do begin
    		;   	for m = nM/6, nM-nM/6-1 do begin

    		;     ealphak_(n,m) = ealphak(n,m) 
    		;     ebetak_(n,m) = ebetak(n,m) 
    		;     ebk_(n,m) = ebk(n,m) 

    		;   	endfor
    		;endfor

			;; While this is not entirely correct due to the 
			;; missing symmetry in the number of left/right or
			;; up/down basis functions (due to requiring a square
			;; matrix) but it's close enough

			;ealpha_	= fft ( ealphak, /inv, /center )
			;ebeta_	= fft ( ebetak, /inv, /center )
			;eb_	= fft ( ebk, /inv, /center )

 			;scale = max ( abs ( [ealpha_[*],ebeta_[*],eb_[*]] ) ) 
			;scalePrl = max ( abs(abs ( [eb_[*]] )) ) 

			;redPID = 4
			;iContour, id = redPID, view_grid = [3,1], dimensions = [1200,300]

			;contour_field, ealpha_,x,y, nLevs, scale, id = redPID, view = 1
			;contour_field, ebeta_, x, y, nLevs, scale, id = redPID, view = 2
			;contour_field, eb_, x, y, nLevs, scalePrl, id = redPID, view = 3

			if keyword_set(full) then begin
			; Old slow version of reconstruction 

			ealpha2_	= complexArr ( nX, nY )
			ebeta2_	= complexArr ( nX, nY )
			eb2_	= complexArr ( nX, nY )

			ealpha2_[*,*]	= 0 
			ebeta2_[*,*]	= 0 
			eb2_[*,*]		= 0 

    		for i = 0, nX-1 do begin
					print, i, nX
    			for j = 0, nY-1 do begin
    		   		for n = 0, nN-1 do begin
    		        	for m = 0, nM-1 do begin

    		                  cexpkxky = dRbFn_bFn(n,i) * xx(n, i) * yy(m, j)
    		                  ealpha2_(i,j) += ealphak(n,m) * cexpkxky
    		                  ;ebeta2_(i,j) += ebetak(n,m) * cexpkxky
    		                  ;eb2_(i,j) += ebk(n,m) * cexpkxky
							  
							  ;eb2_(i,j) += matrix_multiply ( ebk[n,*] , yy[*,j] ) * xx[n,i]

    		        	endfor
    		   	 	endfor
    			endfor
  			endfor

 			scale = max ( abs ( [ealpha2_[*],ebeta2_[*],eb2_[*]] ) )
			scalePrl = max ( abs(abs ( [eb2_[*]] )) ) 

			redPID = 5
			iContour, id = redPID, view_grid = [3,1], dimensions = [1200,300]

			contour_field, ealpha2_,x,y, nLevs, scale, id = redPID, view = 1
			contour_field, ebeta2_, x, y, nLevs, scale, id = redPID, view = 2
			contour_field, eb2_, x, y, nLevs, scalePrl, id = redPID, view = 3
			endif

		endelse


		xAll = [xAll, x2D[*]]
		yAll = [yAll, y2D[*]]
		eAlphaAll = [ eAlphaAll, eAlpha[*]]
		eBetaAll = [ eBetaAll, eBeta[*]]
		eBAll = [ eBAll, eB[*]]
		nAll = [ nAll, n_elements(x2D[*]) ]

		if (size(eAlphaAll_,/dim))[0] eq 0 then begin
			x2D_	= x2D
			eAlphaAll_ = eAlpha
			eBetaAll_ = eBeta
			eBAll_ = eB
		endif else begin
			x2D_	= [ x2D_, x2D ]
			eAlphaAll_ = [ eAlphaAll_, eAlpha ]
			eBetaAll_ = [ eBetaAll_, eBeta ]
			eBAll_ = [ eBAll_, eB ]
		endelse

	endfor

	xAll = xAll[1:*]
	yAll = yAll[1:*]
	eAlphaAll = eAlphaAll[1:*]
	eBetaAll = eBetaAll[1:*]
	eBAll = eBAll[1:*]
	nAll = nAll[1:*]


	; Plot all grid solutions
	; -----------------------

	if size(y,/dim) eq 0 then begin

		iix = n_elements(x)/2
		wrf = freq * 2 * !pi
		wpe = sqrt ( densitySpec[iix,0,0] * e^2 / (me * e0 ) )
		vThe = sqrt ( 2 * tempSpec[iix,0,0]*e / me )
		lambda_de =  sqrt(e0 * tempSpec[iix,0,0]*e/(e^2*densitySpec[iix,0,0]))
		print, 'w/wpe: ', wrf / wpe
		brambillaNumber = 0.25
		print, 'kPar: [1/m]', brambillaNumber * wpe / vThe
		print, 'lambda: [m]', 2*!pi/(brambillaNumber* wpe / vThe)
		print, 'lambda_de: [m]', lambda_de

		; Check with cold plasma theory

		iiMaxE	= where(abs(e_r)^2 eq max(abs(e_r)^2))

		print, 'Er[iiMaxE]: ', e_r[iiMaxE]
		print, 'Jp[iiMaxE]: ', jP_r[iiMaxE]
		print, 'Ja[iiMaxE]: ', jA_r[iiMaxE]

		print, 'Jp/Ja theory: ', (wpe^2/wrf^2)/(1-wpe^2/wrf^2)
		print, 'Jp/Ja actual: ', jP_r[iiMaxE]/jA_r[iiMaxE]

		print, 'Jp/E theory: ', II*wrf*e0*wpe^2/wrf^2
		print, 'Jp/E actual: ', jP_r[iiMaxE]/e_r[iiMaxE]

		print, 'Ja/E theory: ', II*wrf*e0*(1-wpe^2/wrf^2)
		print, 'Ja/E actual: ', jA_r[iiMaxE]/e_r[iiMaxE]


		eRange = max(abs([e_r,e_t,e_z]))
		p_r = plot ( x, e_r, layout=[1,3,1],$
				title='Er',yRange=[-eRange,eRange],ytitle='Er [V/m]',name='Re',window_title='aorsa')
		p_i = plot ( x, imaginary(e_r), color='red',/over,name='Im')
		l = legend(target=[p_r,p_i],position=[0.98,0.95],/norm,font_size=10,horizontal_align='RIGHT')

		p_r = plot ( x, e_t, layout=[1,3,2],/current,$
				title='Et',yRange=[-eRange,eRange],ytitle='Et [V/m]',name='Re')
		p_i = plot ( x, imaginary(e_t), color='red',/over,name='Im')
		l = legend(target=[p_r,p_i],position=[0.98,0.62],/norm,font_size=10,horizontal_align='RIGHT')

		p_r = plot ( x, e_z, layout=[1,3,3],/current,$
				title='Ez',yRange=[-eRange,eRange],ytitle='Ez [V/m]',name='Re')
		p_i = plot ( x, imaginary(e_z), color='red',/over,name='Im')
		l = legend(target=[p_r,p_i],position=[0.98,0.28],/norm,font_size=10,horizontal_align='RIGHT')

		;; Check solution against wave equation

		;;@constants

		;r	= x
		;f	= freq
		;w 	= 2 * !pi * f

		;jP_r_total	= total ( jP_r, 3 )
		;jP_t_total	= total ( jP_t, 3 )
		;jP_z_total	= total ( jP_z, 3 )

		;jP_r_total_	= jP_r_total
		;jP_t_total_	= jP_t_total
		;jP_z_total_	= jP_z_total
		;
		;e_r_ = e_r
		;e_t_ = e_t
		;e_z_ = e_z

		;cdfId = ncdf_open ( 'sigma001.nc', /noWrite ) 

		;	nCdf_varGet, cdfId, 'sigma_re', sigma_re 
		;	nCdf_varGet, cdfId, 'sigma_im', sigma_im

		;ncdf_close, cdfId

		;sigma	= complex ( sigma_re, sigma_im )

		;;restore, '~/code/rsfwc_1D/idl/solutionVals.sav'

		;dlg_deriv1, r, e_r, dydx = de_r, dy2dx2 = d2e_r
		;dlg_deriv1, r, e_t, dydx = de_t, dy2dx2 = d2e_t
		;dlg_deriv1, r, e_z, dydx = de_z, dy2dx2 = d2e_z

		;term1_r = -nPhi * ( nPhi * e_r + II * ( e_t + r * de_t ) ) / r^2
		;term2_r = w^2 / c^2 * ( e_r + II / ( w * e0 ) * jP_r_total )
		;term3_r = -II * w * u0 * jA_r

		;term1_t = 1 / r^2 * ( II * nPhi * e_r - e_t + $
		;	r * ( -II * nPhi * de_r + de_t + r * d2e_t ))
		;term2_t = w^2 / c^2 * ( e_t + II / ( w * e0 ) * jP_t_total )
		;term3_t = -II * w * u0 * jA_t

		;term1_z = (-nPhi^2 * e_z + r * ( de_z + r * d2e_z ) ) / r^2
		;term2_z = w^2 / c^2 * ( e_z + II / ( w * e0 ) * jP_z_total )
		;term3_z = -II * w * u0 * jA_z

		;;term1_r = -nPhi * ( nPhi * e_r + II * r^2 * de_t ) / r^4
		;;term2_r = w^2 / c^2 * ( e_r + II / ( w * e0 ) * jP_r_total )
		;;term3_r = -II * w * u0 * jA_r

		;;term1_t = -1 / r^4 * ( nPhi * ( -2 * II * r ) * e_r $
		;;			+ r^2*e_t $
		;;			- r^2*(-II*nPhi*de_r+r*de_t+r^2*d2e_t))
		;;term2_t = w^2 / c^2 * ( e_t + II / ( w * e0 ) * jP_t_total )
		;;term3_t = -II * w * u0 * jA_t

		;;term1_z = (-nPhi^2*e_z $
		;;	   + r^2*(r*de_z+r^2*d2e_z)) / r^4
		;;term2_z = w^2 / c^2 * ( e_z + II / ( w * e0 ) * jP_z_total )
		;;term3_z = -II * w * u0 * jA_z


		;p = plot ( term1_r, layout=[1,3,1] )
		;p = plot ( imaginary(term1_r), color='red',/over )
		;p = plot ( -term2_r+term3_r, lineStyle='--', thick = 4, transp = 50, /over )
		;p = plot ( imaginary(-term2_r+term3_r), $
		;		color='r',lineStyle='--' , thick = 4, transp = 50, /over )

		;p = plot ( term1_t, layout=[1,3,2], /current )
		;p = plot ( imaginary(term1_t), color='red',/over )
		;p = plot ( -term2_t+term3_t, lineStyle='--', thick = 4, transp = 50, /over )
		;p = plot ( imaginary(-term2_t+term3_t), color='r',lineStyle='--', thick = 4, transp = 50 ,/over )

		;p = plot ( term1_z, layout=[1,3,3], /current )
		;p = plot ( imaginary(term1_z), color='red',/over )
		;p = plot ( -term2_z+term3_z, lineStyle='--', thick = 4, transp = 50, /over )
		;p = plot ( imaginary(-term2_z+term3_z), color='r',lineStyle='--', thick = 4, transp = 50 ,/over )


		nSpec	= n_elements ( jp_r[0,0,*] )

		; jDotE
		; -----

		jDotE	= jp_r[*,*,0] * conj(e_r) $
					+ jp_t[*,*,0] * conj(e_t) $
					+ jp_z[*,*,0] * conj(e_z)

		p = plot (x,jDotE,color='b',thick=3,transp=50,$
				title='J dot E',name='jDote_0',font_size=10,$
				layout=[1,2,1],window_title='aorsa')

		p_array = !NULL
		p_array = [p_array,p]

		for s=1,nSpec-1 do begin

			jDotE	= jp_r[*,*,s] * conj(e_r) $
					+ jp_t[*,*,s] * conj(e_t) $
					+ jp_z[*,*,s] * conj(e_z)

			p = plot ( x, jDotE, /over,name='jDotE_'+strTrim(string(s),2) )
			p_array = [p_array,p]

		endfor

		; jAnt
		; ---

	   	l = legend(target=p_array,position=[0.8,0.9],/norm,font_size=10)

		pr = plot (x,jA_r,color='b',thick=3,transp=50,$
				title='jAnt',name='jAnt_r',font_size=10,$
				layout=[1,2,2],/current)
		pt = plot (x,jA_t,color='r',thick=3,transp=50,$
				name='jAnt_t',/over)
		pz = plot (x,jA_z,color='g',thick=3,transp=50,$
				name='jAnt_z',/over)

	   	l = legend(target=[pr,pt,pz],position=[0.8,0.4],/norm,font_size=10)

		; jP
		; --

		jpRange = max(abs([abs(jp_r),abs(jp_t),abs(jp_z)]))

		s = 0
		p_array = !NULL
		p = plot (x,jp_r[*,0,s],thick=2,transp=50,$
				title='jPr',name='jPr_re_'+strTrim(string(s),2),font_size=10,$
				layout=[1,3,1],yRange=[-jpRange,jpRange],window_title='aorsa')
		p_array = [p_array,p]
		p = plot (x,imaginary(jp_r[*,0,s]),/over,name='jPr_im_'+strTrim(string(s),2),color='r',thick=2)
		p_array = [p_array,p]
		for s=1,nSpec-1 do begin
			p = plot ( x, jp_r[*,0,s],/over,name='jPr_re_'+strTrim(string(s),2),thick=2,lineStyle=s)
			p_array = [p_array,p]
			p = plot ( x, imaginary(jp_r[*,0,s]),/over,name='jPr_im_'+strTrim(string(s),2),thick=2,lineStyle=s,color='r')
			p_array = [p_array,p]
		endfor
	   	l = legend(target=p_array,position=[0.99,0.9],/norm,font_size=10,horizontal_align='RIGHT')

		s = 0
		p_array = !null
		p = plot (x,jp_t[*,0,s],thick=2,transp=50,$
				title='jPt',name='jPt_re_'+strtrim(string(s),2),font_size=10,$
				layout=[1,3,2],yrange=[-jPrange,jPrange],/current)
		p_array = [p_array,p]
		p = plot (x,imaginary(jp_t[*,0,s]),/over,name='jPt_im_'+strtrim(string(s),2),color='r',thick=2)
		p_array = [p_array,p]
		for s=1,nspec-1 do begin
			p = plot ( x, jp_t[*,0,s],/over,name='jPt_re_'+strtrim(string(s),2),thick=2,linestyle=s)
			p_array = [p_array,p]
			p = plot ( x, imaginary(jp_t[*,0,s]),/over,name='jPt_im_'+strtrim(string(s),2),thick=2,linestyle=s,color='r')
			p_array = [p_array,p]
		endfor
	   	l = legend(target=p_array,position=[0.99,0.6],/norm,font_size=10,horizontal_align='right')

		s = 0
		p_array = !null
		p = plot (x,jP_z[*,0,s],thick=2,transp=50,$
				title='jPz',name='jPz_re_'+strtrim(string(s),2),font_size=10,$
				layout=[1,3,3],yrange=[-jPrange,jPrange],/current)
		p_array = [p_array,p]
		p = plot (x,imaginary(jP_z[*,0,s]),/over,name='jPz_im_'+strtrim(string(s),2),color='r',thick=2)
		p_array = [p_array,p]
		for s=1,nspec-1 do begin
			p = plot ( x, jP_z[*,0,s],/over,name='jPz_re_'+strtrim(string(s),2),thick=2,linestyle=s)
			p_array = [p_array,p]
			p = plot ( x, imaginary(jP_z[*,0,s]),/over,name='jPz_im_'+strtrim(string(s),2),thick=2,linestyle=s,color='r')
			p_array = [p_array,p]
		endfor
	   	l = legend(target=p_array,position=[0.99,0.25],/norm,font_size=10,horizontal_align='right')


		; Jp.E vs Jant.E integrals

		jP_r_total	= total ( jP_r, 3 )
		jP_t_total	= total ( jP_t, 3 )
		jP_z_total	= total ( jP_z, 3 )

		jPDotE	= -0.5 * real_part ( conj(jP_r) * e_r $
				+ conj(jP_t) * e_t $
				+ conj(jP_z) * e_z )

		jADotE	= -0.5 * real_part ( conj(jA_r) * e_r $
				+ conj(jA_t) * e_t $
				+ conj(jA_z) * e_z )

		intJPDotE = total ( jPDotE*2*!pi*x )
		intJADotE = total ( jADotE*2*!pi*x )

		print, 'intJPDotE: ', intJPDotE
		print, 'intJADotE: ', intJADotE

	endif else begin


		nLevs = 11
		scale = 40
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		c = contour ( jPAlpha[*,*,0], x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill )
		c = contour ( -jPAlpha[*,*,0], x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
		stop

		; Field contour plot
		; ------------------

		nLevs	= 21 

		if(not keyword_set(scale1)) then $
		scale1 = max ( abs ( [ealpha] ) )
		print, 'Scale1: ', scale1
		if(not keyword_set(scale3)) then $
		scale3 = max ( abs ( [eB] ) )
		print, 'Scale3: ', scale3


		nLevs = 20
		levels	= (fIndGen(nLevs)+1)/(nLevs-1) * scale1 * 0.5
		colors	= 256-(bytScl ( levels, top = 253 )+1)

		window, 0
		device, decomposed = 0
		!p.background = 255
		loadct, 3
		contour, ealphaAll, xAll, yAll, $
			irreg = 1, $
			levels = levels, $
			c_colors = colors, $
			color = 0, /fill
		loadct, 1
		contour, -ealphaAll, xAll, yAll, $
			irreg = 1, $
			levels = levels, $
			c_colors = colors, $
			/over, color = 0, /fill

		levels	= (fIndGen(nLevs)+1)/(nLevs-1) * scale3 * 0.5
		colors	= 256-(bytScl ( levels, top = 253 )+1)

		window, 1
		loadct, 3
		contour, eBAll, xAll, yAll, $
			irreg = 1, $
			levels = levels, $
			c_colors = colors, $
			color = 0, /fill
		loadct, 1
		contour, -eBAll, xAll, yAll, $
			irreg = 1, $
			levels = levels, $
			c_colors = colors, $
			/over, color = 0, /fill


	endelse

	if keyword_set(sav) then begin

		xorig = x
		eAlphaOrig = eAlpha
		eBetaOrig = eBeta
		eBOrig = eB

		save, xorig, eAlphaOrig, eBetaOrig, eBOrig, $
				file = 'soln.sav'
	endif

	stop
end
