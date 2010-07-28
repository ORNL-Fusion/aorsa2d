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

pro plot_solution, oneD = oneD, full = full, $
		scale1 = scale1, scale2 = scale2, scale3_ = scale3_, $
		sav = sav


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

	for ii=0,n_elements(fileList)-1 do begin

		cdfId = ncdf_open ( fileListData[ii], /noWrite ) 
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
			nCdf_varGet, cdfId, 'drBfn_re', dRbFn_bFn_re
			nCdf_varGet, cdfId, 'drBfn_im', dRbFn_bFn_im
			nCdf_varGet, cdfId, 'dzBfn_re', dzbFn_bFn_re
			nCdf_varGet, cdfId, 'dzBfn_im', dzbFn_bFn_im
		ncdf_close, cdfId

		xx	= complex ( xx_re, xx_im )
		yy	= complex ( yy_re, yy_im )

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

		cdfId = ncdf_open ( fileList[ii], /noWrite ) 

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

			nCdf_varGet, cdfId, 'jouleHeating', jouleHeating 

			ealpha	= complex ( ealpha_re, ealpha_im )
			ebeta	= complex ( ebeta_re, ebeta_im )
			eb	= complex ( eb_re, eb_im )

			ealphak	= complex ( ealphak_re, ealphak_im )
			ebetak	= complex ( ebetak_re, ebetak_im )
			ebk	= complex ( ebk_re, ebk_im )

			eR	= complex ( er_re, er_im )
			eTh	= complex ( et_re, et_im )
			eZ	= complex ( ez_re, ez_im )

			jPAlpha = complex ( jAlpha_re, jAlpha_im )
			jPBeta = complex ( jBeta_re, jBeta_im )
			jPB = complex ( jB_re, jB_im )

		ncdf_close, cdfId

		; 1D catch
		; --------

		@load_colors

		if (keyword_set(oneD) ) then begin

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

			print, 'dEdR - ', dEdR_L, dEdR_R, dReAlpha[1], dReAlpha[nX-2]

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

	if keyword_set(oneD) then begin

		;restore, '../smithe_1/soln.sav' 
		;restore, '../tftr_ibw/soln.sav' 
		;restore, '../cmod_ibw_1d/soln.sav'
		;restore, '../d3d_ibw_1d/soln.sav'
		;restore, '../brambilla/soln.sav'
		;restore, '../jaeger_1/soln.sav'
		restore, '../jaeger_2/soln.sav'
		;restore, '../cynthia/soln.sav'

		iPlot, xorig, ealphaorig, thick=2, view_grid=[3,1], /stretch_to_fit, $
				dimensions = [900,300]
		iPlot, xorig, eBetaorig, thick=2, /view_next, /stretch_to_fit
		iPlot, xorig, eBorig, thick=2, /view_next, /stretch_to_fit

		for i=0,n_elements(fileList)-1 do begin

			if i mod 3 eq 0 then color = red
			if i mod 3 eq 1 then color = blue
			if i mod 3 eq 2 then color = green
			if i gt 0 then startII = total(nAll[0:i-1]) else startII=0
			stopII = startII+nAll[i]-1
			print, i, startII, stopII
			iplot, xall[startII:stopII],eAlphaAll[startII:stopII], $
					color=color, /over, sym_index=4, view_number=1
			iplot, xall[startII:stopII],eBetaAll[startII:stopII], $
					color=color, /over, sym_index=4, view_number=2
			iplot, xall[startII:stopII],eBAll[startII:stopII], $
					color=color, /over, sym_index=4, view_number=3
		endfor

	endif else begin

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
