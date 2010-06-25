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
		view_number = view

	if not keyword_set(noLines) then $
	iContour, field, x, y, $
		c_value = levels/2, $
		rgb_indices = colors, $
		/zoom_on_resize, $
		rgb_table = 1, $
		over = id

	iContour, -field, x, y, $
		c_value = levels, $
		rgb_indices = colors, $
		/zoom_on_resize, $
		over = id, $
		/fill, $
		rgb_table = 3

	if not keyword_set(noLines) then $
	iContour, -field, x, y, $
		c_value = levels/2, $
		rgb_indices = colors, $
		/zoom_on_resize, $
		rgb_table = 3, $
		over = id 

end

pro plot_solution, oneD = oneD, full = full, $
		scale1 = scale1, scale2 = scale2, scale3_ = scale3_

	fileList = file_search ( 'solution*.nc' )
	fileListData = file_search ( 'runData*.nc' )

	xAll = 0.0
	nAll = 0
	eAlphaAll = complex (0,0)
	eAlphaAll_ = complex (0,0)
	dReAlphaAll_ = complex (0,0)
	d2ReAlphaAll_ = complex (0,0)
	d3ReAlphaAll_ = complex (0,0)
	d4ReAlphaAll_ = complex (0,0)

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
		nCdf_varGet, cdfId, 'd2rBfn_re', d2RbFn_bFn_re
		nCdf_varGet, cdfId, 'd2rBfn_im', d2RbFn_bFn_im
		nCdf_varGet, cdfId, 'd2zBfn_re', d2zbFn_bFn_re
		nCdf_varGet, cdfId, 'd2zBfn_im', d2zbFn_bFn_im
		nCdf_varGet, cdfId, 'd3rBfn_re', d3rbFn_bFn_re
		nCdf_varGet, cdfId, 'd3rBfn_im', d3rbFn_bFn_im
		nCdf_varGet, cdfId, 'd4rBfn_re', d4rbFn_bFn_re
		nCdf_varGet, cdfId, 'd4rBfn_im', d4rbFn_bFn_im
	ncdf_close, cdfId

	xx	= complex ( xx_re, xx_im )
	yy	= complex ( yy_re, yy_im )

	dRbFn_bFn	= complex ( dRbFn_bFn_re, dRbFn_bFn_im )
	dZbFn_bFn	= complex ( dZbFn_bFn_re, dZbFn_bFn_im )

	d2RbFn_bFn	= complex ( d2RbFn_bFn_re, d2RbFn_bFn_im )
	d2ZbFn_bFn	= complex ( d2ZbFn_bFn_re, d2ZbFn_bFn_im )

	d3RbFn_bFn	= complex ( d3RbFn_bFn_re, d3RbFn_bFn_im )
	d4RbFn_bFn	= complex ( d4RbFn_bFn_re, d4RbFn_bFn_im )

	nX	= n_elements ( xx[0,*] )
	nY	= n_elements ( yy[0,*] )
	nN	= n_elements ( xx[*,0] )
	nM	= n_elements ( yy[*,0] )

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

		ealpha	= complex ( ealpha_re, ealpha_im )
		ebeta	= complex ( ebeta_re, ebeta_im )
		eb	= complex ( eb_re, eb_im )

		ealphak	= complex ( ealphak_re, ealphak_im )
		ebetak	= complex ( ebetak_re, ebetak_im )
		ebk	= complex ( ebk_re, ebk_im )

		eR	= complex ( er_re, er_im )
		eTh	= complex ( et_re, et_im )
		eZ	= complex ( ez_re, ez_im )

	ncdf_close, cdfId

	; 1D catch
	; --------

	@load_colors

	if (keyword_set(oneD) ) then begin

		
		; Field solution

		range = max(abs([eAlpha,eBeta,eB]))
		range = [-range,range]
		iPlot, x, eAlpha, view_grid=[1,3], /stretch_to_fit, yrange=range
		iPlot, x, imaginary(eAlpha), over = 1

		iPlot, x, eBeta, /view_next, /stretch_to_fit, yrange=range
		iPlot, x, imaginary(eBeta), over = 1

		iPlot, x, eB, /view_next, /stretch_to_fit, yrange=range
		iPlot, x, imaginary(eB), over = 1

		; Rotated Field solution

		iPlot, x, eR, view_grid=[1,3], /stretch_to_fit, yrange=range
		iPlot, x, imaginary(eR), over = 1

		iPlot, x, eTh, /view_next, /stretch_to_fit, yrange=range
		iPlot, x, imaginary(eTh), over = 1

		iPlot, x, eZ, /view_next, /stretch_to_fit, yrange=range
		iPlot, x, imaginary(eZ), over = 1


		; Spectrum 

		iPlot, abs(eAlphak)^2, view_grid=[1,3], /stretch_to_fit, /yLog
		iPlot, abs(eBetak)^2, /view_next, /stretch_to_fit, over = 1, /yLog
		iPlot, abs(eBk)^2, /view_next, /stretch_to_fit, over = 1, /yLog


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

	  		                  d2RbFn = d2RbFn_bFn(n,i) * xx(n, i) * yy(m, j)
    		                  d2Realpha(i,j) += ealphak(n,m) * d2RbFn

		  	                  d3RbFn = d3RbFn_bFn(n,i) * xx(n, i) * yy(m, j)
    		                  d3Realpha(i,j) += ealphak(n,m) * d3RbFn

			                  d4RbFn = d4RbFn_bFn(n,i) * xx(n, i) * yy(m, j)
    		                  d4Realpha(i,j) += ealphak(n,m) * d4RbFn
						  
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
		d2EdR2_L = 2*(alpha*eAlpha[0]-eAlpha[1]*(1+alpha)+eAlpha[2])/(alpha*(alpha+1)*h1^2)

		h1 = x[nX-2]-x[nX-3]
		h2 = x[nX-1]-x[nX-2]
		alpha = h2/h1
		dEdR_R = (-1*alpha^2*eAlpha[nX-3]+1*eAlpha[nX-2]*(alpha^2-1)+eAlpha[nX-1])/(h1*(alpha+1)*alpha)
		d2EdR2_R = 2*(alpha*eAlpha[nX-3]-eAlpha[nX-2]*(1+alpha)+eAlpha[nX-1])/(alpha*(alpha+1)*h1^2)

		print, 'dEdR - ', dEdR_L, dEdR_R, dReAlpha[1], dReAlpha[nX-2]
		print, 'd2EdR2 - ', d2EdR2_L, d2EdR2_R, d2ReAlpha[1], d2ReAlpha[nX-2]

		print, 'd3EdR - ', dReAlpha[overlap], dReAlpha[nX-1-overlap]
		print, 'd4EdR - ', d2ReAlpha[overlap], d2ReAlpha[nX-1-overlap]

		print, ''
		print, 'eBet: '
		print, '-----'
		print, eBeta[1],eBeta[nX-2]

		print, ''
		print, 'eB: ', eB[1],eB[nX-2]
		print, '---'



	endif else begin

		; Field contour plot
		; ------------------

		nLevs	= 21 

		if(not keyword_set(scale1)) then $
		scale1 = max ( abs ( [ealpha] ) )
		if(not keyword_set(scale2)) then $
		scale2 = max ( abs ( [ebeta] ) )
		if(not keyword_set(scale3_)) then $
		scale3_ = max ( abs(abs ( [eb[*]] )) ) * 0.2
		print, 'Scale1: ', scale1
		print, 'Scale2: ', scale2
		print, 'Scale3: ', scale3_

		fieldPlot = 2
		iContour, id = fieldPlot, view_grid = [3,1], $
				/zoom_on_resize

		contour_field, ealpha, x, y, nLevs, scale1, $
				id = fieldPlot, view = 1, /noLines
		contour_field, ebeta, x, y, nLevs, scale2, $
				id = fieldPlot, view = 2, /noLines
		contour_field, (real_part(eb)<scale3_)>(-scale3_), x, y, nLevs, scale3_, $
				id = fieldPlot, view = 3, /noLines

		fieldPlot = 5
		iContour, id = fieldPlot, view_grid = [3,1], $
				/zoom_on_resize

		contour_field, abs(ealpha), x, y, nLevs, scale1, $
				id = fieldPlot, view = 1, /noLines
		contour_field, abs(ebeta), x, y, nLevs, scale2, $
				id = fieldPlot, view = 2, /noLines
		contour_field, abs(eb), x, y, nLevs, scale3_, $
				id = fieldPlot, view = 3, /noLines

		scale1 = max ( abs ( er ) )
		scale2 = max ( abs ( eth ) )
		scale3_ = max ( abs(abs ( ez )) ) 
		print, 'Scale1: ', scale1
		print, 'Scale2: ', scale2
		print, 'Scale3: ', scale3_

		fieldPlot = 4
		iContour, id = fieldPlot, view_grid = [3,1], $
				/zoom_on_resize

		contour_field, er, x, y, nLevs, scale1, id = fieldPlot, view = 1
		contour_field, eth, x, y, nLevs, scale2, id = fieldPlot, view = 2
		contour_field, ez, x, y, nLevs, scale3_, id = fieldPlot, view = 3


		; Spectrum contour plot
		; ---------------------

		scale1 = max ( abs ( [ealphak] ) ) * 0.8
		scale2 = max ( abs ( [ebetak] ) ) * 0.8 
		scale3_ = max ( abs ( [ebk] ) ) * 0.8 

		specPID = 3
		iContour, id = specPID, view_grid = [3,1], dimensions = [1200,300]

		contour_field, abs(ealphak),	kx, ky, nLevs, scale1, id = specPID, view = 1, /noLines
		contour_field, abs(ebetak),		kx, ky, nLevs, scale2, id = specPID, view = 2, /noLines
		contour_field, abs(ebk),		kx, ky, nLevs, scale3_, id = specPID, view = 3, /noLines

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
stop
 		scale = max ( abs ( [ealpha2_[*],ebeta2_[*],eb2_[*]] ) )
		scalePrl = max ( abs(abs ( [eb2_[*]] )) ) 

		redPID = 5
		iContour, id = redPID, view_grid = [3,1], dimensions = [1200,300]

		contour_field, ealpha2_,x,y, nLevs, scale, id = redPID, view = 1
		contour_field, ebeta2_, x, y, nLevs, scale, id = redPID, view = 2
		contour_field, eb2_, x, y, nLevs, scalePrl, id = redPID, view = 3
		endif

	endelse


	xAll = [xAll, x[overlap:nX-1-overlap]]
	eAlphaAll = [ eAlphaAll, eAlpha[overlap:nX-1-overlap]]
	nAll = [ nAll, n_elements(x[overlap:nX-1-overlap]) ]
	eAlphaAll_ = [ eAlphaAll_, eAlpha_[overlap:nX-1-overlap]]
	dReAlphaAll_ = [ dReAlphaAll_, dReAlpha[overlap:nX-1-overlap]]
	d2ReAlphaAll_ = [ d2ReAlphaAll_, d2ReAlpha[overlap:nX-1-overlap]]
	d3ReAlphaAll_ = [ d3ReAlphaAll_, d3ReAlpha[overlap:nX-1-overlap]]
	d4ReAlphaAll_ = [ d4ReAlphaAll_, d2ReAlpha[overlap:nX-1-overlap]]

	endfor

	xAll = xAll[1:*]
	eAlphaAll = eAlphaAll[1:*]
	nAll = nAll[1:*]
	eAlphaAll_ = eAlphaAll_[1:*]
	dReAlphaAll_ = dReAlphaAll_[1:*]
	d2ReAlphaAll_ = d2ReAlphaAll_[1:*]
	d3ReAlphaAll_ = d3ReAlphaAll_[1:*]
	d4ReAlphaAll_ = d4ReAlphaAll_[1:*]

	restore, '../smithe_1/soln.sav' 
	iPlot, xorig, ealphaorig, thick=2, color = blue
	iplot, xall,eAlphaAll, color=red, /over, sym_index=4



stop
end
