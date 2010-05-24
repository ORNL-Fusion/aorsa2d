pro contour_field, field, x, y, nLevs, scale, $
		id = id, view = view

	levels	= (fIndGen(nLevs)+1)/(nLevs-1) * scale * 1.1
	colors	= 256-(bytScl ( levels, top = 253 )+1)

	iContour, field, x, y, $
		c_value = levels, $
		rgb_indices = colors, $
		/fill, $
		rgb_table = 1, $
		/zoom_on_resize, $
		xTickFont_size = 26.0, $
		yTickFont_size = 26.0, $
		overplot = id, $
		view_number = view

	iContour, field, x, y, $
		c_value = levels/2, $
		rgb_indices = colors, $
		rgb_table = 1, $
		over = id

	iContour, -field, x, y, $
		c_value = levels, $
		rgb_indices = colors, $
		over = id, $
		/fill, $
		rgb_table = 3
	
	iContour, -field, x, y, $
		c_value = levels/2, $
		rgb_indices = colors, $
		rgb_table = 3, $
		over = id 

end

pro plot_solution, oneD = oneD

	cdfId = ncdf_open ( 'runData.nc', /noWrite ) 
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
		;nCdf_varGet, cdfId, 'kxsav',  kx
		;nCdf_varGet, cdfId, 'kysav', ky 
	ncdf_close, cdfId

	xx	= complex ( xx_re, xx_im )
	yy	= complex ( yy_re, yy_im )

	nX	= n_elements ( xx[0,*] )
	nY	= n_elements ( yy[0,*] )
	nN	= n_elements ( xx[*,0] )
	nM	= n_elements ( yy[*,0] )

	cdfId = ncdf_open ( 'solution.nc', /noWrite ) 

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

		ealpha	= complex ( ealpha_re, ealpha_im )
		ebeta	= complex ( ebeta_re, ebeta_im )
		eb	= complex ( eb_re, eb_im )

		ealphak	= complex ( ealphak_re, ealphak_im )
		ebetak	= complex ( ebetak_re, ebetak_im )
		ebk	= complex ( ebk_re, ebk_im )

	ncdf_close, cdfId

	; 1D catch
	; --------

	if (keyword_set(oneD) ) then begin

		iPlot, x, eAlpha, view_grid=[1,3], /stretch_to_fit
		iPlot, x, imaginary(eAlpha), /over

		iPlot, x, eBeta, /view_next, /stretch_to_fit
		iPlot, x, imaginary(eBeta), /over

		iPlot, x, eB, /view_next, /stretch_to_fit
		iPlot, x, imaginary(eB), /over

	endif else begin

	; Field contour plot
	; ------------------

	nLevs	= 21 
	scale = max ( abs ( [ealpha[*],ebeta[*],eb[*]] ) ) 
	scalePrl = max ( abs(abs ( [eb[*]] )) ) 

	fieldPlot = 2
	iContour, id = fieldPlot, view_grid = [3,1], dimensions = [1200,300]

	contour_field, ealpha, x, y, nLevs, scale, id = fieldPlot, view = 1
	contour_field, ebeta, x, y, nLevs, scale, id = fieldPlot, view = 2
	contour_field, eb, x, y, nLevs, scalePrl, id = fieldPlot, view = 3


	; Spectrum contour plot
	; ---------------------

	scale = mean ( abs ( [ealphak_re[*],ebetak_re[*],ebk_re[*]] ) ) 
	scalePar = max ( abs ( [ebk_re[*]] ) ) 

	specPID = 3
	iContour, id = specPID, view_grid = [3,1], dimensions = [1200,300]

	contour_field, abs(ealphak),	kx, ky, nLevs, scale, id = specPID, view = 1
	contour_field, abs(ebetak),		kx, ky, nLevs, scale, id = specPID, view = 2
	contour_field, abs(ebk),		kx, ky, nLevs, scalePar, id = specPID, view = 3

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

	;; Old slow version of reconstruction 

	;ealpha2_	= complexArr ( nX, nY )
	;ebeta2_	= complexArr ( nX, nY )
	;eb2_	= complexArr ( nX, nY )

	;ealpha2_[*,*]	= 0 
	;ebeta2_[*,*]	= 0 
	;eb2_[*,*]		= 0 

    ;for i = 0, nX-1 do begin
    ;	for j = 0, nY-1 do begin
    ;   		for n = nN/6, nN-nN/6-1 do begin
    ;        	for m = nM/6, nM-nM/6-1 do begin

    ;                  cexpkxky = xx(n, i) * yy(m, j)
    ;                  ealpha2_(i,j) = ealpha2_(i,j) + ealphak(n,m) * cexpkxky
    ;                  ebeta2_(i,j) = ebeta2_(i,j) + ebetak(n,m) * cexpkxky
    ;                  eb2_(i,j) = eb2_(i,j) + ebk(n,m) * cexpkxky

    ;        	endfor
    ;   	 	endfor
    ;	endfor
  	;endfor

 	;scale = max ( abs ( [ealpha2_[*],ebeta2_[*],eb2_[*]] ) )
	;scalePrl = max ( abs(abs ( [eb2_[*]] )) ) 

	;redPID = 5
	;iContour, id = redPID, view_grid = [3,1], dimensions = [1200,300]

	;contour_field, ealpha2_,x,y, nLevs, scale, id = redPID, view = 1
	;contour_field, ebeta2_, x, y, nLevs, scale, id = redPID, view = 2
	;contour_field, eb2_, x, y, nLevs, scalePrl, id = redPID, view = 3

	endelse
stop
end
