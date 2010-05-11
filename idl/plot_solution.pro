pro contour_field, field, x, y, nLevs, scale, $
		initial = initial

	levels	= (fIndGen(nLevs)+1)/(nLevs-1) * scale * 1.1
	colors	= 256-(bytScl ( levels, top = 253 )+1)

	if keyword_set ( initial ) then begin

	iContour, field, x, y, $
		c_value = levels, $
		rgb_indices = colors, $
		/fill, $
		rgb_table = 1, $
		view_grid = [3,2], $
		/stretch_to_fit, $
		/zoom_on_resize, $
		/scale_isotropic

	endif else begin

	iContour, field, x, y, $
		c_value = levels, $
		rgb_indices = colors, $
		/fill, $
		rgb_table = 1, $
		view_grid = [3,2], $
		view_next = 1, $
		/stretch_to_fit, $
		/zoom_on_resize, $
		/scale_isotropic, $
		overPlot = 0 


	endelse

	iContour, field, x, y, $
		c_value = levels/2, $
		rgb_indices = colors, $
		rgb_table = 1, $
		over = 1 

	iContour, -field, x, y, $
		c_value = levels, $
		rgb_indices = colors, $
		over = 1, $
		/fill, $
		rgb_table = 3, $
		/scale_isotropic	
	
	iContour, -field, x, y, $
		c_value = levels/2, $
		rgb_indices = colors, $
		rgb_table = 3, $
		over = 1 

end

pro plot_solution

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
		nCdf_varGet, cdfId, 'kxsav',  kxsav
		nCdf_varGet, cdfId, 'kysav', kysav 
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

	scale = max ( abs ( [ealpha[*],ebeta[*],eb[*]] ) ) 
	scalePrl = max ( abs(abs ( [eb[*]] )) ) 
	nLevs	= 21 

	;contour_field, ealpha, x, y, nLevs, scale, /initial
	;contour_field, ebeta, x, y, nLevs, scale
	;contour_field, eb, x, y, nLevs, scalePrl

	;;contour_field, imaginary(ealpha), x, y, nLevs, scale
	;;contour_field, imaginary(ebeta), x, y, nLevs, scale
	;;contour_field, imaginary(eb), x, y, nLevs, scale
	;
	;contour_field, abs(ealpha), x, y, nLevs, scale
	;contour_field, abs(ebeta), x, y, nLevs, scale
	;contour_field, abs(eb), x, y, nLevs, scalePrl

stop;	
	scale = max ( abs ( [ealphak_re[*],ebetak_re[*],ebk_re[*]] ) ) 
	scalePar = max ( abs ( [ebk_re[*]] ) ) 

	;contour_field, ealphak,	kx, ky, nLevs, scale, /initial
	;contour_field, ebetak,	kx, ky, nLevs, scale
	;contour_field, ebk,		kx, ky, nLevs, scale

	;;contour_field, imaginary(ealphak),	kx, ky, nLevs, scale 
	;;contour_field, imaginary(ebetak),	kx, ky, nLevs, scale 
	;;contour_field, imaginary(ebk),		kx, ky, nLevs, scale 

	;contour_field, abs(ealphak),	kx, ky, nLevs, scale 
	;contour_field, abs(ebetak),	kx, ky, nLevs, scale 
	;contour_field, abs(ebk),		kx, ky, nLevs, scale 

stop

	; Reconstruct the fields using only a specific set of
	; basis vectors
	; ---------------------------------------------------

	ealpha_	= complexArr ( nX, nY )
	ebeta_	= complexArr ( nX, nY )
	eb_	= complexArr ( nX, nY )

    for i = 0, nX-1 do begin
    	for j = 0, nY-1 do begin
       		for n = nN/4, nN-nN/4-1 do begin
            	for m = nM/4, nM-nM/4-1 do begin

                      cexpkxky = xx(n, i) * yy(m, j)
                      ealpha_(i,j) = ealpha_(i,j) + ealphak(n,m) * cexpkxky

            	endfor
       	 	endfor
    	endfor
  	endfor

 	scale = max ( abs ( [ealpha_[*],ebeta_[*],eb_[*]] ) ) 
	contour_field, ealpha_, x, y, nLevs, scale, /initial
	;contour_field, ebeta, x, y, nLevs, scale
	;contour_field, eb, x, y, nLevs, scalePrl


stop
end
