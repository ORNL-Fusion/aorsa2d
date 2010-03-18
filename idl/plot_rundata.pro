pro plot_rundata

	cdfId = ncdf_open ( 'runData.nc', /noWrite ) 
		nCdf_varGet, cdfId, 'x', x 
		nCdf_varGet, cdfId, 'y', y 
		nCdf_varGet, cdfId, 'jy_re', jy_re 
		nCdf_varGet, cdfId, 'jy_im', jy_im 
		nCdf_varGet, cdfId, 'bmod', bmod
		nCdf_varGet, cdfId, 'bx', bx
		nCdf_varGet, cdfId, 'by', by
		nCdf_varGet, cdfId, 'bz', bz

	ncdf_close, cdfId

	cdfId = ncdf_open ( 'amat.nc', /noWrite ) 
		nCdf_varGet, cdfId, 'amat_re', amat_re
		nCdf_varGet, cdfId, 'amat_im', amat_im 
		nCdf_varGet, cdfId, 'nPtsX', nPtsX 
		nCdf_varGet, cdfId, 'nPtsY', nPtsY 
		nCdf_varGet, cdfId, 'nModesX', nModesX 
		nCdf_varGet, cdfId, 'nModesY', nModesY 
		nCdf_varGet, cdfId, 'xx_re', xx_re 
		nCdf_varGet, cdfId, 'xx_im', xx_im
		nCdf_varGet, cdfId, 'yy_re', yy_re 
		nCdf_varGet, cdfId, 'yy_im', yy_im 
		xx	= complex ( xx_re, xx_im )
		yy	= complex ( yy_re, yy_im )
	ncdf_close, cdfId

	contour, jy_im, x, y
stop
	for modeNo=0,n_elements(amat_re[*,0])-1 do begin

		ii	= indGen (nPtsX*nPtsY) * 3
		bfn1_re	= reform ( amat_re[ii+0,modeNo], nPtsY, nPtsX )
		bfn2_re	= reform ( amat_re[ii+1,modeNo], nPtsY, nPtsX  )
		bfn3_re	= reform ( amat_re[ii+2,modeNo], nPtsY, nPtsX  )

		bfn1_im	= reform ( amat_im[ii+0,modeNo], nPtsY, nPtsX  )
		bfn2_im	= reform ( amat_im[ii+1,modeNo], nPtsY, nPtsX  )
		bfn3_im	= reform ( amat_im[ii+2,modeNo], nPtsY, nPtsX  )


		!p.multi = [0,2,2]
		contour, bfn1_re
		contour, bfn1_im
		contour, bfn2_re
		contour, bfn2_im

		stop
	endfor

stop
end
