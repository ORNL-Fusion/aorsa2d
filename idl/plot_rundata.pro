pro plot_rundata

	cdfId = ncdf_open ( 'runData.nc', /noWrite ) 
		nCdf_varGet, cdfId, 'x', x 
		nCdf_varGet, cdfId, 'y', y 
		nCdf_varGet, cdfId, 'jy_re', jy_re 
		nCdf_varGet, cdfId, 'jy_im', jy_im 
	ncdf_close, cdfId

	cdfId = ncdf_open ( 'amat.nc', /noWrite ) 
		nCdf_varGet, cdfId, 'amat_re', amat_re
		nCdf_varGet, cdfId, 'amat_im', amat_im 
	ncdf_close, cdfId

	contour, jy_im, x, y

	for modeNo=0,n_elements(amat_re[*,0])-1 do begin

		ii	= indGen ( n_elements(amat_re[0,*])/3 ) * 3
		bfn1_re	= reform ( amat_re[modeNo,ii], 4, 32 )
		bfn2_re	= reform ( amat_re[modeNo,ii+1], 4, 32 )
		bfn3_re	= reform ( amat_re[modeNo,ii+2], 4, 32 )

		bfn1_im	= reform ( amat_im[modeNo,ii], 4, 32 )
		bfn2_im	= reform ( amat_im[modeNo,ii+1], 4, 32 )
		bfn3_im	= reform ( amat_im[modeNo,ii+2], 4, 32 )


		!p.multi = [0,2,2]
		contour, bfn1_re
		contour, bfn1_im
		contour, bfn2_re
		contour, bfn2_im

		stop
	endfor

stop
end
