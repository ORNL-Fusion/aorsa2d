pro plot_rundata

	cdfId = ncdf_open ( 'runData.nc', /noWrite ) 

	nCdf_varGet, cdfId, 'x', x 
	nCdf_varGet, cdfId, 'y', y 
	nCdf_varGet, cdfId, 'jy_re', jy_re 
	nCdf_varGet, cdfId, 'jy_im', jy_im 

	ncdf_close, cdfId


	contour, jy_im, x, y
stop

end
