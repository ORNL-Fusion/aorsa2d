pro plot_sigma4d

	fileList = file_search ( 'sigma*.nc' )
	fileListData = file_search ( 'runData*.nc' )

	for ii=0,n_elements(fileList)-1 do begin

		cdfId = ncdf_open ( fileList[ii], /noWrite ) 
			nCdf_varGet, cdfId, 'sigma_re', sigma_re
			nCdf_varGet, cdfId, 'sigma_im', sigma_im
		ncdf_close, cdfId

		cdfId = ncdf_open ( fileListData[ii], /noWrite ) 
			nCdf_varGet, cdfId, 'capR', x 
			nCdf_varGet, cdfId, 'y', y 
		ncdf_close, cdfId

		nX = n_elements ( x )
		nY = n_elements ( y )

		sigma = complex ( sigma_re, sigma_im )

		n = 0
		m = 0
		data = reform ( sigma[*,*,n,m,2,2], nX, nY )

	endfor
stop
end
