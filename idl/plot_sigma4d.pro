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
		nN = n_elements ( sigma_re[0,0,*,0,0,0,0] )
		nM = n_elements ( sigma_re[0,0,0,*,0,0,0] )
		nS = n_elements ( sigma_re[0,0,0,0,0,0,*] )

		sigma = complex ( sigma_re, sigma_im )

		n = 0
		m = 0
		data = reform ( sigma[*,*,n,m,2,2,0], nX, nY )


		;for s=0,nS-1 do begin

		;for i=0,nX-1 do begin
		;	for j=0,nY-1 do begin
		;	 	for n=0,nN-1 do begin
		;			for m=0,nM-1 do begin

		;				print, i,j,n,m,s,(sigma[i,j,n,m,*,0,s])[*]

		;			endfor
		;		endfor
		;	endfor
		;endfor

		;endfor


	endfor
stop
end
