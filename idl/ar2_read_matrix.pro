pro ar2_read_matrix

    FileName = 'matrix.nc'

	cdfId = ncdf_open ( FileName, /noWrite ) 
		nCdf_varGet, cdfId, 'reA', reA
		nCdf_varGet, cdfId, 'imA', imA 
    nCdf_close, cdfId	

    A = complex(reA,imA)

    stop

end
