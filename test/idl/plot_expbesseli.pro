pro plot_expBesselI

	cdfId = ncdf_open ( 'test_expBesselI.nc', /noWrite ) 

		ncdf_varget, cdfid, 'res_re', res_re 
		ncdf_varget, cdfid, 'res_im', res_im
		nCdf_varGet, cdfId, 'z_re', z_re 
		nCdf_varGet, cdfId, 'z_im', z_im
		ncdf_varget, cdfid, 'expRes_re', expres_re 
		ncdf_varget, cdfid, 'expRes_im', expres_im
	
	ncdf_close, cdfId

	res = complex ( res_re, res_im )
	expres = complex ( expres_re, expres_im )
	z	= complex ( z_re, z_im )

	iSurface, res, z[*,0], imaginary(z[0,*])
	iSurface, expres, z[*,0], imaginary(z[0,*])

end
