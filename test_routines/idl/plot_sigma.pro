pro plot_sigma

	cdfId = ncdf_open ( 'test_sigma.nc', /noWrite ) 

		ncdf_varget, cdfid, 'sig11_re', sig11_re 
		ncdf_varget, cdfid, 'sig11_im', sig11_im
		ncdf_varget, cdfid, 'sig12_re', sig12_re 
		ncdf_varget, cdfid, 'sig12_im', sig12_im
		ncdf_varget, cdfid, 'sig13_re', sig13_re 
		ncdf_varget, cdfid, 'sig13_im', sig13_im

		ncdf_varget, cdfid, 'sig21_re', sig21_re 
		ncdf_varget, cdfid, 'sig21_im', sig21_im
		ncdf_varget, cdfid, 'sig22_re', sig22_re 
		ncdf_varget, cdfid, 'sig22_im', sig22_im
		ncdf_varget, cdfid, 'sig23_re', sig23_re 
		ncdf_varget, cdfid, 'sig23_im', sig23_im

		ncdf_varget, cdfid, 'sig31_re', sig31_re 
		ncdf_varget, cdfid, 'sig31_im', sig31_im
		ncdf_varget, cdfid, 'sig32_re', sig32_re 
		ncdf_varget, cdfid, 'sig32_im', sig32_im
		ncdf_varget, cdfid, 'sig33_re', sig33_re 
		ncdf_varget, cdfid, 'sig33_im', sig33_im

		ncdf_varget, cdfid, 'r', r 
		ncdf_varget, cdfid, 'k', k 

	ncdf_close, cdfId

	sig11	= dcomplex ( sig11_re, sig11_im )
	sig12	= dcomplex ( sig12_re, sig12_im )
	sig13	= dcomplex ( sig13_re, sig13_im )

	sig21	= dcomplex ( sig21_re, sig21_im )
	sig22	= dcomplex ( sig22_re, sig22_im )
	sig23	= dcomplex ( sig23_re, sig23_im )

	sig31	= dcomplex ( sig31_re, sig31_im )
	sig32	= dcomplex ( sig32_re, sig32_im )
	sig33	= dcomplex ( sig33_re, sig33_im )

	isurface, imaginary(sig11), r, k, view_grid=[3,3], /zoom_on_resize
	isurface, imaginary(sig12), r, k, /view_next
	isurface, imaginary(sig13), r, k, /view_next

	isurface, imaginary(sig21), r, k, /view_next
	isurface, imaginary(sig22), r, k, /view_next
	isurface, imaginary(sig23), r, k, /view_next

	isurface, imaginary(sig31), r, k, /view_next
	isurface, imaginary(sig32), r, k, /view_next
	isurface, imaginary(sig33), r, k, /view_next

	isurface, double(sig11), r, k, view_grid=[3,3], /zoom_on_resize
	isurface, double(sig12), r, k, /view_next
	isurface, double(sig13), r, k, /view_next

	isurface, double(sig21), r, k, /view_next
	isurface, double(sig22), r, k, /view_next
	isurface, double(sig23), r, k, /view_next

	isurface, double(sig31), r, k, /view_next
	isurface, double(sig32), r, k, /view_next
	isurface, double(sig33), r, k, /view_next

stop
end
