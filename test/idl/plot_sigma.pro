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

	sig11	= complex ( sig11_re, sig11_im )
	sig12	= complex ( sig12_re, sig12_im )
	sig13	= complex ( sig13_re, sig13_im )

	sig21	= complex ( sig21_re, sig21_im )
	sig22	= complex ( sig22_re, sig22_im )
	sig23	= complex ( sig23_re, sig23_im )

	sig31	= complex ( sig31_re, sig31_im )
	sig32	= complex ( sig32_re, sig32_im )
	sig33	= complex ( sig33_re, sig33_im )

	isurface, imaginary(sig11), r, k, view_grid=[3,3], /zoom_on_resize
	isurface, imaginary(sig12), r, k, /view_next
	isurface, imaginary(sig13), r, k, /view_next

	isurface, imaginary(sig21), r, k, /view_next
	isurface, imaginary(sig22), r, k, /view_next
	isurface, imaginary(sig23), r, k, /view_next

	isurface, imaginary(sig31), r, k, /view_next
	isurface, imaginary(sig32), r, k, /view_next
	isurface, imaginary(sig33), r, k, /view_next

	isurface, sig11, r, k, view_grid=[3,3], /zoom_on_resize
	isurface, sig12, r, k, /view_next
	isurface, sig13, r, k, /view_next

	isurface, sig21, r, k, /view_next
	isurface, sig22, r, k, /view_next
	isurface, sig23, r, k, /view_next

	isurface, sig31, r, k, /view_next
	isurface, sig32, r, k, /view_next
	isurface, sig33, r, k, /view_next

stop
end
