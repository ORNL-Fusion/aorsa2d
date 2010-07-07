pro plot_spline

	cdfId = ncdf_open ( 'test_spline.nc', /noWrite ) 

		ncdf_varget, cdfid, 'x', x 
		ncdf_varget, cdfid, 'y', y 
		ncdf_varget, cdfid, 'x2', x2
		ncdf_varget, cdfid, 's', s 

	ncdf_close, cdfId

	iPlot, x, y, sym_index = 4
	iPlot, x2, s, sym_index = 6, thick = 2, /over

end
