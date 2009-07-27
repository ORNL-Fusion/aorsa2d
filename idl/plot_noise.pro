pro plot_noise

	cdfId	= ncdf_open ( 'output/noise.nc', /noWrite )
	glob	= ncdf_inquire ( cdfId )
	
	ncdf_varGet, cdfId, 'upar', upar 
	ncdf_varGet, cdfId, 'upar0', upar0
	ncdf_varGet, cdfId, 'dfdupar', dfdupar 
	ncdf_varGet, cdfId, 'dfduper', dfduper 

	!p.multi = [0,1,2]

	plot, upar, dfdupar
	oplot, [upar0,upar0], [min(dfdupar),max(dfdupar)]
	plot, upar, dfduper
	oplot, [upar0,upar0], [min(dfduper),max(dfduper)]
	stop

end
