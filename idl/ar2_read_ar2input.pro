pro ar2_read_ar2input, fileName, rLim=rLim, zLim=zLim, LimMask=LimMask

	cdfId = ncdf_open ( fileName, /noWrite ) 
	nCdf_varGet, cdfid, 'rMin', rMin
	nCdf_varGet, cdfid, 'rMax', rMax
	nCdf_varGet, cdfid, 'zMin', zMin
	nCdf_varGet, cdfid, 'zMax', zMax

	nCdf_varGet, cdfid, 'r', r 
	nCdf_varGet, cdfid, 'z', z
	nCdf_varGet, cdfid, 'br', br
	nCdf_varGet, cdfid, 'bz', bt 
	nCdf_varGet, cdfid, 'bz', bz

	nCdf_varGet, cdfid, 'AtomicZ', AtomicZ 
	nCdf_varGet, cdfid, 'amu', amu 

	nCdf_varGet, cdfid, 'Density_m3', Density_m3
	nCdf_varGet, cdfid, 'Temp_eV', Temp_eV

	nCdf_varGet, cdfid, 'LimMask', mask_lim 
	nCdf_varGet, cdfid, 'Lim_r', rlim 
	nCdf_varGet, cdfid, 'Lim_z', zlim 
	ncdf_close, cdfId

	cdfId = ncdf_open ( 'runData001.nc', /noWrite ) 
		nCdf_varGet, cdfid, 'LimMask', LimMask
	ncdf_close, cdfId


end
