pro plot_f_rho

	if strCmp ( getEnv ( 'MACHINE' ), 'jaguar' ) then scratchDir = getEnv ( 'SYSTEM_USERDIR' )
	if strCmp ( getEnv ( 'MACHINE' ), 'franklin' ) then scratchDir = getEnv ( 'SCRATCH' )
	if strCmp ( getEnv ( 'MACHINE' ), 'dlghp' ) then scratchDir = '/home/dg6/scratch' 

	filename	= 'output/p_f_rho.nc'

	cdfId	= ncdf_open ( fileName, /noWrite )
	glob	= ncdf_inquire ( cdfId )
	ncdf_varGet, cdfId, 'f_rho_vv', f_rho_vv
	ncdf_varGet, cdfId, 'dfduPer', dfduPer
	ncdf_varGet, cdfId, 'dfduPar', dfduPar
	ncdf_varGet, cdfId, 'f_all', f_all

	ncdf_close, cdfId



stop

end

