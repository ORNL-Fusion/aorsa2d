pro plot_f

	fileName	= 'output/p_f.nc'
	cdfId	= ncdf_open ( fileName, /noWrite )
	glob	= ncdf_inquire ( cdfId )
	ncdf_varGet, cdfId, 'dfduper', dfduper
	ncdf_varGet, cdfId, 'dfdupar', dfdupar
	ncdf_varGet, cdfId, 'vPer', vper
	ncdf_varGet, cdfId, 'vPar', vpar

	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'nR' ), name, nR 
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'nz' ), name, nz 
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'nuper' ), name, nuper 
	ncdf_dimInq, cdfId, ncdf_dimId ( cdfId, 'nupar' ), name, nupar 

	ncdf_close, cdfId

	device, decomposed = 0
	loadct, 13, file='davect.tbl'
	levels	= (fIndGen(21)-10)/10 * 1e1
	colors	= bytScl ( levels, top = 253 ) + 1

	!p.multi = [0,8,8]
	iStride	= nR/8
	jStride	= nz/8

	window, 0, ySize = 800
	for j=0,nz-1,jStride do begin
	for i=0,nR-1,iStride do begin

	contour, transpose ( reform ( dfduper[i,j,*,*] ) ), $
			vPar, vPer, $
			levels = levels, $
			c_colors = colors, $
			/fill, charSize = 0.01
	endfor
	endfor
	window, 1, ySize = 800
	for j=0,nz-1,jStride do begin
	for i=0,nR-1,iStride do begin

	contour, transpose ( reform ( dfdupar[i,j,*,*] ) ), $
			vPar, vPer, $
			levels = levels, $
			c_colors = colors, $
			/fill, charSize = 0.01
	endfor
	endfor
stop

end
