pro plot_aorsa_3d, $
		range = range

	eqDskFileName = 'eqdsk'

	eqdsk	= readGEQDSK ( eqDskFileName )
	cdfId = ncdf_open ( 'output/E3D.nc', /noWrite ) 
	
	ncdf_varGet, cdfId, 'eMod',  eMod
	ncdf_varGet, cdfId, 'R', R 
	ncdf_varGet, cdfId, 'z', z_
	ncdf_varGet, cdfId, 'phi', phi 

	ncdf_close, cdfId


	!p.multi = [0,2,2]
	window, 0, xSize = 600, ySize = 800
	device, decomposed = 0
	!p.backGround = 255

	eRange = 4e4
	eLevs	= fIndGen ( 21 ) / 20 * eRange
	eCols	= 255 - ( bytScl ( eLevs, top = 254 ) + 1 )

	for i = 0, 3 do begin

		loadct, 3, /sil
		contour, eMod[*,*,i*n_elements(eMod[0,0,*])/4], R, z_, $
			   /iso, lev = eLevs, $
			   c_col = eCols, /fill, $
			   color = 0

		loadct, 12
		oPlot, eqdsk.rbbbs, eqdsk.zbbbs, color = 1*16-1
		oPlot, eqdsk.rLim, eqdsk.zLim, color = 8*16-1

	endfor

	stop
end
