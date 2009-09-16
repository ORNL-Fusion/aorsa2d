pro over_sample_boundary,  r, z, newR, newz

	newR	= r[0]
	newZ	= z[0]

	for i=0,n_elements(r)-2 do begin

		m	= ( z[i+1]-z[i] ) $
				/ ( r[i+1] - r[i] )
		b	= z[i] - m * r[i]

		d	= sqrt ( ( r[i+1] - r[i] )^2 $
				+ ( z[i+1] - z[i] )^2 )

		nExtra	= 40
		dStep	= (r[i+1] - r[i]) / nExtra

		for j = 0, nExtra - 1 do begin

			if dStep ne 0 then begin
				newR	= [ newR, r[i] + dStep*j ]
				newZ	= [ newZ, m * (r[i] + dStep*j) + b ]
			endif 

		endfor

	endfor

	;newR	= [r, newR]
	;newz	= [z, newz]	
	
end 


pro create_antenna_current

	e	= 1.602e-19

	eqdsk_fileName	 = '~/data/eqdsk/g128797.00400_nstx_efit2'
	eqdsk	= readgeqdsk ( eqdsk_fileName )

	cdfId = ncdf_open ( '/home/dg6/scratch/aorsa2d/nstx/bench_solps_256x256_nPhi_-18/output/plotData.nc', /noWrite ) 
	ncdf_varGet, cdfId, 'rho', rho 
	ncdf_varget, cdfId, 'capR', capR 
	ncdf_varget, cdfId, 'zLoc', zLoc 
	nCdf_varGet, cdfId, 'mask', mask
	nCdf_varGet, cdfId, 'janty', janty_
	nCdf_varGet, cdfId, 'jantx', jantx_
	nCdf_varGet, cdfId, 'pscale', pscale
	ncdf_close, cdfId


	; test the div of jant from old aorsa file

	divJ_old	= fltArr ( size ( janty_, /dim ) )
	dR	= capR[1]-capR[0]
	dz	= zLoc[1]-zloc[0]

	for i=1,n_elements(janty_[*,0])-2 do begin
		for j=1,n_elements(janty_[0,*])-2 do begin
			
			divJ_old[i,j]	= ( janty_[i,j+1]-janty_[i,j-1] ) / dz $
					+ ( jantx_[i+1,j]-jantx_[i-1,j] ) / dR

		endfor
	endfor

	window, 1, xSize = 600, ySize = 600
	device, decomposed = 0
	loadct, 3
	!p.background = 255 
	levels	= fIndGen(21)*10
	colors	= 255 - ( bytScl ( levels, top = 253 ) + 1 )
	contour, janty_, capr, zloc, $
			yrange=[-0.5,0.6], $
			xrange=[1.4,1.58], $
			color=0, $
			/fill, $
			levels = levels, $
			c_colors = colors, /noData

	xStep	= capR[1] - capR[0]
	yStep 	= zLoc[1] - zLoc[0]

	for i=0,n_elements(janty_[*,0])-1 do begin
		for j=0,n_elements(janty_[0,*])-1 do begin

			if janty_[i,j] gt 0.1 then begin

			plots,[capR[i]-xStep/2.0, $
				capR[i]+xStep/2.0, $
				capR[i]+xStep/2.0, $
				capR[i]-xStep/2.0, $
				capR[i]-xStep/2.0], $
				[zLoc[j]-yStep/2.0, $
				zLoc[j]-yStep/2.0, $
				zLoc[j]+yStep/2.0, $
				zLoc[j]+yStep/2.0, $
				zLoc[j]-yStep/2.0], $
				color = 255 - bytScl ( janty_[i,j], max = 100, min = 0 ), $
				/data
	
			polyFill,[capR[i]-xStep/2.0, $
				capR[i]+xStep/2.0, $
				capR[i]+xStep/2.0, $
				capR[i]-xStep/2.0, $
				capR[i]-xStep/2.0], $
				[zLoc[j]-yStep/2.0, $
				zLoc[j]-yStep/2.0, $
				zLoc[j]+yStep/2.0, $
				zLoc[j]+yStep/2.0, $
				zLoc[j]-yStep/2.0], $
				color = 255-janty_[i,j], $
				/data

				endif
				
		endfor
	endfor
	loadct, 12, /sil
	oplot, eqdsk.rlim, eqdsk.zlim, $
			psym = -4, $
			color = 0
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
			color = 8*16-1



stop
	
	window, 0, xSize = 1600, ySize = 1200
	!p.multi = [0,3,2]
	!p.charSize = 2.0
	device, decomposed = 0
	loadct, 12, /sil
	!p.background = 255

	plot, eqdsk.rlim, eqdsk.zlim, $
			psym = -4, $
			color = 0
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
			color = 8*16-1

	oPlot, (eqdsk.rLim)[13:23], (eqdsk.zLim)[13:23], $
			thick = 3, $
			color = 4*16-1
	rShift	= 0.03
	oPlot, (eqdsk.rLim)[13:23]+rShift, (eqdsk.zLim)[13:23], $
			thick = 3, $
			color = 12*16-1

	topDistance	= $
			sqrt ( ((eqdsk.rLim)[12] - ((eqdsk.rLim)[13]+rShift))^2 $
				+ ((eqdsk.zLim)[12] - (eqdsk.zLim)[13])^2 )
	topDirx	=  ((eqdsk.rLim)[12] - ((eqdsk.rLim)[13]+rShift)) / topDistance
	topDiry	=  ((eqdsk.zLim)[12] - (eqdsk.zLim)[13]) / topDistance

	print, 'top:', topDirX, topDirY, sqrt ( topDirX^2 + topDirY^2 )

	topPtX	= topDistance / 2.0  * topDirX + ((eqdsk.rLim)[13] + rShift)
	topPtY	= topDistance / 2.0  * topDirY + (eqdsk.zLim)[13]
	
	plots, [(eqdsk.rLim)[13]+rShift,topPtX], [(eqdsk.zLim)[13],topPtY], $
			color = 5*16-1, $
			thick = 4

	botDistance	= $
			sqrt ( ((eqdsk.rLim)[24] - ((eqdsk.rLim)[23]+rShift))^2 $
				+ ((eqdsk.zLim)[24] - (eqdsk.zLim)[23])^2 )
	botDirx	=  ((eqdsk.rLim)[24] - ((eqdsk.rLim)[23]+rShift)) / botDistance
	botDiry	=  ((eqdsk.zLim)[24] - (eqdsk.zLim)[23]) / botDistance

	print, 'bot:', botDirX, botDirY, sqrt ( botDirX^2 + botDirY^2 )

	botPtX	= botDistance / 2.0  * botDirX + ((eqdsk.rLim)[23] + rShift)
	botPtY	= botDistance / 2.0  * botDirY + (eqdsk.zLim)[23]
	
	plots, [(eqdsk.rLim)[23]+rShift,botPtX], [(eqdsk.zLim)[23],botPtY], $
			color = 5*16-1, $
			thick = 4

	antR	= [ topPtX, (eqdsk.rLim)[13:23], botPtX ]
	antz	= [ topPtY, (eqdsk.zLim)[13:23], botPtY ]

	eqdsk.rLim[13:23]	= eqdsk.rLim[13:23] + rShift

	plot, eqdsk.rlim, eqdsk.zlim, $
			psym = -4, $
			color = 0
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
			color = 8*16-1

	oPlot, antR, antZ, $
			color = 8*16-1, $
			thick = 4

	plot, (eqdsk.rlim)[12:24], (eqdsk.zlim)[12:24], $
			psym = -4, $
			color = 0, $
			xRange = [1.2, 1.9]
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
			color = 8*16-1

;	oversample ant coords

	over_sample_boundary, antR, antz, newAntR, newAntz

	iiSort 	= sort ( newAntz )

	newAntR	= newAntR[iiSort]
	newAntZ	= newAntZ[iiSort]

	oPlot, newantR, newantZ, $
			color = 8*16-1, $
			thick = 4

;	create jX and jY along the line

	nAnt	= n_elements ( newantR )

	antJX_line	= fltArr ( nAnt )
	antJY_line	= fltArr ( nAnt )

	for i=0,nAnt-2 do begin

		distance = sqrt ( (newAntR[i] - newAntR[i+1])^2 + (newAntZ[i] - newAntZ[i+1])^2 )
		jXDir	= (newAntR[i] - newAntR[i+1]) / distance
		jYDir	= -(newAntz[i] - newAntz[i+1]) / distance

		antJX_line[i]	= jXDir
		antJY_line[i]	= jYDir

		if distance eq 0 then stop

	endfor

	antJX_line[nAnt-1]	= antJX_line[nAnt-2]
	antJY_line[nAnt-1]	= antJY_line[nAnt-2]

;	normalise to Am^-1 for 1 Amp

	antThick = 0.03
	antJX_line	= antJX_line / antThick ;/ (2.0*!pi*newAntR)
	antJY_line	= antJY_line / antThick ;/ (2.0*!pi*newAntR)

	; test the div of jant new 

	divJ_new	= fltArr ( size ( jAnty_, /dim ) )
	dR	= xGrid[1]-xGrid[0]
	dz	= yGrid[1]-yGrid[0]

	for i=1,n_elements(jAnty[*,0])-2 do begin
		for j=1,n_elements(jAnty[0,*])-2 do begin
			
			divJ_new[i,j]	= ( janty[i,j+1]-janty[i,j-1] ) / dz $
					+ ( jantx[i+1,j]-jantx[i-1,j] ) / dR

		endfor
	endfor


	for i=1,nX-2 do begin
		for j=1,nY-2 do begin

			if jAntCnt[i,j] gt 0 then begin

				plots, 	[xGrid[i]-xStep/2.0, $
						xGrid[i]+xStep/2.0, $
						xGrid[i]+xStep/2.0, $
						xGrid[i]-xStep/2.0, $
						xGrid[i]-xStep/2.0], $
						[yGrid[j]-yStep/2.0, $
						yGrid[j]-yStep/2.0, $
						yGrid[j]+yStep/2.0, $
						yGrid[j]+yStep/2.0, $
						yGrid[j]-yStep/2.0], $
						color = 0, /data

			endif

		endfor
	endfor	

	loadct, 3, /sil
	for i=1,nXDiv-2 do begin
		for j=1,nYDiv-2 do begin

			if jAntDiv[i,j] gt 0 then begin

				plots, 	[xGrid[i*xSam-(xSam/2-1)], $
						xGrid[i*xSam+(xSam/2-1)], $
						xGrid[i*xSam+(xSam/2-1)], $
						xGrid[i*xSam-(xSam/2-1)], $
						xGrid[i*xSam-(xSam/2-1)]], $
						[yGrid[j*ySam-(ySam/2-1)], $
						yGrid[j*ySam-(ySam/2-1)], $
						yGrid[j*ySam+(ySam/2-1)], $
						yGrid[j*ySam+(ySam/2-1)], $
						yGrid[j*ySam-(ySam/2-1)]], $
						color = 231-(bytScl ( abs ( jAntDiv[i,j] ), max=10, min = 0, top=230 ) + 1), $
						/data, $
						thick = 3.0
				print, ''
				print, xGrid[i*xSam], yGrid[j*ySam]
				print, jAntDiv[i,j]

				jyTop	= total ( jAntY[i*xSam-(xSam/2-1):i*xSam+(xSam/2-1),j*ySam-(ySam/2-1)] )
				jyBot	= total ( jAntY[i*xSam-(xSam/2-1):i*xSam+(xSam/2-1),j*ySam+(ySam/2-1)] )
				jyL	= total ( jAntY[i*xSam-(xSam/2-1),j*ySam-(ySam/2-1):j*ySam+(ySam/2-1)] )
				jyR	= total ( jAntY[i*xSam+(xSam/2-1),j*ySam-(ySam/2-1):j*ySam+(ySam/2-1)] )

				jxTop	= total ( jAntX[i*xSam-(xSam/2-1):i*xSam+(xSam/2-1),j*ySam-(ySam/2-1)] )
				jxBot	= total ( jAntX[i*xSam-(xSam/2-1):i*xSam+(xSam/2-1),j*ySam+(ySam/2-1)] )
				jxL	= total ( jAntX[i*xSam-(xSam/2-1),j*ySam-(ySam/2-1):j*ySam+(ySam/2-1)] )
				jxR	= total ( jAntX[i*xSam+(xSam/2-1),j*ySam-(ySam/2-1):j*ySam+(ySam/2-1)] )


				print, jyTop, jyBot, jyL, jyR
				print, jxTop, jxBot, jxL, jxR

			endif

		endfor
	endfor	

	loadct, 12, /sil
	contour, janty_, capr, zloc, $
			color = 0, $
			xRange = [1.2, 1.9], $
			yRange = [-0.6, 0.6] 
	oPlot, eqdsk.rlim, eqdsk.zlim, $
			psym = -4, $
			color = 0
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
			color = 8*16-1

	plot, xGrid, jAnty[*,n_elements(janty[0,*])/2-3], $
			color = 0
	plot, capR, jAnty_[*,n_elements(janty_[0,*])/2-3], $
			color = 0


;	save modified rLim/zLim boundary in eqdsk file

	openw, lun, eqdsk_fileName+'.dlgMod', /get_lun

	f1  = '(6a8,3i4)'
	f2  = '(5e16.9)'
	f3  = '(2i5)'

	printf, lun, format = f1, eqdsk.case_, eqdsk.idum, eqdsk.nW, eqdsk.nH
	printf, lun, format = f2, eqdsk.rdim, eqdsk.zdim, eqdsk.rcentr, eqdsk.rleft, eqdsk.zmid
	printf, lun, format = f2, eqdsk.rmaxis, eqdsk.zmaxis, eqdsk.simag, eqdsk.sibry, eqdsk.bcentr
	printf, lun, format = f2, eqdsk.current, eqdsk.simag, eqdsk.xdum, eqdsk.rmaxis, eqdsk.xdum
	printf, lun, format = f2, eqdsk.zmaxis, eqdsk.xdum, eqdsk.sibry, eqdsk.xdum, eqdsk.xdum
	printf, lun, format = f2, eqdsk.fpol
	printf, lun, format = f2, eqdsk.pres
	printf, lun, format = f2, eqdsk.ffprim 
	printf, lun, format = f2, eqdsk.pprime 
	printf, lun, format = f2, eqdsk.psizr 
	printf, lun, format = f2, eqdsk.qpsi
	printf, lun, format = f3, eqdsk.nbbbs, eqdsk.limitr
	printf, lun, format = f2, [eqdsk.rbbbs,eqdsk.zbbbs]
	printf, lun, format = f2, [eqdsk.rlim,eqdsk.zlim]

	close, lun


;	write netCDF data file containing antenna current to be 
;	read by aorsa

	outFileName	= '/home/dg6/data/aorsa_ant/dlg_nstx_ant.nc'
	nc_id	= nCdf_create ( outFileName, /clobber )
	nCdf_control, nc_id, /fill
	
	nR_id	= nCdf_dimDef ( nc_id, 'nR', nX )
	nz_id	= nCdf_dimDef ( nc_id, 'nz', nY )
	scalar_id	= nCdf_dimDef ( nc_id, 'scalar', 1 )
	
	R_id = nCdf_varDef ( nc_id, 'R_binCenters', [ nR_id ], /float )
	z_id = nCdf_varDef ( nc_id, 'z_binCenters', [ nz_id ], /float )
	jantx_id = nCdf_varDef ( nc_id, 'jantx', [nR_id, nz_id], /float )
	janty_id = nCdf_varDef ( nc_id, 'janty', [nR_id, nz_id], /float )

	nCdf_control, nc_id, /enDef
	
	nCdf_varPut, nc_id, R_id, xGrid
	nCdf_varPut, nc_id, z_id, yGrid 
	nCdf_varPut, nc_id, jantx_id, jAntX 
	nCdf_varPut, nc_id, janty_id, jAntY 

	nCdf_close, nc_id



	stop

end
