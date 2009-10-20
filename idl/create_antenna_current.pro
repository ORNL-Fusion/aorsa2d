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

	newR	= [r, newR]
	newz	= [z, newz]	
	
end 

pro is_inside, r, z, rBry, zBry, $
		mask = mask

	mask   = intArr ( size ( r, /dim ) )
	nR	= n_elements ( r[*,0] )
	nz	= n_elements ( r[0,*] )
	
	for i = 0, nR - 1 do begin
	    for j = 0, nZ - 1 do begin 
	
	        q1_  = n_elements ( where ( ( r[i,j] - rBry gt 0 ) and ( z[i,j] - zBry gt 0 ), q1 ) )
	        q2_  = n_elements ( where ( ( r[i,j] - rBry gt 0 ) and ( z[i,j] - zBry le 0 ), q2 ) )
	        q3_  = n_elements ( where ( ( r[i,j] - rBry le 0 ) and ( z[i,j] - zBry gt 0 ), q3 ) )
	        q4_  = n_elements ( where ( ( r[i,j] - rBry le 0 ) and ( z[i,j] - zBry le 0 ), q4 ) )
	
	        if ( q1 gt 0 ) and ( q2 gt 0 ) and ( q3 gt 0 ) and ( q4 gt 0 ) then $
	            mask[i,j]  = 1
	               
	    endfor
	endfor

end 

pro create_antenna_current, $
		gridMatch = gridMatch

	e	= 1.602e-19

	;eqdsk_fileName	 = '~/data/nstx/eqdsk/g120742.00275.EFIT02.mds.uncorrected.qscale_1.00000'
	;eqdsk_fileName	 = '~/data/nstx/eqdsk/g120745.00277.EFIT02.mds.uncorrected.qscale_1.00000'
	eqdsk_fileName	 = '~/data/nstx/eqdsk/g120740.00275.EFIT02.mds.uncorrected.qscale_1.00000'


	eqdsk	= readgeqdsk ( eqdsk_fileName )

	if keyword_set ( gridMatch ) then begin

		cdfId = ncdf_open ( '/home/dg6/scratch/aorsa2d/nstx/bench/output/plotData.nc', /noWrite ) 
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

	endif

	window, 0, xSize = 800, ySize = 600
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
	
	;	first step is to just use a straight vertical antenna connecting
	;	the top and bottom points, do overwrite this for the time being

	;antR	= antR * 0 + topPtX*1.04
	;antz	= (antz<0.15)>(-0.15)

	eqdsk.rLim[13:23]	= eqdsk.rLim[13:23] + rShift
	;iiShift	= where ( eqdsk.rLim gt 0.6, iiShiftCnt )
	;eqdsk.rLim[iiShift]	= 1.8

	;	custom limiter

	eqdsk.rLim[*]	= 0.18
   	eqdsk.zLim[*]	= -1.6
	eqdsk.rlim	= [ 0.18, 1.4, 1.74, 1.74, 1.4, 0.18, 0.18 ]
	eqdsk.zLim	= [ -1.6, -1.6, -0.3, 0.3, 1.6, 1.6, -1.6 ]
	eqdsk.limitr	= n_elements ( eqdsk.rlim )

	;	grow custom limiter

	grow_rlim	= [ 0.0, 1.56, 1.92, 1.92, 1.56, 0.0, 0.0 ]
	grow_zLim	= [ -1.75, -1.75, -0.32, 0.32, 1.75, 1.75, -1.75 ]

	plot, eqdsk.rlim, eqdsk.zlim, $
			psym = -4, $
			color = 0
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
			color = 8*16-1

	oplot, grow_rlim, grow_zlim, $
			psym = -4, $
			color = 12*16-1 
	
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

	antCnt = 0
	for i=0,nAnt-2 do begin

		distance = sqrt ( (newAntR[i] - newAntR[i+1])^2 + (newAntZ[i] - newAntZ[i+1])^2 )
		jXDir	= (newAntR[i] - newAntR[i+1]) / distance
		jYDir	= -(newAntz[i] - newAntz[i+1]) / distance

		if distance gt 0 then begin

			if size(antJX_line,/ty) eq 0 then antJX_line = jXDir $
					else antJX_line = [antJX_line,jXDir]
			if size(antJY_line,/ty) eq 0 then antJY_line = jYDir $
					else antJY_line = [antJY_line,jYDir]

			antCnt++

		endif

	endfor

	antJX_line	= [antJX_line, antJX_line[antCnt-1]]
	antJY_line	= [antJY_line, antJY_line[antCnt-1]]
	nAnt	= antCnt+1

;	put the antenna line on a grid

	nX	= 120 
	nY	= 120 

	rMax	= max ( eqdsk.r )
	rMin	= -0.2
	xRange	= rMax - rMin
	dX	= xRange / nX
	antGrid_x	= fIndGen ( nX ) * dX + rMin + dX / 2.0
	yRange	= max ( eqdsk.z ) - min ( eqdsk.z )
	antGrid_y	= fIndGen ( nY )/(nY-1) * yRange + min ( eqdsk.z )

	antGrid_x2D	= rebin ( antGrid_x, nX, nY )
	antGrid_y2D	= transpose ( rebin ( antGrid_y, nY, nX ) )

	if keyword_set ( gridMatch ) then begin

		nX	= n_elements ( capR )
		nY	= n_elements ( zLoc )
		xRange = max ( capR ) - min ( capR )
		yRange = max ( zLoc ) - min ( zLoc )
		antGrid_x	= capR
		antGrid_y	= zLoc

	endif

	antJX_grid	= fltArr ( nX, nY )
	antJY_grid	= fltArr ( nX, nY )

	;; this only works for the vertical line current 

	;antXii	=  ( newAntR[0] - min(antGrid_x) ) / xRange * nX
	;antYii_top	=  ( max(newAntZ) - min(antGrid_y) ) / yRange * nY
	;antYii_bot	=  ( min(newAntZ) - min(antGrid_y) ) / yRange * nY

	;	define connecting path for the antenna, 
	;	at the moment it requires right angle connectors,
	;	but this will be changed later

	xAnt	= 1.7
	yAnt1	= -1.4
	yAnt2	=  1.4
	iiFeedLength	= where ( antGrid_x ge xAnt, iiFeedCnt )
	iiAntLength		= where ( antGrid_y ge yAnt1 and antGrid_y le yAnt2, iiAntCnt )
	bottomFeederX	= nX - indGen ( iiFeedCnt ) -1
	bottomFeederY	= bottomFeederX*0+min(iiAntLength)
	antLengthX		= intArr ( iiAntCnt ) + bottomFeederX[iiFeedCnt-1]
	antLengthY		= indGen ( iiAntCnt ) + min ( iiAntLength )
	topFeederX		= reverse ( bottomFeederX ) 
	topFeederY		= bottomFeederX*0+max(iiAntLength)
	
	;ant_grid_path_X	= [ bottomFeederX, antLengthX[1:*], topFeederX[1:*] ]
	;ant_grid_path_Y	= [ bottomFeederY, antLengthY[1:*], topFeederY[1:*] ]

	;trace_ant_path, ant_grid_path_X, ant_grid_path_Y, nX, nY, $
	;		JxInterp = antJX_grid, JyInterp = antJY_grid, $
	;		dX = antGrid_x[1]-antGrid_x[0], $
	;		dY = antGrid_y[1]-antGrid_y[0]

;	overwrite with simple form to start

	dz	= antGrid_x[1] - antGrid_x[0]
	dR	= antGrid_y[1] - antGrid_y[0]

	antJX_grid[*]	= 0
	antJY_grid[*]	= 0
	;antJX_grid[bottomFeederX[0:iiFeedCnt-2],bottomFeederY[0:iiFeedCnt-2]]	= -1.0*dR
	;antJY_grid[antLengthX[1:iiAntCnt-2],antLengthY[1:iiAntCnt-2]]	= 1.0*dz
	;antJY_grid[antLengthX,antLengthY]	= 1.0*dz


	iiAntYRange	= where ( antGrid_y gt yAnt1 and antGrid_y lt yAnt2 )
	iiTop	= max ( iiAntYRange, min = iiBot )
	sigX	= 0.01
	sigY	= 0.2
	antJY_grid	= exp ( -( $
			( antGrid_x2D-xAnt )^2 / sigX + ( antGrid_y2D-antGrid_y[iiTop] )^2 / sigY ) )
	antJY_grid	+= exp ( -( $
			( antGrid_x2D-xAnt )^2 / sigX + ( antGrid_y2D-antGrid_y[iiBot] )^2 / sigY ) )

	antJY_grid[*,iiAntYRange]	= exp ( -( antGrid_x[*,iiAntYRange]-xAnt )^2 / sigX )

	antJY_grid	= antJY_grid / max ( antJY_grid )

	;	align direction of J to be perp to b

	ant_byTop = interpolate ( eqdsk.bz, ( xAnt - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nw-1), $
    				( 0.4 - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nh-1) )
	ant_byBot = interpolate ( eqdsk.bz, ( xAnt - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nw-1), $
    				( -0.4 - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nh-1) )
	ant_bPhiTop = interpolate ( eqdsk.bPhi, ( xAnt - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nw-1), $
    				( 0.4 - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nh-1) )
	ant_bPhiBot = interpolate ( eqdsk.bPhi, ( xAnt - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nw-1), $
    				( -0.4 - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nh-1) )
	ant_by = interpolate ( eqdsk.bz, ( xAnt - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nw-1), $
    				( 0.0 - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nh-1) )
	ant_bPhi = interpolate ( eqdsk.bPhi, ( xAnt - eqdsk.rleft ) / eqdsk.rdim * (eqdsk.nw-1), $
    				( 0.0 - min ( eqdsk.z ) ) / eqdsk.zdim * (eqdsk.nh-1) )
			
	bAntTop_y	= ant_byTop / sqrt ( ant_byTop^2 + ant_bPhiTop^2 )
	bAntTop_phi	= ant_bPhiTop / sqrt ( ant_byTop^2 + ant_bPhiTop^2 )
	bAntBot_y	= ant_byBot / sqrt ( ant_byBot^2 + ant_bPhiBot^2 )
	bAntBot_phi	= ant_bPhiBot / sqrt ( ant_byBot^2 + ant_bPhiBot^2 )
	bAnt_y	= ant_by / sqrt ( ant_by^2 + ant_bPhi^2 )
	bAnt_phi	= ant_bPhi / sqrt ( ant_by^2 + ant_bPhi^2 )

	mean_bY	= mean ( [ bAntTop_y, bAntBot_y, bAnt_y ] )
	mean_bPhi	= mean ( [ bAntTop_phi, bAntBot_phi, bAnt_phi ] )

	antJY_grid	= antJY_grid * mean_bPhi
	antJZ_grid	= antJY_grid * mean_by 

stop
	;iiZeroOut	= where ( antJY_grid lt max(antJY_grid)/100.0*5 )

	;antJY_grid[iiZeroOut]	= 0.0

	;antJX_grid[topFeederX[1:*],topFeederY[1:*]]	= 1.0*dR

	;iiNewAntR	= (newAntR-min(capR))/(max(capR)-min(capR))*nX
	;iiNewAntz	= (newAntz-min(zLoc))/(max(zLoc)-min(zLoc))*nY

	;for i=0,n_elements(newAntR)-1 do begin

	;		antJY_grid[iiNewAntR[i],iiNewAntz[i]]	= 1.0

	;endfor


;	antJX_grid	= smooth ( antJX_grid, 4 )
;	antJY_grid	= smooth ( antJY_grid, 4 )
;
;	antJX_bak	= antJX_grid
;	antJY_bak	= antJY_grid
;
	;for i=0,nX-1 do antJX_grid[i,*]	= ts_smooth ( (antJX_grid[i,*])[*], 4 )
	;for j=0,nY-1 do antJY_grid[*,j]	= ts_smooth ( antJY_grid[*,j], 4 )

	;iiBadX	= where ( antJX_grid ne antJX_grid )
	;antJX_grid[iiBadX]	= antJX_bak[iiBadX]
	;iiBadY	= where ( antJY_grid ne antJY_grid )
	;antJY_grid[iiBadY]	= antJY_bak[iiBadY]


;	normalise to Am^-1 for 1 Amp

	antThick = 0.03
	antJX_line	= antJX_line / antThick ;/ (2.0*!pi*newAntR)
	antJY_line	= antJY_line / antThick ;/ (2.0*!pi*newAntR)

	; test the div of jant new 

	divJ_new	= fltArr ( nX, nY )
	for i=1,nX-2 do begin
		for j=1,nY-2 do begin
			
			divJ_new[i,j]	= ( antJY_grid[i,j+1]-antJY_grid[i,j-1] ) / dz $
					+ ( antJX_grid[i+1,j]-antJX_grid[i-1,j] ) / dR

		endfor
	endfor


	;	create the smooth wall sigma

	over_sample_boundary, (eqdsk.rlim)[*], (eqdsk.zlim)[*], rlim_os, zlim_os
	is_inside, antGrid_x2D, antGrid_y2D, rlim_os, zlim_os, $
			mask = mask_innerWall

	over_sample_boundary, grow_rlim, grow_zlim, grow_rlim_os, grow_zlim_os
	is_inside, antGrid_x2D, antGrid_y2D, grow_rlim_os, grow_zlim_os, $
			mask = mask_outerWall

	smoothSigma	= fltArr ( nX, nY )
	wallSigX	= 0.04
	wallSigY	= 0.02
	
	for i=0,nX-1 do begin
		for j=0,nY-1 do begin

			if mask_outerWall[i,j] ne 1 then begin

				smoothSigma	+= exp ( -( $
					( antGrid_x2D-antGrid_x2D[i,j] )^2 / wallSigX $
					+ ( antGrid_y2D-antGrid_y2D[i,j] )^2 / wallSigY ) )
	
			endif

		endfor
	endfor

	smoothSigma	= smoothSigma / max ( smoothSigma )

	loadct, 12, /sil
	window, 3, ySize = 900
	!p.multi = [0,2,3]
	contour, sqrt(antJX_grid^2+antJY_grid^2), antGrid_x, antGrid_y, $
			color = 0, $
			xRange = [1.2, 1.9], $
			yRange = [-0.6, 0.6] 
	oPlot, eqdsk.rlim, eqdsk.zlim, $
			psym = -4, $
			color = 0
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
			color = 8*16-1

	contour, divJ_new, antGrid_x, antGrid_y, $
			color = 0, $
			xRange = [1.3, 1.7], $
			yRange = [-0.8, 0.8] 
	oPlot, eqdsk.rlim, eqdsk.zlim, $
			psym = -4, $
			color = 0
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
			color = 8*16-1
	veloVect, antJX_grid, antJY_grid, antGrid_x, antGrid_y,$
		   	/over, color = 12*16-1, length= 1.5
	
	veloVect, antJX_grid, antJY_grid, antGrid_x, antGrid_y,$
		   	color = 12*16-1, length= 2.5

	plot, antjy_grid[nX/4*3,*],color=0, psym = -4
	plot, antjy_grid[*,nY/2],color=0, psym = -4
	contour, antjy_grid, antGrid_x, antGrid_y, $
			color = 0, nlev=30
	oPlot, eqdsk.rlim, eqdsk.zlim, $
			psym = -4, $
			color = 0
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
			color = 8*16-1


	window, 9, ySize = 1200
	!p.multi = 0
	loadct, 3
	contour, smoothsigma, antgrid_x, antgrid_y, $
		color = 0, levels = (findgen(10)+1)/10, $
		c_labels = fIndGen(10)+1, /iso, $
		/fill

	loadct, 12
	oplot, eqdsk.rlim, eqdsk.zlim, $
			psym = -4, $
			color = 8*16-1, $
			thick = 2
	oPlot, eqdsk.rbbbs, eqdsk.zbbbs, $
			color = 8*16-1, $
			thick = 2

	oplot, grow_rlim, grow_zlim, $
			psym = -4, $
			color = 12*16-1, $
			thick = 2
	
	!p.multi = 0

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
	jantz_id = nCdf_varDef ( nc_id, 'jantz', [nR_id, nz_id], /float )
	limSigma_id = nCdf_varDef ( nc_id, 'limiter_sigma', [nR_id, nz_id], /float )

	nCdf_control, nc_id, /enDef
	
	nCdf_varPut, nc_id, R_id, antGrid_x
	nCdf_varPut, nc_id, z_id, antGrid_y 
	nCdf_varPut, nc_id, jantx_id, antJX_grid 
	nCdf_varPut, nc_id, janty_id, antJY_grid
	nCdf_varPut, nc_id, jantz_id, antJZ_grid
	nCdf_varPut, nc_id, limSigma_id, smoothSigma

	nCdf_close, nc_id

	stop

end
