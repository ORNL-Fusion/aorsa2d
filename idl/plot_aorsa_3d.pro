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

	;	get cartesian coords

	nR	=	n_elements (R)
	nz_	=  	n_elements (z_)
	nPhi	= n_elements (phi)

	R3D	= rebin ( R, nR, nz_, nPhi )
	z3D	= transpose ( rebin ( z_, nz_, nPhi, nR ), [2,0,1] )
	phi3D	= transpose ( rebin ( phi, nPhi, nR, nz_ ), [1,2,0] )

	x	= R3D * cos ( phi3D )
	y	= R3D * sin ( phi3D )
	z	= z3D

	nPts	= n_elements ( x[*] )

	xRange	= 1.7
	yRange	= 1.7
	zRange	= 2.5

	nX	= 101
	nY	= 101
	nz	= 101
	cartData	= fltArr ( nX, nY, nZ )
	cartCnt		= intArr ( nX, nY, nZ )

	cartX	= ( fIndGen ( nX ) - nX / 2 ) / ( nX / 2 ) * xRange
	cartY	= ( fIndGen ( nY ) - nY / 2 ) / ( nY / 2 ) * yRange 
	cartZ	= ( fIndGen ( nZ ) - nZ / 2 ) / ( nZ / 2 ) * zRange 


	;	plot toroidal slice

	;	create an interpolated limiter boundary

	newR	= (eqdsk.rLim)[0]
	newZ	= (eqdsk.zLim)[0]

	for i=0,n_elements(eqdsk.rLim)-2 do begin

		;	get slope

		m	= ( (eqdsk.zLim)[i+1]-(eqdsk.zLim)[i] ) $
				/ ( (eqdsk.rLim)[i+1] - (eqdsk.rLim)[i] )
		b	= (eqdsk.zLim)[i] - m * (eqdsk.rLim)[i]

		;	distance

		d	= sqrt ( ( (eqdsk.rLim)[i+1] - (eqdsk.rLim)[i] )^2 $
				+ ( (eqdsk.zLim)[i+1] - (eqdsk.zLim)[i] )^2 )

		;dMin	= ( capR[0] - capR[1] ) / 2.0

		;if d gt abs(dMin) then begin

			nExtra	= 10;fix ( d / abs(dMin) )
			dStep	= ((eqdsk.rLim)[i+1] - (eqdsk.rLim)[i]) / nExtra

			for j = 0, nExtra - 1 do begin

				if dStep ne 0 then begin
					newR	= [ newR, (eqdsk.rLim)[i] + dStep*j ]
					newZ	= [ newZ, m * ((eqdsk.rLim)[i] + dStep*j) + b ]
				endif

			endfor

		;endif

	endfor

	newRbbbs	= (eqdsk.rbbbs)[0]
	newZbbbs	= (eqdsk.zbbbs)[0]

	for i=0,n_elements(eqdsk.rbbbs)-2 do begin

		;	get slope

		m	= ( (eqdsk.zbbbs)[i+1]-(eqdsk.zbbbs)[i] ) $
				/ ( (eqdsk.rbbbs)[i+1] - (eqdsk.rbbbs)[i] )
		b	= (eqdsk.zbbbs)[i] - m * (eqdsk.rbbbs)[i]

		;	distance

		d	= sqrt ( ( (eqdsk.rbbbs)[i+1] - (eqdsk.rbbbs)[i] )^2 $
				+ ( (eqdsk.zbbbs)[i+1] - (eqdsk.zbbbs)[i] )^2 )

		;dMin	= ( capR[0] - capR[1] ) / 2.0

		;if d gt abs(dMin) then begin

			nExtra	= 10;fix ( d / abs(dMin) )
			dStep	= ((eqdsk.rbbbs)[i+1] - (eqdsk.rbbbs)[i]) / nExtra

			for j = 0, nExtra - 1 do begin

				if dStep ne 0 then begin
					newRbbbs	= [ newRbbbs, (eqdsk.rbbbs)[i] + dStep*j ]
					newZbbbs	= [ newZbbbs, m * ((eqdsk.rbbbs)[i] + dStep*j) + b ]
				endif

			endfor

		;endif

	endfor


	;loadct, 3, /silent
	set_plot, 'ps'
	!p.multi = [0,2,1]
	nLevs = 21
	if(not keyword_set(range)) then range = max ( eMod[*,100,*] ) / 2.0

	b3Phi	=fIndGen (360) *!dtor 	
	b3R		=fltArr(360)+min(eqdsk.rlim)
	b3X		= b3R * cos ( b3Phi )
	b3Y		= b3R * sin ( b3Phi )

	levels	= fIndGen ( nLevs ) / nLevs * range 
	colors	= 255 - ( bytScl ( levels, top = 253 ) + 1 )

	iiRight	= where ( newR gt 1.2 )
	useZ	= newZ[iiRight]
	useR	= newR[iiRight]

	strapR	= fltArr(12) + 1.56
	sep	= 19.4e-2/(2.0*!pi*1.56)*360*!dtor
	strapPhi	= fIndGen(12)*sep-6.5*sep
	strapX	= strapR * cos(strapPhi)
	strapY	= strapR * sin(strapPhi)

	zspan =128 
	for i=170,500 do begin
			print, i
		device, fileName = 'output/movie/'+string(i,format='(i3.3)')+'aorsa_dlg_3d.eps', $
			/color, $
			bits_per_pixel = 8, $
			/encap, $
			xSize = 8.7*2, $
			ySize = 8, $
			/inc


		;	find nearest index in both bbbs and emod
		iiMod	= (((newZbbbs)[i]-min(z_))/(max(z_)-min(z_))*n_elements(z_))
		
		iiClose	= (where ( abs(useZ-newZbbbs[i]) eq min ( abs(useZ - newZbbbs[i]) ) ))[0]
		;zii	= n_elements(eMod[0,*,0])/2-zspan/2 + zspan/(6-1)*i
		loadct, 3, /sil
		polar_contour, transpose(reform(eMod[*,iiMod,*])), phi, R, $
				levels = levels, $
			   	c_colors = colors, $
				/fill, $
				title = 'z = '+ string ( newZbbbs[i] ), $
				color = 255, $
				xRange = [-1.4, 1.6], $
				yRange = [-1.4,1.4], $
				xSty = 1, ySty = 1
		loadct, 12, /sil

	b1Phi	=fIndGen (360) *!dtor 	
	b1R		=fltArr(360)+newrbbbs[i]
	b1X		= b1R * cos ( b1Phi )
	b1Y		= b1R * sin ( b1Phi )

	b2Phi	=fIndGen (360) *!dtor 	
	b2R		=fltArr(360)+useR[iiClose]
	b2X		= b2R * cos ( b2Phi )
	b2Y		= b2R * sin ( b2Phi )

		lineR=[0,useR[iiClose]]
		linePhi	= [phi[(i-170) mod 200],phi[(i-170) mod 200]]
		lineX = lineR * cos ( linePhi )
		lineY = lineR * sin ( linePhi )
		plots,lineX, lineY, $
				color = 14*16-1, $
				thick = 40, $
				lineStyle = 0
	
		loadct, 3, /sil
		polar_contour, transpose(reform(eMod[*,iiMod,*])), phi, R, $
				levels = levels, $
			   	c_colors = colors, $
				title = 'z = '+ string ( newZbbbs[i] ), $
				color = 255, $
				xRange = [-1.4, 1.6], $
				yRange = [-1.4,1.4], $
				xSty = 1, ySty = 1, /over
		loadct, 12, /sil


		plots, b1X, b1Y, $
				color = 8 * 16 - 1, $
				thick = 6, $
				lineStyle = 0
		plots, b2X, b2Y, $
				color = 1 * 16 - 1, $
				thick = 6, $
				lineStyle = 0
		plots, b3X, b3Y, $
				color = 1 * 16 - 1, $
				thick = 6, $
				lineStyle = 0
		plots, strapX, strapY, $
				psym = 4, $
				color = 12 * 16 - 1, $
				thick = 6
	loadct, 3, /sil
		contour, eMod[*,*,(i-170) mod 200], R,z_, $
				color = 255, $
				xRange = [0,1.7], $
				yRange = [-1.8, 1.8], $
				levels = levels, $
				c_colors = colors, $
				/fill
	loadct, 12, /sil
	plots, findGen(20)/20*(max(eqdsk.rlim)-min(eqdsk.rlim))+min(eqdsk.rlim),fltArr(20)+(newzbbbs)[i], $
			linestyle = 0, $
			color = 14*16-1, $
			thick = 40 
		loadct, 3, /sil
		contour, eMod[*,*,(i-170) mod 200], R,z_, $
				color = 255, $
				xRange = [0,1.7], $
				yRange = [-1.8, 1.8], $
				levels = levels, $
				c_colors = colors, $
				/over
	loadct, 12, /sil
	
	plots, eqdsk.rbbbs, eqdsk.zbbbs, $
			color = 8*16-1, $
			thick = 6
	plots, eqdsk.rlim, eqdsk.zlim, $
			color = 0, $
			thick = 6
device, /close
	endfor


	;device, /close

	;for i=0,nR-1 do begin
	;	print, i, nR-1
	;	for j=0,nz_-1 do begin
	;		for k=0,nPhi-1 do begin
	;				
	;			iIndex	=  ( x[i,j,k] - min ( cartX ) ) / ( max(cartX)-min(cartX) ) * nX
	;			jIndex	=  ( y[i,j,k] - min ( cartY ) ) / ( max(cartY)-min(cartY) ) * nY
	;			kIndex	=  ( z[i,j,k] - min ( cartZ ) ) / ( max(cartZ)-min(cartZ) ) * nZ

	;			cartData[iIndex,jIndex,kIndex]	+= eMod[i,j,k]
	;			cartCnt[iIndex,jIndex,kIndex]	++

	;		endfor
	;	endfor
	;endfor

	;iiCnt	= where ( cartCnt gt 0 )
	;cartData[iiCnt]	= cartData[iiCnt] / cartCnt[iiCnt]

	stop
	restore, 'cartData.sav'

	;	use object graphics

	loadct, 3, rgb_table = orange
	orange	= reverse ( orange )

	data	= byte ( ( bytScl ( cartData, top = 255, max = 1e4, min = 1e0	) ) )

	myVol	= obj_new ( 'IDLgrVolume', data, interpolate = 0 )
	xAxis	= obj_new ( 'IDLgrAxis', 0 )	
	yAxis	= obj_new ( 'IDLgrAxis', 1 )	
	zAxis	= obj_new ( 'IDLgrAxis', 2 )	

	myVol -> setProperty, xCoord_conv = [-xRange, xRange/nX*2.0], $
			yCoord_conv = [-yRange, yRange / nY*2.0], $
			zCoord_conv = [-zRange, zRange / nZ*2.0	]

	myVol -> SetProperty, ZERO_OPACITY_SKIP=1
	myVol -> SetProperty, ZBUFFER=1

	myVol -> getProperty, xRange = xr, yRange = yr, zRange = zr
	xAxis -> setProperty, range = xr * xRange / nX 
	yAxis -> setProperty, range = yr * yRange / nY 
	zAxis -> setProperty, range = zr * zRange / nZ 

	myVol -> setProperty, rgb_table0 = orange

	;power = 0.5
	;opac = reverse ( fIndGen(256)^power/(255.0^power)*255 )
	;opac[0:100]=50
	;opac[101:200]=20
	;opac[201:220]=5
	;opac[221:*]=1
	;opac[255] = 0
	opac = findGen(256)*0+50
	opac[0:50]=0
	myVol -> SetProperty, OPACITY_TABLE0 = opac
	;myVol -> setProperty, composite_function = 1

	myModel	= obj_new ( 'IDLgrModel' )
	myView	= obj_new ( 'IDLgrView' )
	windowSize = [ 400, 400 ]
	myWin	= obj_new ( 'IDLgrWindow', dim = windowSize, retain = 1, qual = 2 )

	;myModel -> add, xAxis
	;myModel -> add, yAxis	
	;myModel -> add, zAxis

	myView -> add, myModel

	;	create limiter polyLine

	limPhi1	= !pi*0/6 
	limX1	= eqdsk.rlim * cos ( limPhi1 )
	limY1	= eqdsk.rlim * sin ( limPhi1 )
	limZ1	= eqdsk.zlim
	myLim1	= obj_new ( 'IDLgrPolyline', [ limx1, limy1, limz1 ] )
	myLim1 -> setProperty, color = transpose ( orange[50,*] )
	;myModel -> add, myLim1

	limPhi2	= !pi*1/6
	limX2	= eqdsk.rlim * cos ( limPhi2 )
	limY2	= eqdsk.rlim * sin ( limPhi2 )
	limZ2	= eqdsk.zlim
	myLim2	= obj_new ( 'IDLgrPolyline', [ limx2, limy2, limz2 ] )
	myLim2 -> setProperty, color = transpose ( orange[50,*] )
	;myModel -> add, myLim2

	limPhi3	= !pi*2/6
	limX3	= eqdsk.rlim * cos ( limPhi3 )
	limY3	= eqdsk.rlim * sin ( limPhi3 )
	limZ3	= eqdsk.zlim
	myLim3	= obj_new ( 'IDLgrPolyline', [ limx3, limy3, limz3 ] )
	myLim3 -> setProperty, color = transpose ( orange[50,*] )
	;myModel -> add, myLim3

	limPhi4	= !pi*3/6
	limX4	= eqdsk.rlim * cos ( limPhi4 )
	limY4	= eqdsk.rlim * sin ( limPhi4 )
	limZ4	= eqdsk.zlim
	myLim4	= obj_new ( 'IDLgrPolyline', [ limx4, limy4, limz4 ] )
	myLim4 -> setProperty, color = transpose ( orange[50,*] )
	;myModel -> add, myLim4

	nPoly = n_elements ( limX1 )-1
	for i=0,nPoly do begin
		if size(vessel,/type) eq 0 then $
			vessel	= [limx1[i],limy1[i],limz1[i]] $
		else $
			vessel	= [ [vessel], [limx1[i],limy1[i],limz1[i]] ]
		vessel	= [ [vessel], [limx2[i],limy2[i],limz2[i]] ]
		vessel	= [ [vessel], [limx3[i],limy3[i],limz3[i]] ]
		vessel	= [ [vessel], [limx4[i],limy4[i],limz4[i]] ]
		
	endfor

	for i=0,nPoly-1 do begin

		if size(polygons,/type) eq 0 then $
			polygons = [4,i*4 + [0,1,5,4]] $
		else $
			polygons	= [ polygons, [4,i*4 + [0,1,5,4]] ]

		polygons	= [ polygons, [4,i*4 + [1,2,6,5]] ]
		polygons	= [ polygons, [4,i*4 + [2,3,7,6]] ]

	endfor

	loadct, 1, rgb_table = blue
	myAnt1	= obj_new ( 'IDLgrPolygon', vessel, polygons = polygons, style = 2 )
	myAnt1 -> setProperty, color = transpose(blue[200,*])
	myModel -> add, myAnt1

	;	take a chunck out of the lcfs for antenna

	wd	= 19.4e-2
	circ	= 2.0 * !Pi * max ( eqdsk.rbbbs )
	phiSep	= wd / circ * 2.0 * !pi
	rShift	= 1.50 - max ( eqdsk.rbbbs )

	for i=0,11 do begin
		iiAntBbbs	= where ( eqdsk.rbbbs gt 1.38 )
		antR	= (eqdsk.rbbbs)[iiAntBbbs]+rShift
		antz	= (eqdsk.zbbbs)[iiAntBbbs]
		antP	= antR * 0.0 + i*phiSep -phiSep*6

		antX	= antR * cos ( antP )
		antY	= antR * sin ( antP )

		myAntA	= obj_new ( 'IDLgrPolyline', transpose ( [ [antx], [anty], [antz] ] ) )
		myAntA -> setProperty, color = transpose ( orange[255,*] ), thick = 10, alpha_channel=0.5 
		myModel -> add, myAntA
	endfor


	;	add contours

	;contour, data[*,*,nZ/2], cartx, carty, $
	;		path_xy = contours, $
	;		path_info = contour_info, $
	;		levels = fIndGen(21)/21*255, $
	;		/path_data_coords

	;for i=0,n_elements(contour_info)-1 do begin
	;	this = contour_info[i]
	;	print, this.n, this.value, this.type
	;	if this.type ne 0 then begin
	;		cX	= contours[0,this.offset:this.offset+this.n-1]
	;		cY	= contours[1,this.offset:this.offset+this.n-1]
	;		cZ	= cX * 0 + cartz[nZ/2]
	;		cData	= transpose ( [ [[cX]],[[cY]],[[cZ]] ] )
	;		myAntA	= obj_new ( 'IDLgrPolyline', cData )
	;		myModel -> add, myAntA
	;	endif

	;endfor

	for j=0,19 do begin

	contour, reform(eMod[*,*,j*10]), R, z_, $
			path_xy = contours, $
			path_info = contour_info, $
			/path_data_coords, $
			color = 0, $
			levels	= (fIndGen (20)+1) / 20 * 4e4

	for i=0,n_elements(contour_info)-1 do begin
		this = contour_info[i]
		if this.type ne 0 then begin
			cX	= contours[0,this.offset:this.offset+this.n-1]*cos(phi[j*10])
			cY	= contours[0,this.offset:this.offset+this.n-1]*sin(phi[j*10])
			cZ	= contours[1,this.offset:this.offset+this.n-1]
			cData	= transpose ( [ [[cX]],[[cY]],[[cZ]] ] )
			myAntA	= obj_new ( 'IDLgrPolyline', cData )
			myAntA -> setProperty, color = orange[50];[255-(bytScl(this.value,top=253,max=1e3,min=0)+1)]
			myAntA -> setProperty, alpha_channel = (bytScl(this.value,top=253,max=1e4,min=0)/253)>0.5
			;myModel -> add, myAntA
		endif

	endfor

	endfor

	myLight = obj_new ( 'IDLgrLight', type = 1, location = [3,0,0] )
	myModel -> add, myLight
	myLight = obj_new ( 'IDLgrLight', type = 1, location = [-1,1,1] )
	;myModel -> add, myLight


	myModel -> add, myVol
	myModel -> rotate, [0,1,0], 90 
	myModel -> rotate, [1,0,0], 0 
	myModel -> rotate, [0,0,1], 90 
	set_view, myView, myWin
	myWin -> draw, myView
	
	for i =0, 360 do begin
	;myModel -> rotate, [0,1,0], 1
	xLoc	= 3 * cos ( i * !dtor )
	yLoc	= 3 * sin ( i * !dtor )
	myLight -> setProperty, location = [xLoc,yLoc,0]
	myWin -> draw, myView
	endfor

	myWin -> getProperty, resolution = screenResolution

	oClip	= obj_new ( 'IDLgrClipboard', quality = 2, $
		   graphics_tree = myView, $
		  	resolution = screenResolution, $
		   	dimensions = windowSize	)
	fileName	= 'testA.ps'
	oClip -> draw, fileName = fileName, vector=0, post = 0 
	oClip -> setProperty, graphics_tree = obj_new ()
	obj_destroy, oClip
	stop
xObjView, myModel, renderer = 1
end
