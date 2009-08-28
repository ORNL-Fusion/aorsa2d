pro plot_aorsa_3d

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

	myLight = obj_new ( 'IDLgrLight', type = 2 )
	myModel -> add, myLight

	myModel -> add, myVol
	myModel -> rotate, [0,1,0], 90 
	myModel -> rotate, [1,0,0], 0 
	myModel -> rotate, [0,0,1], 90 
	set_view, myView, myWin
	myWin -> draw, myView

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
