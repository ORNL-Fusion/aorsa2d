function ar2_read_vorpal, runFolderName=runFolderName, oneD = oneD, freq = freq

	@constants
	plotB0 = 0
	plotDensity = 0
	plotE = 1
	plotSource = 0
	plotjP = 1

    RunId = 'output/test'
	if keyword_set(runFolderName) then RunId = expand_path(runFolderName)+'/output/test'
    if not keyword_set(freq) then freq = 53.0e6
    files_edgeE = file_search(RunId+'_edgeE_*')
    files_rfSource = file_search(RunId+'_rfSource_*')
    files_Jp0 = file_search(RunId+'_linearJelectron_*')
    files_Jp1 = file_search(RunId+'_linearJion1_*')
    files_globals = file_search(RunId+'_Globals_*')

    nF = n_elements(files_edgeE)
    if n_elements(files_globals) ne nF then stop

    d=h5_parse(RunId+'_universe_1.h5',/read)
    g = d.globalGridGlobal
    nCells = g.vsNumCells._data
    lBounds = g.vsLowerBounds._data
    uBounds = g.vsUpperBounds._data

    dGrid = (uBounds - lBounds)/nCells
    xGrid = fIndGen(nCells[0]+1)*dGrid[0]+lBounds[0] 
	if not oneD then begin
    yGrid = fIndGen(nCells[1]+1)*dGrid[1]+lBounds[1] 
    zGrid = fIndGen(nCells[2]+1)*dGrid[2]+lBounds[2] 
	endif

    files_B0 = file_search(RunId+'_B0_*')
    b=h5_parse(files_B0[0],/read)
    b0 = b.b0._data

    files_density_e = file_search(RunId+'_density0electron_*')
    _n0=h5_parse(files_density_e[0],/read)
    n0 = _n0.density0electron._data

    files_density_e = file_search(RunId+'_density0ion1_*')
    _n1=h5_parse(files_density_e[0],/read)
    n1 = _n1.density0ion1._data


    for f=0,nF-1 do begin 
        data = h5_parse(files_edgeE[f],/read_data)
        result = data.edgeE._data
 
        if nX eq !null then begin
			if not oneD then begin
            nX = n_elements(result[0,0,0,*])
			endif else begin
            nX = n_elements(result[0,*])
			endelse
            time = dblArr(nF)

            eEdgeX = fltArr(nX,nF)
            eEdgeY = fltArr(nX,nF)
            eEdgeZ = fltArr(nX,nF)

            rfSourceX = fltArr(nX,nF)
            rfSourceY = fltArr(nX,nF)
            rfSourceZ = fltArr(nX,nF)

            jP0X = fltArr(nX,nF)
            jP0Y = fltArr(nX,nF)
            jP0Z = fltArr(nX,nF)

	        jP1X = fltArr(nX,nF)
            jP1Y = fltArr(nX,nF)
            jP1Z = fltArr(nX,nF)
	
			if not oneD then begin
            nY = n_elements(result[0,0,*,0])
            nZ = n_elements(result[0,*,0,0])
            ySlice = nY/2-1
            zSlice = nZ/2-1
			endif
        endif
    
        time[f] = data.time.vstime._data
		if not oneD then begin
        eEdgeX[*,f] = result[0,zSlice,ySlice,*]
        eEdgeY[*,f] = result[1,zSlice,ySlice,*]
        eEdgeZ[*,f] = result[2,zSlice,ySlice,*]
		endif else begin
        eEdgeX[*,f] = (result[0,*])[*]
        eEdgeY[*,f] = (result[1,*])[*]
        eEdgeZ[*,f] = (result[2,*])[*]
		endelse

        data = h5_parse(files_rfSource[f],/read_data)
        result = data.rfSource._data

		if not oneD then begin
        rfSourceX[*,f] = result[0,zSlice,ySlice,*]
        rfSourceY[*,f] = result[1,zSlice,ySlice,*]
        rfSourceZ[*,f] = result[2,zSlice,ySlice,*]
		endif else begin
        rfSourceX[*,f] = (result[0,*])[*]
        rfSourceY[*,f] = (result[1,*])[*]
        rfSourceZ[*,f] = (result[2,*])[*]
		endelse

        data = h5_parse(files_Jp0[f],/read_data)
        result = data.linearJelectron._data

		if not oneD then begin
        Jp0X[*,f] = result[0,zSlice,ySlice,*]
        Jp0Y[*,f] = result[1,zSlice,ySlice,*]
        Jp0Z[*,f] = result[2,zSlice,ySlice,*]
		endif else begin
        Jp0X[*,f] = (result[0,*])[*]
        Jp0Y[*,f] = (result[1,*])[*]
        Jp0Z[*,f] = (result[2,*])[*]
		endelse

        data = h5_parse(files_Jp1[f],/read_data)
        result = data.linearJion1._data

		if not oneD then begin
        Jp1X[*,f] = result[0,zSlice,ySlice,*]
        Jp1Y[*,f] = result[1,zSlice,ySlice,*]
        Jp1Z[*,f] = result[2,zSlice,ySlice,*]
		endif else begin
        Jp1X[*,f] = (result[0,*])[*]
        Jp1Y[*,f] = (result[1,*])[*]
        Jp1Z[*,f] = (result[2,*])[*]
		endelse

    endfor

	if not oneD then begin
		n0 = n0[0,zSlice,ySlice,*]
		n1 = n1[0,zSlice,ySlice,*]
		b0 = b0[*,zSlice,ySlice,*]
	endif

	if plotDensity then begin
		p=plot(xGrid,n0,title='density')
		p=plot(xGrid,n1,/over,color='r')
	endif
	if plotB0 then begin
    	p=plot(xGrid,b0[0,*])
    	p=plot(xGrid,b0[1,*],/over,color='red')
    	p=plot(xGrid,b0[2,*],/over,color='blue')
		p.save, "v_b0.png", resolution=300
	endif

    SortII = sort(time)
    time = time[SortII]
    for i=0,nX-1 do begin

        eEdgeX[i,*] = eEdgeX[i,SortII]
        eEdgeY[i,*] = eEdgeY[i,SortII]
        eEdgeZ[i,*] = eEdgeZ[i,SortII]

        rfSourceX[i,*] = rfSourceX[i,SortII]
        rfSourceY[i,*] = rfSourceY[i,SortII]
        rfSourceZ[i,*] = rfSourceZ[i,SortII]

        jP0X[i,*] = jP0X[i,SortII]
        jP0Y[i,*] = jP0Y[i,SortII]
        jP0Z[i,*] = jP0Z[i,SortII]

        jP1X[i,*] = jP1X[i,SortII]
        jP1Y[i,*] = jP1Y[i,SortII]
        jP1Z[i,*] = jP1Z[i,SortII]
 
    endfor

    ; Create the freq domain version

    eEdgeX_freq = complexArr(nX)
    eEdgeY_freq = complexArr(nX)
    eEdgeZ_freq = complexArr(nX)

    rfSourceX_freq = complexArr(nX)
    rfSourceY_freq = complexArr(nX)
    rfSourceZ_freq = complexArr(nX)

    jP0X_freq = complexArr(nX)
    jP0Y_freq = complexArr(nX)
    jP0Z_freq = complexArr(nX)

    jP1X_freq = complexArr(nX)
    jP1Y_freq = complexArr(nX)
    jP1Z_freq = complexArr(nX)

    iiSteadyState = where(time gt 0);3.0*max(time)/4.0)
    for i=0,nX-1 do begin
        eEdgeX_freq[i] = ar2_time_to_freq(eEdgeX[i,iiSteadyState],time[iiSteadyState],freq,i=i)
        eEdgeY_freq[i] = ar2_time_to_freq(eEdgeY[i,iiSteadyState],time[iiSteadyState],freq,i=i)
        eEdgeZ_freq[i] = ar2_time_to_freq(eEdgeZ[i,iiSteadyState],time[iiSteadyState],freq,i=i)

        rfSourceX_freq[i] = ar2_time_to_freq(rfSourceX[i,iiSteadyState],time[iiSteadyState],freq,i=i)
        rfSourceY_freq[i] = ar2_time_to_freq(rfSourceY[i,iiSteadyState],time[iiSteadyState],freq,i=i)
        rfSourceZ_freq[i] = ar2_time_to_freq(rfSourceZ[i,iiSteadyState],time[iiSteadyState],freq,i=i);,/PlotSpec, NStop=126)

        jP0X_freq[i] = e0*ar2_time_to_freq(jP0X[i,iiSteadyState],time[iiSteadyState],freq,i=i)
        jP0Y_freq[i] = e0*ar2_time_to_freq(jP0Y[i,iiSteadyState],time[iiSteadyState],freq,i=i);,/PlotSpec, NStop=120)
        jP0Z_freq[i] = e0*ar2_time_to_freq(jP0Z[i,iiSteadyState],time[iiSteadyState],freq,i=i)

        jP1X_freq[i] = e0*ar2_time_to_freq(jP1X[i,iiSteadyState],time[iiSteadyState],freq,i=i)
        jP1Y_freq[i] = e0*ar2_time_to_freq(jP1Y[i,iiSteadyState],time[iiSteadyState],freq,i=i);,/PlotSpec, NStop=120)
        jP1Z_freq[i] = e0*ar2_time_to_freq(jP1Z[i,iiSteadyState],time[iiSteadyState],freq,i=i)
 
    endfor

	jPx_freq = jP0X_freq + jP1X_freq 
	jPy_freq = jP0Y_freq + jP1Y_freq 
	jPz_freq = jP0Z_freq + jP1Z_freq 

    ;; There is a factor of pi mismatch in the amplitude of the E and Jp relative to 
    ;; Ja when compared with AORSA. Not sure which (AORSA or VORPAL) is correct, so 
    ;; just applying a factor of pi correction here. It is ONLY on the E and Jp, i.e.,
    ;; the Ja compare correctly.

    ;CorrectionFactor = 1.0/!pi

    ;eEdgeX_freq = eEdgeX_freq * CorrectionFactor
    ;eEdgeY_freq = eEdgeY_freq * CorrectionFactor
    ;eEdgeZ_freq = eEdgeZ_freq * CorrectionFactor

    ;jPX_freq = jPX_freq * CorrectionFactor
    ;jPY_freq = jPY_freq * CorrectionFactor
    ;jPZ_freq = jPZ_freq * CorrectionFactor

    ;; There is also a phase difference that I'm as yet unable to account for when comparing 
    ;; E and Jp relative to Ja when comparing AORSA and VORPAL. Again, just correcting for 
    ;; it here until I figure our where the cause is.

    ;eEdgeX_freq = -conj(eEdgeX_freq) 
    ;eEdgeY_freq = -conj(eEdgeY_freq) 
    ;eEdgeZ_freq = -conj(eEdgeZ_freq) 

    ;jPX_freq = -conj(jPX_freq) 
    ;jPY_freq = -conj(jPY_freq) 
    ;jPZ_freq = -conj(jPZ_freq) 

	if plotE then begin
	xRange = [1.0,3.0]
    range=1.0
    p=plot(xGrid,eEdgeX_freq,layout=[1,3,1],xRange=xRange,title='Vorpal eEdge',$
            yRange=[-range,range], yTitle='Ex [V/m]')
    p=plot(xGrid,imaginary(eEdgeX_freq),/over,color='r')
    range=0.15
    p=plot(xGrid,eEdgeY_freq,layout=[1,3,2],/current,xRange=xRange,$
            yRange=[-range,range], yTitle='Ey [V/m]')
    p=plot(xGrid,imaginary(eEdgeY_freq),/over,color='r')
    range=0.4
    p=plot(xGrid,eEdgeZ_freq,layout=[1,3,3],/current,xRange=xRange,$
            yRange=[-range,range], yTitle='Ez [V/m]')
    p=plot(xGrid,imaginary(eEdgeZ_freq),/over,color='r')
	p.save, "v_eEdge.png", resolution=300
	endif

	if plotSource then begin
    p=plot(xGrid,rfSourceX_freq,layout=[1,3,1])
    p=plot(xGrid,imaginary(rfSourceX_freq),/over,color='r')
    p=plot(xGrid,rfSourceY_freq,layout=[1,3,2],/current)
    p=plot(xGrid,imaginary(rfSourceY_freq),/over,color='r')
    p=plot(xGrid,rfSourceZ_freq,layout=[1,3,3],/current)
    p=plot(xGrid,imaginary(rfSourceZ_freq),/over,color='r')
	p.save, "v_rfSource.png", resolution=300
	endif

	if plotjP then begin
    p=plot(xGrid,jPX_freq,layout=[1,3,1])
    p=plot(xGrid,imaginary(jPX_freq),/over,color='r')
    p=plot(xGrid,jPY_freq,layout=[1,3,2],/current)
    p=plot(xGrid,imaginary(jPY_freq),/over,color='r')
    p=plot(xGrid,jPZ_freq,layout=[1,3,3],/current)
    p=plot(xGrid,imaginary(jPZ_freq),/over,color='r')
	p.save, "v_jP.png", resolution=300
	endif

    solution = { $
            r: xGrid, $
            x: xGrid, $
			z: xGrid*0, $
            jP_r: jPx_freq, $
            jP_t: jPy_freq, $
            jP_z: jPz_freq, $
            E_r: eEdgeX_freq, $
            E_t: eEdgeY_freq, $
            E_z: eEdgeZ_freq, $
            jA_r: rfSourceX_freq, $
            jA_t: rfSourceY_freq, $
            jA_z: rfSourceZ_freq }
            
    return, solution
end
