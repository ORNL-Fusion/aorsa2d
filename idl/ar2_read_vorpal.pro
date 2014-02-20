pro ar2_read_vorpal

    RunId = 'Benchmark03'
    freq = 53.0e6
    files_edgeE = file_search(RunId+'_edgeE_*')
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
    yGrid = fIndGen(nCells[1]+1)*dGrid[1]+lBounds[1] 
    zGrid = fIndGen(nCells[2]+1)*dGrid[2]+lBounds[2] 

    b=h5_parse(RunId+'_B0_1.h5',/read)
    b0 = b.b0._data


    for f=0,nF-1 do begin 
        data = h5_parse(files_edgeE[f],/read_data)
        result = data.edgeE._data
 
        if nX eq !null then begin
            nX = n_elements(result[0,0,0,*])
            nY = n_elements(result[0,0,*,0])
            nZ = n_elements(result[0,*,0,0])
            eEdgeX = dblArr(nX,nF)
            eEdgeY = dblArr(nX,nF)
            eEdgeZ = dblArr(nX,nF)
            time = dblArr(nF)
            ySlice = nY/2-1
            zSlice = nZ/2-1
        endif
    
        time[f] = data.time.vstime._data
        eEdgeX[*,f] = result[0,zSlice,ySlice,*]
        eEdgeY[*,f] = result[1,zSlice,ySlice,*]
        eEdgeZ[*,f] = result[2,zSlice,ySlice,*]
        
    endfor

    p=plot(xGrid,b0[0,zSlice,ySlice,*])
    p=plot(xGrid,b0[1,zSlice,ySlice,*],/over,color='red')
    p=plot(xGrid,b0[2,zSlice,ySlice,*],/over,color='blue')

    SortII = sort(time)
    time = time[SortII]
    for i=0,nX-1 do begin
        eEdgeX[i,*] = eEdgeX[i,SortII]
        eEdgeY[i,*] = eEdgeY[i,SortII]
        eEdgeZ[i,*] = eEdgeZ[i,SortII]
    endfor

    ; Create the freq domain version

    eEdgeX_freq = complexArr(nX)
    eEdgeY_freq = complexArr(nX)
    eEdgeZ_freq = complexArr(nX)

    iiSteadyState = where(time gt 3.0*max(time)/4.0)
    for i=0,nX-1 do begin
        eEdgeX_freq[i] = ar2_time_to_freq(eEdgeX[i,iiSteadyState],time[iiSteadyState],freq)
        eEdgeY_freq[i] = ar2_time_to_freq(eEdgeY[i,iiSteadyState],time[iiSteadyState],freq)
        eEdgeZ_freq[i] = ar2_time_to_freq(eEdgeZ[i,iiSteadyState],time[iiSteadyState],freq)
    endfor

    p=plot(xGrid,eEdgeX_freq,layout=[1,3,1])
    p=plot(xGrid,imaginary(eEdgeX_freq),/over,color='r')
    p=plot(xGrid,eEdgeY_freq,layout=[1,3,2],/current)
    p=plot(xGrid,imaginary(eEdgeY_freq),/over,color='r')
    p=plot(xGrid,eEdgeZ_freq,layout=[1,3,3],/current)
    p=plot(xGrid,imaginary(eEdgeZ_freq),/over,color='r')



stop
end
