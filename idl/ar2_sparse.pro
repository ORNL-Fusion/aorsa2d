function testSparse, nx, ny, xMin, xMax, yMin, yMax, x2d=x, y2d=y

    xRng = (xMax-xMin)
    xGrid = fIndGen(nx)/(nx-1)*xRng+xMin 

    yRng = (yMax-yMin)
    yGrid = fIndGen(ny)/(ny-1)*yRng+yMin 

    x = rebin(xGrid,nx,ny)
    y = transpose(rebin(yGrid,ny,nx))

    sigx = 0.1
    sigy = 0.1
    x0 = 0
    y0 = 0
    A = 0.1

    return, A * exp( -( (x-x0)^2/(2*sigx^2) + (y-y0)^2/(2*sigy^2) ) ) 

end

function testSparse2, nx, ny, xMin, xMax, yMin, yMax, x2d=x, y2d=y

    referenceSolution = expand_path("~/scratch/aorsa2d/sparse-template-128x128")
    refSolutionRaw = ar2_read_solution(referenceSolution,1) 

    xRng = (xMax-xMin)
    xGrid = fIndGen(nx)/(nx-1)*xRng+xMin 

    yRng = (yMax-yMin)
    yGrid = fIndGen(ny)/(ny-1)*yRng+yMin 

    x = rebin(xGrid,nx,ny)
    y = transpose(rebin(yGrid,ny,nx))

    ii = (x - refSolutionRaw.r[0]) / (refSolutionRaw.r[-1]-refSolutionRaw.r[0]) * (n_elements(refSolutionRaw.r)-1)
    jj = (y - refSolutionRaw.z[0]) / (refSolutionRaw.z[-1]-refSolutionRaw.z[0]) * (n_elements(refSolutionRaw.z)-1)

    refSolution = interpolate(refSolutionRaw.e_t,ii,jj)

    return, refSolution 

end

pro ar2_sparse

@constants

readOnly = 1
testFunction = 0
testFunction2 = 0

if testFunction then begin
    scale = 0.1 
endif else begin
    scale = 10.0
endelse

; Test application of sparse grids to AORSA

d = 2
l = 6
lMin = 2

nX = [256,128,128,64,64,32,32,16,16, 8,  8,  4,  4]-1
nY = [  4,  4,  8, 8,16,16,32,32,64,64,128,128,256]-1
i = alog2(nX+1)
j = alog2(nY+1)
sgn = [+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1]

nEval = n_elements(nX)

template = expand_path("~/scratch/aorsa2d/sparse-template")
referenceSolution = expand_path("~/scratch/aorsa2d/sparse-template-128x128")

cd, current = rootDir

runAR2Cmd = 'mpirun -n 24 ~/code/aorsa2d/xaorsa2d'

runList = []

; Setup and run AORSA instances

for run=0,nEval-1 do begin

    ; Stage run
    
    thisRun = "sparse-test-"+string(run,format='(i3.3)')
    runList = [runList,thisRun]

    if not readOnly then begin

        print, 'Staging run ... '+thisRun
        file_copy, template, thisRun, /recursive
        
        ; Set the grid resolution for this run
        
        cd, thisRun
        
        this_nR = nX[run] 
        sedCmd = "sed -i .in-bak 's/nRAll(1) = 16/nRAll(1) = "+string(this_nR,format='(i3.3)')+"/' aorsa2d.in"
        spawn, sedCmd
        this_nZ = nY[run] 
        sedCmd = "sed -i .in-bak 's/nZAll(1) = 128/nZAll(1) = "+string(this_nZ,format='(i3.3)')+"/' aorsa2d.in"
        spawn, sedCmd
        
        spawn, runAR2Cmd
        
        cd, rootDir

    endif

endfor

; Build the desired grid

nRGrid = 256 
nZGrid = 256

if testFunction then begin

    rMin = -1.0      
    rMax = +1.0
    rRng = (rMax-rMin)
    rGrid = fIndGen(nRGrid)/(nRGrid-1)*rRng+rMin 

    zMin = -1.0      
    zMax = +1.0
    zRng = (zMax-zMin)
    zGrid = fIndGen(nZGrid)/(nZGrid-1)*zRng+zMin 

    r2D = rebin(rGrid,nRGrid,nZGrid)
    z2D = transpose(rebin(zGrid,nZGrid,nRGrid))

endif else begin

    thisSolution = ar2_read_solution(runList[0],1) 
        
    rMin = thisSolution.r[0]      
    rMax = thisSolution.r[-1]
    rRng = (rMax-rMin)
    rGrid = fIndGen(nRGrid)/(nRGrid-1)*rRng+rMin 

    zMin = thisSolution.z[0]      
    zMax = thisSolution.z[-1]
    zRng = (zMax-zMin)
    zGrid = fIndGen(nZGrid)/(nZGrid-1)*zRng+zMin 

    r2D = rebin(rGrid,nRGrid,nZGrid)
    z2D = transpose(rebin(zGrid,nZGrid,nRGrid))

endelse

; Evaluate reference solution

    if testFunction then begin
       
        refSolution = testSparse(nRGrid,nZGrid,rMin,rMax,zMin,zMax) 

    endif else begin

        refSolutionRaw = ar2_read_solution(referenceSolution,1) 
    
        ;scaleFac = max(abs(refSolutionRaw.e_r))
        ;print, 'ref scalefac: ',scaleFac

        ;refSolutionRaw.e_r = refSolutionRaw.e_r / scaleFac

        ii = (r2D - refSolutionRaw.r[0]) / (refSolutionRaw.r[-1]-refSolutionRaw.r[0]) * (n_elements(refSolutionRaw.r)-1)
        jj = (z2D - refSolutionRaw.z[0]) / (refSolutionRaw.z[-1]-refSolutionRaw.z[0]) * (n_elements(refSolutionRaw.z)-1)

        refSolution = interpolate(refSolutionRaw.e_t,ii,jj)

    endelse

; Evaluate the function on the sparse grid

f = complexArr(nRGrid,nZGrid)

for thisEval=0,nEval-1 do begin

    if testFunction then begin

        this_nx = nX[thisEval]
        this_ny = nY[thisEval]       

        this_f = testSparse(this_nx,this_ny,rMin,rMax,zMin,zMax,x2d=this_r2D,y2d=this_z2D) 
        this_r = (this_r2D[*,0])[*]
        this_z = (this_z2D[0,*])[*]

    endif else if testFunction2 then begin

        this_nx = nX[thisEval]
        this_ny = nY[thisEval]       

        this_f = testSparse2(this_nx,this_ny,rMin,rMax,zMin,zMax,x2d=this_r2D,y2d=this_z2D) 
        this_r = (this_r2D[*,0])[*]
        this_z = (this_z2D[0,*])[*]

    endif else begin

        thisSolution = ar2_read_solution(runList[thisEval],1) 

        ;scaleFac = max(abs(thisSolution.e_r))
        ;print, 'this scalefac: ',scaleFac

        ;thisSolution.e_r = thisSolution.e_r / scaleFac

        this_f = thisSolution.e_t
        this_r = thisSolution.r
        this_z = thisSolution.z

        print, runList[thisEval]

    endelse

    ; Interpolate to full grid

    this_nR = n_elements(this_f[*,0])
    this_nZ = n_elements(this_f[0,*])

    ii = (r2D - rGrid[0]) / (rMax-rMin) * (this_nR-1)
    jj = (z2D - zGrid[0]) / (zMax-zMin) * (this_nZ-1)

    thisSolutionFullGrid = interpolate(this_f,ii,jj)

    f = f + sgn[thisEval] * thisSolutionFullGrid

    nLevs = 10

	thisField = this_f[*,*]
    title = 'f (eval)'
	levels = fIndGen(nLevs)/(nLevs-1)*scale
	colors = reverse(bytScl(levels, top=253)+1)
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
    aspect = 1.0
    wTitle = string(thisEval)+', '+string(sgn[thisEval])+', '+string(i[thisEval])+', '+string(j[thisEval])
    if testFunction then begin
	    s = surface ( PlotField, this_r, this_z, style=1, aspect_ratio=aspect, title=title, layout=[3,2,1], $
           window_title=wTitle )
    endif else begin
	    c = contour ( PlotField, this_r, this_z, c_value=levels, rgb_indices=colors, rgb_table=3, $
            aspect_ratio=aspect, dimensions=dimensions, title=title, layout=[3,2,1], c_label_show=0, window_title=wTitle, /fill )
	    c = contour ( -PlotField, this_r, this_z, c_value=levels, rgb_indices=colors, rgb_table=1, /over, c_label_show=0, /fill  )
    endelse

	thisField = thisSolutionFullGrid[*,*]
    title = 'f (interp)'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
    if testFunction then begin 
        s = surface ( PlotField, rGrid, zGrid, style=1, aspect_ratio=aspect, title=title, layout=[3,2,2],/current )
    endif else begin
	    c = contour ( PlotField, rGrid, zGrid, c_value=levels, rgb_indices=colors, rgb_table=3, $
            aspect_ratio=aspect, dimensions=dimensions, title=title, xRange=xRange, yRange=yRange, $
            layout=[3,2,2], /current, c_label_show=0, /fill  )
	    c = contour ( -PlotField, rGrid, zGrid, c_value=levels, rgb_indices=colors, rgb_table=1, /over, c_label_show=0, /fill  )
    endelse

    p=plot(this_r,this_f[*,this_nZ/2],layout=[3,2,5],/current)
    p=plot(rGrid,thisField[*,nZGrid/2],layout=[3,2,5],/current,/over,color='r')

    thisField = f[*,*]
    title = 'f (sparse)'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
    if testFunction then begin
        s = surface ( PlotField, rGrid, zGrid, style=1, aspect_ratio=aspect, title=title, layout=[3,2,3],/current )
    endif else begin
	    c = contour ( PlotField, rGrid, zGrid, c_value=levels, rgb_indices=colors, rgb_table=3, $
            aspect_ratio=aspect, dimensions=dimensions, title=title, xRange=xRange, yRange=yRange, $
            layout=[3,2,3], /current, c_label_show=0, /fill  )
	    c = contour ( -PlotField, rGrid, zGrid, c_value=levels, rgb_indices=colors, rgb_table=1, /over, c_label_show=0, /fill  )
    endelse

    thisField = refSolution-f 
    title = 'fRef - fSparse'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
    if testFunction then begin
        s = surface ( PlotField, rGrid, zGrid, style=1, aspect_ratio=aspect, title=title, layout=[3,2,4],/current )
    endif else begin
	    c = contour ( PlotField, rGrid, zGrid, c_value=levels, rgb_indices=colors, rgb_table=3, $
            aspect_ratio=aspect, dimensions=dimensions, title=title, xRange=xRange, yRange=yRange, $
            layout=[3,2,4], /current, c_label_show=0, /fill  )
	    c = contour ( -PlotField, rGrid, zGrid, c_value=levels, rgb_indices=colors, rgb_table=1, /over, c_label_show=0, /fill )
    endelse

    thisField = refSolution
    title = 'fRef'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
    if testFunction then begin
        s = surface ( PlotField, rGrid, zGrid, style=1, aspect_ratio=aspect, title=title, layout=[3,2,6],/current )
    endif else begin
	    c = contour ( PlotField, rGrid, zGrid, c_value=levels, rgb_indices=colors, rgb_table=3, $
            aspect_ratio=aspect, dimensions=dimensions, title=title, xRange=xRange, yRange=yRange, $
            layout=[3,2,6], /current, c_label_show=0, /fill  )
	    c = contour ( -PlotField, rGrid, zGrid, c_value=levels, rgb_indices=colors, rgb_table=1, /over, c_label_show=0, /fill )
    endelse

endfor

stop
end
