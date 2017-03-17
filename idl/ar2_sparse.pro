pro ar2_sparse

@constants

readOnly = 1

scale = 0.1 

; Test application of sparse grids to AORSA

d = 2
l = 6
lMin = 2

nX = [256,128,128,64,64,32,32,16,16, 8,  8,  4,  4]
nY = [  4,  4,  8, 8,16,16,32,32,64,64,128,128,256]
sgn = -[+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1,-1,+1]

nRuns = n_elements(nX)

stop
template = expand_path("~/scratch/aorsa2d/sparse-template")
referenceSolution = expand_path("~/scratch/aorsa2d/sparse-template-128x128")

cd, current = rootDir

runAR2Cmd = 'mpirun -n 24 ~/code/aorsa2d/xaorsa2d'

runList = []

; Setup and run AORSA instances

for run=0,nRuns-1 do begin

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

; Read solutions

nRGrid = 64 
nZGrid = 64 

eAlpha = complexArr(nRGrid,nZGrid)
eB = complexArr(nRGrid,nZGrid)

e_r = complexArr(nRGrid,nZGrid)
e_z = complexArr(nRGrid,nZGrid)

for run=0,nRuns-1 do begin

    thisSolution = ar2_read_solution(runList[run],1) 

    if size(RGrid,/type) eq 0 then begin
        
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

    endif

    eAlphak = thisSolution.eAlphak
    eBk = thisSolution.eBk

    nM = n_elements(ealphak[*,0])
    nN = n_elements(ealphak[0,*])

    mMin = -nM/2
    mMax =  nM/2

    if (nM mod 2) eq 0 then mMin = mMin+1

    nMin = -nN/2
    nMax =  nN/2

    if (nN mod 2) eq 0 then nMin = nMin+1

    ;eAlpha[*] = 0
    ;eB[*] = 0

    normFacX = 2*!pi/rRng
    normFacY = 2*!pi/zRng

    InterpFromGrid = 1

    if InterpFromGrid then begin

        ii = (r2D - rGrid[0]) / (rMax-rMin) * (n_elements(thisSolution.r)-1)
        jj = (z2D - zGrid[0]) / (zMax-zMin) * (n_elements(thisSolution.z)-1)

        thisSolutionFullGrid_r = interpolate(thisSolution.e_r,ii,jj)
        thisSolutionFullGrid_z = interpolate(thisSolution.e_z,ii,jj)

        e_r += sgn[run] * thisSolutionFullGrid_r
        e_z += sgn[run] * thisSolutionFullGrid_z

        eAlpha = e_r
        eB = e_z

    endif else begin

        for i=0,nRGrid-1 do begin
            for j=0,nZGrid-1 do begin

                ;rNorm = (rGrid[i]-rMin)*normFacX
                ;zNorm = (zGrid[j]-zMin)*normFacY

                ; WHY is this not nRGrid-1 in AORSA?

                rNorm = i * 2 * !pi / ( nRGrid )
                zNorm = j * 2 * !pi / ( nZGrid )

                for m=mMin,mMax do begin
                    for n=nMin,nMax do begin

                        xx = exp( _ii * m * rNorm )
                        yy = exp( _ii * n * zNorm )

                        eAlpha[i,j] += sgn[run] * eAlphak[m-mMin,n-nMin] * xx * yy 
                        eB[i,j] += sgn[run] * eBk[m-mMin,n-nMin] * xx * yy 

                    endfor
                endfor

            endfor
        endfor

    endelse

    print, runList[run]
    print, 'M: ', mMin, mMax
    print, 'N: ', nMin, nMax

        nLevs = 10


		thisField = ealpha[*,*]
        title = 'E_alpha'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
        aspect = 1.0
		c = contour ( PlotField, rGrid, zGrid, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
            aspect_ratio=aspect, dimensions=dimensions, title=title, xRange=xRange, yRange=yRange )
		c = contour ( -PlotField, rGrid, zGrid, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )

		thisField = eB[*,*]
        title = 'E_b'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
        aspect = 1.0
		c = contour ( PlotField, rGrid, zGrid, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
            aspect_ratio=aspect, dimensions=dimensions, title=title, xRange=xRange, yRange=yRange )
		c = contour ( -PlotField, rGrid, zGrid, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )

		thisField = thisSolutionFullGrid_z[*,*]
        title = 'this E_b'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
        aspect = 1.0
		c = contour ( PlotField, rGrid, zGrid, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
            aspect_ratio=aspect, dimensions=dimensions, title=title, xRange=xRange, yRange=yRange )
		c = contour ( -PlotField, rGrid, zGrid, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )


endfor

refSolution = ar2_read_solution(referenceSolution,1) 

thisField = refSolution.e_z
title = 'Reference E_b'
levels = fIndGen(nLevs)/(nLevs-1)*scale
colors = reverse(bytScl(levels, top=253)+1)
PlotField = (real_part(thisField)<max(levels))>min(-levels)
aspect = 1.0
c = contour ( PlotField, refSolution.r, refSolution.z, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
    aspect_ratio=aspect, dimensions=dimensions, title=title, xRange=xRange, yRange=yRange )
c = contour ( -PlotField, refSolution.r, refSolution.z, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )

stop
end
