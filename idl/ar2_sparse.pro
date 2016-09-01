pro ar2_sparse

; Test application of sparse grids to AORSA

d = 2
l = 5
lMin = 2

levels = 2^(indgen(l)+lMin)

; Get max resolution to run at

if l mod 2 gt 0 then begin
    maxRes = levels[l/2]^2
endif else begin
    maxRes = levels[l/2] * levels[l/2-1]
endelse

stop
template = expand_path("~/scratch/aorsa2d/sparse-template")

cd, current = rootDir

runAR2Cmd = 'mpirun -n 24 ~/code/aorsa2d/xaorsa2d'

; Setup and run AORSA instances

for i=0,l-1 do begin

    for j=0,l-1 do begin

        if levels[i] * levels[j] le maxRes then begin

            ; Stage run

            thisRun = "sparse-test-"+string(i,format='(i3.3)')+'-'+string(j,format='(i3.3)')
            print, 'Staging run ... '+thisRun
            file_copy, template, thisRun, /recursive

            ; Set the grid resolution for this run

            cd, thisRun
            
            this_nR = levels[i] 
            sedCmd = "sed -i .in-bak 's/nRAll(1) = 16/nRAll(1) = "+string(this_nR,format='(i3.3)')+"/' aorsa2d.in"
            spawn, sedCmd
            this_nZ = levels[j] 
            sedCmd = "sed -i .in-bak 's/nZAll(1) = 128/nZAll(1) = "+string(this_nZ,format='(i3.3)')+"/' aorsa2d.in"
            spawn, sedCmd

            spawn, runAR2Cmd

            cd, rootDir

        endif

    endfor

endfor


; Read solutions

stop
end
