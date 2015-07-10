global_nuOmg = r2d*0

bdry_nuOmg = 1.0
PercentageWidth = 2
nIterations = 10
dR = 0.05
dZ = 0.05

width = fix(nR*PercentageWidth/100)
for n=0,nIterations-1 do begin

    ii = where(r2D gt rMax-dR)
    global_nuOmg[ii] = 1
    ii = where(r2D lt rMin+dR)
    global_nuOmg[ii] = 1
    
    ii = where(z2D gt zMax-dz)
    global_nuOmg[ii] = 1
    ii = where(z2D lt zMin+dz)
    global_nuOmg[ii] = 1

    global_nuOmg = smooth(global_nuOmg,width,/edge_truncate) * bdry_nuOmg

endfor

for s=0,nSpec-1 do begin
    nuOmg[*,*,s] = global_nuOmg>0.01
endfor


