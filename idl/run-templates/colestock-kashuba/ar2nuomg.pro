global_nuOmg = r2D*0

; Create a smooth absorbing layer at the boundary edge

MinNuOmg = 0.025

for s=0,nSpec-1 do begin
    nuOmg[*,*,s] = global_nuOmg>MinNuOmg
endfor


