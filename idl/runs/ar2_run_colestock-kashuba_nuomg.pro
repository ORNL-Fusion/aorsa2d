global_nuOmg = r2D*0

; Create a smooth absorbing layer at the boundary edge

MinNuOmg = 0.0025
absorbingNuOmg = 1.0
absorbing_left_r = 1.60 ; Left domain edge
absorbing_right_r = 1.72
cos_arg_2D = -(r2D - absorbing_right_r)/(absorbing_right_r-absorbing_left_r)*!pi;-!pi
iiNuOmgSet1 = where(r2D gt absorbing_left_r,iiNuOmgSetCnt) 
global_nuOmg[iiNuOmgSet1] = (cos(cos_arg_2D[iiNuOmgSet1])+1)*0.5*absorbingNuOmg+MinNuOmg
iiNuOmgSet2 = where(r2D gt absorbing_right_r,iiNuOmgSetCnt) 
global_nuOmg[iiNuOmgSet2] = absorbingNuOmg+MinNuOmg

absorbing_left_r = 0.92 ; Left domain edge
absorbing_right_r = 1.04
cos_arg_2D = -(r2D - absorbing_right_r)/(absorbing_right_r-absorbing_left_r)*!pi-!pi
iiNuOmgSet1 = where(r2D lt absorbing_right_r,iiNuOmgSetCnt) 
global_nuOmg[iiNuOmgSet1] = (cos(cos_arg_2D[iiNuOmgSet1])+1)*0.5*absorbingNuOmg+MinNuOmg
iiNuOmgSet2 = where(r2D lt absorbing_left_r,iiNuOmgSetCnt) 
global_nuOmg[iiNuOmgSet2] = absorbingNuOmg+MinNuOmg

for s=0,nSpec-1 do begin
    nuOmg[*,*,s] = global_nuOmg>MinNuOmg
endfor


