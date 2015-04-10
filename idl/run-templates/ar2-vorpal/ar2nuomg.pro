global_nuOmg = r2D*0

; Create a smooth absorbing layer at the boundary edge

MinNuOmg = 0.00

absorbingNuOmg = 10.0

left_absorber = 1
right_absorber = 0

_rMin = 1.1 
_rMax = 1.7 

if left_absorber then begin
	cos_arg_2D = -(r2D - _rMax)/(_rMax-_rMin)*!pi-!pi
	iiNuOmgSet1 = where(r2D lt _rMax,iiNuOmgSetCnt) 
	global_nuOmg[iiNuOmgSet1] = (cos(cos_arg_2D[iiNuOmgSet1])+1)*0.5*absorbingNuOmg
	iiNuOmgSet2 = where(r2D lt _rMin,iiNuOmgSetCnt) 
	global_nuOmg[iiNuOmgSet2] = absorbingNuOmg
endif

if right_absorber then begin
	cos_arg_2D = -(r2D - _rMax)/(_rMax-_rMin)*!pi
	iiNuOmgSet1 = where(r2D gt _rMin,iiNuOmgSetCnt) 
	global_nuOmg[iiNuOmgSet1] = (cos(cos_arg_2D[iiNuOmgSet1])+1)*0.5*absorbingNuOmg
	iiNuOmgSet2 = where(r2D gt _rMax,iiNuOmgSetCnt) 
	global_nuOmg[iiNuOmgSet2] = absorbingNuOmg
endif


for s=0,nSpec-1 do begin
    nuOmg[*,*,s] = global_nuOmg>MinNuOmg
endfor

