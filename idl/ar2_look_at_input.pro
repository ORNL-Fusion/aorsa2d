pro ar2_look_at_input, FileName=FileName, Compare=Compare

if keyword_set(FileName) then InputFile = FileName else InputFile = 'ar2Input.nc'

cdfId = ncdf_open ( InputFile, /noWrite ) 
	nCdf_varGet, cdfId, 'AtomicZ', _Z
	nCdf_varGet, cdfId, 'amu', amu
	nCdf_varGet, cdfId, 'r', r 
	nCdf_varGet, cdfId, 'z', z 
	nCdf_varGet, cdfId, 'br', br 
	nCdf_varGet, cdfId, 'bt', bt
	nCdf_varGet, cdfId, 'bz', bz 
	nCdf_varGet, cdfId, 'Density_m3', Density_m3
	nCdf_varGet, cdfId, 'Temp_eV', Temp_eV
ncdf_close, cdfId

nR = n_elements(r)
nZ = n_elements(z)
nS = n_elements(amu)

print, 'nR: ', nR
print, 'nZ: ', nZ
print, 'nS: ', nS

if keyword_set(compare) then begin
	cdfId = ncdf_open ( Compare, /noWrite ) 
		nCdf_varGet, cdfId, 'AtomicZ', _Z_c
		nCdf_varGet, cdfId, 'amu', amu_c
		nCdf_varGet, cdfId, 'r', r_c 
		nCdf_varGet, cdfId, 'z', z_c 
		nCdf_varGet, cdfId, 'br', br_c 
		nCdf_varGet, cdfId, 'bt', bt_c
		nCdf_varGet, cdfId, 'bz', bz_c 
		nCdf_varGet, cdfId, 'Density_m3', Density_m3_c
		nCdf_varGet, cdfId, 'Temp_eV', Temp_eV_c
	ncdf_close, cdfId
	
	nR_c = n_elements(r_c)
	nZ_c = n_elements(z_c)
	nS_c = n_elements(amu_c)

	print, 'nR_c: ', nR_c
	print, 'nZ_c: ', nZ_c
	print, 'nS_c: ', nS_c

	i_c = nZ_c/2
endif

i = nZ/2

yRange = [0,1.2]*median(Density_m3[*,i,0])
p=plot(r,Density_m3[*,i,0],title='Density_m3',yRange=yRange)
for s=1,nS-1 do begin
		p=plot(r,Density_m3[*,i,s],/over)
endfor
if keyword_set(compare) then begin
	for s=0,nS_c-1 do begin
			p=plot(r_c,Density_m3_c[*,i_c,s],/over,thick=2,transp=50,color='b')
	endfor
endif

p=plot(r,Temp_eV[*,i,0],title='Temp_eV')
for s=1,nS-1 do begin
		p=plot(r,Temp_eV[*,i,s],/over)
endfor
if keyword_set(compare) then begin
	for s=0,nS_c-1 do begin
			p=plot(r_c,Temp_eV_c[*,i_c,s],/over,thick=2,transp=50,color='b')
	endfor
endif

p=plot(r,bt[*,i],title='b')
p=plot(r,br[*,i],/over)
p=plot(r,bz[*,i],/over)

if keyword_set(compare) then begin
	p=plot(r_c,bt_c[*,i_c],/over,thick=2,transp=50,color='b')
	p=plot(r_c,br_c[*,i_c],/over,thick=2,transp=50,color='b')
	p=plot(r_c,bz_c[*,i_c],/over,thick=2,transp=50,color='b')
endif


stop

end
