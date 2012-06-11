function prof1, x, a
	return, a[0] + ( a[1] - a[0] ) * ( 1.0 - x^a[3] )^a[2]
end

pro ar2_create_flux_profiles, nSpec, nn, tt, nR, nZ, PsiNorm, Mask_bbb, d_bbb, Density_m3, Temp_eV, $
		DensityMin = DensityMin, TempMin = TempMin

	@constants

	; exp decay profile from Lamelle et al., Nucl. Fusion, 46 (2006) 432-443
	; see pages 436-437
	; Update: had to modify this since Lamelles goes way too low.

	l_sol = 0.03
	lambda = 1.0/l_sol

	for s=0,nspec-1 do begin

		for i=0,nR-1 do begin
			for j=0,nZ-1 do begin

				if Mask_bbb[i,j] eq 1 then begin

					density_m3[i,j,s] = prof1(sqrt(PsiNorm[i,j]),nn[*,s])
					temp_eV[i,j,s] = prof1(sqrt(PsiNorm[i,j]),tt[*,s])

				endif else begin

					eta = d_bbb[i,j]*lambda
					density_m3[i,j,s] = nn[0,s] * exp (-eta)
					temp_eV[i,j,s] = tt[0,s] * exp (-eta)

				endelse

			endfor
		endfor

	endfor

	Density_m3 = Density_m3>DensityMin
	Temp_eV = Temp_eV>TempMin

   	s_ne = surface(density_m3[*,*,0], layout=[3,1,1])
   	s_te = surface(temp_eV[*,*,0], layout=[3,1,2], /current)

	p=plot(Density_m3[*,nZ/2,0])
	for s=1,nSpec-1 do begin
		p=plot(Density_m3[*,nZ/2,s],/over)
	endfor

	p=plot(Temp_eV[*,nZ/2,0])
	for s=1,nSpec-1 do begin
		p=plot(Temp_eV[*,nZ/2,s],/over)
	endfor


	; Sanity checking ...

	iiBad = where(density_m3 lt 1e5,iiBadCnt)
	if iiBadCnt gt 0 then begin
		print, 'ERROR: density failed sanity check'
		stop
	endif

	iiBad = where(temp_eV lt 0,iiBadCnt)
	if iiBadCnt gt 0 then begin
		print, 'ERROR: temp failed sanity check'
		stop
	endif

	iiBad = where(temp_eV*e lt 0,iiBadCnt)
	if iiBadCnt gt 0 then begin
		print, 'ERROR: temp*q failed sanity check'
		stop
	endif

end
