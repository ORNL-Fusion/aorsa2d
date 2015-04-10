function prof1, x, a
	return, a[0] + ( a[1] - a[0] ) * ( 1.0 - x^a[3] )^a[2]
end

pro ar2_create_flux_profiles, nSpec, nn, tt, nR, nZ, PsiNorm, Mask_bbb, d_bbb, Density_m3, Temp_eV, $
		DensityMin = DensityMin, TempMin = TempMin, $
        NumericProfiles = NumericProfiles, NumericData_n_m3 = NumericData_n_m3, $
        NumericData_T_eV = NumericData_T_eV, $ 
		r2d=r2d, z2d=z2d

	@constants

    if keyword_set(NumericProfiles) then numeric = NumericProfiles else numeric = 0

	; exp decay profile from Lamelle et al., Nucl. Fusion, 46 (2006) 432-443
	; see pages 436-437
	; Update: had to modify this since Lamelles goes way too low.

	for s=0,nspec-1 do begin

		for i=0,nR-1 do begin
			for j=0,nZ-1 do begin

				if Mask_bbb[i,j] eq 1 then begin

                    if numeric then begin

                        rho = sqrt(PsiNorm[i,j])

                        thisRho = numericData_n_m3[*,0]
                        thisn = numericData_n_m3[*,s+1]
                        density_m3[i,j,s] = interpol(thisn,thisRho,rho, /spline) 

                        thisRho = numericData_T_eV[*,0]
                        thisT = numericData_T_eV[*,s+1]
                        temp_eV[i,j,s] = interpol(thisT,thisRho,rho, /spline) 

                    endif else begin
					    density_m3[i,j,s] = prof1(sqrt(PsiNorm[i,j]),nn[*,s])
					    temp_eV[i,j,s] = prof1(sqrt(PsiNorm[i,j]),tt[*,s])
                    endelse

				endif

			endfor
		endfor

		for i=0,nR-1 do begin
			for j=0,nZ-1 do begin

				if Mask_bbb[i,j] eq 0 then begin

                	; find closest point inside the mask
					_r = r2d[i,j]
					_z = z2d[i,j]
					_d = sqrt((r2d-_r)^2+(z2d-_z)^2)
					_d = _d + ((-(mask_bbb-1))*1e6) ; just add a large number outside of lcfs
					_ii = where(abs(_d) eq min(abs(_d)), iiCnt) 	
					_ii = _ii[0]
					_nTmp = density_m3[*,*,s]
					_tTmp = temp_eV[*,*,s]
					; find the n and t values at that closest point
					_n = _nTmp[_ii]
					_t = _tTmp[_ii]
					;_mask = 1+(-(mask_bbb-1))*1e30
					;_n = min(_nTmp*_mask)
					;_t = min(_tTmp*_mask)

	                l_sol = 0.02
	                n_ant = densityMin[s]
	                t_ant = tempMin

					;eta = -d_bbb[i,j]/l_sol
                    eta = -(psiNorm[i,j]-1)/l_sol
					density_m3[i,j,s] = n_ant + (_n - n_ant) * exp (eta)
					temp_eV[i,j,s] = t_ant + (_t - t_ant) * exp (eta)
					;density_m3[i,j,s] = _n * exp (eta) > _n/maxDecay
					;temp_eV[i,j,s] = _t * exp (eta) > _t/maxDecay

				endif

			endfor
		endfor

	endfor

    ;; Run simple LCFS to wall smoothing algorithm

    ;nIter = 45 
    ;PercentWidth =  2.0
    ;nWidth = fix(nR*PercentWidth/100.0)

    ;density_m3_ = density_m3
    ;temp_eV_ = temp_eV
    ;iiMask = where(mask_bbb eq 1) 

    ;for s=0,nSpec-1 do begin

    ;    original = density_m3_[*,*,s]
    ;    for nIter = 0,nIter-1 do begin
    ;           density_m3[*,*,s] = density_m3[*,*,s]>DensityMin
    ;           temp = density_m3[*,*,s]
    ;           temp[iiMask] = original[iiMask]
    ;           density_m3[*,*,s] = temp
    ;           density_m3[*,*,s] = smooth(density_m3[*,*,s],nWidth,/edge_truncate)
    ;    endfor

    ;    original = temp_eV_[*,*,s]
    ;    for nIter = 0,nIter-1 do begin
    ;           temp_eV[*,*,s] = temp_eV[*,*,s]>TempMin
    ;           temp = temp_eV[*,*,s]
    ;           temp[iiMask] = original[iiMask]
    ;           temp_eV[*,*,s] = temp
    ;           temp_eV[*,*,s] = smooth(temp_eV[*,*,s],nWidth,/edge_truncate)
    ;    endfor
 
    ;endfor

	;Density_m3 = Density_m3>DensityMin
	;Temp_eV = Temp_eV>TempMin

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
