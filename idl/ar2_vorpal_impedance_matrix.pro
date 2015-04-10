function ar2_vorpal_impedance_matrix, rSrc, rSmp, ar2RunDir=_ar2RunDir

	if keyword_set(_ar2RunDir) then ar2RunDir = _ar2RunDir else ar2RunDir = './'

	ar2DataFiles = file_search(ar2RunDir+'/output/runData*.nc')
	ar2SolnFiles = file_search(ar2RunDir+'/output/solution*.nc')

	nRHS = n_elements(ar2SolnFiles)

	; Determine dimensionality for interpolations
	twoD = 0
    a = ar2_read_solution (ar2RunDir, 1)
	if size(a.x2d,/n_dim) gt 1 then twoD = 1

	for rhs=0,NRHS-1 do begin

		print, 'RHS: ', rhs+1
        a = ar2_read_solution (ar2RunDir, rhs+1)

		; Interpolate the aorsa E field values to the xSrc and xSmp locations

		if not twoD then begin

			; source plane
			srcEr_re = interpol(real_part(a.E_r),a.r,rSrc,/spline)
			srcEt_re = interpol(real_part(a.E_t),a.r,rSrc,/spline)
			srcEz_re = interpol(real_part(a.E_z),a.r,rSrc,/spline)

			srcEr_im = interpol(imaginary(a.E_r),a.r,rSrc,/spline)
			srcEt_im = interpol(imaginary(a.E_t),a.r,rSrc,/spline)
			srcEz_im = interpol(imaginary(a.E_z),a.r,rSrc,/spline)

			srcEr = complex(srcEr_re,srcEr_im)
			srcEt = complex(srcEt_re,srcEt_im)
			srcEz = complex(srcEz_re,srcEz_im)
	
			; sample plane
			smpEr_re = interpol(real_part(a.E_r),a.r,rsmp,/spline)
			smpEt_re = interpol(real_part(a.E_t),a.r,rsmp,/spline)
			smpEz_re = interpol(real_part(a.E_z),a.r,rsmp,/spline)

			smpEr_im = interpol(imaginary(a.E_r),a.r,rsmp,/spline)
			smpEt_im = interpol(imaginary(a.E_t),a.r,rsmp,/spline)
			smpEz_im = interpol(imaginary(a.E_z),a.r,rsmp,/spline)

			smpEr = complex(smpEr_re,smpEr_im)
			smpEt = complex(smpEt_re,smpEt_im)
			smpEz = complex(smpEz_re,smpEz_im)

			ratio_r = srcEr/smpEr
			ratio_t = srcEt/smpEt
			ratio_z = srcEz/smpEz

			print, ratio_r, ratio_t, ratio_z

		endif

		if twoD then begin

		endif

		

   	endfor 

end
