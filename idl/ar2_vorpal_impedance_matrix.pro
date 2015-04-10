function ar2_vorpal_impedance_matrix, rSrc, rSmp, ar2RunDir=_ar2RunDir

	if keyword_set(_ar2RunDir) then ar2RunDir = _ar2RunDir else ar2RunDir = './'

	ar2DataFiles = file_search(ar2RunDir+'/output/runData*.nc')
	ar2SolnFiles = file_search(ar2RunDir+'/output/solution*.nc')

	nRHS = n_elements(ar2SolnFiles)
	nPts = n_elements(rSrc)

	; Determine dimensionality for interpolations
	twoD = 0
    a = ar2_read_solution (ar2RunDir, 1)
	if size(a.x2d,/n_dim) gt 1 then twoD = 1

	srcM = complexArr(nPts*3,nRHS)
	smpM = complexArr(nPts*3,nRHS)

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

			srcM[*,rhs] = [$
					complex(srcEr_re,srcEr_im), $
					complex(srcEt_re,srcEt_im), $
					complex(srcEz_re,srcEz_im) ]
	
			; sample plane
			smpEr_re = interpol(real_part(a.E_r),a.r,rsmp,/spline)
			smpEt_re = interpol(real_part(a.E_t),a.r,rsmp,/spline)
			smpEz_re = interpol(real_part(a.E_z),a.r,rsmp,/spline)

			smpEr_im = interpol(imaginary(a.E_r),a.r,rsmp,/spline)
			smpEt_im = interpol(imaginary(a.E_t),a.r,rsmp,/spline)
			smpEz_im = interpol(imaginary(a.E_z),a.r,rsmp,/spline)

			smpM[*,rhs] = [$
					complex(smpEr_re,smpEr_im), $
					complex(smpEt_re,smpEt_im), $
					complex(smpEz_re,smpEz_im) ]

		endif

		if twoD then begin

		endif

   	endfor 

	; Calculate S matrix

	la_svd, smpM, w, u, v, /double
	print, w
	iiSingular = where(w/max(w) lt 1e-5,iiSingularCnt)
	print, iiSingular
	_w = 1/w
	_w[iiSingular] = 0
	_w = diag_matrix(_w) 

	inv = v ## _w ## TRANSPOSE( u )

	S = srcM ## inv

	k = 31.5
	dx = abs(rSmp-rSrc)
	_i = complex(0,1)
	print, exp(-_i*k*dx)

stop

end
