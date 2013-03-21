pro ar2_fit_sources

	SourceLocationsFile = 'AR2SourceLocations.nc'
	RunDataFile = 'runData001.nc'

	cdfId = ncdf_open ( SourceLocationsFile, /noWrite ) 

		nCdf_varGet, cdfId, 'cs_r', CurrentSource_r
		nCdf_varGet, cdfId, 'cs_z', CurrentSource_z
		nCdf_varGet, cdfId, 'component_ident', CurrentSource_ComponentID

	ncdf_close, cdfId

	NRHS = n_elements(CurrentSource_r)

	cdfId = ncdf_open ( RunDataFile, /noWrite ) 

		nCdf_varGet, cdfId, 'capR', r
		nCdf_varGet, cdfId, 'y', z
		nCdf_varGet, cdfId, 'jr_re', jr_re
		nCdf_varGet, cdfId, 'jr_im', jr_im
		nCdf_varGet, cdfId, 'jt_re', jt_re
		nCdf_varGet, cdfId, 'jt_im', jt_im
		nCdf_varGet, cdfId, 'jz_re', jz_re
		nCdf_varGet, cdfId, 'jz_im', jz_im

	ncdf_close, cdfId

	NRHS_check = n_elements(jr_re[0,0,*])

	if NRHS ne NRHS_check then begin
			print, 'Error: NRHS does not match'
			stop
	endif

	jr_2D = complex(jr_re,jr_im)
	jt_2D = complex(jt_re,jt_im)
	jz_2D = complex(jz_re,jz_im)

	; Fit to data at these locations

	nr_Fit = 1
	nt_Fit = 50
	nz_Fit = 20
	
	n_Fit = nr_Fit * nt_Fit * nz_Fit

	r_Fit_1D = fltArr(nr_Fit)
	t_Fit_1D = fltArr(nt_Fit)
	z_Fit_1D = fltArr(nz_Fit)

	r_Fit = fltArr(n_Fit)
	t_Fit = fltArr(n_Fit)
	z_Fit = fltArr(n_Fit)

	jr_vorpal = complexarr(n_Fit)
	jt_vorpal = complexarr(n_Fit)
	jz_vorpal = complexarr(n_Fit)

	rmin = 7.0

	tmin = -20
	tmax = +20

	zmin = -1.0
	zmax = +1.0

	r_Fit_1D[*] = rmin
	t_Fit_1D = tmin+fIndGen(nt_Fit)*(tmax-tmin)/(nt_Fit-1)		
	z_Fit_1D = zmin+fIndGen(nz_Fit)*(zmax-zmin)/(nz_Fit-1)		

	c = 0

	for i=0,nr_Fit-1 do begin
		for j=0,nt_Fit-1 do begin
			for k=0,nz_Fit-1 do begin

				r_Fit[c] = rmin
				t_Fit[c] = tmin+j*(tmax-tmin)/(nt_Fit-1)		
				z_Fit[c] = zmin+k*(zmax-zmin)/(nz_Fit-1)		

				kz = 3.0
				kt = 1.0
				jr_vorpal[c] = exp(complex(0,1)*$
						(kz*z_Fit[c]+kt*r_Fit[c]*t_Fit[c]))

				c++
			endfor
		endfor
	endfor

	t_Fit = t_Fit*!dtor
	t_Fit_1D = t_Fit_1D*!dtor

	nr = n_elements(jr_2D[*,0,0])
	nz = n_elements(jr_2D[0,*,0])

	n_nPhi = 20 
	nPhi_stretch = 10
	nPhi = IndGen(n_nPhi)*nPhi_stretch-n_nPhi*nPhi_stretch/2

	n_Basis = NRHS * n_nPhi / 6

	amat = complexarr(n_Fit,n_Basis)

	c = 0
	for p=0,n_nPhi-1 do begin
		for rhs=0,NRHS-1 do begin
			if CurrentSource_ComponentID[rhs] eq 0 then begin
				for i=0,n_Fit-1 do begin
					; interpolate 2D
					r_index = (r_Fit[i]-min(r))/(max(r)-min(r))*nr
					z_index = (z_Fit[i]-min(z))/(max(z)-min(z))*nz
					tmp_2D = interpolate(jr_2D[*,*,rhs],[r_index],[z_index])
					; expand to 3D for t angle and nPhi
					tmp_3D  = tmp_2D*exp(complex(0,1)*nPhi[p]*t_Fit[i])
					amat[i,c] = tmp_3D
				endfor
				;c3=contour(transpose(reform(amat[*,c],nz_Fit,nt_Fit)),$
				;	t_Fit_1D*!radeg,z_Fit_1D,n_levels=21,/fill,$
				;	xtitle='theta [degrees]', ytitle='z [m]', title='a basis fn')

				;stop
				c++
			endif
		endfor
	endfor

	print, n_Basis, c

	coeffs = LA_LEAST_SQUARES(transpose(amat),jr_vorpal, status=stat)

	input = transpose(reform(jr_vorpal,nz_Fit,nt_Fit))
	c1=contour(input,t_Fit_1D*!radeg,z_Fit_1D,n_levels=21,/fill)

	fit = transpose(reform(transpose(amat)##coeffs,nz_Fit,nt_Fit))
	c2=contour(fit,t_Fit_1D*!radeg,z_Fit_1D,n_levels=21,/fill)

	slice=3
	p=plot(t_fit_1d*!radeg,input[*,slice])
	p=plot(t_fit_1d*!radeg,fit[*,slice],color='red',/over)
	p=plot(t_fit_1d*!radeg,imaginary(input[*,slice]),/over,linestyle='dash')
	p=plot(t_fit_1d*!radeg,imaginary(fit[*,slice]),color='red',/over,linestyle='dash')

	c3=contour(imaginary(input),t_Fit_1D*!radeg,z_Fit_1D,n_levels=21,/fill,rgb_table=3)
	c4=contour(imaginary(fit),t_Fit_1D*!radeg,z_Fit_1D,n_levels=21,/fill,rgb_table=3)


stop

end
