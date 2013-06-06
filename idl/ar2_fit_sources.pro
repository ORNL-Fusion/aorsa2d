pro ar2_fit_sources, $
	ThisComponentID=ThisComponentID, $
    RunDataFiles = RunDataFiles, $
	CoeffsOut2D = CoeffsOut2D, $ ; nPhi x nSources
    TwoDim = TwoDim

	; 0 = r
	; 1 = t
	; 2 = z

	if not keyword_set(ThisComponentID) then ThisComponentID=0

	SourceLocationsFile = 'ar2SourceLocations.nc'
	if not keyword_set(RunDataFiles) then ThisRunDataFiles = file_search('output/runData*.nc')

	cdfId = ncdf_open ( SourceLocationsFile, /noWrite ) 

		nCdf_varGet, cdfId, 'cs_r', CurrentSource_r
		nCdf_varGet, cdfId, 'cs_z', CurrentSource_z
		nCdf_varGet, cdfId, 'component_ident', CurrentSource_ComponentID

	ncdf_close, cdfId

	NRHS = n_elements(CurrentSource_r)

    NRHS_check = n_elements(ThisRunDataFiles)
	if NRHS ne NRHS_check then begin
			print, 'Error: NRHS does not match'
			stop
	endif


    for rhs=0,NRHS-1 do begin

	    cdfId = ncdf_open ( ThisRunDataFiles[rhs], /noWrite ) 

	    	nCdf_varGet, cdfId, 'capR', r
	    	nCdf_varGet, cdfId, 'y', z
            ncdf_VarGet, cdfId, 'nPhi', nPhi
	    	nCdf_varGet, cdfId, 'jr_re', jr_re
	    	nCdf_varGet, cdfId, 'jr_im', jr_im
	    	nCdf_varGet, cdfId, 'jt_re', jt_re
	    	nCdf_varGet, cdfId, 'jt_im', jt_im
	    	nCdf_varGet, cdfId, 'jz_re', jz_re
	    	nCdf_varGet, cdfId, 'jz_im', jz_im

	    ncdf_close, cdfId

        nR = n_elements(r)
        nZ = n_elements(z)

        if(jr_2D eq !NULL) then begin
	        jr_2D = ComplexArr(nR,nZ,NRHS)
	        jt_2D = ComplexArr(nR,nZ,NRHS)
	        jz_2D = ComplexArr(nR,nZ,NRHS)
            RHS_nPhi = FltArr(NRHS)
        endif

        jr_2D[*,*,rhs] = complex(jr_re,jr_im)
        jt_2D[*,*,rhs] = complex(jt_re,jt_im)
        jz_2D[*,*,rhs] = complex(jz_re,jz_im)
        RHS_nPhi[rhs] = nPhi

    endfor

	nSources = n_elements(CurrentSource_r)/3

; ---------------------------------------
; Create a test Vorpal dataset to fit to.
;

	nr_data = 64d0
	nt_data = 1d0
	nz_data = 32d0
	
	n_data = nr_data * nt_data * nz_data

	r_data_1D = fltArr(nr_data)
	t_data_1D = fltArr(nt_data)
	z_data_1D = fltArr(nz_data)

	r_data = fltArr(n_data)
	t_data = fltArr(n_data)
	z_data = fltArr(n_data)

	jr_data = complexarr(n_data)
	jt_data = complexarr(n_data)
	jz_data = complexarr(n_data)

    rmin = 1.7d0
    rmax = 2.1d0

	tmin = 0d0 ; degrees
	tmax = 0d0

	zmin = -0.2
	zmax = +0.2

	if nr_data eq 1 then r_data_1D[*] = rmin else $
	    r_data_1D = rmin+fIndGen(nr_data)*(rmax-rmin)/(nr_data-1d0)		
	if nt_data eq 1 then t_data_1D[*] = tmin else $
	    t_data_1D = tmin+fIndGen(nt_data)*(tmax-tmin)/(nt_data-1d0)		
	if nz_data eq 1 then z_data_1D[*] = zmin else $
	    z_data_1D = zmin+fIndGen(nz_data)*(zmax-zmin)/(nz_data-1d0)		

    r_data = transpose(rebin(r_data_1D,nr_data,nt_data,nz_data),[0,1,2])
    t_data = transpose(rebin(t_data_1D,nt_data,nz_data,nr_data),[1,2,0])
    z_data = transpose(rebin(z_data_1D,nz_data,nt_data,nr_data),[2,1,0])

	;c = 0

	;for i=0,nr_data-1 do begin
	;	for j=0,nt_data-1 do begin
	;		for k=0,nz_data-1 do begin

	;			r_data[c] = r_data_1D[i]
	;			t_data[c] = t_data_1D[j]	
	;			z_data[c] = z_data_1D[k]

	;			c++
	;		endfor
	;	endfor
	;endfor

	t_data = t_data*!dtor
	t_data_1D = t_data_1D*!dtor

	;kz = 4.0
	;kt = 5.5
	;jr_data = exp(complex(0,1)*(kz*z_data+kt*r_data*t_data))

    antR = 2.0
    antZ = 0.0
    sigR = 0.03
    sigZ = 0.1
    jr_data = 1*exp ( -( (r_data-antR)^2/sigR^2 $
                            + (z_data-antZ)^2/sigZ^2 ) )

    jt_data = jr_data*0
    jz_data = jt_data*0


;
; End of creating test provided jP data.
; --------------------------------------


    ; Build fit martix A

    if keyword_set(TwoDim) then begin

        n_nPhi = 1
        nPhi = RHS_nPhi[0]

    endif else begin

	    n_nPhi = 21 
	    nPhi_stretch = 4 
	    nPhi = IndGen(n_nPhi)*nPhi_stretch-(n_nPhi/2)*nPhi_stretch

    endelse

	n_Basis = NRHS * n_nPhi / 3

	amat = complexarr(n_data,n_Basis)

	c = 0
	for p=0,n_nPhi-1 do begin
		for rhs=0,NRHS-1 do begin
			if CurrentSource_ComponentID[rhs] eq ThisComponentID then begin

				; interpolate 2D
				r_index = (r_data-min(r))/(max(r)-min(r))*nR
				z_index = (z_data-min(z))/(max(z)-min(z))*nZ
				tmp_2D = interpolate(jr_2D[*,*,rhs],[r_index],[z_index],cubic=-0.5)

				for i=0,n_data-1 do begin

					; expand to 3D for t angle and nPhi
                    if keyword_set(TwoDim) then begin
					    amat[i,c] = tmp_2D[i]
                    endif else begin
					    tmp_3D  = tmp_2D[i]*exp(complex(0,1)*nPhi[p]*t_data[i])
					    amat[i,c] = tmp_3D
                    endelse

				endfor
				c++
			endif
		endfor
	endfor

	
	;coeff_reform_test = intarr(n_nPhi*NRHS/3)
	;c=0
	;for p=0,n_nPhi-1 do begin
	;	for rhs=0,NRHS-1,3 do begin
	;			coeff_reform_test[c] = rhs
	;			c++
	;	endfor
	;endfor

	;amat = [[amat],[complex(0,1)*amat]]

	print, n_Basis, c

	help, amat
	b = jr_data[*]
	coeffs = LA_LEAST_SQUARES(transpose(amat),b, status=stat,method=3)

    if keyword_set(TwoDim) then begin

        x = r_data_1D
        y = z_data_1D

        nLevs = 20
        ScaleFac = 0.3
        scale = max(abs([jr_data]))*ScaleFac
        dimensions = [600,600]

        ThisField = reform(jr_data,nR_data,nZ_data)

		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)

        title = 'Re(Data)'
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
            aspect_ratio=1.0, layout=[2,2,1], dimensions=dimensions, title=title )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
        p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )

        title = 'Im(Data)'
		PlotField = (imaginary(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
            aspect_ratio=1.0, layout=[2,2,2], /current, title=title )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
        p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )

        ThisField = reform(transpose(amat)##coeffs,nR_data,nZ_data)

        title = 'Re(Fit)'
 		PlotField = (real_part(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
            aspect_ratio=1.0, layout=[2,2,3], /current, title=title )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
        p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )

        title = 'Im(Fit)'
		PlotField = (imaginary(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
            aspect_ratio=1.0, layout=[2,2,4], /current, title=title )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
        p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )

        Coeffs2D = coeffs

        stop

    endif else begin

	    input = transpose(reform(jr_data,nz_data,nt_data))
	    c1=contour(input,t_data_1D*!radeg,z_data_1D,n_levels=21,/fill)

	    fit = transpose(reform(transpose(amat)##coeffs,nz_data,nt_data))
	    c2=contour(fit,t_data_1D*!radeg,z_data_1D,n_levels=21,/fill)

	    zSlice=3
	    p=plot(t_data_1d*!radeg,input[*,zSlice])
	    p=plot(t_data_1d*!radeg,fit[*,zSlice],color='red',/over)
	    p=plot(t_data_1d*!radeg,imaginary(input[*,zSlice]),/over,linestyle='dash')
	    p=plot(t_data_1d*!radeg,imaginary(fit[*,zSlice]),color='red',/over,linestyle='dash')

	    tSlice=3
	    p=plot(z_data_1d,input[tSlice,*])
	    p=plot(z_data_1d,fit[tSlice,*],color='red',/over)
	    p=plot(z_data_1d,imaginary(input[tSlice,*]),/over,linestyle='dash')
	    p=plot(z_data_1d,imaginary(fit[tSlice,*]),color='red',/over,linestyle='dash')

	    c3=contour(imaginary(input),t_data_1D*!radeg,z_data_1D,n_levels=21,/fill,rgb_table=3)
	    c4=contour(imaginary(fit),t_data_1D*!radeg,z_data_1D,n_levels=21,/fill,rgb_table=3)

	    coeffs_2D = transpose(reform(coeffs,nSources,n_nPhi))

	    c5=contour(coeffs_2D,nPhi,CurrentSource_z[IndGen(nSources)*3],n_levels=21,/fill)

    endelse

	CoeffsOut2D = Coeffs2D ; (n_nPhi,nSources)

end
