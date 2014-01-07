pro ar2_fit_sources, $
    FitThis = FitThis, $
    WithTheseFiles = WithTheseFiles, $
	CoeffsOut = CoeffsOut, $ ; [nSourceFiles]
    TwoDim = TwoDim, $
    NRHS = NRHS, $
    FitLayer = FitLayer

	; 0 = r
	; 1 = t
	; 2 = z

	DataFiles = file_search(WithTheseFiles+'/output/solution*.nc')
    sFitMe = ar2_read_solution (WithTheseFiles, 0)
    nSpecies = n_elements(sFitMe.jP_r[0,0,*])

	if(not keyword_set(NRHS))then NRHS = n_elements(DataFiles)

    ; Determine which points lay in the fit layer

    iiFitThese = where(sFitMe.r ge FitLayer[0] and sFitMe.r le FitLayer[1], iiFitTheseCnt)
    if iiFitTheseCnt le 0 then begin
        print, 'Error: no points to fit within FitLayer'
        stop
    endif

    N = iiFitTheseCnt

    ; Build fit martix A and solve for each component seperately.
    ; Could do each species seperately too I guess.

	amat = complexarr(N,NRHS)

    ComponentArray = [0,1,2]

for spec = 0, nSpecies-1 do begin

    foreach component, ComponentArray do begin

        c = 0

	    for rhs=0,NRHS-1 do begin

            s = ar2_read_solution (WithTheseFiles, rhs)

            if component eq 0 then thisJp = s.jP_r[iiFitThese,0,spec]
            if component eq 1 then thisJp = s.jP_t[iiFitThese,0,spec]
            if component eq 2 then thisJp = s.jP_z[iiFitThese,0,spec]

            aMat[*,rhs] = thisJp[*]

	    endfor

	    help, amat

        if component eq 0 then b = sFitMe.jP_r[iiFitThese,0,spec]
        if component eq 1 then b = sFitMe.jP_t[iiFitThese,0,spec]
        if component eq 2 then b = sFitMe.jP_z[iiFitThese,0,spec]

	    coeffs = LA_LEAST_SQUARES(transpose(amat),b, status=stat,method=3)
        
        if component eq 0 then Coeffs_R = coeffs
        if component eq 1 then Coeffs_T = coeffs
        if component eq 2 then Coeffs_Z = coeffs

    endforeach

endfor

    ; Create the "perFileList" coefficient list

    CoeffsOut = ComplexArr(NRHS)

    CoeffsOut[ComponentIndicies[*,0]] = Coeffs_R 
    CoeffsOut[ComponentIndicies[*,1]] = Coeffs_T 
    CoeffsOut[ComponentIndicies[*,2]] = Coeffs_Z 

    p=plot(Coeffs_R,layout=[1,4,1],dimension=[600,600],title='Coeffs R')
    p=plot(imaginary(Coeffs_R),/over,color='b')
    p=plot(Coeffs_T,layout=[1,4,2],/current,title='Coeffs T')
    p=plot(imaginary(Coeffs_T),/over,color='b')
    p=plot(Coeffs_Z,layout=[1,4,3],/current,title='Coeffs Z')
    p=plot(imaginary(Coeffs_Z),/over,color='b')
    p=plot(CoeffsOut,layout=[1,4,4],/current,title='Coeffs All')
    p=plot(imaginary(CoeffsOut),/over,color='b')


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

        ThisField = FitR_2D 

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

end
