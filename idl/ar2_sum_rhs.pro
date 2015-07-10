pro ar2_sum_rhs


    ar2InputFile = 'ar2Input.nc'
	SourceLocationsFile = 'ar2SourceLocations.nc'
	RunDataFiles = file_search('output/runData*.nc')
    SolutionFiles = file_search('output/solution*.nc')

 	ar2_read_ar2input, ar2InputFile, RunDataFiles[0], $
			rLim=rLim,zLim=zLim,LimMask=LimMask
   
    NRHS_check = n_elements(RunDataFiles)
    assert, NRHS_check > 0, 'Oops, No runFiles?'

	SourceLocationsFile = 'ar2SourceLocations.nc'

	cdfId = ncdf_open ( SourceLocationsFile, /noWrite ) 

		nCdf_varGet, cdfId, 'cs_r', CurrentSource_r
		nCdf_varGet, cdfId, 'cs_z', CurrentSource_z
		nCdf_varGet, cdfId, 'component_ident', CurrentSource_ComponentID

	ncdf_close, cdfId

    ar2_fit_sources, CoeffsOut = Coeffs, RunDataFiles = RunDataFiles, /TwoDim

    NRHS = n_elements(Coeffs)

    assert, NRHS eq NRHS_check, 'NRHS ne nRunDataFiles'

    ; Load the fields for all RHSs

    RHS_nPhi = IntArr(NRHS)

    for rhs=0,NRHS-1 do begin

	    cdfId = ncdf_open ( SolutionFiles[rhs], /noWrite ) 

	    	nCdf_varGet, cdfId, 'r', r
	    	nCdf_varGet, cdfId, 'z', z
            ncdf_VarGet, cdfId, 'nPhi', nPhi

	    	nCdf_varGet, cdfId, 'er_re', er_re
	    	nCdf_varGet, cdfId, 'er_im', er_im
	    	nCdf_varGet, cdfId, 'et_re', et_re
	    	nCdf_varGet, cdfId, 'et_im', et_im
	    	nCdf_varGet, cdfId, 'ez_re', ez_re
	    	nCdf_varGet, cdfId, 'ez_im', ez_im

	    	nCdf_varGet, cdfId, 'jP_r_re', jPr_re
	    	nCdf_varGet, cdfId, 'jP_r_im', jPr_im
	    	nCdf_varGet, cdfId, 'jP_t_re', jPt_re
	    	nCdf_varGet, cdfId, 'jP_t_im', jPt_im
	    	nCdf_varGet, cdfId, 'jP_z_re', jPz_re
	    	nCdf_varGet, cdfId, 'jP_z_im', jPz_im


	    ncdf_close, cdfId

        nR = n_elements(r)
        nZ = n_elements(z)

        RHS_nPhi[rhs] = nPhi

        if er_2D eq !null then begin
                
            er_2D = coeffs[rhs] * complex(er_re,er_im)
            et_2D = coeffs[rhs] * complex(et_re,et_im)
            ez_2D = coeffs[rhs] * complex(ez_re,ez_im)

            jPr_2D = coeffs[rhs] * complex(jPr_re,jPr_im)
            jPt_2D = coeffs[rhs] * complex(jPt_re,jPt_im)
            jPz_2D = coeffs[rhs] * complex(jPz_re,jPz_im)

        endif else begin

            er_2D = er_2D + coeffs[rhs] * complex(er_re,er_im)
            et_2D = et_2D + coeffs[rhs] * complex(et_re,et_im)
            ez_2D = ez_2D + coeffs[rhs] * complex(ez_re,ez_im)

            print, max(abs(er_2D)), rhs

            jPr_2D = jPr_2D + coeffs[rhs] * complex(jPr_re,jPr_im)
            jPt_2D = jPt_2D + coeffs[rhs] * complex(jPt_re,jPt_im)
            jPz_2D = jPz_2D + coeffs[rhs] * complex(jPz_re,jPz_im)

        endelse

        RHS_nPhi[rhs] = nPhi

    endfor

    x = r
    y = z

    nLevs = 20
    ScaleFac = 0.3
    scale = max(abs([er_2D,et_2D,ez_2D]))*ScaleFac
    dimensions = [600,600]


	levels = fIndGen(nLevs)/(nLevs-1)*scale
	colors = reverse(bytScl(levels, top=253)+1)

    ThisField = er_2D
    title = 'Re(er_2D)'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
	c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
        aspect_ratio=1.0, layout=[2,3,1], dimensions=dimensions, title=title )
	c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
    p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )
    p = plot ( rLim, zLim, /over )

    ThisField = et_2D
    title = 'Re(et_2D)'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
	c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
        aspect_ratio=1.0, layout=[2,3,3], /current, title=title )
	c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
    p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )
    p = plot ( rLim, zLim, /over )

    ThisField = ez_2D
    title = 'Re(ez_2D)'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
	c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
        aspect_ratio=1.0, layout=[2,3,5], /current, title=title )
	c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
    p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )
    p = plot ( rLim, zLim, /over )

    ThisField = Abs(er_2D)
    title = 'Abs(er_2D)'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
	c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
        aspect_ratio=1.0, layout=[2,3,2], /current, title=title )
	c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
    p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )
    p = plot ( rLim, zLim, /over )

    ThisField = Abs(et_2D)
    title = 'Abs(et_2D)'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
	c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
        aspect_ratio=1.0, layout=[2,3,4], /current, title=title )
	c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
    p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )
    p = plot ( rLim, zLim, /over )

    ThisField = Abs(ez_2D)
    title = 'Abs(ez_2D)'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
	c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
        aspect_ratio=1.0, layout=[2,3,6], /current, title=title )
	c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
    p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )
    p = plot ( rLim, zLim, /over )


    ThisSpecies = 0
    ScaleFac = 0.3
    scale = max(abs([jPr_2D[*,*,ThisSpecies],jPt_2D[*,*,ThisSpecies],jPz_2D[*,*,ThisSpecies]]))*ScaleFac

	levels = fIndGen(nLevs)/(nLevs-1)*scale
	colors = reverse(bytScl(levels, top=253)+1)

    ThisField = jPr_2D[*,*,ThisSpecies]
    title = 'Re(jPr_2D)'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
	c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
        aspect_ratio=1.0, layout=[2,3,1], dimensions=dimensions, title=title )
	c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
    p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )
    p = plot ( rLim, zLim, /over )

    ThisField = jPt_2D[*,*,ThisSpecies]
    title = 'Re(jPt_2D)'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
	c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
        aspect_ratio=1.0, layout=[2,3,3], /current, title=title )
	c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
    p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )
    p = plot ( rLim, zLim, /over )

    ThisField = jPz_2D[*,*,ThisSpecies]
    title = 'Re(jPz_2D)'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
	c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
        aspect_ratio=1.0, layout=[2,3,5], /current, title=title )
	c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
    p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )
    p = plot ( rLim, zLim, /over )

    ThisField = Abs(jPr_2D[*,*,ThisSpecies])
    title = 'Abs(jPr_2D)'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
	c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
        aspect_ratio=1.0, layout=[2,3,2], /current, title=title )
	c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
    p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )
    p = plot ( rLim, zLim, /over )

    ThisField = Abs(jPt_2D[*,*,ThisSpecies])
    title = 'Abs(jPt_2D)'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
	c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
        aspect_ratio=1.0, layout=[2,3,4], /current, title=title )
	c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
    p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )
    p = plot ( rLim, zLim, /over )

    ThisField = Abs(jPz_2D[*,*,ThisSpecies])
    title = 'Abs(jPz_2D)'
	PlotField = (real_part(thisField)<max(levels))>min(-levels)
	c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
        aspect_ratio=1.0, layout=[2,3,6], /current, title=title )
	c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
    p = plot( CurrentSource_r, CurrentSource_z, /over, symbol = "s", lineStyle='none' )
    p = plot ( rLim, zLim, /over )

end
