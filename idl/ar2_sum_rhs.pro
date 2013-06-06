pro ar2_sum_rhs

	SourceLocationsFile = 'ar2SourceLocations.nc'
	RunDataFiles = file_search('output/runData*.nc')
    SolutionFiles = file_search('output/solution*.nc')
    
    NRHS_check = n_elements(RunDataFiles)

    ar2_fit_sources, CoeffsOut2D = Coeffs, RunDataFiles = RunDataFiles, /TwoDim

    NRHS = n_elements(Coeffs)

    ; Load the fields for all RHSs

    for rhs=0,NRHS-1 do begin

	    cdfId = ncdf_open ( RunDataFiles[rhs], /noWrite ) 

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



end
