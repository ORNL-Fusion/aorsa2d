pro ar2_create_current_source_location_file

	; Set the R,Z locations (cylindrical coordinates)
	; for the RHS current sources. Each current source
	; location means 3 RHS vectors as we need to account
	; for each component of the source vector. 

    patch2 = 1

	nRSources = 5 ; x 3 for all components, real & imag
    nZSources = 15 
    nSourcesTotal = nRSources * nZSources
	NRHS = nSourcesTotal * 3

    if patch2 eq 1 then begin

        nRSources2 = 3 
        nZSources2 = 20 
        NRHS2 = nRSources2 * nZSources2 * 3

        NRHS = NRHS + NRHS2
        nSourcesTotal = nSourcesTotal + nRSources2*nZSources2

    endif

	cs_r = fltArr(NRHS)
	cs_z = fltArr(NRHS)
	component_ident = intArr(NRHS)

	ii=0

	rMin = 1.8
    rMax = 1.85
	zMin = -0.17
	zMax = 0.17

	for i=0,nRSources-1 do begin
        for j = 0,nZSources-1 do begin
		    for c=0,2 do begin
                if nRSources gt 1 then $
                    cs_r[ii] = rMin+(rMax-rMin)/(nRSources-1)*i else cs_r[ii] = rMin
		    	if nZSources gt 1 then $
                    cs_z[ii] = zMin+(zMax-zMin)/(nZSources-1)*j else cs_z[ii] = zMin
		    	component_ident[ii] = c
		    	ii++
		    endfor
        endfor
	endfor

    if patch2 eq 1 then begin

	    rMin = 1.98
        rMax = 2.02
	    zMin = -0.17
	    zMax = 0.17

	    for i=0,nRSources2-1 do begin
            for j = 0,nZSources2-1 do begin
	    	    for c=0,2 do begin
                    if nRSources2 gt 1 then $
	    	    	    cs_r[ii] = rMin+(rMax-rMin)/(nRSources2-1)*i else cs_r[ii] = rMin
                    if nZSources2 gt 1 then $
	    	    	    cs_z[ii] = zMin+(zMax-zMin)/(nZSources2-1)*j else cs_z[ii] = zMin
	    	    	component_ident[ii] = c
	    	    	ii++
	    	    endfor
            endfor
	    endfor

    endif


	if ii ne NRHS then message, 'ERROR'
    
    iiBad = where(cs_r ne cs_r,iiBadCnt)
    if iiBadCnt gt 0 then message, 'ERROR'

    iiBad = where(cs_z ne cs_z,iiBadCnt)
    if iiBadCnt gt 0 then message, 'ERROR' 

	nc_id = nCdf_create ('ar2SourceLocations.nc', /clobber )

	nCdf_control, nc_id, /fill

	nSources_id = nCdf_dimDef ( nc_id, 'nSources', nSourcesTotal ) 
	if nSources_id lt 0 then stop
	nRHS_id = nCdf_dimDef ( nc_id, 'NRHS', NRHS )
	if nRHS_id lt 0 then stop

	cs_r_id = nCdf_varDef ( nc_id, 'cs_r', nRHS_id, /float )
	if cs_r_id lt 0 then stop

	cs_z_id = nCdf_varDef ( nc_id, 'cs_z', nRHS_id, /float )
	if cs_z_id lt 0 then stop

	component_ident_id = nCdf_varDef ( nc_id, 'component_ident', nRHS_id, /short )
	if component_ident_id lt 0 then stop

	nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, cs_r_id, cs_r
	nCdf_varPut, nc_id, cs_z_id, cs_z
	nCdf_varPut, nc_id, component_ident_id, component_ident

	nCdf_close, nc_id

	print, 'Success'

    p = plot( cs_r, cs_z, symbol = "s", lineStyle='none' )

end
