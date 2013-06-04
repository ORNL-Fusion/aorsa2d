pro ar2_create_current_source_location_file

	; Set the R,Z locations (cylindrical coordinates)
	; for the RHS current sources. Each current source
	; location means 3 RHS vectors as we need to account
	; for each component of the source vector. 

	nRSources = 10 ; x 3 x 2 for all components, real & imag
    nZSources = 30
	NRHS = nRSources * nZSources * 3

	cs_r = fltArr(NRHS)
	cs_z = fltArr(NRHS)
	component_ident = intArr(NRHS)

	ii=0

	;rLoc = 2.075
	;zMin = -0.2
	;zMax = 0.2

	rMin = 1.8
    rMax = 1.9
	zMin = -0.17
	zMax = 0.17

	for i=0,nRSources-1 do begin
        for j = 0,nZSources-1 do begin
		    for c=0,2 do begin
		    	cs_r[ii] = rMin+(rMax-rMin)/(nRSources-1)*i
		    	cs_z[ii] = zMin+(zMax-zMin)/(nZSources-1)*i
		    	component_ident[ii] = c
		    	ii++
		    endfor
        endfor
	endfor
	if ii ne NRHS then stop

	nc_id = nCdf_create ('AR2SourceLocations.nc', /clobber )

	nCdf_control, nc_id, /fill
	
	nSources_id = nCdf_dimDef ( nc_id, 'nSources', nRSources * nZSources ) 
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

end
