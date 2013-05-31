pro ar2_create_current_source_location_file

	; Set the R,Z locations (cylindrical coordinates)
	; for the RHS current sources. Each current source
	; location means 3 RHS vectors as we need to account
	; for each component of the source vector. 

	nSources = 30 ; x 3 x 2 for all components, real & imag
	NRHS = nSources * 3

	cs_r = fltArr(NRHS)
	cs_z = fltArr(NRHS)
	component_ident = intArr(NRHS)

	j=0

	;rLoc = 2.075
	;zMin = -0.2
	;zMax = 0.2

	rLoc = 2.2
	zMin = -0.15
	zMax = 0.15

	for i=0,nSources-1 do begin
		for c=0,2 do begin
			cs_r[j] = rLoc
			cs_z[j] = zMin+(zMax-zMin)/(nSources-1)*i
			component_ident[j] = c
			j++
		endfor
	endfor
	if j ne NRHS then stop

	nc_id = nCdf_create ('AR2SourceLocations.nc', /clobber )

	nCdf_control, nc_id, /fill
	
	nSources_id = nCdf_dimDef ( nc_id, 'nSources', nSources ) 
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
