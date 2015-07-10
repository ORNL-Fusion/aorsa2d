pro ar2_create_current_source_location_file

	; Set the R,Z locations (cylindrical coordinates)
	; for the RHS current sources. Each current source
	; location means 3 RHS vectors as we need to account
	; for each component of the source vector. 

    ; Right simple
    ;SingleCurrentAtEnd = 1
    ; Left simple
    SingleCurrentAtEnd = 0

    patch2 = 0

	nRSources = 1 ; x 3 for all components, real & imag
    nZSources = 1 
    nSourcesTotal = nRSources * nZSources
	NRHS = nSourcesTotal * 3

    if patch2 eq 1 then begin

        nRSources2 = 3 
        nZSources2 = 20 
        NRHS2 = nRSources2 * nZSources2 * 3

        NRHS = NRHS + NRHS2
        nSourcesTotal = nSourcesTotal + nRSources2*nZSources2

    endif

    if SingleCurrentAtEnd eq 1 then begin

        NRHS = NRHS + 1

    endif

	cs_r = fltArr(NRHS)
	cs_z = fltArr(NRHS)
	component_ident = intArr(NRHS)
    realFac = fltArr(NRHS)
    imagFac = fltArr(NRHS)

	ii=0

    ; Right simple
	rMin = 1.85
    rMax = 2.00
    ; Left simple
    rMin = 2.5
    rMax = 2.5

	zMin = 0.3
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

    realFac[*] = 1
    imagFac[*] = 0

    ;NRHS = NRHS+n_elements(cs_r)
    ;ii = ii + n_elements(cs_r)
    ;cs_r = [cs_r,cs_r]
    ;cs_z = [cs_z,cs_z]
    ;component_ident = [component_ident,component_ident]
    ;realFac = [realFac,realFac*0]
    ;imagFac = [imagFac,imagFac+1]

    if SingleCurrentAtEnd eq 1 then begin

        cs_r[-1] = 2.5
        cs_z[-1] = 0.0
        component_ident[-1] = 2 ; z-direction     
        realFac[-1] = 1
        imagFac[-1] = 0
        ii++
    endif

	if ii ne NRHS then message, 'ERROR'
    
    iiBad = where(cs_r ne cs_r,iiBadCnt)
    if iiBadCnt gt 0 then message, 'ERROR'

    iiBad = where(cs_z ne cs_z,iiBadCnt)
    if iiBadCnt gt 0 then message, 'ERROR' 


    ; Add the standard single RHS current source as the last element


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

	realFac_id = nCdf_varDef ( nc_id, 'realFac', nRHS_id, /float )
	if realFac_id lt 0 then stop
	imagFac_id = nCdf_varDef ( nc_id, 'imagFac', nRHS_id, /float )
	if imagFac_id lt 0 then stop


	nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, cs_r_id, cs_r
	nCdf_varPut, nc_id, cs_z_id, cs_z
	nCdf_varPut, nc_id, component_ident_id, component_ident
	nCdf_varPut, nc_id, realFac_id, realFac
	nCdf_varPut, nc_id, imagFac_id, imagFac

	nCdf_close, nc_id

	print, 'Success'

    doPlot = 0
    if doPlot then p = plot( cs_r, cs_z, symbol = "s", lineStyle='none' )

end
