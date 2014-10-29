pro ar2_create_kj_input

    runFolderName = './'
    RHS = 1

    ar2 = ar2_read_namelist(runFolderName = runFolderName)
    data = ar2_read_rundata(runFolderName, RHS)
    s = ar2_read_solution(runFolderName, RHS)

    nR = n_elements(data.r)
    nS = n_elements(data.densitySpec[0,0,*])

	nc_id = nCdf_create ('ar2_kj.nc', /clobber )

	nCdf_control, nc_id, /fill
	
	nr_id = nCdf_dimDef ( nc_id, 'nR', nR )
	scalar_id = nCdf_dimDef ( nc_id, 'scalar', 1 )
	nSpec_id = nCdf_dimDef ( nc_id, 'nSpec', nS )

	freq_id = nCdf_varDef ( nc_id, 'freq', scalar_id, /float )
	r_id = nCdf_varDef ( nc_id, 'r', nr_id, /float )
	z_id = nCdf_varDef ( nc_id, 'z', nr_id, /float )

	B0_r_id = nCdf_varDef ( nc_id, 'B0_r', nr_id, /float )
	B0_p_id = nCdf_varDef ( nc_id, 'B0_p', nr_id, /float )
	B0_z_id = nCdf_varDef ( nc_id, 'B0_z', nr_id, /float )

	e_r_re_id = nCdf_varDef ( nc_id, 'e_r_re', nr_id, /float )
	e_r_im_id = nCdf_varDef ( nc_id, 'e_r_im', nr_id, /float )
	e_p_re_id = nCdf_varDef ( nc_id, 'e_p_re', nr_id, /float )
	e_p_im_id = nCdf_varDef ( nc_id, 'e_p_im', nr_id, /float )
	e_z_re_id = nCdf_varDef ( nc_id, 'e_z_re', nr_id, /float )
	e_z_im_id = nCdf_varDef ( nc_id, 'e_z_im', nr_id, /float )
    
	b_r_re_id = nCdf_varDef ( nc_id, 'b_r_re', nr_id, /float )
	b_r_im_id = nCdf_varDef ( nc_id, 'b_r_im', nr_id, /float )
	b_p_re_id = nCdf_varDef ( nc_id, 'b_p_re', nr_id, /float )
	b_p_im_id = nCdf_varDef ( nc_id, 'b_p_im', nr_id, /float )
	b_z_re_id = nCdf_varDef ( nc_id, 'b_z_re', nr_id, /float )
	b_z_im_id = nCdf_varDef ( nc_id, 'b_z_im', nr_id, /float )

	jP_r_re_id = nCdf_varDef ( nc_id, 'jP_r_re', nr_id, /float )
	jP_r_im_id = nCdf_varDef ( nc_id, 'jP_r_im', nr_id, /float )
	jP_p_re_id = nCdf_varDef ( nc_id, 'jP_p_re', nr_id, /float )
	jP_p_im_id = nCdf_varDef ( nc_id, 'jP_p_im', nr_id, /float )
	jP_z_re_id = nCdf_varDef ( nc_id, 'jP_z_re', nr_id, /float )
	jP_z_im_id = nCdf_varDef ( nc_id, 'jP_z_im', nr_id, /float )

	jA_r_re_id = nCdf_varDef ( nc_id, 'jA_r_re', nr_id, /float )
	jA_r_im_id = nCdf_varDef ( nc_id, 'jA_r_im', nr_id, /float )
	jA_p_re_id = nCdf_varDef ( nc_id, 'jA_p_re', nr_id, /float )
	jA_p_im_id = nCdf_varDef ( nc_id, 'jA_p_im', nr_id, /float )
	jA_z_re_id = nCdf_varDef ( nc_id, 'jA_z_re', nr_id, /float )
	jA_z_im_id = nCdf_varDef ( nc_id, 'jA_z_im', nr_id, /float )

	jP_r_re_spec_id = nCdf_varDef ( nc_id, 'jP_r_re_spec', [nr_id,nSpec_id], /float )
	jP_r_im_spec_id = nCdf_varDef ( nc_id, 'jP_r_im_spec', [nr_id,nSpec_id], /float )
	jP_p_re_spec_id = nCdf_varDef ( nc_id, 'jP_p_re_spec', [nr_id,nSpec_id], /float )
	jP_p_im_spec_id = nCdf_varDef ( nc_id, 'jP_p_im_spec', [nr_id,nSpec_id], /float )
	jP_z_re_spec_id = nCdf_varDef ( nc_id, 'jP_z_re_spec', [nr_id,nSpec_id], /float )
	jP_z_im_spec_id = nCdf_varDef ( nc_id, 'jP_z_im_spec', [nr_id,nSpec_id], /float )

	Density_id = nCdf_varDef ( nc_id, 'density_m3', [nr_id,nSpec_id], /float )

	nCdf_control, nc_id, /enDef

	nCdf_varPut, nc_id, freq_id, ar2['freq']

	nCdf_varPut, nc_id, r_id, data.r

	nCdf_varPut, nc_id, z_id, data.z

	nCdf_varPut, nc_id, B0_r_id, data.br
	nCdf_varPut, nc_id, B0_p_id, data.bt
	nCdf_varPut, nc_id, B0_z_id, data.bz

	nCdf_varPut, nc_id, e_r_re_id, real_part(s.e_r) 
	nCdf_varPut, nc_id, e_r_im_id, imaginary(s.e_r) 
	nCdf_varPut, nc_id, e_p_re_id, real_part(s.e_t) 
	nCdf_varPut, nc_id, e_p_im_id, imaginary(s.e_t) 
	nCdf_varPut, nc_id, e_z_re_id, real_part(s.e_z) 
	nCdf_varPut, nc_id, e_z_im_id, imaginary(s.e_z) 

	nCdf_varPut, nc_id, b_r_re_id, real_part(s.e_r*0) 
	nCdf_varPut, nc_id, b_r_im_id, imaginary(s.e_r*0) 
	nCdf_varPut, nc_id, b_p_re_id, real_part(s.e_t*0) 
	nCdf_varPut, nc_id, b_p_im_id, imaginary(s.e_t*0) 
	nCdf_varPut, nc_id, b_z_re_id, real_part(s.e_z*0) 
	nCdf_varPut, nc_id, b_z_im_id, imaginary(s.e_z*0) 

	nCdf_varPut, nc_id, jP_r_re_id, real_part(total(s.jP_r,3))
	nCdf_varPut, nc_id, jP_r_im_id, imaginary(total(s.jP_r,3)) 
	nCdf_varPut, nc_id, jP_p_re_id, real_part(total(s.jP_t,3)) 
	nCdf_varPut, nc_id, jP_p_im_id, imaginary(total(s.jP_t,3)) 
	nCdf_varPut, nc_id, jP_z_re_id, real_part(total(s.jP_z,3)) 
	nCdf_varPut, nc_id, jP_z_im_id, imaginary(total(s.jP_z,3)) 

	nCdf_varPut, nc_id, jA_r_re_id, real_part(s.jA_r) 
	nCdf_varPut, nc_id, jA_r_im_id, imaginary(s.jA_r) 
	nCdf_varPut, nc_id, jA_p_re_id, real_part(s.jA_t)
	nCdf_varPut, nc_id, jA_p_im_id, imaginary(s.jA_t)
	nCdf_varPut, nc_id, jA_z_re_id, real_part(s.jA_z)
	nCdf_varPut, nc_id, jA_z_im_id, imaginary(s.jA_z)

    SpecShift = nS-1 
	nCdf_varPut, nc_id, jP_r_re_spec_id, real_part(shift(reform(s.jP_r),[0,SpecShift]))
	nCdf_varPut, nc_id, jP_r_im_spec_id, imaginary(shift(reform(s.jP_r),[0,SpecShift])) 
	nCdf_varPut, nc_id, jP_p_re_spec_id, real_part(shift(reform(s.jP_t),[0,SpecShift]))
	nCdf_varPut, nc_id, jP_p_im_spec_id, imaginary(shift(reform(s.jP_t),[0,SpecShift])) 
	nCdf_varPut, nc_id, jP_z_re_spec_id, real_part(shift(reform(s.jP_z),[0,SpecShift]))
	nCdf_varPut, nc_id, jP_z_im_spec_id, imaginary(shift(reform(s.jP_z),[0,SpecShift])) 

    TmpDensity = FltArr(nR,nS)
    for s=0,nS-1 do begin
        TmpDensity[*,s] = data.densitySpec[*,0,s]
    endfor

    nCdf_varPut, nc_id, Density_id, TmpDensity

	nCdf_close, nc_id

    stop
end
