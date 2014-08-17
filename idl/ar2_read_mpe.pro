function ar2_read_mpe, FileName

	File = expand_path(FileName)

	cdfId = ncdf_open ( File, /noWrite ) 
		nCdf_varGet, cdfId, 'r', r 

		nCdf_varGet, cdfId, 'E_r_re', Er_re 
		nCdf_varGet, cdfId, 'E_r_im', Er_im 
		nCdf_varGet, cdfId, 'E_p_re', Et_re 
		nCdf_varGet, cdfId, 'E_p_im', Et_im 
		nCdf_varGet, cdfId, 'E_z_re', Ez_re 
		nCdf_varGet, cdfId, 'E_z_im', Ez_im 

		nCdf_varGet, cdfId, 'jP_r_re', jPr_re 
		nCdf_varGet, cdfId, 'jP_r_im', jPr_im 
		nCdf_varGet, cdfId, 'jP_p_re', jPt_re 
		nCdf_varGet, cdfId, 'jP_p_im', jPt_im 
		nCdf_varGet, cdfId, 'jP_z_re', jPz_re 
		nCdf_varGet, cdfId, 'jP_z_im', jPz_im 
	
	ncdf_close, cdfId

    E_r = complex(Er_re,Er_im)
    E_t = complex(Et_re,Et_im)
    E_z = complex(Ez_re,Ez_im)

    jP_r = complex(jPr_re,jPr_im)
    jP_t = complex(jPt_re,jPt_im)
    jP_z = complex(jPz_re,jPz_im)

    solution = { $
                r: r, $
                e_r: e_r, $
                e_t: e_t, $
                e_z: e_z, $
                jP_r: jP_r, $
                jP_t: jP_t, $
                jP_z: jP_z }

    return, solution

end



