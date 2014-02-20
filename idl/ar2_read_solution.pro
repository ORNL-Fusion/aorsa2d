function ar2_read_solution, runFolderName, RHS

	; This is just to ensure we get consistent files, i.e., not
	; relying on the order of the list of a file read.
	SolutionFiles = file_search(runFolderName+'/output/solution*.nc')
	SolutionStr = StrMid(file_basename(SolutionFiles[0]),0,18)
	RHS_Str = string(RHS,format='(i6.6)')
    SolutionFile = File_DirName(SolutionFiles[0])+'/'+SolutionStr+RHS_STr+'.nc'
	;print, SolutionFile

	cdfId = ncdf_open ( SolutionFile, /noWrite ) 

		nCdf_varGet, cdfId, 'r', r
		nCdf_varGet, cdfId, 'z', z
		nCdf_varGet, cdfId, 'nPhi', nPhi 
		nCdf_varGet, cdfId, 'freqcy', freq 
	
		nCdf_varGet, cdfId, 'ealpha_re', ealpha_re 
		ealpha_re = temporary(ealpha_re[*,*])
		nCdf_varGet, cdfId, 'ebeta_re', ebeta_re 
		ebeta_re = temporary(ebeta_re[*,*])
		nCdf_varGet, cdfId, 'eB_re', eB_re 
		eb_re = temporary(eb_re[*,*])

		nCdf_varGet, cdfId, 'ealpha_im', ealpha_im 
		ealpha_im = temporary(ealpha_im[*,*])
		nCdf_varGet, cdfId, 'ebeta_im', ebeta_im 
		ebeta_im = temporary(ebeta_im[*,*])
		nCdf_varGet, cdfId, 'eB_im', eB_im 
		eb_im = temporary(eb_im[*,*])

		nCdf_varGet, cdfId, 'ealphak_re', ealphak_re 
		ealphak_re = temporary(ealphak_re[*,*])
		nCdf_varGet, cdfId, 'ebetak_re', ebetak_re 
		ebetak_re = temporary(ebetak_re[*,*])
		nCdf_varGet, cdfId, 'eBk_re', eBk_re 
		ebk_re = temporary(ebk_re[*,*])

		nCdf_varGet, cdfId, 'ealphak_im', ealphak_im 
		ealphak_im = temporary(ealphak_im[*,*])
		nCdf_varGet, cdfId, 'ebetak_im', ebetak_im 
		ebetak_im = temporary(ebetak_im[*,*])
		nCdf_varGet, cdfId, 'eBk_im', eBk_im 
		ebk_im = temporary(ebk_im[*,*])

		nCdf_varGet, cdfId, 'er_re', er_re 
		er_re = temporary(er_re[*,*])
		nCdf_varGet, cdfId, 'et_re', et_re 
		et_re = temporary(et_re[*,*])
		nCdf_varGet, cdfId, 'ez_re', ez_re 
		ez_re = temporary(ez_re[*,*])

		nCdf_varGet, cdfId, 'er_im', er_im 
		er_im = temporary(er_im[*,*])
		nCdf_varGet, cdfId, 'et_im', et_im 
		et_im = temporary(et_im[*,*])
		nCdf_varGet, cdfId, 'ez_im', ez_im 
		ez_im = temporary(ez_im[*,*])

		nCdf_varGet, cdfId, 'jalpha_re', jalpha_re 
		nCdf_varGet, cdfId, 'jbeta_re', jbeta_re 
		nCdf_varGet, cdfId, 'jB_re', jB_re 

		nCdf_varGet, cdfId, 'jalpha_im', jalpha_im 
		nCdf_varGet, cdfId, 'jbeta_im', jbeta_im 
		nCdf_varGet, cdfId, 'jB_im', jB_im 

		nCdf_varGet, cdfId, 'jP_r_re', jPr_re 
		nCdf_varGet, cdfId, 'jP_t_re', jPt_re 
		nCdf_varGet, cdfId, 'jP_z_re', jPz_re 

		nCdf_varGet, cdfId, 'jP_r_im', jPr_im 
		nCdf_varGet, cdfId, 'jP_t_im', jPt_im 
		nCdf_varGet, cdfId, 'jP_z_im', jPz_im 

		nCdf_varGet, cdfId, 'jouleHeating', jouleHeating 

		ealpha	= complex ( ealpha_re, ealpha_im )
		ebeta	= complex ( ebeta_re, ebeta_im )
		eb	= complex ( eb_re, eb_im )

		ealphak	= complex ( ealphak_re, ealphak_im )
		ebetak	= complex ( ebetak_re, ebetak_im )
		ebk	= complex ( ebk_re, ebk_im )

		e_r	= complex ( er_re, er_im )
		e_t	= complex ( et_re, et_im )
		e_z	= complex ( ez_re, ez_im )

		jP_r	= complex ( jPr_re, jPr_im )
		jP_t	= complex ( jPt_re, jPt_im )
		jP_z	= complex ( jPz_re, jPz_im )

		jPAlpha = complex ( jAlpha_re, jAlpha_im )
		jPBeta = complex ( jBeta_re, jBeta_im )
		jPB = complex ( jB_re, jB_im )

	ncdf_close, cdfId

	x = r
	y = z
	nx = n_elements(r)
	ny = n_elements(z)
	x2D	= rebin ( r, nX, nY )
	if nY gt 1 then begin
		y2D = transpose(rebin ( z, nY, nX ))
	endif else begin
		y2D = fltArr(1)+z
	endelse

    solution = { $
                r: r, $
                z: z, $
                x: x, $
                y: y, $
                x2d: x2d, $
                y2d: y2d, $
                jPAlpha: jPAlpha, $
                jPBeta: jPBeta, $
                jPB: jPB, $
                jP_r: jP_r, $
                jP_t: jP_t, $
                jP_z: jP_z, $
                e_r: e_r, $
                e_t: e_t, $
                e_z: e_z }


    return, solution

end



