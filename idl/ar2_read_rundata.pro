function ar2_read_rundata, runFolderName, RHS

	;RunDataFiles = file_search(runFolderName+'/output/solution*.nc')
    ;RunDataFile = RunDataFiles[RHS-1]

	; This is just to ensure we get consistent files, i.e., not
	; relying on the order of the list of a file read.
	RunDataFiles = file_search(runFolderName+'/output/runData*.nc')
	RunDataStr = StrMid(file_basename(RunDataFiles[0]),0,17)
	RHS_Str = string(RHS,format='(i6.6)')
    RunDataFile = File_DirName(RunDataFiles[0])+'/'+RunDataStr+RHS_STr+'.nc'
	;print, RunDataFile


	;ar2_read_ar2input, ar2InputFile, RunDataFile, $
	;		rLim=rLim,zLim=zLim,LimMask=LimMask

	cdfId = ncdf_open ( RunDataFile, /noWrite ) 
		nCdf_varGet, cdfId, 'nPhi', nPhi 
		nCdf_varGet, cdfId, 'freq', freq 
		nCdf_varGet, cdfId, 'capR', x 
		nCdf_varGet, cdfId, 'y', y 
		nCdf_varGet, cdfId, 'jr_re', jr_re 
		nCdf_varGet, cdfId, 'jr_im', jr_im 
		nCdf_varGet, cdfId, 'jt_re', jt_re 
		nCdf_varGet, cdfId, 'jt_im', jt_im 
		nCdf_varGet, cdfId, 'jz_re', jz_re 
		nCdf_varGet, cdfId, 'jz_im', jz_im 
		nCdf_varGet, cdfId, 'bmod', bmod
		nCdf_varGet, cdfId, 'brU', brU
		nCdf_varGet, cdfId, 'btU', btU
		nCdf_varGet, cdfId, 'bzU', bzU
		nCdf_varGet, cdfId, 'densitySpec', densitySpec
		nCdf_varGet, cdfId, 'tempSpec', tempSpec
		nCdf_varGet, cdfId, 'xx_re', xx_re 
		nCdf_varGet, cdfId, 'xx_im', xx_im
		nCdf_varGet, cdfId, 'yy_re', yy_re 
		nCdf_varGet, cdfId, 'yy_im', yy_im
		nCdf_varGet, cdfId, 'drBfn_re', dRbFn_bFn_re
		nCdf_varGet, cdfId, 'drBfn_im', dRbFn_bFn_im
		nCdf_varGet, cdfId, 'dzBfn_re', dzbFn_bFn_re
		nCdf_varGet, cdfId, 'dzBfn_im', dzbFn_bFn_im
		nCdf_varGet, cdfId, 'nuOmg', nuOmg
		nCdf_varGet, cdfId, 'LimMask', LimMask
	ncdf_close, cdfId

	xx	= complex ( xx_re, xx_im )
	yy	= complex ( yy_re, yy_im )

    ;nuOmg = nuOmg[*,*,ThisSpec]

	jA_r = complex ( jr_re[*,*], jr_im[*,*] )
	jA_t = complex ( jt_re[*,*], jt_im[*,*] )
	jA_z = complex ( jz_re[*,*], jz_im[*,*] )

	dRbFn_bFn	= complex ( dRbFn_bFn_re, dRbFn_bFn_im )
	dZbFn_bFn	= complex ( dZbFn_bFn_re, dZbFn_bFn_im )

	nN	= n_elements ( xx[*,0] )
	nM	= n_elements ( yy[*,0] )

    runData = { $
        r : x, $
        z : y, $
        jA_r : complex(jr_re,jr_im), $ 
        jA_t : complex(jt_re,jt_im), $ 
        jA_z : complex(jz_re,jz_im), $
        densitySpec : densitySpec, $
        br : brU * bmod, $
        bt : btU * bmod, $
        bz : bzU * bmod, $
	   	nPhi : nPhi, $
	   	freq : freq }

    return, runData
end
