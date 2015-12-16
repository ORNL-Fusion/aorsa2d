function ar2_read_solution, runFolderName, RHS
   
    @constants

	;; This is just to ensure we get consistent files, i.e., not
	;; relying on the order of the list of a file read.
	;SolutionFiles = file_search(runFolderName+'/output/solution*.nc')
	;SolutionStr = StrMid(file_basename(SolutionFiles[0]),0,18)
	;RHS_Str = string(RHS,format='(i6.6)')
    ;SolutionFile = File_DirName(SolutionFiles[0])+'/'+SolutionStr+RHS_STr+'.nc'

	ar2Input = ar2_read_namelist( RunFolderName = RunFolderName)
	ThisGridNo = 1
	GridNoStr = string(ThisGridNo,format='(i3.3)')
	ThisNPhi = ar2Input['nPhi']
    ThisRHS = RHS
	nPhiStr = string(ThisNPhi,format='(i+4.3)')
	rhsStr = string(ThisRHS,format='(i6.6)')

	SolutionFile = expand_path(RunFolderName)+'/output/solution_'+GridNoStr+'_'+nPhiStr+'_'+rhsStr+'.nc'
	RunDataFile = expand_path(RunFolderName)+'/output/runData_'+GridNoStr+'_'+nPhiStr+'_'+rhsStr+'.nc'

	;print, SolutionFile

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
        nCdf_varGet, cdfId, 'nZ_1D', nz_1D
	ncdf_close, cdfId


	cdfId = ncdf_open ( SolutionFile, /noWrite ) 

		nCdf_varGet, cdfId, 'r', r
		nCdf_varGet, cdfId, 'z', z
		nCdf_varGet, cdfId, 'nPhi', nPhi 
		nCdf_varGet, cdfId, 'freq', freq 
	
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

    if nY gt 1 then begin
            print, 'NOT CALCULATING B1 ... yet to implement in 2-D'
    endif else begin

        ; Calculate the H vector & the Poynting vector

         
        h_r = complexArr(nX)
        h_t = complexArr(nX)
        h_z = complexArr(nX)

        k_z = nZ_1D*0 ; this is NOT right, fix for non-zero nZ

        dr = r[1]-r[0]
        for i=2,nX-3 do begin

            dEz_dr = (1.0/12.0*e_z[i-2] - 2.0/3.0*e_z[i-1]$
                    +2.0/3.0*e_z[i+1] - 1.0/12.0*e_z[i+2])/dr
            drEt_dr = (1.0/12.0*r[i-2]*e_t[i-2] - 2.0/3.0*r[i-1]*e_t[i-1]$
                    +2.0/3.0*r[i+1]*e_t[i+1] - 1.0/12.0*r[i+2]*e_t[i+2])/dr

            h_r[i] = -II*k_z*e_t[i] + II*nPhi*e_z[i]/r[i]
            h_t[i] = II*k_z*e_r[i] - dEz_dr 
            h_z[i] = (-II*nPhi*e_r[i] + drEt_dr )/r[i]

        endfor

        wrf = 2 * !pi * freq
        
        h_r = h_r / (II*wRF*u0)
        h_t = h_t / (II*wRF*u0)
        h_z = h_z / (II*wRF*u0)

        b_r = u0 * h_r
        b_t = u0 * h_t
        b_z = u0 * h_z

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
                e_z: e_z, $
                b1_r: b_r, $
                b1_t: b_t, $
                b1_z: b_z, $
                jA_r: complex(jr_re,jr_im), $
                jA_t: complex(jt_re,jt_im), $
                jA_z: complex(jz_re,jz_im) }

    return, solution

end



