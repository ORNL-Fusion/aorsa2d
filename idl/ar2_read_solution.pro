function ar2_read_solution, runFolderName, RHS
   
    @dlg_constants

	ar2Input = ar2_read_namelist( RunFolderName = RunFolderName)
	ThisGridNo = 1
	GridNoStr = string(ThisGridNo,format='(i3.3)')
	ThisNPhi = ar2Input['nPhi']
    ThisRHS = RHS
	nPhiStr = string(ThisNPhi,format='(i+7.6)')
	rhsStr = string(ThisRHS,format='(i6.6)')

	SolutionFile = expand_path(RunFolderName)+'/output/solution_'+GridNoStr+'_'+nPhiStr+'_'+rhsStr+'.nc'
	RunDataFile = expand_path(RunFolderName)+'/output/runData_'+GridNoStr+'_'+nPhiStr+'_'+rhsStr+'.nc'

    arR = dlg_read_netcdf(RunDataFile)

    arS = dlg_read_netcdf(SolutionFile)

    return, arS

    ;r = arS['r']
    ;z = arS['z']

    ;E_r = arS['E_r']
    ;E_t = arS['E_t']
    ;E_z = arS['E_z']

    ;nPhi = arS['nPhi']
    ;freq = arS['freq']

    ;kz_1d = arR['kz_1d']

	;x = r
	;y = z
	;nx = n_elements(r)
	;ny = n_elements(z)
	;x2D	= rebin ( r, nX, nY )
	;if nY gt 1 then begin
	;	y2D = transpose(rebin ( z, nY, nX ))
	;endif else begin
	;	y2D = fltArr(1)+z
	;endelse

    ;if nY gt 1 then begin
    ;        print, 'NOT CALCULATING B1 ... yet to implement in 2-D'
    ;endif else begin

    ;    ; Calculate the H vector & the Poynting vector

    ;     
    ;    h_r = complexArr(nX)
    ;    h_t = complexArr(nX)
    ;    h_z = complexArr(nX)

    ;    k_z = kz_1d; this is NOT right, fix for non-zero nZ

    ;    dr = r[1]-r[0]
    ;    for i=2,nX-3 do begin

    ;        dEz_dr = (1.0/12.0*e_z[i-2] - 2.0/3.0*e_z[i-1]$
    ;                +2.0/3.0*e_z[i+1] - 1.0/12.0*e_z[i+2])/dr
    ;        drEt_dr = (1.0/12.0*r[i-2]*e_t[i-2] - 2.0/3.0*r[i-1]*e_t[i-1]$
    ;                +2.0/3.0*r[i+1]*e_t[i+1] - 1.0/12.0*r[i+2]*e_t[i+2])/dr

    ;        h_r[i] = -_ii*k_z*e_t[i] + _ii*nPhi*e_z[i]/r[i]
    ;        h_t[i] = _ii*k_z*e_r[i] - dEz_dr 
    ;        h_z[i] = (-_ii*nPhi*e_r[i] + drEt_dr )/r[i]

    ;    endfor

    ;    wrf = 2 * !pi * freq
    ;    
    ;    h_r = h_r / (_ii*wRF*_u0)
    ;    h_t = h_t / (_ii*wRF*_u0)
    ;    h_z = h_z / (_ii*wRF*_u0)

    ;    b_r = _u0 * h_r
    ;    b_t = _u0 * h_t
    ;    b_z = _u0 * h_z

    ;endelse

    ;return, arS

    ;;solution = { $
    ;;            r: r, $
    ;;            z: z, $
    ;;            x: x, $
    ;;            y: y, $
    ;;            x2d: x2d, $
    ;;            y2d: y2d, $
    ;;            jPAlp: jPAlpha, $
    ;;            jPBet: jPBeta, $
    ;;            jPPrl: jPB, $
    ;;            jP_r: jP_r, $
    ;;            jP_t: jP_t, $
    ;;            jP_z: jP_z, $
    ;;            e_r: e_r, $
    ;;            e_t: e_t, $
    ;;            e_z: e_z, $
    ;;            ealpk: ealphak, $
    ;;            ebetk: ebetak, $
    ;;            eprlk: ebk, $
    ;;            ealp: ealpha, $
    ;;            ebet: ebeta, $
    ;;            eprl: eb, $
    ;;            b1_r: b_r, $
    ;;            b1_t: b_t, $
    ;;            b1_z: b_z, $
    ;;            jA_r: complex(jr_re,jr_im), $
    ;;            jA_t: complex(jt_re,jt_im), $
    ;;            jA_z: complex(jz_re,jz_im), $
    ;;            sig: sig, $
    ;;            kr : kr, $
    ;;            kz : kz }
    ;;
    ;;return, solution

end



