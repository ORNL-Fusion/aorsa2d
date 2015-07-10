pro ar2_iterate_right_left

    cd, current = WorkingDir
    leftName = '~/scratch/aorsa2d/ar2_vorpal/left_simple'
    rightName = '~/scratch/rsfwc_1d/ar2_vorpal/right_simple_'

    LeftFitLayer = [1.7,1.9]
    RightFitLayer = [1.9,2.1]

    nSubCycles = 3 

for MM = 0,10 do begin ; Cyclic MPE loop


for NN = 0, nSubCycles-1 do begin ; full iteration loop

    right_ThisDir =  RightName+StrCompress(string(NN),/rem)
    right_NextDir =  RightName+StrCompress(string(NN+1),/rem)
    right_PrevDir =  RightName+StrCompress(string(NN-1),/rem)

    cd, right_ThisDir
    if MM eq 0 or NN gt 0 then begin
            print, 'Running RSFWC-1D ...'
            print, right_ThisDir
            print, 'MM: ', MM
            print, 'NN: ', NN
            spawn, 'idl run_rsfwc.pro'
    endif

    DirExists = file_test(right_NextDir,/directory)
    if not DirExists then begin
        CreateNextDirectory = 'cp -r '+right_ThisDir+' '+right_NextDir
        print, 'Running ... '+CreateNextDirectory
        spawn, CreateNextDirectory 
    endif   
    cd, right_NextDir

    right_s = rsfwc_read_solution (right_ThisDir)
    ar2_fit_sources, FitThis = Right_s, WithTheseFiles = leftName, $
			FitLayer = RightFitLayer, $
			CoeffsOut_r = Coeffs_r, $
			CoeffsOut_t = Coeffs_t, $
			CoeffsOut_z = Coeffs_z

	NRHS = n_elements(Coeffs_r[*,0])
	nSpec = n_elements(Coeffs_r[0,*])

    This_s = ar2_read_solution (LeftName, 1)
    nL = n_elements(This_s.jP_r[*,0,0])

	This_jP_r = complexArr(nL) 
	This_jP_t = complexArr(nL)
	This_jP_z = complexArr(nL)

	This_jA_r = complexArr(nL) 
	This_jA_t = complexArr(nL)
	This_jA_z = complexArr(nL)

	This_E_r = complexArr(nL) 
	This_E_t = complexArr(nL)
	This_E_z = complexArr(nL)

	for rhs = 0, NRHS-1 do begin
    	This_s = ar2_read_solution (LeftName, rhs+1)
    	This_r = ar2_read_runData (LeftName, rhs+1)

        ; Remember to sum over the species :)
        basis_r = (total(This_s.jP_r,3))[*]
        basis_t = (total(This_s.jP_t,3))[*]
        basis_z = (total(This_s.jP_z,3))[*]

		This_jP_r[*] += Coeffs_r[rhs] * basis_r
		This_jP_t[*] += Coeffs_t[rhs] * basis_t
		This_jP_z[*] += Coeffs_z[rhs] * basis_z

        basis_r = (This_s.E_r)[*]
        basis_t = (This_s.E_t)[*]
        basis_z = (This_s.E_z)[*]

		This_E_r[*] += Coeffs_r[rhs] * basis_r
		This_E_t[*] += Coeffs_t[rhs] * basis_t
		This_E_z[*] += Coeffs_z[rhs] * basis_z

        basis_r = (This_r.jA_r)[*]
        basis_t = (This_r.jA_t)[*]
        basis_z = (This_r.jA_z)[*]

		This_jA_r[*] += Coeffs_r[rhs] * basis_r
		This_jA_t[*] += Coeffs_t[rhs] * basis_t
		This_jA_z[*] += Coeffs_z[rhs] * basis_z

	endfor

    ; Retart the iteration by scaling the magnitude of 
    ; what AORSA feeds back ...

    ScaleFac = 0.1
    This_jP_r = This_jP_r * ScaleFac 
    This_jP_t = This_jP_t * ScaleFac 
    This_jP_z = This_jP_z * ScaleFac 

    DoPlots = 0
    if DoPlots then begin
    p=plot(This_s.r,This_E_r,layout=[1,3,1], window_title='AORSA Fit E')
    p=plot(This_s.r,imaginary(This_E_r),/over,color='r')
    p=plot(This_s.r,This_E_t,layout=[1,3,2],/current)
    p=plot(This_s.r,imaginary(This_E_t),/over,color='r')
    p=plot(This_s.r,This_E_z,layout=[1,3,3],/current)
    p=plot(This_s.r,imaginary(This_E_z),/over,color='r')

    p=plot(This_s.r,This_jP_r,layout=[1,3,1], window_title='AORSA Fit jP')
    p=plot(This_s.r,imaginary(This_jP_r),/over,color='r')
    p=plot(This_s.r,This_jP_t,layout=[1,3,2],/current)
    p=plot(This_s.r,imaginary(This_jP_t),/over,color='r')
    p=plot(This_s.r,This_jP_z,layout=[1,3,3],/current)
    p=plot(This_s.r,imaginary(This_jP_z),/over,color='r')

    ;p=plot(This_s.r,This_jA_z,layout=[1,3,3],/current,color='b')
    ;p=plot(This_s.r,imaginary(This_jA_z),/over,color='b',linestyle='--')

    p=plot(right_s.r,right_s.jP_z,/over,thick=2)
    p=plot(right_s.r,imaginary(right_s.jP_z),/over,thick=2,color='r')

    stop
    endif


    ; Create netCdf file for RSFWC to read
    ; ------------------------------------

    r_rsfwc = right_s.r
    r_rsfwc_ = right_s.r_

    RSFWC_replace = intArr(n_elements(r_rsfwc))
    RSFWC_replace_ = intArr(n_elements(r_rsfwc_))

    iiRSFWC = where(r_rsfwc ge LeftFitLayer[0] and r_rsfwc le LeftFitLayer[1],iiCntRSFWC)
    iiRSFWC_ = where(r_rsfwc_ ge LeftFitLayer[0] and r_rsfwc_ le LeftFitLayer[1],iiCntRSFWC_)

    RSFWC_replace[iiRSFWC] = 1
    RSFWC_replace_[iiRSFWC_] = 1

    Jp_r_RSFWC = complexArr(n_elements(r_rsfwc))
    Jp_r_RSFWC_ = complexArr(n_elements(r_rsfwc_))
    Jp_t_RSFWC = complexArr(n_elements(r_rsfwc))
    Jp_t_RSFWC_ = complexArr(n_elements(r_rsfwc_))
    Jp_z_RSFWC = complexArr(n_elements(r_rsfwc))
    Jp_z_RSFWC_ = complexArr(n_elements(r_rsfwc_))

	jP_r_re =  interpol(real_part(This_jP_r),This_s.r, r_rsfwc,/spline)
	jP_r_re_ = interpol(real_part(This_jP_r),This_s.r,r_rsfwc_,/spline)
	jP_r_im =  interpol(imaginary(This_jP_r),This_s.r, r_rsfwc,/spline)
	jP_r_im_ = interpol(imaginary(This_jP_r),This_s.r,r_rsfwc_,/spline)

	jP_t_re =  interpol(real_part(This_jP_t),This_s.r, r_rsfwc,/spline)
	jP_t_re_ = interpol(real_part(This_jP_t),This_s.r,r_rsfwc_,/spline)
	jP_t_im =  interpol(imaginary(This_jP_t),This_s.r, r_rsfwc,/spline)
	jP_t_im_ = interpol(imaginary(This_jP_t),This_s.r,r_rsfwc_,/spline)

	jP_z_re =  interpol(real_part(This_jP_z),This_s.r, r_rsfwc,/spline)
	jP_z_re_ = interpol(real_part(This_jP_z),This_s.r,r_rsfwc_,/spline)
	jP_z_im =  interpol(imaginary(This_jP_z),This_s.r, r_rsfwc,/spline)
	jP_z_im_ = interpol(imaginary(This_jP_z),This_s.r,r_rsfwc_,/spline)

    jP_r_RSFWC  = complex(jP_r_re, jP_r_im)
    jP_r_RSFWC_ = complex(jP_r_re_,jP_r_im_)

    jP_t_RSFWC  = complex(jP_t_re, jP_t_im)
    jP_t_RSFWC_ = complex(jP_t_re_,jP_t_im_)

    jP_z_RSFWC  = complex(jP_z_re, jP_z_im)
    jP_z_RSFWC_ = complex(jP_z_re_,jP_z_im_)

    ncId = ncdf_create('aorsa_to_rsfwc.nc', /clobber)
    ncdf_control, ncId, /fill
    nrId = ncdf_dimdef(ncId, 'nR', n_elements(right_s.r))
    nrId_ = ncdf_dimdef(ncId, 'nR_', n_elements(right_s.r_))

    replace_id = ncdf_vardef(ncId, 'replace',nrId,/short)
    replace_id_ = ncdf_vardef(ncId, 'replace_',nrId_,/short)

    jP_r_re_id = ncdf_vardef(ncId,'jP_r_re',nrId,/float)
    jP_t_re_id = ncdf_vardef(ncId,'jP_p_re',nrId,/float)
    jP_z_re_id = ncdf_vardef(ncId,'jP_z_re',nrId,/float)
    jP_r_im_id = ncdf_vardef(ncId,'jP_r_im',nrId,/float)
    jP_t_im_id = ncdf_vardef(ncId,'jP_p_im',nrId,/float)
    jP_z_im_id = ncdf_vardef(ncId,'jP_z_im',nrId,/float)

    jP_r_re_id_ = ncdf_vardef(ncId,'jP_r_re_',nrId_,/float)
    jP_t_re_id_ = ncdf_vardef(ncId,'jP_p_re_',nrId_,/float)
    jP_z_re_id_ = ncdf_vardef(ncId,'jP_z_re_',nrId_,/float)
    jP_r_im_id_ = ncdf_vardef(ncId,'jP_r_im_',nrId_,/float)
    jP_t_im_id_ = ncdf_vardef(ncId,'jP_p_im_',nrId_,/float)
    jP_z_im_id_ = ncdf_vardef(ncId,'jP_z_im_',nrId_,/float)

    ncdf_control, ncId, /endef

    ncdf_varput, ncId, replace_id, RSFWC_replace 
    ncdf_varput, ncId, replace_id_, RSFWC_replace_

    ncdf_varput, ncId, jP_r_re_id, real_part(jP_r_RSFWC)
    ncdf_varput, ncId, jP_t_re_id, real_part(jP_t_RSFWC)
    ncdf_varput, ncId, jP_z_re_id, real_part(jP_z_RSFWC)
    ncdf_varput, ncId, jP_r_im_id, imaginary(jP_r_RSFWC)
    ncdf_varput, ncId, jP_t_im_id, imaginary(jP_t_RSFWC)
    ncdf_varput, ncId, jP_z_im_id, imaginary(jP_z_RSFWC)

    ncdf_varput, ncId, jP_r_re_id_, real_part(jP_r_RSFWC_)
    ncdf_varput, ncId, jP_t_re_id_, real_part(jP_t_RSFWC_)
    ncdf_varput, ncId, jP_z_re_id_, real_part(jP_z_RSFWC_)
    ncdf_varput, ncId, jP_r_im_id_, imaginary(jP_r_RSFWC_)
    ncdf_varput, ncId, jP_t_im_id_, imaginary(jP_t_RSFWC_)
    ncdf_varput, ncId, jP_z_im_id_, imaginary(jP_z_RSFWC_)

    ncdf_close, ncId

    ; Create ASCII file for Vorpal to read
    ; ------------------------------------

	VorpalFileName = 'aorsaToVorpal_E.txt'
	openw, lun, VorpalFileName, /get_lun

	printf, lun, 'nX: '+string(n_elements(this_s.r), format='(i4.4)')
	printf, lun, 'X, Re(Ex)[V/m], Im(Ex)[V/m], Re(Ey)[V/m], Im(Ey)[V/m], Re(Ez)[V/m], Im(Ez)[V/m]'

	this_E_x = this_E_r
	this_E_y = this_E_z
	this_E_z = -this_E_t

	for i=0,n_elements(this_s.r)-1 do begin
		printf, lun, this_s.r[i], $
			real_part(this_E_x), imaginary(this_E_x), $
			real_part(this_E_y), imaginary(this_E_y), $
			real_part(this_E_z), imaginary(this_E_z), $
			format='(7(f10.3,1x)'
	endfor
 	close, lun


endfor ; full iteration loop

cd, WorkingDir
MPE_FileName = 'restart_mpe_'+StrCompress(string(MM),/rem)+'.nc'

rsfwc_read_iterations, nSubCycles, MPE_FileName = MPE_FileName

Copy_MPE_FileIntoIteration0 = 'cp '+MPE_FileName+' '+RightName+'0/rsfwc_1d_r0.nc'
spawn, Copy_MPE_FileIntoIteration0 

endfor

end
