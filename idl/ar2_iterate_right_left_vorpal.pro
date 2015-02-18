pro ar2_iterate_right_left_vorpal, titan=titan

    cd, current = WorkingDir

	if keyword_set(titan)then begin
    	leftName = '/lustre/atlas/proj-shared/fus048/dg6/ar2_vorpal/left_simple'
    	rightName = '/lustre/atlas/proj-shared/fus048/dg6/vorpal/right_simple_'
	endif else begin
    	leftName = '/Users/dg6/scratch/aorsa2d/ar2_vorpal/left_simple'
    	rightName = '/Users/dg6/scratch/vorpal/right_simple_'
	endelse

    LeftFitLayer = [1.7,1.9]
    RightFitLayer = [1.9,2.1]

    nSubCycles = 2 

for MM = 0,0 do begin ; Cyclic MPE loop


for NN = 1, nSubCycles-1 do begin ; full iteration loop

    right_ThisDir =  RightName+StrCompress(string(NN),/rem)
    right_NextDir =  RightName+StrCompress(string(NN+1),/rem)
    right_PrevDir =  RightName+StrCompress(string(NN-1),/rem)

    cd, right_ThisDir
    if MM eq 0 or NN gt 0 then begin
            print, 'Running VORPAL ...'
            print, right_ThisDir
            print, 'MM: ', MM
            print, 'NN: ', NN
			
			if keyword_set(titan)then begin
            	spawn, './vorpal-titan-interactive.pbs'
			endif else begin
            	;spawn, './vorpal-osx.sh'
			endelse
    endif

    DirExists = file_test(right_NextDir,/directory)
    if not DirExists then begin
        CreateNextDirectory = 'cp -r '+right_ThisDir+' '+right_NextDir
        print, 'Running ... '+CreateNextDirectory
        spawn, CreateNextDirectory 
    endif   
    cd, right_NextDir

    ;right_s = rsfwc_read_solution (right_ThisDir)
	freq = 53.0e6

    ; Use the MPE solution instead of an initial Vorpal run
    ; for the first step in non-zero MPE cycles
    if MM gt 0 and NN eq 0 then begin
        right_s = ar2_read_mpe (right_ThisDir+'/restart_mpe.nc') 
        stop
    endif else begin
        right_s = ar2_read_vorpal (runFolderName = right_ThisDir, freq = freq )
    endelse

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

    ;ScaleFac = 0.1
    ;This_jP_r = This_jP_r * ScaleFac 
    ;This_jP_t = This_jP_t * ScaleFac 
    ;This_jP_z = This_jP_z * ScaleFac 

    DoPlots = 1
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

    ;p=plot(right_s.r,right_s.jP_z,/over,thick=2)
    ;p=plot(right_s.r,imaginary(right_s.jP_z),/over,thick=2,color='r')

    R = right_s
    L = This_s

    iiV = where(R.r gt RightFitLayer[0]) 
    iiA = where(L.r lt RightFitLayer[1])

    FullSolutionFile = expand_path('~/scratch/vorpal/full_simple/vorpal.sav')
    restore, FullSolutionFile

    p=plot(R.r[iiV],R.e_R[iiV],color='b',thick=2,layout=[1,3,1],title='L/R E field')
    p=plot(R.r[iiV],imaginary(R.e_R[iiV]),color='b',thick=1,lineStyle='-',/over)
    p=plot(L.r[iiA],This_e_R[iiA],/over,color='g',thick=2)
    p=plot(L.r[iiA],imaginary(This_e_R[iiA]),color='g',thick=1,lineStyle='-',/over)

    p=plot(v.r,v.e_r,/over)
    p=plot(v.r,imaginary(v.e_r),/over,lineStyle='-')

    p=plot(R.r[iiV],R.e_T[iiV],color='b',thick=2,layout=[1,3,2],/current)
    p=plot(R.r[iiV],imaginary(R.e_T[iiV]),color='b',thick=1,lineStyle='-',/over)
    p=plot(L.r[iiA],This_e_T[iiA],/over,color='g',thick=2)
    p=plot(L.r[iiA],imaginary(This_e_T[iiA]),color='g',thick=1,lineStyle='-',/over)

    p=plot(v.r,v.e_t,/over)
    p=plot(v.r,imaginary(v.e_t),/over,lineStyle='-')

    p=plot(R.r[iiV],R.e_Z[iiV],color='b',thick=2,layout=[1,3,3],/current)
    p=plot(R.r[iiV],imaginary(R.e_Z[iiV]),color='b',thick=1,lineStyle='-',/over)
    p=plot(L.r[iiA],This_e_Z[iiA],/over,color='g',thick=2)
    p=plot(L.r[iiA],imaginary(This_e_Z[iiA]),color='g',thick=1,lineStyle='-',/over)

    p=plot(v.r,v.e_z,/over)
    p=plot(v.r,imaginary(v.e_z),/over,lineStyle='-')


    p=plot(R.r[iiV],R.jP_R[iiV],color='b',thick=2,layout=[1,3,1], title='L/R jP')
    p=plot(R.r[iiV],imaginary(R.jP_R[iiV]),color='b',thick=1,lineStyle='-',/over)
    p=plot(L.r[iiA],This_jP_R[iiA],/over,color='g',thick=2)
    p=plot(L.r[iiA],imaginary(This_jP_R[iiA]),color='g',thick=1,lineStyle='-',/over)

    p=plot(R.r[iiV],R.jP_T[iiV],color='b',thick=2,layout=[1,3,2], /current)
    p=plot(R.r[iiV],imaginary(R.jP_T[iiV]),color='b',thick=1,lineStyle='-',/over)
    p=plot(L.r[iiA],This_jP_T[iiA],/over,color='g',thick=2)
    p=plot(L.r[iiA],imaginary(This_jP_T[iiA]),color='g',thick=1,lineStyle='-',/over)

    p=plot(R.r[iiV],R.jP_Z[iiV],color='b',thick=2,layout=[1,3,3], /current)
    p=plot(R.r[iiV],imaginary(R.jP_Z[iiV]),color='b',thick=1,lineStyle='-',/over)
    p=plot(L.r[iiA],This_jP_Z[iiA],/over,color='g',thick=2)
    p=plot(L.r[iiA],imaginary(This_jP_Z[iiA]),color='g',thick=1,lineStyle='-',/over)




stop
    endif


    ; Create ASCII file for Vorpal to read
    ; ------------------------------------

	VorpalFileName = 'aorsaToVorpal_E.txt'
	openw, lun, VorpalFileName, /get_lun

	printf, lun, 'nX: '+string(n_elements(this_s.r), format='(i4.4)')
	printf, lun, 'X, Re(Ex)[V/m], Im(Ex)[V/m], Re(Ey)[V/m], Im(Ey)[V/m], Re(Ez)[V/m], Im(Ez)[V/m]'

	this_E_x = this_Jp_r
	this_E_y = this_Jp_t
	this_E_z = this_Jp_z

	for i=0,n_elements(this_s.r)-1 do begin
		printf, lun, this_s.r[i], $
			real_part(this_E_x[i]), imaginary(this_E_x[i]), $
			real_part(this_E_y[i]), imaginary(this_E_y[i]), $
			real_part(this_E_z[i]), imaginary(this_E_z[i]), $
			format='(7(f13.6,1x))'
	endfor
 	close, lun

endfor ; full iteration loop

cd, WorkingDir
MPE_FileName = 'restart_mpe_'+StrCompress(string(MM),/rem)+'.nc'

rsfwc_read_iterations, nSubCycles, RightName, $
		MPE_FileName = MPE_FileName, /vorpal, freq = freq

Copy_MPE_FileIntoIteration0 = 'cp '+MPE_FileName+' '+RightName+'0/rsfwc_1d_r0.nc'
Copy_MPE_FileIntoIteration0 = 'cp '+MPE_FileName+' '+RightName+'0/restart_mpe.nc'

spawn, Copy_MPE_FileIntoIteration0 

endfor

end
