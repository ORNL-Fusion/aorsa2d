pro ar2_iterate_right_left

    leftName = 'left_simple'
    rightName = 'right_simple'

    LeftFitLayer = [1.85,1.9]
    RightFitLayer = [2.0,2.05]

    driverSideName = rightName

    ; Get the initial driver side jP

	driverSideFiles = file_search(driverSideName+'/output/solution*.nc')
    NRHS_driverSide = n_elements(driverSideFiles)

    ; Assume the last element of the driver side is the standard
    ; AORSA driver (single component), and the rest are the unit response
    ; drivers.

    print, 'Driver RHS number: ', NRHS_driverSide
    Driver_s = ar2_read_solution (driverSideName, NRHS_driverSide)

	IteratedRight_s = Driver_s

for NN = 0, 4 do begin ; full iteration loop

    ar2_fit_sources, FitThis = IteratedRight_s, WithTheseFiles = leftName, $
			FitLayer = RightFitLayer, $
			CoeffsOut_r = Coeffs_r, $
			CoeffsOut_t = Coeffs_t, $
			CoeffsOut_z = Coeffs_z

	NRHS = n_elements(Coeffs_r[*,0])
	nSpec = n_elements(Coeffs_r[0,*])

    This_s = ar2_read_solution (LeftName, 1)

	This_jP_r = This_s.jP_r*0
	This_jP_t = This_s.jP_t*0
	This_jP_z = This_s.jP_z*0

	for ss = 0, nSpec-1 do begin
		for rhs = 0, NRHS-1 do begin
    		This_s = ar2_read_solution (LeftName, rhs+1)
			This_jP_r[*,*,ss] += Coeffs_r[rhs,ss] * This_s.jP_r[*,*,ss]
			This_jP_t[*,*,ss] += Coeffs_t[rhs,ss] * This_s.jP_t[*,*,ss]
			This_jP_z[*,*,ss] += Coeffs_z[rhs,ss] * This_s.jP_z[*,*,ss]
		endfor
	endfor

    IteratedLeft_s = ar2_read_solution (LeftName, 1)
	IteratedLeft_s.jP_r = This_jP_r
	IteratedLeft_s.jP_t = This_jP_t
	IteratedLeft_s.jP_z = This_jP_z

	iiLeftFit = where(This_s.r lt 999);RightFitLayer[1])
	p=plot(This_s.r,real_part(This_jP_r[iiLeftFit,0,0]),layout=[1,3,1],title="Iteration: "+string(NN,format='(i2.2)'))
	p=plot(This_s.r,imaginary(This_jP_r[iiLeftFit,0,0]),/over,color="red")
	p=plot(Driver_s.r,real_part(Driver_s.jP_r[*,0,0]),/over,thick=3)
	p=plot(Driver_s.r,imaginary(Driver_s.jP_r[*,0,0]),/over,thick=3,color="red")

	p=plot(This_s.r,real_part(This_jP_t[iiLeftFit,0,0]),layout=[1,3,2],/current)
	p=plot(This_s.r,imaginary(This_jP_t[iiLeftFit,0,0]),/over,color="red")
	p=plot(Driver_s.r,real_part(Driver_s.jP_t[*,0,0]),/over,thick=3)
	p=plot(Driver_s.r,imaginary(Driver_s.jP_t[*,0,0]),/over,thick=3,color="red")

    p=plot(This_s.r,real_part(This_jP_z[iiLeftFit,0,0]),layout=[1,3,3],/current)
	p=plot(This_s.r,imaginary(This_jP_z[iiLeftFit,0,0]),/over,color="red")
	p=plot(Driver_s.r,real_part(Driver_s.jP_z[*,0,0]),/over,thick=3)
	p=plot(Driver_s.r,imaginary(Driver_s.jP_z[*,0,0]),/over,thick=3,color="red")

    ar2_fit_sources, FitThis = IteratedLeft_s, WithTheseFiles = RightName, $
			FitLayer = LeftFitLayer, $
			CoeffsOut_r = Coeffs_r, $
			CoeffsOut_t = Coeffs_t, $
			CoeffsOut_z = Coeffs_z, $
			NRHS = NRHS

	NRHS = n_elements(Coeffs_r[*,0])
	nSpec = n_elements(Coeffs_r[0,*])

    This_s = ar2_read_solution (RightName, 1)

	This_jP_r = This_s.jP_r*0
	This_jP_t = This_s.jP_t*0
	This_jP_z = This_s.jP_z*0

	for ss = 0, nSpec-1 do begin
		for rhs = 0, NRHS-1 do begin
    		This_s = ar2_read_solution (RightName, rhs+1)
			This_jP_r[*,*,ss] += Coeffs_r[rhs,ss] * This_s.jP_r[*,*,ss]
			This_jP_t[*,*,ss] += Coeffs_t[rhs,ss] * This_s.jP_t[*,*,ss]
			This_jP_z[*,*,ss] += Coeffs_z[rhs,ss] * This_s.jP_z[*,*,ss]
		endfor
	endfor

	IteratedRight_s.jP_r = This_jP_r
	IteratedRight_s.jP_t = This_jP_t
	IteratedRight_s.jP_z = This_jP_z

	iiRightFit = where(This_s.r gt -999);LeftFitLayer[0])
	p=plot(This_s.r,real_part(This_jP_r[iiRightFit,0,0]),layout=[1,3,1])
	p=plot(This_s.r,imaginary(This_jP_r[iiRightFit,0,0]),/over,color="red")
	p=plot(IteratedLeft_s.r,real_part(IteratedLeft_s.jP_r[iiLeftFit,0,0]),/over,thick=3)
	p=plot(IteratedLeft_s.r,imaginary(IteratedLeft_s.jP_r[iiLeftFit,0,0]),/over,thick=3,color="red")

	p=plot(This_s.r,real_part(This_jP_t[iiRightFit,0,0]),layout=[1,3,2],/current)
	p=plot(This_s.r,imaginary(This_jP_t[iiRightFit,0,0]),/over,color="red")
	p=plot(IteratedLeft_s.r,real_part(IteratedLeft_s.jP_t[iiLeftFit,0,0]),/over,thick=3)
	p=plot(IteratedLeft_s.r,imaginary(IteratedLeft_s.jP_t[iiLeftFit,0,0]),/over,thick=3,color="red")

    p=plot(This_s.r,real_part(This_jP_z[iiRightFit,0,0]),layout=[1,3,3],/current)
	p=plot(This_s.r,imaginary(This_jP_z[iiRightFit,0,0]),/over,color="red")
	p=plot(IteratedLeft_s.r,real_part(IteratedLeft_s.jP_z[iiLeftFit,0,0]),/over,thick=3)
	p=plot(IteratedLeft_s.r,imaginary(IteratedLeft_s.jP_z[iiLeftFit,0,0]),/over,thick=3,color="red")
    
    stop

endfor ; full iteration loop

end
