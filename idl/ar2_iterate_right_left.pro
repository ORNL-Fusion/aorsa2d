pro ar2_iterate_right_left

    leftName = 'left_simple'
    rightName = 'right_simple'

    LeftFitLayer = [1.85,1.90]
    RightFitLayer = [2.00,2.05]

    driverSideName = rightName

    ; Get the initial driver side jP

	driverSideFiles = file_search(driverSideName+'/output/solution*.nc')
    NRHS_driverSide = n_elements(driverSideFiles)

    ; Assume the last element of the driver side is the standard
    ; AORSA driver (single component), and the rest are the unit response
    ; drivers.

    print, 'Driver RHS number: ', NRHS_driverSide
    s = ar2_read_solution (driverSideName, NRHS_driverSide)

    ar2_fit_sources, FitThis = s, WithTheseFiles = leftName, FitLayer = LeftFitLayer

    stop

end
