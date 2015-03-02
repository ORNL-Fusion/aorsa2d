
	freq = 500.0d6
	nphi = -56

	bField_eqdsk = 1
	eqdskFileName = 'g200201.00000.dlgMod'
	bTorFactor = 1.0

	flux_profiles = 1
    numeric_flux_profiles = 1

    FileName = 'densityT.txt'
    nRho = file_lines(FileName)
    openr, lun, FileName, /get_lun
    data = replicate( {rho: 0.0, n: 0.0, ni1: 0.0, ni2: 0.0, ni3: 0.0, ni4: 0.0} ,nRho)
    readf, lun, data    
    free_lun, lun

	atomicZ	= [-1,1,1,1,6]
	amu = [me/mi,2,1,2,12]
    nS = n_elements(amu)
    NumericData_n_m3 = fltArr(nRho,nS+1)
    NumericData_n_m3[*,0] = data[*].rho
    NumericData_n_m3[*,1] = data[*].n
    NumericData_n_m3[*,2] = data[*].ni1
    NumericData_n_m3[*,3] = data[*].ni2
    NumericData_n_m3[*,4] = data[*].ni3
    NumericData_n_m3[*,5] = data[*].ni4

    FileName = 'tempT.txt'
    nRho = file_lines(FileName)
    openr, lun, FileName, /get_lun
    data = replicate( {rho: 0.0, t: 0.0, ti1: 0.0, ti2: 0.0, ti3: 0.0, ti4: 0.0} ,nRho)
    readf, lun, data    
    free_lun, lun

    NumericData_T_eV = fltArr(nRho,nS+1)
    NumericData_T_eV[*,0] = data[*].rho
    NumericData_T_eV[*,1] = data[*].t*1e3
    NumericData_T_eV[*,2] = data[*].ti1 * 1e3
    NumericData_T_eV[*,3] = data[*].ti2 * 1e3
    NumericData_T_eV[*,4] = data[*].ti3 * 1e3
    NumericData_T_eV[*,5] = data[*].ti4 * 1e3

    nn = fltArr(4,nS) 
    tt = nn

	DensityMin = 0.2e17
	TempMin = 1e-3

	; Grid
	nR = 500
	nZ = 500
	rMin = 0.8
	rMax = 2.55
	zMin = -1.40
	zMax = +1.35

