
	freq = 30.0d6
	nphi = -21

	bField_eqdsk = 1
	eqdskFileName = 'NSTX_130608A03_trxpl_0.410_plasma_state.geq'
	bTorFactor = 1.0

	flux_profiles = 1
   	numeric_flux_profiles = 1
	fred_namelist_input = 1

	DensityMin = [2e18,2e18,5e9,1e16,5e9,6e14]*0.25 ; closed
	;DensityMin = [2e18,2e18,5e9,1e16,5e9,6e14]*2 ; med open
	;DensityMin = [2e18,2e18,5e9,1e16,5e9,6e14]*4 ; full open

	TempMin = 1e-3

	; Grid
	nR = 300
	nZ = 300
	rMin = 0.05
	rMax = 1.7
	zMin = -1.80
	zMax = +1.80

