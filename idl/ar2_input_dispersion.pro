pro ar2_input_dispersion, wrf, amu, AtomicZ, nn, nSpec, nR, nZ, Density_m3, bMag, $
		r2D, nuOmg, kPar2D, kPerSq_F=kPerSq_F, kPerSq_S=kPerSq_S

	@dlg_constants

	StixP = 1d0
	StixL = 1d0
	StixR = 1d0

	for s=0,nSpec-1 do begin

        mass = amu[s] * _amu * dComplex(1,0);nuOmg[*,*,s])

		;	these are in si not cgs units (i think)
		
		wp	= sqrt ( Density_m3[*,*,s] * (AtomicZ[s]*_e)^2 / ( mass * _e0 ) )
		wc	=  AtomicZ[s]*_e * bMag / ( mass )

		StixP	= StixP + ( 1d0 - ( wp^2/wrf^2 ) )
		StixL	= StixL + ( 1d0 - ( wp^2 / ( wrf * ( wrf - wc ) ) ) )
		StixR	= StixR + ( 1d0 - ( wp^2 / ( wrf * ( wrf + wc ) ) ) )

	endfor

	StixS	= 0.5d0 * ( StixR + StixL )

	nPar2D	= kPar2D * _c / wrf

	kPerSq_F    = -(nPar2D^2 - stixR)*(nPar2D^2 - stixL) / (nPar2D^2 - stixS ) * wrf^2/_c^2
	kPerSq_S    = -(nPar2D^2 - stixS) * stixP / stixS * wrf^2/_c^2

end
