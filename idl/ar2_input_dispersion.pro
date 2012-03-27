pro ar2_input_dispersion, wrf, amu, AtomicZ, nn, nPhi, nSpec, nR, nZ, Density_m3, bMag, $
		r2D

	@constants

	StixP = 1d0
	StixL = 1d0
	StixR = 1d0

	resonances = fltarr(nR,nZ,nSpec)

	for s=0,nSpec-1 do begin

		;	calculate cutoffs etc

		;	these are in si not cgs units (i think)
		
		wp	= sqrt ( Density_m3[*,*,s] * (AtomicZ[s]*e)^2 / ( amu[s]*mi * e0 ) )
		wc	=  AtomicZ[s]*e * bMag / ( amu[s]*mi )

		for i=0,nR-1 do begin
			for j=0,nZ-1 do begin
				resonances[i,j,s] = 1.0 / ( wc[i,j] mod wrf )
			endfor
		endfor

		StixP	= StixP + ( 1d0 - ( wp^2/wrf^2 ) )
		StixL	= StixL + ( 1d0 - ( wp^2 / ( wrf * ( wrf - wc ) ) ) )
		StixR	= StixR + ( 1d0 - ( wp^2 / ( wrf * ( wrf + wc ) ) ) )

	endfor

	StixS	= 0.5d0 * ( StixR + StixL )

	kPar2D	= nPhi / r2D
	nPar2D	= kPar2D * c / wrf

	kPerSq_F    = complex(-(nPar2D^2 - stixR)*(nPar2D^2 - stixL) / (nPar2D^2 - stixS ) * wrf^2/c^2,nPar2D*0)
	kPerSq_S    = complex(-(nPar2D^2 - stixS) * stixP / stixS * wrf^2/c^2,nPar2D*0)

	zSlice = nZ/2
	kPerSq_1 = kPerSq_F[*,zSlice]
	kPerSq_2 = kPerSq_S[*,zSlice]

	p=plot(r2D[*,zSlice],sqrt(kPerSq_F[*,zSlice]),color='b')
	p=plot(r2D[*,zSlice],sqrt(kPerSq_S[*,zSlice]),/over,color='r')

end
