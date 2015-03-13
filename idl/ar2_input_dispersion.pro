pro ar2_input_dispersion, wrf, amu, AtomicZ, nn, nPhi, nSpec, nR, nZ, Density_m3, bMag, $
		r2D, resonances = resonances, $
		IonIonHybrid_res_freq=IonIonHybrid_res_freq, $
		IonIonHybrid_cut_freq=IonIonHybrid_cut_freq, $
		Spec1=Spec1,Spec2=Spec2, $
		kPerSq_F=kPerSq_F, kPerSq_S=kPerSq_S, $
		StixP=StixP,StixL=StixL,StixR=StixR,StixS=StixS

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
				resonances[i,j,s] = (wc[i,j] mod wrf)/wrf
			endfor
		endfor

		StixP	= StixP + ( 1d0 - ( wp^2/wrf^2 ) )
		StixL	= StixL + ( 1d0 - ( wp^2 / ( wrf * ( wrf - wc ) ) ) )
		StixR	= StixR + ( 1d0 - ( wp^2 / ( wrf * ( wrf + wc ) ) ) )

	endfor

	if keyword_set(Spec1) and keyword_set(Spec2) and nSpec gt 2 then begin
		;print, 'Calculating Ion-Ion Hybrid Resonance Freq.'
		; Ion-Ion Hybrid Freq (pg. 248 Brambilla)
		nuSpec1 	= Density_m3[*,*,Spec1] / Density_m3[*,*,0]
		nuSpec2 	= Density_m3[*,*,Spec2] / Density_m3[*,*,0]
		IonIonHybrid_res = $
				(nuSpec1*AtomicZ[Spec1]*(AtomicZ[Spec2]/amu[Spec2])+nuSpec2*AtomicZ[Spec2]*(AtomicZ[Spec1]/amu[Spec1])) $
				/ (nuSpec1*AtomicZ[Spec1]*(amu[Spec2]/AtomicZ[Spec2])+nuSpec2*AtomicZ[Spec2]*(amu[Spec1]/AtomicZ[Spec1]))
		IonIonHybrid_cut = nuSpec1*AtomicZ[Spec1]*(AtomicZ[Spec2]/amu[Spec2])+nuSpec2*AtomicZ[Spec2]*(AtomicZ[Spec1]/amu[Spec1])
		wc_H2 = (e*bMag/mi)^2
		IonIonHybrid_res_freq = sqrt(IonIonHybrid_res * wc_H2)
		IonIonHybrid_cut_freq = IonIonHybrid_cut * sqrt(wc_H2)
	endif

	StixS	= 0.5d0 * ( StixR + StixL )

	kPar2D	= nPhi / r2D
	nPar2D	= kPar2D * c / wrf

	kPerSq_F    = complex(-(nPar2D^2 - stixR)*(nPar2D^2 - stixL) / (nPar2D^2 - stixS ) * wrf^2/c^2,nPar2D*0)
	kPerSq_S    = complex(-(nPar2D^2 - stixS) * stixP / stixS * wrf^2/c^2,nPar2D*0)

	zSlice = nZ/2
	kPerSq_1 = kPerSq_F[*,zSlice]
	kPerSq_2 = kPerSq_S[*,zSlice]

	;p=plot(r2D[*,zSlice],sqrt(kPerSq_F[*,zSlice]),color='b',/ylog)
	;p=plot(r2D[*,zSlice],sqrt(kPerSq_S[*,zSlice]),/over,color='r',/ylog)

end
