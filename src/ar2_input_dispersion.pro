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

		StixP	= StixP - ( wp^2/wrf^2 )
		StixL	= StixL - ( wp^2 / ( wrf * ( wrf - wc ) ) )
		StixR	= StixR - ( wp^2 / ( wrf * ( wrf + wc ) ) )

	endfor

	StixS	= 0.5d0 * ( StixR + StixL )
	StixD	= 0.5d0 * ( StixR - StixL )
	StixA	= StixS

	; without upshift

	kPar2D	= nPhi / r2D
	nPar2D	= kPar2D * c / wrf

	StixB	= -1.0 * ( StixR * StixL + StixP * StixS - nPar2D^2 * ( StixP + StixS ) )
	StixC	= StixP * ( nPar2D^2 - StixR ) * ( nPar2D^2 - StixL )

	b24ac_	= Stixb^2 - 4.0 * Stixa * StixC

	kPerpSq_2D_1 =  ( -StixB + sqrt ( complex( b24ac_, b24ac_*0) ) ) / ( 2d0 * StixA )  * wrf / c
	kPerpSq_2D_2 =  ( -StixB - sqrt ( complex( b24ac_, b24ac_*0) ) ) / ( 2d0 * StixA )  * wrf / c

	zSlice = nZ/2
	kPerSq_1 = kPerpSq_2D_1[*,zSlice]
	kPerSq_2 = kPerpSq_2D_1[*,zSlice]

	p=plot(r2D[*,zSlice],r2D[*,zSlice],/noData,yRange=[0,1500])

	kPer = complexArr(nR,4)
   	kPer[*,0] = +sqrt(kPerSq_1)
	kPer[*,1] = -sqrt(kPerSq_1)	
   	kPer[*,2] = +sqrt(kPerSq_2)
	kPer[*,3] = -sqrt(kPerSq_2)	

	kRPlot	= real_part ( kPer )
	iiNeg	= where ( kRPlot lt 0, iiNegCnt )
	if iiNegCnt gt 0 then kRPlot[iiNeg] = -1.0 * kRPlot[iiNeg] 

	dlg_rootPlotter, r2D[*,zSlice], kRPlot

	kRPlot	= imaginary ( kPer )
	iiNeg	= where ( kRPlot lt 0, iiNegCnt )
	if iiNegCnt gt 0 then kRPlot[iiNeg]	= -1.0 * kRPlot[iiNeg] 

	dlg_rootPlotter, r2D[*,zSlice], kRPlot, /imag

end
