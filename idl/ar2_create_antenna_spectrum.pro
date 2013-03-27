pro ar2_create_antenna_spectrum, $
		nPhiArray, $
		Coeffs = Coeffs

	; Create an antenna spectrum
	; --------------------------

	; NSTX
	; ----

	nStrap	= 12
	xlt		= 7.62e-2 ; strap width
	wd		= 19.40e-2 ; center to center separation
	phase_deg	= -30.0
	rAnt	= 1.57

	;; ITER
	;; ----

	;nStrap	= 4
	;xlt		= 17.0e-2
	;wd		= 38.5e-2
	;phase_deg	= -90.0
	;rAnt	= 8.165

	;; C-Mod
	;; -----

	;nStrap = 4
	;xlt = 7.5e-2
	;wd = 17.3e-2
	;phase_deg = [0.0,180.0,180.0,0.0]
	;rAnt = 0.9

	gap		= wd-xlt
	nPts	= 1024+1
	x		= fIndGen(nPts)/(nPts-1)*2*!pi*rAnt-!pi*rAnt
	amp		= fltArr ( nPts )
	phs		= fltArr ( nPts )
   
	offset	= 0;-(nStrap*(wd+xlt)-xlt)/2.0

	nstxOffset	= 0
	for i=0,nStrap-1 do begin

		minVal	= i*(gap+xlt)+offset
		maxVal	= minVal + gap

		if n_elements(phase_deg) gt 1 then begin
			phase = phase_deg[i]*!dtor
		endif else begin
			phase = phase_deg*i*!dtor
		endelse
		
		if i eq 6 and abs(phase) eq 180*!dtor then begin 
			;print, 'NSTX +/- 180 phase ? ... if not the check this'
			;print, (abs(phase) mod 2*!pi) * !radeg
			;nstxOffset = -180.0*!dtor
		endif

		phase = phase + nstxOffSet

		iiStrap	= where ( x gt minVal and x le maxVal, iiCnt )
		
		amp[iiStrap]	= 1
		phs[iiStrap]	= phase

	endfor	

	j	= amp * cos ( phs ) + complex(0,1) * sin ( phs )

	osf		= 10
	jPad	= complexArr ( nPts * osf )
	jPad[0:nPts-1]	= j

	dx	= abs(x[1]-x[0])
	k	= (fIndGen ( nPts*osf )-(nPts*osf+1)/2) / ( nPts * dx * osf ) * 2 * !pi
	nPhi_ant	= k * rAnt

	ant_spec	= fft(jPad, -1, /center)

	p = plot ( x, amp, layout = [2,2,1], title = 'Amplitude' )
	p = plot ( x, phs, layout = [2,2,2], title = 'Phase', /current )
	p = plot ( x, j, layout = [2,2,3], title = 'Real part', /current )
	p = plot ( x, imaginary(j), /current, layout = [2,2,4], title = 'Imag part')
	p = plot ( nPhi_ant, abs(ant_spec)^2, xrange=[-50,50], $
			title = 'Antenna Spectrum', xTitle = 'nPhi', yTitle = 'arb. units', $
			thick = 2)


	Coeffs = ComplexArr(n_elements(nPhiArray))
	for nphi=0,n_elements(nPhiArray)-1 do begin

		iiNPhi	= where ( nPhi_ant gt nPhi - 0.5 and nPhi_ant lt nPhi + 0.5, iiCnt )
		;Coeffs[nphi] = abs(mean ( ant_spec[iiNPhi] ))
		Coeffs[nphi] = mean ( ant_spec[iiNPhi] )

	endfor

	stop

end
