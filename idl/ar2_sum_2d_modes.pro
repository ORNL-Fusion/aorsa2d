pro ar2_sum_2d_modes

	fileList = file_search ( 'output/plotData*.nc' )
	fileListJp	= file_search ( 'output/jpStix*.nc' )

	nFiles	= n_elements ( fileList )

	cdfId = ncdf_open ( fileList[0], /noWrite ) 

		ncdf_varget, cdfId, 'capR', R 
		ncdf_varget, cdfId, 'zLoc', z 

	ncdf_close, cdfId

	nR	= n_elements ( R )
	nz	= n_elements ( z )

	eR2D_all	= complexArr ( nR, nz, nFiles )
	eT2D_all	= complexArr ( nR, nz, nFiles )
	eZ2D_all	= complexArr ( nR, nz, nFiles )

	nPhi_all	= intArr ( nFiles )
	;nPhi_ant_all	= complexArr ( nFiles )

	for f=0,nFiles-1 do begin

		cdfId = ncdf_open ( fileList[f], /noWrite ) 

			ncdf_varget, cdfId, 'nPhi', nPhi
			ncdf_varGet, cdfId, 'wdote_rz', wdote_rz
			ncdf_varget, cdfId, 'wdoti1_rz', wdoti1_rz
			ncdf_varget, cdfId, 'wdoti2_rz', wdoti2_rz
			nCdf_varGet, cdfId, 'pscale', pscale
			nCdf_varGet, cdfId, 'eAlpha_real', eAlpha_real
			nCdf_varGet, cdfId, 'eAlpha_imag', eAlpha_imag
			nCdf_varGet, cdfId, 'eBeta_real', eBeta_real
			nCdf_varGet, cdfId, 'eBeta_imag', eBeta_imag
			nCdf_varGet, cdfId, 'eb_real', eb_real
			nCdf_varGet, cdfId, 'eb_imag', eb_imag
			nCdf_varGet, cdfId, 'jeDotE', jeDotE
			nCdf_varGet, cdfId, 'jtDotE', jtDotE
			nCdf_varGet, cdfId, 'j1DotE', j1DotE

			nCdf_varGet, cdfId, 'janty', jAntY

			nCdf_varGet, cdfId, 'xkperp_real', xkperp_real
			nCdf_varGet, cdfId, 'xkperp_imag', xkperp_imag
			nCdf_varGet, cdfId, 'pscale', pscale
			print, 'pscale :', pscale, nPhi

		ncdf_close, cdfId

		cdfId = ncdf_open ( fileListJp[f], /noWrite ) 

			nCdf_varGet, cdfId, 'jpAlpha_re', jPAlpha_re
			nCdf_varGet, cdfId, 'jpAlpha_im', jPAlpha_im
			nCdf_varGet, cdfId, 'jpBeta_re', jPBeta_re
			nCdf_varGet, cdfId, 'jpBeta_im', jPBeta_im
			nCdf_varGet, cdfId, 'jpB_re', jPB_re
			nCdf_varGet, cdfId, 'jpB_im', jPB_im

			nCdf_varGet, cdfId, 'ex_re', eR_re
			nCdf_varGet, cdfId, 'ex_im', eR_im
			nCdf_varGet, cdfId, 'ey_re', eZ_re
			nCdf_varGet, cdfId, 'ey_im', eZ_im
			nCdf_varGet, cdfId, 'ez_re', eT_re
			nCdf_varGet, cdfId, 'ez_im', eT_im

			nCdf_varGet, cdfId, 'bx_re', bR_re
			nCdf_varGet, cdfId, 'bx_im', bR_im
			nCdf_varGet, cdfId, 'by_re', bZ_re
			nCdf_varGet, cdfId, 'by_im', bZ_im
			nCdf_varGet, cdfId, 'bz_re', bT_re
			nCdf_varGet, cdfId, 'bz_im', bT_im

		ncdf_close, cdfId


		nPhi_all[f]	= nPhi

		eAlpha2D_all[*,*,f]	= complex ( eAlpha_real, eAlpha_imag ) / sqrt(pscale)
		eBeta2D_all[*,*,f]	= complex ( eBeta_real, eBeta_imag ) / sqrt(pscale)
		eB2D_all[*,*,f]	= complex ( eB_real, eB_imag ) / sqrt(pscale)

		jPAlpha2D_all[*,*,f]	= complex ( jPAlpha_re, jPAlpha_im ) / sqrt(pscale)
		jPBeta2D_all[*,*,f]	= complex ( jPBeta_re, jPBeta_im ) / sqrt(pscale)
		jPB2D_all[*,*,f]	= complex ( jPB_re, jPB_im ) / sqrt(pscale)

		eR2D_all[*,*,f]	= complex ( eR_re, eR_im ) / sqrt(pscale)
		eT2D_all[*,*,f]	= complex ( eT_re, eT_im ) / sqrt(pscale)
		eZ2D_all[*,*,f]	= complex ( eZ_re, eZ_im ) / sqrt(pscale)

		bR2D_all[*,*,f]	= complex ( bR_re, bR_im ) / sqrt(pscale)
		bT2D_all[*,*,f]	= complex ( bT_re, bT_im ) / sqrt(pscale)
		bZ2D_all[*,*,f]	= complex ( bZ_re, bZ_im ) / sqrt(pscale)

		kPerp_all[*,*,f]	= complex ( xkPerp_real, xkPerp_imag )

		jAntY	= complex ( jAntY*0, jAntY )

		wdote_all[*,*,f]	= wdote_rz / pscale
		wdoti1_all[*,*,f]	= wdoti1_rz / pscale

		;if nPhi eq 0 then eAlpha2D_all[*,*,f] = 0

		;; Get antenna scaling for this nPhi

		;iiNPhi	= where ( nPhi_ant gt nPhi - 0.5 and nPhi_ant lt nPhi + 0.5, iiCnt )
		;nPhi_ant_all[f]	= abs(mean ( ant_spec[iiNPhi] ))

	endfor

	ar2_create_antenna_spectrum ( nPhiArray = nPhi_all, Coeffs = nPhi_ant_all )

	wdot_all	= wdote_all + wdoti1_all


	; Sort according to nPhi
	; ----------------------

	iiSort	= sort ( nPhi_all )
	nPhi_all		= nPhi_all[iiSort]

	eAlpha2D_all	= eAlpha2D_all[*,*,iiSort]
	eBeta2D_all		= eBeta2D_all[*,*,iiSort]
	eB2D_all		= eB2D_all[*,*,iiSort]

	wdot_all		= wdot_all[*,*,iiSort]
	nPhi_ant_all	= nPhi_ant_all[iiSort]

	jPAlpha2D_all	= jPAlpha2D_all[*,*,iiSort]
	jPBeta2D_all	= jPBeta2D_all[*,*,iiSort]
	jPB2D_all		= jPB2D_all[*,*,iiSort]

	eR2D_all		= eR2D_all[*,*,iiSort]
	eT2D_all		= eT2D_all[*,*,iiSort]
	eZ2D_all		= eZ2D_all[*,*,iiSort]

	bR2D_all		= bR2D_all[*,*,iiSort]
	bT2D_all		= bT2D_all[*,*,iiSort]
	bZ2D_all		= bZ2D_all[*,*,iiSort]
	

	;; Median filter the spectrum to remove cavity mode spikes
	;; -------------------------------------------------------
	;width	= 6

	;absE	= sqrt ( abs(eAlpha2d_all)^2 + abs(eBeta2d_all)^2 + abs(eB2d_all)^2 )
	;meanE	= mean(mean(absE,dim=1),dim=1)
	;medianE	= median(meanE,width,dim=1,/even)

	;meanWDot	= mean(mean(wdot_all,dim=1),dim=1)
	;medianWDot	= median(meanWDot,width,dim=1,/even)

	;for f=0,nFiles-1 do begin

	;	eAlpha2D_all[*,*,f]	= eAlpha2D_all[*,*,f] / meanE[f] * medianE[f]
	;	eBeta2D_all[*,*,f]	= eBeta2D_all[*,*,f] / meanE[f] * medianE[f]
	;	eB2D_all[*,*,f]	= eB2D_all[*,*,f] / meanE[f] * medianE[f]

	;	wdot_all[*,*,f]	= wdot_all[*,*,f] / meanWDot[f] * medianWDot[f]

	;endfor


	; E*.J power 
	; ----------

	eDotJ	= 0.5 * real_part ( $
				conj ( eAlpha2D_all ) * jPAlpha2D_all + $
				conj ( eBeta2D_all ) * jPBeta2D_all + $
				conj ( eB2D_all ) * jPB2D_all )

	p = plot ( total ( total ( eDotJ, 1 ), 1 ), $
			title = 'eDotJ spectrum' )


	;; Apply linear -ve or +ve spectral weighting
	;; ------------------------------------------

	;totalFac	= total ( nPhi_ant_all )
	;slope	= -2
	;weight	= (fIndGen(nFiles)+1)*abs(slope)
	;if slope lt 0 then weight = reverse(weight)
	;
	;nPhi_ant_all_orig	= nPhi_ant_all
	;nPhi_ant_all	= nPhi_ant_all * weight
	;nPhi_ant_all	= nPhi_ant_all / total ( nPhi_ant_all ) * totalFac


	meanModeSpectrumAlpha	= mean(mean(abs(eAlpha2d_all),dim=1),dim=1)
	meanModeSpectrumBeta	= mean(mean(abs(eBeta2d_all),dim=1),dim=1)
	meanModeSpectrumB		= mean(mean(abs(eB2d_all),dim=1),dim=1)


	meanModeSpectrumWdot	= mean(mean(abs(wdot_all),dim=1),dim=1)

	p = plot ( nPhi_all, meanModeSpectrumAlpha, title = 'Mean mode spectrum' )
	p = plot ( nPhi_all, meanModeSpectrumBeta, /over )
	p = plot ( nPhi_all, meanModeSpectrumB, /over )

	p = plot ( nPhi_all, nPhi_ant_all, thick = 2, title = 'Antenna Spectrum', $
		   xTitle = 'nPhi', yTitle = 'arb. units')

   	; Smooth out the spectrum a bit

	meanModeSpectrum = mean(mean(sqrt(abs(eAlpha2d_all)^2+abs(eBeta2d_all)^2+abs(eB2d_all)^2),dim=1),dim=1)

	for i=0,n_elements(meanModeSpectrum)-1 do begin

		sf = sqrt(meanModeSpectrum[i])/meanModeSpectrum[i]

		eAlpha2D_all[*,*,i]		= eAlpha2D_all[*,*,i]*sf
		eBeta2D_all[*,*,i]		= eBeta2D_all[*,*,i]*sf
		eB2D_all[*,*,i]			= eB2D_all[*,*,i]*sf

		wdot_all[*,*,i]			= wdot_all[*,*,i]*sf

		jPAlpha2D_all[*,*,i]	= jPAlpha2D_all[*,*,i]*sf
		jPBeta2D_all[*,*,i]		= jPBeta2D_all[*,*,i]*sf
		jPB2D_all[*,*,i]		= jPB2D_all[*,*,i]*sf

		eR2D_all[*,*,i]	= eR2D_all[*,*,i]*sf	
		eT2D_all[*,*,i]	= eT2D_all[*,*,i]*sf	
		eZ2D_all[*,*,i]	= eZ2D_all[*,*,i]*sf	
                                         
		bR2D_all[*,*,i]	= bR2D_all[*,*,i]*sf	
		bT2D_all[*,*,i]	= bT2D_all[*,*,i]*sf	
		bZ2D_all[*,*,i]	= bZ2D_all[*,*,i]*sf	

	endfor
	

	; Smooth 2D field data prior to Poyntin vector calc
	; -------------------------------------------------

	pSlice = where ( nPhi_all eq -22 )

	;smoothing = 5

	;eR2D_all[pSlice]	= smooth ( eR2D_all[pSlice], smoothing )
	;eT2D_all[pSlice]	= smooth ( eT2D_all[pSlice], smoothing )
	;eZ2D_all[pSlice]	= smooth ( eZ2D_all[pSlice], smoothing )

	;bR2D_all[pSlice]	= smooth ( bR2D_all[pSlice], smoothing )
	;bT2D_all[pSlice]	= smooth ( bT2D_all[pSlice], smoothing )
	;bZ2D_all[pSlice]	= smooth ( bZ2D_all[pSlice], smoothing )


	; 2D Poynting vector
	; ------------------

	u0	= 1.256e-6

	SpR2D	= 0.5 / u0 * real_part ( $
			conj ( eT2D_all[*,*,pSlice] ) * bZ2D_all[*,*,pSlice] $
			- conj ( eZ2D_all[*,*,pSlice] ) * bT2D_all[*,*,pSlice] )
	SpT2D	= -0.5 / u0 * real_part ( $
			conj ( eR2D_all[*,*,pSlice] ) * bZ2D_all[*,*,pSlice] $
			- conj ( eZ2D_all[*,*,pSlice] ) * bR2D_all[*,*,pSlice] )
	SpZ2D	= 0.5 / u0 * real_part ( $
			conj ( eR2D_all[*,*,pSlice] ) * bT2D_all[*,*,pSlice] $
			- conj ( eT2D_all[*,*,pSlice] ) * bR2D_all[*,*,pSlice] )



	; Create the phi grid
	; -------------------

	nP	= 360 
	dPhi	= 360.0 / (nP-1)
	phi	= fIndGen ( nP ) * dPhi * !dtor
	
	eAlpha3D	= complexArr ( nR, nz, nP )
	eBeta3D	= complexArr ( nR, nz, nP )
	eB3D	= complexArr ( nR, nz, nP )

	jPAlpha3D	= complexArr ( nR, nz, nP )
	jPBeta3D	= complexArr ( nR, nz, nP )
	jPB3D	= complexArr ( nR, nz, nP )

	eR3D	= complexArr ( nR, nz, nP )
	eT3D	= complexArr ( nR, nz, nP )
	eZ3D	= complexArr ( nR, nz, nP )

	bR3D	= complexArr ( nR, nz, nP )
	bT3D	= complexArr ( nR, nz, nP )
	bZ3D	= complexArr ( nR, nz, nP )

	;wdot3D	= fltArr ( nR, nz, nP )
	wdot3D	= fltArr ( nR, nz, nP )

	II	= complex ( 0, 1 )

	print, 'Summing fields ...'

	for f=0,nFiles-1 do begin
		for p=0,nP-1 do begin

			expFac	= complex(0,1) * nPhi_all[f] * phi[p]

			eAlpha3D[*,*,p]	+= nPhi_ant_all[f] * eAlpha2D_all[*,*,f] * exp ( expFac ) 
			eBeta3D[*,*,p]	+= nPhi_ant_all[f] * eBeta2D_all[*,*,f] * exp ( expFac ) 
			eB3D[*,*,p]		+= nPhi_ant_all[f] * eB2D_all[*,*,f] * exp ( expFac ) 

			jPAlpha3D[*,*,p]	+= nPhi_ant_all[f] * jPAlpha2D_all[*,*,f] * exp ( expFac ) 
			jPBeta3D[*,*,p]		+= nPhi_ant_all[f] * jPBeta2D_all[*,*,f] * exp ( expFac ) 
			jPB3D[*,*,p]		+= nPhi_ant_all[f] * jPB2D_all[*,*,f] * exp ( expFac ) 

			eR3D[*,*,p]	+= nPhi_ant_all[f] * eR2D_all[*,*,f] * exp ( expFac ) 
			eT3D[*,*,p]	+= nPhi_ant_all[f] * eT2D_all[*,*,f] * exp ( expFac ) 
			eZ3D[*,*,p]	+= nPhi_ant_all[f] * eZ2D_all[*,*,f] * exp ( expFac ) 

			bR3D[*,*,p]	+= nPhi_ant_all[f] * bR2D_all[*,*,f] * exp ( expFac ) 
			bT3D[*,*,p]	+= nPhi_ant_all[f] * bT2D_all[*,*,f] * exp ( expFac ) 
			bZ3D[*,*,p]	+= nPhi_ant_all[f] * bZ2D_all[*,*,f] * exp ( expFac ) 

			iiBad = where ( eR2D_all[*,*,f] ne eR2D_all[*,*,f], iiBadCnt )
			if iiBadCnt gt 0 then begin
				print, 'NaNs in the e field solution, arrgghh.'
				stop
			endif

		endfor
	endfor

	print, 'DONE'


	; 3D E*.J
	; -------

	eDotJ3D	= 0.5 * real_part ( $
			conj ( eAlpha3D ) * jPAlpha3D + $
			conj ( eBeta3D ) * jPBeta3D + $
			conj ( eB3D ) * jPB3D )


	; 3D Poynting vector
	; ------------------

	SpR	= 0.5 / u0 * real_part ( $
			conj ( eT3D ) * bZ3D - conj ( eZ3D ) * bT3D )
	SpT	= -0.5 / u0 * real_part ( $
			conj ( eR3D ) * bZ3D - conj ( eZ3D ) * bR3D )
	SpZ	= 0.5 / u0 * real_part ( $
			conj ( eR3D ) * bT3D - conj ( eT3D ) * bR3D )

	;SpScale	= max ( [ SpR[*], SpT[*], SpZ[*] ] )

	;SpR	= SpR / SpScale * 100
	;SpT	= SpT / SpScale * 100
	;SpZ	= SpZ / SpScale * 100


	;; Plot a 2D slice from the 3D field reconstruction
	;; ------------------------------------------------

	;slice = 0
	;scaleFac = 0.2 
	;loadct, 13, file = 'davect.tbl', rgb_table = rgb_table
	;loadct, 3, rgb_table = rgb_table
	;rgb_table = reverse ( rgb_table )

	;nLevs	= 20 
	;scale	= max ( abs([eAlpha3D[*,*,slice],eBeta3D[*,*,slice]]) ) * scaleFac
	;levels	= (fIndGen ( nLevs )/(nLevs-1)-0.5) * scale  * 2 * 1.1
	;colors	= bytScl ( levels, top = 253 ) + 1

	;p = contour ( (abs(eAlpha3D[*,*,slice])<scale)>(-scale), R, z, $
	;		c_value = levels, $
	;		rgb_indices = colors, $
	;		rgb_table = rgb_table, $
	;		title = 'eAlpha', $
	;		aspect_ratio = 1.0, $
	;		/fill )


	; Smooth Poynting vector but remember the periodic 
	; direction must be treated specially
	; ------------------------------------------------

	print, 'Smoothing Poynting vector for plotting ...'
	smoothWidth = 10

	SpR	= smooth ( [ [[SpR]],[[SpR]],[[SpR]] ], smoothWidth, /edge )
	SpT	= smooth ( [ [[SpT]],[[SpT]],[[SpT]] ], smoothWidth, /edge )
	SpZ	= smooth ( [ [[SpZ]],[[SpZ]],[[SpZ]] ], smoothWidth, /edge )

	SpR	= SpR[*,*,nP:nP+nP-1]
	SpT	= SpT[*,*,nP:nP+nP-1]
	SpZ	= SpZ[*,*,nP:nP+nP-1]

	;; Rotate Poynting vector to Cartesian
	;; for plotting in VisIt. This is only because
	;; the vector rotation in VisIt seems broken.
	;; -------------------------------------------

	;print, 'Rotating Poynting vector to Cartesian ...'
	;for p=0,nP-1 do begin

	;	rotMat	= [ [ cos ( phi[p] ), sin ( phi[p] ), 0 ], $
	;				[-sin ( phi[p] ), cos ( phi[p] ), 0 ], $
	;				[ 0, 0, 1 ] ]
	;	rotMat	= transpose ( rotMat )

	;	for j=0,nZ-1 do begin
	;		for i=0,nR-1 do begin

	;			cylVec	= [ [SpR[i,j,p]],[SpT[i,j,p]],[SpZ[i,j,p]] ]
	;			carVec	= rotMat ## cylVec

	;			SpR[i,j,p]	= carVec[0]
	;			SpT[i,j,p]	= carVec[1]
	;			SpZ[i,j,p]	= carVec[2]

	;		endfor
	;	endfor
	;endfor


	; Create 3D coordinate variables
	; ------------------------------

	r3D	= rebin ( R, nR, nz, nP )
	z3D	= transpose ( rebin ( z, nz, nR, nP ), [1,0,2] )
	p3D	= transpose ( rebin ( phi, nP, nR, nz ), [1,2,0] )


	; Get E field structure along a field line
	; ----------------------------------------

	eqDskFileName	= 'eqdsk'
	eqdsk	= readGEQDSK ( eqDskFileName, bInterpS = bInterpS, $
			fieldLineIn=[0.479,10.33,1.06], fieldLineOut = fieldLineOut, $
		   	use_dlg_bField = use_dlg_bField	)

	if keyword_set ( use_dlg_bField ) then begin

		print, 'Using dlg_bField.nc'
		cdfId = ncdf_open ( 'dlg_bField.nc', /noWrite ) 

			ncdf_varget, cdfId, 'R', R_dlg_bField 
			ncdf_varget, cdfId, 'z', z_dlg_bField 
			ncdf_varget, cdfId, 'bR', bR_dlg_bField 
			ncdf_varget, cdfId, 'bPhi', bPhi_dlg_bField 
			ncdf_varget, cdfId, 'bz', bz_dlg_bField 

		ncdf_close, cdfId

    	bInterpS    = { bR : bR_dlg_bField, $
    	            rleft : r_dlg_bField[0], $
    	            rdim : r_dlg_bField[-1]-r_dlg_bField[0], $
    	            nW : n_elements(r_dlg_bField), $
    	            z : z_dlg_bField, $
    	            zdim : z_dlg_bField[-1]-z_dlg_bField[0], $
    	            nH : n_elements(z_dlg_bField), $
    	            bPhi : bPhi_dlg_bField, $
    	            bz : bz_dlg_bField }   

	endif
	;
	;rIntCoord	= ( fieldLineOut[*,0] - R[0] ) / (R[-1]-R[0]) * (nR-1)
	;tIntCoord	= ( fieldLineOut[*,1] - p[0] ) / (phi[-1]-phi[0]) * (np-1)
	;zIntCoord	= ( fieldLineOut[*,2] - z[0] ) / (z[-1]-z[0]) * (nz-1)

	;eAlphaFieldLine	= interpolate ( eAlpha3D, rIntCoord, zIntCoord, tIntCoord, cubic = -0.5 )
	;eBetaFieldLine	= interpolate ( eBeta3D, rIntCoord, zIntCoord, tIntCoord, cubic = -0.5 )
	;eBFieldLine	= interpolate ( eB3D, rIntCoord, zIntCoord, tIntCoord, cubic = -0.5 )


	;; Look at field line time dependence for standing wave
	;; ----------------------------------------------------

	;freq	= 30e6
	;w		= 2*!pi*freq
	;dt		= 1/freq/25

	;for t=0,1000 do begin $
	;&	plot, real_part(eAlphaFieldLine)*cos(w*dt*t) $
	;		+ imaginary(eAlphaFieldLine)*sin(w*dt*t), $
	;	   	yRange = [-0.015, 0.015] $
	;&	oplot, abs(eAlphaFieldLine), thick = 2 $
	;&	oplot, real_part(eAlphaFieldLine), thick = 2 $
	;&	wait, 0.04 $
	;&endfor


	;; Create a list of starting locations for VisIt streamlines
	;; ---------------------------------------------------------

	;openw, lun, 'streamlineLocs.txt', /get_lun

	;SpMag	= sqrt ( SpR^2 + SpT^2 + SpZ^2 )
	;for p=0,nP-1 do begin

	;	thisSpMag	= SpMag[*,*,p]

	;	thisR3D		= R3D[*,*,p]
	;	thisZ3D		= z3D[*,*,p]

	;	iiZero	= where ( thisz3D lt -0.5 or thisz3D gt 0.5, iiCnt )
	;	thisSpMag[iiZero] = 0
	;	iiMax	= where ( thisSpMag eq max ( thisSpMag ), iiCnt ) 

	;	rMax	= thisR3D[iiMax]
	;	zMax	= thisz3D[iiMax]
	;	pMax	= phi[p]

	;	xMax	= rMax * cos ( pMax )
	;	yMax	= rMax * sin ( pMax )

	;	printf, lun, xMax, yMax, zMax, format = '(3(f9.6,1x))'

	;endfor

	;close, lun

	;; Testing VisIt plug-in
	;; ---------------------

	;SpR[*]	= 0
	;SpT[*]	= 0
	;SpZ[*]	= 0

	;ii	= where ( p3D lt !pi/4 and z3d gt 0 and r3d lt 0.4, iiCnt )
	;SpZ[ii]	= 1


	; Normalise to power absorbed
	; ---------------------------

	dr	= R[1]-R[0]
	dz	= z[1]-z[0]
	int_eDotJ	= total ( eDotJ3D * r3D * dr * dz * dPhi )
	
	print, 'Total wDot: ', int_eDotJ

	powerRF	= 1.8d6; [W]
	pScale_	= powerRF / int_eDotJ

	eAlpha3D	= eAlpha3D * sqrt ( pScale_ )
	eBeta3D	= eBeta3D * sqrt ( pScale_ )
	eB3D	= eB3D * sqrt ( pScale_ )

	SpR	= SpR * sqrt ( pScale_ )
	SpT	= SpT * sqrt ( pScale_ )
	SpZ	= SpZ * sqrt ( pScale_ )

	eDotJ3D	= eDotJ3D / pScale_


	; Create an antenna structure for visit
	; -------------------------------------

	ant3D	= fltArr ( nR, nz, nP )

	antShift = 0.05
	rMin	= (rAnt+antShift)-0.025
	rMax	= (rAnt+antShift)+0.025
	zMin	= -0.25
	zMax	= +0.25

	c2cPhi	= wd / (2*!pi*rAnt) * 2 * !pi; [rad]
	strapPhi	= xlt / (2*!pi*rAnt) * 2 * !pi; [rad]

	phiOffset	= -12 * c2cPhi / 2.0

	for i=0, nStrap-1 do begin

		pMin	= phiOffset+i*c2cPhi
		pMax	= phiOffset+i*c2cPhi+strapPhi
		if pMin lt 0*!dtor then pMin += 360*!dtor
		if pMax lt 0*!dtor then pMax += 360*!dtor
		if pMax gt 360*!dtor then pMax -= 360*!dtor

		iiStrap	= where ( r3D gt rMin and r3D lt rMax $
					and z3D gt zMin and z3D lt zMax $
					and p3D gt pMin and p3D lt pMax, iiCnt )

		print, pMin*!radeg, pMax*!radeg, iiCnt
		if iiCnt gt 0 then ant3D[iiStrap] = 1

	endfor


	; Cut out section of plot based on LCFS
	; -------------------------------------

    bHere   = interpB ( bInterpS, (r3D[*,*,0])[*], (z3D[*,*,0])[*] )

	bEqR3D	= fltArr ( nR, nz, nP )
	bEqT3D	= fltArr ( nR, nz, nP )
	bEqZ3D	= fltArr ( nR, nz, nP )

	for p=0,nP-1 do begin
		bEqR3D[*,*,p]	= reform((bHere.bR)[*],nR,nZ)
		bEqT3D[*,*,p]	= -reform((bHere.bPhi)[*],nR,nZ)
		bEqZ3D[*,*,p]	= reform((bHere.bZ)[*],nR,nZ)
	endfor

	overSample_boundary, (eqdsk.rbbbs)[*], (eqdsk.zbbbs)[*], rLCFS, zLCFS	
	overSample_boundary, (eqdsk.rlim)[*], (eqdsk.zlim)[*], rLIM, zLIM	

	maskRho	= intArr ( nR, nz )
	for i=0,nz-1 do begin
		maskRho[*,i] = is_in_shape ( R, fltArr(nR)+z[i], rLCFS, zLCFS )
	endfor

	maskLim	= intArr ( nR, nz )
	for i=0,nz-1 do begin
		maskLim[*,i] = is_in_shape ( R, fltArr(nR)+z[i], rLIM, zLIM )
	endfor

	maskRho3D	= rebin ( maskRho, nR, nz, nP )
	maskLim3D	= rebin ( maskLim, nR, nz, nP )

	iiZero	= where ( r3D gt 0.3 or $
			z3d gt 1.5 or $
			z3d lt -1.5, iiCnt ) 
			;(p3D gt 2*!pi/6 and p3D lt 2*!pi/6*5), iiCnt )

	maskLim3D	-= 1
	maskLim3D	= abs ( maskLim3D )
	maskLim3D[iiZero]	= 0

	;R3D	= rebin ( R, nR, nz, nP )
	;z3D	= transpose ( rebin ( z, nz, nR, nP ), [1,0,2] )
	;p3D	= transpose ( rebin ( phi, np, nR, nz ), [1,2,0] )

	;iiZero	= where ( mask3D eq 0 and p3D gt !pi and p3d lt 2*!pi, iiCnt )
	;eAlpha3D[iiZero] = 0

	;; Smooth the eMag in the p direction (periodic)
	;; ----------------------------------

	eMag	= sqrt ( abs(eAlpha3D)^2 + abs(eBeta3D)^2 + abs (eB3D)^2 )

	eMagOrig = eMag

	;for i=0,nR-1 do begin
	;	for j=0,nZ-1 do begin
	;		for k=0,nP-1 do begin

	;			nSmooth = 10
	;			kkArr = indGen(nSmooth)-nSmooth/2+k
	;			kkL = where(kkArr gt nP-1,kkLCnt)
	;			if kkLCnt gt 0 then kkArr[kkL]=kkArr[kkL]-nP

	;			tmp = 0
	;			for s=0,nSmooth-1 do begin
	;				tmp = tmp + eMagOrig[i,j,kkArr[s]]
	;			endfor

	;			eMag[i,j,k] = tmp / nSmooth	

	;		endfor
	;	endfor
	;endfor



	; Write netCDF file for visit (r,t,z)
	; -----------------------------------
	;
	; REMEMBER: A visit database plugin has been 
	; written for this specific format, so if you 
	; change it you will need to correct and recompile
	; the avt*FileFormat.C file for the plugin.
	;

	print, 'Writing netCdf file ...'	

	nTimeSlices	= 1 
	freq	= 30e6
	w		= 2*!pi*freq
	dt		= 1/freq/nTimeSlices

	for t=0,nTimeSlices-1 do begin

		e3D	= real_part ( eAlpha3D ) * cos ( w*t*dt) $
				+ imaginary ( eAlpha3D ) * sin ( w*t*dt )

		outFileName = '3d_fields_'+string(t,format='(i3.3)')+'.nca'

		nc_id	= nCdf_create ( outFileName, /clobber, /netcdf4_format )
		nCdf_control, nc_id, /nofill
		
		nR_id	= nCdf_dimDef ( nc_id, 'nR', nR )
		nz_id	= nCdf_dimDef ( nc_id, 'nz', nz )
		np_id	= nCdf_dimDef ( nc_id, 'nP', nP )

		r_id = ncdf_vardef ( nc_id, 'r', [nr_id], /float )
		z_id = ncdf_vardef ( nc_id, 'z', [nz_id], /float )
		p_id = ncdf_vardef ( nc_id, 'p', [np_id], /float )

		eAbs_id = ncdf_vardef ( nc_id, 'eAbs', [nr_id, np_id, nz_id], /float )
		eTime_id = ncdf_vardef ( nc_id, 'eTime', [nr_id, np_id, nz_id], /float )
		ant_id = ncdf_vardef ( nc_id, 'ant', [nr_id, np_id, nz_id], /float )
		lim_id = ncdf_vardef ( nc_id, 'lim', [nr_id, np_id, nz_id], /float )
		bbbs_id = ncdf_vardef ( nc_id, 'bbbs', [nr_id, np_id, nz_id], /float )
		eDotJ_id = ncdf_vardef ( nc_id, 'eDotJ', [nr_id, np_id, nz_id], /float )

		spR_id = ncdf_vardef ( nc_id, 'SpR', [nr_id, np_id, nz_id], /float )
		spT_id = ncdf_vardef ( nc_id, 'SpT', [nr_id, np_id, nz_id], /float )
		spZ_id = ncdf_vardef ( nc_id, 'SpZ', [nr_id, np_id, nz_id], /float )

		bEqR_id = ncdf_vardef ( nc_id, 'bEqR', [nr_id, np_id, nz_id], /float )
		bEqT_id = ncdf_vardef ( nc_id, 'bEqT', [nr_id, np_id, nz_id], /float )
		bEqZ_id = ncdf_vardef ( nc_id, 'bEqZ', [nr_id, np_id, nz_id], /float )


		nCdf_control, nc_id, /enDef
		
		nCdf_varPut, nc_id, r_id, R
		nCdf_varPut, nc_id, z_id, z 
		nCdf_varPut, nc_id, p_id, phi 
		nCdf_varPut, nc_id, eAbs_id, transpose ( eMag, [0,2,1] )
		nCdf_varPut, nc_id, eTime_id, transpose ( e3D, [0,2,1] )
		nCdf_varPut, nc_id, ant_id, transpose ( ant3D, [0,2,1] )
		nCdf_varPut, nc_id, lim_id, transpose ( maskLim3D, [0,2,1] )
		nCdf_varPut, nc_id, bbbs_id, transpose ( maskRho3D, [0,2,1] )
		nCdf_varPut, nc_id, eDotJ_id, transpose ( eDotJ3D, [0,2,1] )

		nCdf_varPut, nc_id, spR_id, transpose ( SpR, [0,2,1] )
		nCdf_varPut, nc_id, spT_id, transpose ( SpT, [0,2,1] )
		nCdf_varPut, nc_id, spZ_id, transpose ( SpZ, [0,2,1] )

		nCdf_varPut, nc_id, bEqR_id, transpose ( bEqR3D, [0,2,1] )
		nCdf_varPut, nc_id, bEqT_id, transpose ( bEqT3D, [0,2,1] )
		nCdf_varPut, nc_id, bEqZ_id, transpose ( bEqZ3D, [0,2,1] )

		nCdf_close, nc_id

	endfor

	print, 'DONE'
stop	
end
