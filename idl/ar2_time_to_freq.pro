function ar2_time_to_freq, s, t, freq

    nT = n_elements(t)
    if n_elements(s) ne nT then stop

	hanWindow = hanning (nT, alpha=0.5 )

	dt = t[1]-t[0]

	;; Test code
	;ampRe = 1000.0
	;ampIm = -300.0
	;amp = sqrt(ampRe^2+ampIm^2)
	;phs = atan(ampIm,ampRe)
	;s = (amp * exp (-II*(wrf*t+phs)))

	sFFT = fft ( s[*]*hanWindow, /center )
	freqAxis = (fIndGen(nT)-(nT/2)) / (nT*dt)
	freqNorm = freqAxis / freq

	;p=plot(freqNorm,abs(sFFT)^2);,xRange=[0,5])
	;p=plot(freqNorm,imaginary(sFFT),color='blue');,xRange=[0,5])

	; Positive (right) frequency 
	iiAntFreq = where(abs(freqNorm-1) eq min(abs(freqNorm-1)),iiAntFreqCnt)
	rpL = real_part ( sFFT[iiAntFreq[0]] )
	ipL = imaginary ( sFFT[iiAntFreq[0]] )
	; Negative (left) frequency 
	iiAntFreq = where(abs(freqNorm+1) eq min(abs(freqNorm+1)),iiAntFreqCnt)
	rpR = real_part ( sFFT[iiAntFreq[0]] )
	ipR = -imaginary ( sFFT[iiAntFreq[0]] )

	;print, 'Re: ', rpL, rpR, rpL+rpR
	;print, 'Im: ', ipL, ipR, ipL+ipR

	return, complex ( rpL+rpR, ipL+ipR )

end
