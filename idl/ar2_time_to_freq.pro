function ar2_time_to_freq, s, t, freq, i=i


	;tMin = t[0]
	;tMax = t[-1]
	;t_pad = fIndGen(nT*PadFac)/(nT*PadFac-1)*(tMax-tMin)+tMin
	;s_pad = interpol(s,t,t_pad,/spline)
	;dt_pad = t_pad[1]-t_pad[0]

	; Re-sample to an integer number of periods of the 
	; drive freq

	; Test code

	nPeriod = 15
    nT = n_elements(t)
	t = fIndGen(nT)/(nT-1)*(nPeriod/freq)
	;t = t[0:-2]
	sRe = 1000.0
	sIm = -300.0
	amp = sqrt(sRe^2+sIm^2)
	phs = atan(sIm,sRe)
	wrf = freq*2*!pi
	II	= complex(0,1)
	s = amp * cos (wrf*t+phs)

	tRange = t[-1]-t[0]
	nPeriod = floor(tRange * freq)
    nT = n_elements(t)
	ThisT = fIndGen(nT)/(nT-1)*(nPeriod/freq)
	ThisT = ThisT[0:-2]+t[0]
	ThisS = interpol(s,t,ThisT,/spline)

    nT = n_elements(ThisT)
    if n_elements(ThisS) ne nT then stop
	hanWindow = hanning (nT, alpha=0.5 )
	hanFac = 2.0
	PadFac = 10
	dt = ThisT[1]-ThisT[0]
	s_pad = FltArr(nT*PadFac)
	s_pad[0:nT-1] = ThisS*hanWindow

	sFFT = fft ( s_pad[*], /center )
	freqAxis = (fIndGen(nT*PadFac)-(nT*PadFac/2)) / (nT*PadFac*dt)
	freqNorm = freqAxis / freq

	; Positive (right) frequency 
	iiAntFreq = where(abs(freqNorm-1) eq min(abs(freqNorm-1)),iiAntFreqCnt)
	rpL = real_part ( sFFT[iiAntFreq[0]] )
	ipL = imaginary ( sFFT[iiAntFreq[0]] )

	xRange = [-5,+5]*PadFac+iiAntFreq[0]
	p=plot(real_part(sFFT),xRange=xRange,layout=[2,2,1])
	p=plot([0,0]+iiAntFreq[0],[-1,1]*max(abs(real_part(sFFT))),color='b',/over)
	p=plot(imaginary(sFFT),color='blue',xRange=xRange,layou=[2,2,2],/current)
	p=plot([0,0]+iiAntFreq[0],[-1,1]*max(abs(imaginary(sFFT))),color='b',/over)


	; Negative (left) frequency 
	iiAntFreq = where(abs(freqNorm+1) eq min(abs(freqNorm+1)),iiAntFreqCnt)
	rpR = real_part ( sFFT[iiAntFreq[0]] )
	ipR = -imaginary ( sFFT[iiAntFreq[0]] )

	xRange = [-5,+5]*PadFac+iiAntFreq[0]
	p=plot(real_part(sFFT),xRange=xRange,layout=[2,2,3],/current)
	p=plot([0,0]+iiAntFreq[0],[-1,1]*max(abs(real_part(sFFT))),color='r',/over)
	p=plot(imaginary(sFFT),color='blue',xRange=xRange,layou=[2,2,4],/current)
	p=plot([0,0]+iiAntFreq[0],[-1,1]*max(abs(imaginary(sFFT))),color='r',/over)

	rpL = rpL * padFac * hanFac
	rpR = rpR * padFac * hanFac

	ipL = ipL * padFac * hanFac
	ipR = ipR * padFac * hanFac


	print, 'Re: ', rpL, rpR, rpL+rpR
	print, 'Im: ', ipL, ipR, ipL+ipR

stop
if keyword_set(i) then begin
	if i eq 190 then stop
endif
	return, complex ( rpL+rpR, ipL+ipR )
end
