function ar2_time_to_freq, s, t, freq, i=i, PlotSpec=PlotSpec, NStop = NStop

    ThisNStop = 0
    if keyword_set(NStop)then ThisNStop=NStop

	;t = t-t[0]
    nT = n_elements(t)
	II	= complex(0,1)

	;; Test code
	;; ---------
	;;nPeriod = 8 
    ;;nT = n_elements(t)
	;;t = fIndGen(nT)/(nT-1)*(nPeriod/freq)
	;;t = t[0:-2]
	;sRe = 1000.0
	;sIm = -300.0
	;amp = sqrt(sRe^2+sIm^2)
	;phs = atan(sIm,sRe)
	;wrf = freq*2*!pi
	;;s = amp * cos (wrf*t+phs)
	;s = real_part(amp * exp (-II*(wrf*t+phs)))

	; Extract the LAST exact period with specified
	; initial phase to calibrate to AORSA, i.e., try
	; setting the start point of the FFT window such that
	; we preserve the phase relative to the initial 
	; driver phase.

	; Find the number of the last complete cycle

	PeriodFractionOffset = -0.75; Should always be -ve
	period = 1.0/freq
	t_period = floor(t/period)
	last_period = max(t_period)-1 ; starting at 0
	StartTime = period*last_period + PeriodFractionOffset*period 

    print, 'LastPeriod: ', last_period
    print, 'StartTime: ', StartTime

	;tRange = t[-1]-t[0]
	;nPeriod = floor(tRange * freq)
	nPeriod = 1
	;print, 'nPeriod: ', nPeriod
    nT = n_elements(t)
	ThisT = fIndGen(nT)/(nT-1)*(nPeriod/freq)
	ThisT = ThisT[0:-2]
	ThisT = ThisT + StartTime
	ThisS = interpol(s,t,ThisT,/spline)

    nT = n_elements(ThisT)
    if n_elements(ThisS) ne nT then stop
	hanWindow = hanning (nT, alpha=0.5 )
	hanFac = 1
	PadFac = 1
	dt = ThisT[1]-ThisT[0]
	s_pad = FltArr(nT*PadFac)
	s_pad[0:nT-1] = ThisS;*hanWindow

	sFFT = fft ( s_pad[*], /center )
	freqAxis = (fIndGen(nT*PadFac)-(nT*PadFac/2)) / (nT*PadFac*dt)
	freqNorm = freqAxis / freq

	; Positive (right) frequency 
	iiAntFreq = where(abs(freqNorm-1) eq min(abs(freqNorm-1)),iiAntFreqCnt)
	rpL = real_part ( sFFT[iiAntFreq[0]] )
	ipL = imaginary ( sFFT[iiAntFreq[0]] )

PlotThis = 0
if keyword_set(i) then begin
    if i eq ThisNStop and keyword_set(PlotSpec) then PlotThis = 1
endif

if PlotThis then begin

	xRange = [-5,+5]*PadFac+iiAntFreq[0]
	p=plot(real_part(sFFT),xRange=xRange,layout=[2,2,1])
	p=plot([0,0]+iiAntFreq[0],[-1,1]*max(abs(real_part(sFFT))),color='b',/over)
	p=plot(imaginary(sFFT),color='blue',xRange=xRange,layou=[2,2,2],/current)
	p=plot([0,0]+iiAntFreq[0],[-1,1]*max(abs(imaginary(sFFT))),color='b',/over)

endif

	; Negative (left) frequency 
	iiAntFreq = where(abs(freqNorm+1) eq min(abs(freqNorm+1)),iiAntFreqCnt)
	rpR = real_part ( sFFT[iiAntFreq[0]] )
	ipR = -imaginary ( sFFT[iiAntFreq[0]] )

if PlotThis then begin

	xRange = [-5,+5]*PadFac+iiAntFreq[0]
	p=plot(real_part(sFFT),xRange=xRange,layout=[2,2,3],/current)
	p=plot([0,0]+iiAntFreq[0],[-1,1]*max(abs(real_part(sFFT))),color='r',/over)
	p=plot(imaginary(sFFT),color='blue',xRange=xRange,layou=[2,2,4],/current)
	p=plot([0,0]+iiAntFreq[0],[-1,1]*max(abs(imaginary(sFFT))),color='r',/over)

endif

	rpL = rpL * padFac * hanFac
	rpR = rpR * padFac * hanFac

	ipL = ipL * padFac * hanFac
	ipR = ipR * padFac * hanFac

	print, 'Re: ', rpL, rpR, rpL+rpR
	print, 'Im: ', ipL, ipR, ipL+ipR

	rp = rpR;+rpR
	ip = ipR;+ipR

	_amp = sqrt(rp^2+ip^2)
	_phs = atan(ip,rp)

	;print, amp, phs*!radeg
	;print, _amp, _phs*!radeg


if keyword_set(i) and keyword_set(PlotSpec) and i eq ThisNStop then stop

	answer = complex ( rp, ip ) 
	; To perhaps correct for an e^-iwt sign error I think i need
	; to conj the answer. At least that makes vorpal consistent
	; with aorsa and rsfwc_1d
	return, conj(answer)*2
end
