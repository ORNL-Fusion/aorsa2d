pro test_chebyshev_k

	nPts	= 20000
	xBasis = ( fIndGen ( nPts+1 ) / nPts - 0.5 ) * 2 
	xLeft = 0.0
	xRight = 2.0
	xRange = xRight-xLeft
	x	= (xBasis + 1)/2*xRange+xLeft
	n	= 50
	chebT	= cos ( n * acos ( xBasis ) )

	!p.multi = [0,1,2]
	plot, x, chebT	

	nFFT	= 100.0
	fftStep	= nPts / (nFFT-1)
	fftWidth	= nPts / 10.0 
	fftWidth	+= (fftWidth mod 2 )
	hannWin	= hanning ( fftWidth )
	osf	= 10

	fft2d	= fltArr ( nFFT, fftWidth*osf )
	tmp	= complexArr ( fftWidth * osf )

	for i=0,nFFT-1 do begin
			iiLow	= i * fftStep-fftWidth/2
			iiHig	= i * fftStep+fftWidth/2-1
			if iiLow gt 0 and iiHig lt nPts-2 then begin
				tmp2	= chebT[iiLow:iiHig]*hannWin
				tmp[fftWidth*osf/2-fftWidth/2:fftWidth*osf/2+fftWidth/2-1] = $
						tmp2
				fft2d[i,*]	= abs(fft(tmp))^2
			endif
	endfor

	dx	= x[1]-x[0]
	k	= (fIndGen(fftWidth*osf)) / ( fftWidth * osf * dx ) * 2 * !pi

	fft2d_pm	= fft2d
	nqII	= n_elements(fft2d[0,*])/2
	fft2d_pm[*,nqii:*] = fft2d[*,0:nqii-1]
	fft2d_pm[*,0:nqii-1] = fft2d[*,nqii:*]
	k_pm = k - k[nqii]

	levels	= 9.0^findgen(10)*1e-9
	contour, fft2d_pm, fIndGen(nFFT)/(nFFT-1)*xRange+xLeft, k_pm, /fill, $
			yRange = [-6*n,6*n], $
			xRange = [xLeft*0.8,xRight*1.2], $
			levels = levels

	j	= fIndGen(nPts+1)
	oPlot, x[1:nPts], n/(sin(!pi * j[1:nPts]/nPts))^0.5 * 2 / xRange, thick = 1
	;oPlot, x, -n/cos(xBasis*!Pi/2)^0.5/xrange^0.5, thick = 1

	stop
end
