pro ar2_read_matrix

    FileName = 'matrix.nc'

	cdfId = ncdf_open ( FileName, /noWrite ) 
		nCdf_varGet, cdfId, 'reA', reA
		nCdf_varGet, cdfId, 'imA', imA 
    nCdf_close, cdfId	

    A = dcomplex(reA,imA)

    nPts = n_elements(A[*,0])
    i1 = lIndGen(nPts/3)*3
    i2 = i1+1
    i3 = i2+1

    A[*,i1] = fft(A[*,i1],-1,dimension=2);,/center)
    A[*,i2] = fft(A[*,i2],-1,dimension=2);,/center)
    A[*,i3] = fft(A[*,i3],-1,dimension=2);,/center)

    evals = la_eigenproblem(A,double=1)

    x = signum(real_part(evals))*alog10(abs(real_part(evals)))
    y = signum(imaginary(evals))*alog10(abs(imaginary(evals)))

    ;x1 = real_part(evals)^2-imaginary(evals)^2
    ;y1 = 2*real_part(evals)*imaginary(evals)

    x1 = real_part(evals)
    y1 = imaginary(evals)

    p=plot(x1,y1,symbol='D',lineStyle='none',xTitle='real',yTitle='imag',title='l',$
        layout=[1,2,1],aspect_ratio=1.0)
    p=plot(x,y,symbol='D',lineStyle='none',xTitle='real',yTitle='imag',$
        title='sgn(l)*alog10(abs(l))',layout=[1,2,2],/current,aspect_ratio=1.0)

    stop

end
