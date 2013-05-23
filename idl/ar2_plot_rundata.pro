pro ar2_plot_rundata

	fileList = file_search ( 'runData*.nc' )
    FileListIndex = 0

	cdfId = ncdf_open ( fileList[FileListIndex], /noWrite ) 
		nCdf_varGet, cdfId, 'capR', x 
		nCdf_varGet, cdfId, 'y', y 
		nCdf_varGet, cdfId, 'jr_re', jr_re 
		nCdf_varGet, cdfId, 'jr_im', jr_im 
		nCdf_varGet, cdfId, 'jt_re', jt_re 
		nCdf_varGet, cdfId, 'jt_im', jt_im 
		nCdf_varGet, cdfId, 'jz_re', jy_re 
		nCdf_varGet, cdfId, 'jz_im', jy_im 
		nCdf_varGet, cdfId, 'bmod', bmod
		nCdf_varGet, cdfId, 'brU', bxn
		nCdf_varGet, cdfId, 'btU', bzn
		nCdf_varGet, cdfId, 'bzU', byn

		nCdf_varGet, cdfId, 'densitySpec', densitySpec
		nCdf_varGet, cdfId, 'tempSpec', tempSpec
		nCdf_varGet, cdfId, 'nuOmg', nuOmgSpec
		nCdf_varGet, cdfId, 'omgc', omgc 
		nCdf_varGet, cdfId, 'omgp2', omgp2
	
		;nCdf_varGet, cdfId, 'xx_re', xx_re 
		;nCdf_varGet, cdfId, 'xx_im', xx_im
		;nCdf_varGet, cdfId, 'yy_re', yy_re 
		;nCdf_varGet, cdfId, 'yy_im', yy_im

		nCdf_varGet, cdfId, 'drUrr', drUrr
		nCdf_varGet, cdfId, 'drUrt', drUrt
		nCdf_varGet, cdfId, 'drUrz', drUrz

		nCdf_varGet, cdfId, 'drUtr', drUtr
		nCdf_varGet, cdfId, 'drUtt', drUtt
		nCdf_varGet, cdfId, 'drUtz', drUtz

		nCdf_varGet, cdfId, 'drUzr', drUzr
		nCdf_varGet, cdfId, 'drUzt', drUzt
		nCdf_varGet, cdfId, 'drUzz', drUzz

		nCdf_varGet, cdfId, 'dzUrr', dzUrr
		nCdf_varGet, cdfId, 'dzUrt', dzUrt
		nCdf_varGet, cdfId, 'dzUrz', dzUrz

		nCdf_varGet, cdfId, 'dzUtr', dzUtr
		nCdf_varGet, cdfId, 'dzUtt', dzUtt
		nCdf_varGet, cdfId, 'dzUtz', dzUtz

		nCdf_varGet, cdfId, 'dzUzr', dzUzr
		nCdf_varGet, cdfId, 'dzUzt', dzUzt
		nCdf_varGet, cdfId, 'dzUzz', dzUzz

		nCdf_varGet, cdfId, 'Urr', Urr
		nCdf_varGet, cdfId, 'Urt', Urt
		nCdf_varGet, cdfId, 'Urz', Urz

		nCdf_varGet, cdfId, 'Utr', Utr
		nCdf_varGet, cdfId, 'Utt', Utt
		nCdf_varGet, cdfId, 'Utz', Utz

		nCdf_varGet, cdfId, 'Uzr', Uzr
		nCdf_varGet, cdfId, 'Uzt', Uzt
		nCdf_varGet, cdfId, 'Uzz', Uzz

	ncdf_close, cdfId

	;xx	= complex ( xx_re, xx_im )
	;yy	= complex ( yy_re, yy_im )

	bx	= bxn * bmod
	by	= byn * bmod
	bz	= bzn * bmod

	nR	= n_elements ( x )
	nz	= n_elements ( y )
	nSpec	= n_elements ( densitySpec[0,0,*] )

	p=plot(x, bMod[*,nz/2.0], $
			title = 'z=0 equil B field')
	p=plot(x, bx[*,nz/2.0], /over)
	p=plot(x, by[*,nz/2.0], /over)
	p=plot(x, bz[*,nz/2.0], /over)

	p=plot(x, densitySpec[*,nz/2.0,0], title = 'z=0 Density',yrange=range,/ylog)
	for i=1,nSpec-1 do p=plot(x, densitySpec[*,nz/2,i], /over)
    
	p=plot(x, tempSpec[*,nz/2,0], title = 'z=0 Temp')
	for i=1,nSpec-1 do p=plot(x, tempSpec[*,nz/2,i], /over)

	p=plot(x, nuOmgSpec[*,nz/2,0], title = 'z=0 nuOmg')
	for i=1,nSpec-1 do p=plot(x, nuOmgSpec[*,nz/2,i], /over)


stop
end
