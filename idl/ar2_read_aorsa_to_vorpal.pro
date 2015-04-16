pro ar2_read_aorsa_to_vorpal

	fileName = 'aorsaToVorpal_E.txt'
	n = file_lines(fileName)
	nHeader = 2
	n = n-nHeader

	data = replicate( {x:0.0,$
			ex_re:0.0, ex_im:0.0, $
			ey_re:0.0, ey_im:0.0, $
			ez_re:0.0, ez_im:0.0 }, n)

	openr, lun, fileName, /get_lun
	skip_lun, lun, nHeader, /lines
	readf, lun, data	
	close, lun

	p=plot(data.x,data.ex_re)
	p=plot(data.x,data.ex_im,color='r',/over)
stop
end
