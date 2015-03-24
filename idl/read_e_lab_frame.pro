function read_e_lab_frame, fileName

	openr, lun, fileName, /get_lun
	nX = 400L
	nY = 400L
	skip_lun, lun, 1, /lines
	datar=replicate({r:0.0},nX)
	readf, lun, datar
	dataz=replicate({z:0.0},nY)
	readf, lun, dataz

	r = datar.r
	z = dataz.z

	n=nX*nY
	data1=replicate({er:0.0,ei:0.0},n)	
	readf, lun, data1
	data2=replicate({er:0.0,ei:0.0},n)	
	readf, lun, data2
	data3=replicate({er:0.0,ei:0.0},n)	
	readf, lun, data3

	er = reform(complex(data1.er,data1.ei),nX,nY)
	et = reform(complex(data2.er,data2.ei),nX,nY)
	ez = reform(complex(data3.er,data3.ei),nX,nY)

	return, {r:r,z:z,er:er,et:et,ez:ez}

end
