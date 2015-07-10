function read_vtk, fileName

	openr, lun, fileName, /get_lun
	nX = 500L
	nY = 320L
	skip_lun, lun, 10, /lines
	data=replicate({e:0.0},nX*nY)
	readf, lun, data
	skip_lun, lun, 2, /lines

	rMin = 0.947
	rMax = 2.298
	zMin = -1.2
	zMax = +1.2

	r = fIndGen(nX)/(nX-1)*(rMax-rMin)+rMin
	z = fIndGen(nY)/(nY-1)*(zMax-zMin)+zMin

	er = reform(data.e,nX,nY)
	et = er*0
	ez = er*0

	return, {r:r,z:z,er:er,et:et,ez:ez}

end
