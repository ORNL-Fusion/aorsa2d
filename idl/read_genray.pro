function read_genray, fileName

	openr, lun, fileName, /get_lun

	nLines = file_lines(fileName)

	data = replicate({distance:0.0,r:0.0,z:0.0,nPar:0.0,nPer:0.0,power:0.0,br:0.0,bz:0.0,bt:0.0},nLines)	
	readf, lun, data

	close, lun
	free_lun, lun	

	return, data
end


