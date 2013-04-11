pro ar2_read_namelist, ar2Input=ar2Input

	FileName = 'aorsa2d.in'

	openr, lun, FileName, /get_lun
	array = ''
	line = ''
	while not eof(lun) do begin
		readf, lun, line
	   	array = [array,line]	
	endwhile
	free_lun, lun
	
	ar2Input = hash ( $
			"AR2InputFileName", '', $
			"useAR2Input", 0, $
		    "nPhi",0, $
			"AR2SourceLocationsFileName", '', $
			"useAR2SourceLocationsFile", 0)

	foreach value, ar2Input, entryName do begin

		for i=0,n_elements(array)-1 do begin

			ThisLine = array[i]

			if StrMatch(ThisLine,'*'+entryName+'*',/fold_case) then begin

					if size(ar2Input[entryName],/type) eq 7 then begin
						pos = StRegEx(ThisLine,"'.*'",length=len)
						ThisLine = strtrim(strmid(ThisLine,pos+1,len-2),2)
						ar2Input[entryName] = ThisLine	
					endif

					if size(ar2Input[entryName],/type) eq 2 then begin
						pos = StRegEx(ThisLine,"=.*,",length=len)
						ThisLine = strtrim(strmid(ThisLine,pos+1,len-2),2)
						if strcmp(ThisLine,'.true.') then begin
								ThisInt = 1
						endif else if strcmp(ThisLine, '.false.') then begin
								ThisInt = 0
						endif else begin
								ThisInt = fix(ThisLine)	
						endelse
						ar2Input[entryName] = ThisInt
					endif
			endif

		endfor

	endforeach
end
