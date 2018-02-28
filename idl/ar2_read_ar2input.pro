function ar2_read_ar2input, runFolderName, $
        fileName = _fileName

    if keyword_set(_fileName) then begin
        ar2InFileName = _fileName 
    endif else begin
        ar2InFileName = runFolderName + '/input/ar2Input.nc'
    endelse

    ar2 = dlg_read_netcdf(ar2InFileName)

    return, ar2

end
