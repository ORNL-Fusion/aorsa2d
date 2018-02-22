function ar2_read_ar2input, runFolderName, $
		rLim=rLim, zLim=zLim, LimMask=LimMask, $
        rlcfs=rlcfs, zlcfs=zlcfs, fileName = _fileName

    if keyword_set(_fileName) then begin
        ar2InFileName = _fileName 
    endif else begin
        ar2InFileName = runFolderName + '/input/ar2Input.nc'
    endelse

	cdfId = ncdf_open ( ar2InFileName, /noWrite ) 
	nCdf_varGet, cdfid, 'rMin', rMin
	nCdf_varGet, cdfid, 'rMax', rMax
	nCdf_varGet, cdfid, 'zMin', zMin
	nCdf_varGet, cdfid, 'zMax', zMax

	nCdf_varGet, cdfid, 'r', r 
	nCdf_varGet, cdfid, 'z', z
	nCdf_varGet, cdfid, 'br', br
	nCdf_varGet, cdfid, 'bt', bt 
	nCdf_varGet, cdfid, 'bz', bz

	nCdf_varGet, cdfid, 'AtomicZ', AtomicZ 
	nCdf_varGet, cdfid, 'amu', amu 

	nCdf_varGet, cdfid, 'Density_m3', Density_m3
	nCdf_varGet, cdfid, 'Temp_eV', Temp_eV
	nCdf_varGet, cdfid, 'nuOmg', nuOmg

	nCdf_varGet, cdfid, 'LimMask', mask_lim 
	nCdf_varGet, cdfid, 'Lim_r', rlim 
	nCdf_varGet, cdfid, 'Lim_z', zlim 

    if ncdf_varid(cdfid,'rlcfs') lt 0 then begin
        print, 'WARNING: arrInput.nc file missing lcfs varible'
        print, 'Setting default'
        rlcfs = [0,5,5,0,0]
        zlcfs = [-2,-2,2,2,-2]
    endif else begin 
	    nCdf_varGet, cdfid, 'rlcfs', rlcfs 
	    nCdf_varGet, cdfid, 'zlcfs', zlcfs
    endelse

    if ncdf_varid(cdfid,'kPerSq_F') lt 0 then begin
        print, 'WARNING: arrInput.nc file missing kPerSq varible'
        print, 'Setting default'
        kPerSq_F = mask_lim*0 
        kPerSq_S = mask_lim*0 
    endif else begin
	    nCdf_varGet, cdfid, 'kPerSq_F', kPerSq_F 
	    nCdf_varGet, cdfid, 'kPerSq_S', kPerSq_S 
    endelse
	
	ncdf_close, cdfId

	;cdfId = ncdf_open ( runDataFileName, /noWrite ) 
	;	nCdf_varGet, cdfid, 'LimMask', LimMask
	;ncdf_close, cdfId

    return, ar2 = { $
            rMin: rMin, $
            rMax: rMax, $
            zMin: zMin, $
            zMax: zMax, $
            r: r, $
            z: z, $
            br: br, $
            bt: bt, $
            bz: bz, $
            AtomicZ: AtomicZ, $
            amu: amu, $
            Density_m3: Density_m3, $
            Temp_eV: Temp_eV, $
            nuOmg: nuOmg, $
            LimMask: mask_lim, $
            Lim_r: rlim, $
            Lim_z: zlim, $
            rlcfs: rlcfs, $
            zlcfs: zlcfs, $
            kPerSq_F: kPerSq_F, $
            kPerSq_S: kPerSq_S }

end
