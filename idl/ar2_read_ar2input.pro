pro ar2_read_ar2input, ar2InFileName, $
		rLim=rLim, zLim=zLim, LimMask=LimMask, ar2 = ar2, $
        rlcfs=rlcfs, zlcfs=zlcfs

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

	nCdf_varGet, cdfid, 'rlcfs', rlcfs 
	nCdf_varGet, cdfid, 'zlcfs', zlcfs

	nCdf_varGet, cdfid, 'kPerSq_F', kPerSq_F 
	nCdf_varGet, cdfid, 'kPerSq_S', kPerSq_S 
	
	ncdf_close, cdfId

	;cdfId = ncdf_open ( runDataFileName, /noWrite ) 
	;	nCdf_varGet, cdfid, 'LimMask', LimMask
	;ncdf_close, cdfId

    ar2 = { $
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
