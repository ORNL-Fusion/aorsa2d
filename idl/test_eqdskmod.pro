pro test_eqdskMod

fName1	= '~/data/nstx/eqdsk/g120740.00275.EFIT02.mds.uncorrected.qscale_1.00000'
fName2	= '~/data/nstx/eqdsk/g120740.00275.EFIT02.mds.uncorrected.qscale_1.00000.dlgLimFix'

eqdsk1	= readgeqdsk ( fName1 )
eqdsk2	= readgeqdsk ( fName2 )

window, 0, xSize = 1200, ySize = 900
device, decomposed = 0
!p.multi = [0,2,2]

veloVect, eqdsk1.br, eqdsk1.bz, eqdsk1.r, eqdsk1.z, length = 150 
loadct, 12
oPlot, eqdsk1.rbbbs, eqdsk1.zbbbs, color = 12*16-1, thick = 3
oPlot, eqdsk1.rlim, eqdsk1.zlim, color = 8*16-1, thick = 3



;	modify the polidal field for use with an extended limiter

	over_sample_boundary, eqdsk2.rlim[*], eqdsk2.zlim[*], rlim_os, zlim_os
	mask = is_in_shape ( eqdsk2.r2d, eqdsk2.z2d, rlim_os, zlim_os )
			

	iiOutside	= where ( mask ne 1 )
	eqdsk2.br[iiOutside]	= 0
	eqdsk2.bz[iiOutside]	= 0

	brNew	= fltArr ( size ( eqdsk2.br, /dim ) )
	bzNew	= fltArr ( size ( eqdsk2.br, /dim ) )

	for i = 0, n_elements ( brNew[*,0] ) - 1 do begin
		for j = 0, n_elements ( brNew[0,*] ) - 1 do begin
	
			if mask[i,j] ne 1 then begin

				distance	= sqrt ( (eqdsk1.r[i] - eqdsk1.r2d)^2 + (eqdsk2.z[j] - eqdsk2.z2d)^2 )
				iiNZ	= where ( sqrt ( eqdsk2.br^2 + eqdsk2.bz^2 ) gt 0 )
				iiUse	= where ( distance[iiNZ] eq min ( distance[iiNZ] ) )
				ratio	= eqdsk2.bPhi[i,j] / eqdsk2.bPhi[iiNZ[iiUse]]
				if eqdsk2.r[i] lt 0.4 then ratio = 1
				brNew[i,j]	= eqdsk2.br[iiNZ[iiUse]] * ratio
				bzNew[i,j]	= eqdsk2.bz[iiNZ[iiUse]] * ratio 

			endif

		endfor
	endfor

	eqdsk2.br	+= brNew
	eqdsk2.bz	+= bzNew

	;eqdsk2.br	= smooth ( eqdsk2.br, 5, /edge )
	;eqdsk2.bz	= smooth ( eqdsk2.bz, 5, /edge )

	nX	= 1024 
	nY	= 1024 

	;cubPar = -0.5 
	;brSave	= conGrid ( eqdsk2.br, nX, nY, cubic = cubPar, /minu )
	;bzSave	= conGrid ( eqdsk2.bz, nX, nY, cubic = cubPar, /minu  )
	;bPhiSave	= conGrid ( eqdsk2.bPhi, nX, nY, cubic = cubPar, /minu )
	;rSave	= conGrid ( eqdsk2.r, nX, cubic = cubPar, /minu  )
	;zSave	= conGrid ( eqdsk2.z, nY, cubic = cubPar, /minu  )

	rSave	= eqdsk2.rLeft + fIndGen ( nX ) * ( eqdsk2.rDim / ( nX - 1 ) )
	zSave	= min(eqdsk2.z) + fIndGen ( nY ) * ( eqdsk2.zDim / ( nY - 1 ) )

	rSave2D	= rebin ( rSave, nX, nY )
	zSave2D	= transpose ( rebin ( zSave, nY, nX ) )

	brSave = interpolate ( eqdsk2.br, ( rSave2D - eqdsk2.rleft ) / eqdsk2.rdim * (eqdsk2.nw-1), $
   			( zSave2D + eqdsk2.zmaxis - min ( eqdsk2.z ) ) / eqdsk2.zdim * (eqdsk2.nh-1) )
	bzSave = interpolate ( eqdsk2.bz, ( rSave2D - eqdsk2.rleft ) / eqdsk2.rdim * (eqdsk2.nw-1), $
   			( zSave2D + eqdsk2.zmaxis - min ( eqdsk2.z ) ) / eqdsk2.zdim * (eqdsk2.nh-1) )
	bPhiSave = interpolate ( eqdsk2.bPhi, ( rSave2D - eqdsk2.rleft ) / eqdsk2.rdim * (eqdsk2.nw-1), $
   			( zSave2D + eqdsk2.zmaxis - min ( eqdsk2.z ) ) / eqdsk2.zdim * (eqdsk2.nh-1) )


;	write this in a netCDF file for use in aorsa

	outFileName	= '/home/dg6/data/aorsa_ant/dlg_bField.nc'
	nc_id	= nCdf_create ( outFileName, /clobber )
	nCdf_control, nc_id, /fill
	
	nR_id	= nCdf_dimDef ( nc_id, 'nR', nX )
	nz_id	= nCdf_dimDef ( nc_id, 'nz', nY )
	scalar_id	= nCdf_dimDef ( nc_id, 'scalar', 1 )
	
	R_id = nCdf_varDef ( nc_id, 'R', [ nR_id ], /float )
	z_id = nCdf_varDef ( nc_id, 'z', [ nz_id ], /float )
	bR_id = nCdf_varDef ( nc_id, 'bR', [nR_id, nz_id], /float )
	bz_id = nCdf_varDef ( nc_id, 'bz', [nR_id, nz_id], /float )
	bPhi_id = nCdf_varDef ( nc_id, 'bPhi', [nR_id, nz_id], /float )

	nCdf_control, nc_id, /enDef
	
	nCdf_varPut, nc_id, R_id, rSave 
	nCdf_varPut, nc_id, z_id, zSave 
	nCdf_varPut, nc_id, bR_id, brSave 
	nCdf_varPut, nc_id, bz_id, bzSave 
	nCdf_varPut, nc_id, bPhi_id, bPhiSave 

	nCdf_close, nc_id



loadct, 0
veloVect, eqdsk2.br, eqdsk2.bz, eqdsk2.r, eqdsk2.z, length = 5 
loadct, 12
oPlot, eqdsk2.rbbbs, eqdsk2.zbbbs, color = 12*16-1, thick = 3
oPlot, eqdsk2.rlim, eqdsk2.zlim, color = 8*16-1, thick = 3

loadct, 0
veloVect, brNew, bzNew, eqdsk2.r, eqdsk2.z, length = 10 
loadct, 12
oPlot, eqdsk2.rbbbs, eqdsk2.zbbbs, color = 12*16-1, thick = 3
oPlot, eqdsk2.rlim, eqdsk2.zlim, color = 8*16-1, thick = 3

loadct, 0
veloVect, brSave, bzSave, rSave, zSave, length = 10 
loadct, 12
oPlot, eqdsk2.rbbbs, eqdsk2.zbbbs, color = 12*16-1, thick = 3
oPlot, eqdsk2.rlim, eqdsk2.zlim, color = 8*16-1, thick = 3

	cdfId = ncdf_open ( '/home/dg6/scratch/aorsa2d/nstx/bench/output/mchoi_dlg.nc', /noWrite ) 
	nCdf_varGet, cdfId, 'rho_pla', rho_pla 
	nCdf_varGet, cdfId, 'rho_ant', rho_ant
	nCdf_varGet, cdfId, 'antJ_x', antJ_x
	nCdf_varGet, cdfId, 'antJ_y', antJ_y
	nCdf_varGet, cdfId, 'antJ_z', antJ_z
	nCdf_varGet, cdfId, 'plaJ_y', plaJ_x
	nCdf_varGet, cdfId, 'plaJ_y', plaJ_y
	nCdf_varGet, cdfId, 'plaJ_y', plaJ_z
	nCdf_varGet, cdfId, 'antOmega', antOmega
	nCdf_varGet, cdfId, 'nPhi', nPhi 
	nCdf_varGet, cdfId, 'bxn', bxn 
	nCdf_varGet, cdfId, 'byn', byn 
	nCdf_varGet, cdfId, 'bzn', bzn 
	nCdf_varGet, cdfId, 'R', choiR 
	nCdf_varGet, cdfId, 'z', choiz 


	ncdf_close, cdfId

	choiR2D	= rebin ( choiR, n_elements(choiR), n_elements(choiz) )
	choiz2D	= transpose(rebin ( choiz, n_elements(choiz), n_elements(choiR) ))

	brSaveA = interpolate ( eqdsk2.br, ( choiR2D - eqdsk2.rleft ) / eqdsk2.rdim * (eqdsk2.nw-1), $
   			( choiz2D + eqdsk2.zmaxis - min ( eqdsk2.z ) ) / eqdsk2.zdim * (eqdsk2.nh-1) )
	bzSaveA = interpolate ( eqdsk2.bz, ( choiR2D - eqdsk2.rleft ) / eqdsk2.rdim * (eqdsk2.nw-1), $
   			( choiz2D + eqdsk2.zmaxis - min ( eqdsk2.z ) ) / eqdsk2.zdim * (eqdsk2.nh-1) )
	bPhiSaveA = interpolate ( eqdsk2.bPhi, ( choiR2D - eqdsk2.rleft ) / eqdsk2.rdim * (eqdsk2.nw-1), $
   			( choiz2D + eqdsk2.zmaxis - min ( eqdsk2.z ) ) / eqdsk2.zdim * (eqdsk2.nh-1) )




	window, 7
	bModSave	= sqrt ( brSave^2 + bzSave^2 + bPhiSave^2 )
	bModSaveA	= sqrt ( brSaveA^2 + bzSaveA^2 + bPhiSaveA^2 )

	loadct, 12
	plot, choiR, bxn[*,42]
	oplot, choiR, (brSaveA/bModSaveA)[*,42], color = 12*16-1
	;oplot, rSave, (bRSave/bModSave)[*,42] 
	plot, choiR, byn[*,42]
	oplot, choiR, (bzSaveA/bModSaveA)[*,42], color = 12*16-1
	;oplot, rSave, (bzSave/bModSave)[*,42] 
	plot, choiR, bzn[*,42], yRange = [0.9,1.0]
	oplot, choiR, -(bphiSaveA/bModSaveA)[*,42], color = 12*16-1 
	;oplot, rSave, -(bPhiSave/bModSave)[*,42] 


!p.multi = 0
stop
end
