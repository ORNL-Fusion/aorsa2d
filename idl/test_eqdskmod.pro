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

	eqdsk2.br	= smooth ( eqdsk2.br, 5, /edge )
	eqdsk2.bz	= smooth ( eqdsk2.bz, 5, /edge )

	nX	= 512 
	nY	= 512 
	brSave	= conGrid ( eqdsk2.br, nX, nY, cubic = -0.5 )
	bzSave	= conGrid ( eqdsk2.bz, nX, nY, cubic = -0.5 )
	bPhiSave	= conGrid ( eqdsk2.bPhi, nX, nY, cubic = -0.5 )
	rSave	= conGrid ( eqdsk2.r, nX, cubic = -0.5 )
	zSave	= conGrid ( eqdsk2.z, nY, cubic = -0.5 )

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




!p.multi = 0
stop
end
