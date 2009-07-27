pro plot_ql_nc

	if strCmp ( getEnv ( 'MACHINE' ), 'jaguar' ) then scratchDir = getEnv ( 'SYSTEM_USERDIR' )
	if strCmp ( getEnv ( 'MACHINE' ), 'franklin' ) then scratchDir = getEnv ( 'SCRATCH' )
	if strCmp ( getEnv ( 'MACHINE' ), 'dlghp' ) then scratchDir = '/home/dg6/scratch' 

	fileName	= 'output/p_ql.nc'

	cdfId	= ncdf_open ( fileName, /noWrite )
	glob	= ncdf_inquire ( cdfId )
	nCdf_varGet, cdfId, 'ql_b', b 
	;!nCdf_varGet, cdfId, 'ql_c', c 
	;!nCdf_varGet, cdfId, 'ql_e', e 
	;!nCdf_varGet, cdfId, 'ql_f', f 
;	nCdf_varGet, cdfId, 'maxV', maxV
;	nCdf_varGet, cdfId, 'minV', minV
	nCdf_varGet, cdfId, 'R', R
	nCdf_varGet, cdfId, 'z', z
	nCdf_varGet, cdfId, 'vPer', vPer
	nCdf_varGet, cdfId, 'vPar', vPar 
	nCdf_varGet, cdfId, 'vc_mks', vc_mks
	nCdf_varGet, cdfId, 'pscale', pscale 

	nR	= n_elements(b[*,0,0,0])
	nz	= n_elements(b[0,*,0,0])
	nvPer	= n_elements(b[0,0,*,0])
	nvPar	= n_elements(b[0,0,0,*])

	mi	= 2.0 * 1.673e-27
	e	= 1.602e-19
	print, 'max energy (keV): ', max(vPer)^2*0.5*mi/e/1e3

	;vPar4D	= transpose(rebin ( vPar, nvPar, nvPer, nZ, nR ))
	;vPer4D	= transpose(rebin ( vPer, nvPer, nvPar, nZ, nR ),[3,2,0,1]) 
	;vPer2D	= rebin ( vPer, nvPer, nvPar )
	;vPar2D	= transpose ( rebin ( vPar, nvPar, nvPer ) )
	;th4D	= atan ( vPer4D, vPar4D )

	;v4D	= sqrt ( vPar4D^2 + vPer4D^2 )

	;b 	= b / v4D^2
	;c	= c / v4D
	;e	= e / (v4D*sin(th4D))
	;f	= f / sin(th4D)

	;b_cyl	= 0.5 * ( b + f + (-b + f ) * cos ( 2d0 * th4D ) + ( c + e ) * sin ( 2d0 * th4D ) )
	;c_cyl	= 0.5 * ( -c + e + ( c + e ) * cos ( 2d0 * th4D ) + ( b -f ) * sin ( 2d0 * th4D ) ) 
	;e_cyl	= 0.5 * ( c - e + ( c + e ) * cos ( 2d0 * th4D ) + ( b - f ) * sin ( 2d0 * th4D ) )
	;f_cyl	= 0.5 * ( b + f + ( b - f ) * cos ( 2d0 * th4D ) - ( c + e ) * sin ( 2d0 * th4D ) )

	;b_cyl	= ( f * v4D * vPar4D^2 + vPer4D * ( -1d0 * ( c + e * v4D^2 ) * vPar4D + b * v4D * vPer4D ) ) / v4D^4 
	;c_cyl	= (-e * v4D^2 * vPar4D^2 + vPer4D * ( ( b - f ) * v4D * vPar4D + c * vPer4D ) ) / v4D^4
	;e_cyl	= (b * v4D * vPar4D^2 + vPer4D * ( ( c * vPar4D + v4D * ( e * v4D * vPar4D + f * vPer4D ) ) ) ) / v4D^4
	;f_cyl	= (b * v4D * vPar4D^2 + vPer4D * ( c * vPar4D + v4D * ( e * v4D * vPar4D + f * vPer4D ) ) ) / v4D^4 

	;b_cyl2D	= total ( total ( b_cyl, 1), 1 )
	;c_cyl2D	= total ( total ( c_cyl, 1), 1 )
	;e_cyl2D	= total ( total ( e_cyl, 1), 1 )
	;f_cyl2D	= total ( total ( f_cyl, 1), 1 )

	dvPerp	= abs(vPer[0]-vPer[1])
	dvPar	= abs(vPar[0]-vPar[1])

	;divD_per1	= dlg_pDeriv ( vPer2D * b_cyl2D, 1, dvPerp ) / vPer2D 
	;divD_per2	= dlg_pDeriv ( e_cyl2D, 2, dvPar )

	b2D	= total ( total ( b, 1), 1 )
	;c2D	= total ( total ( c, 1), 1 )
	;e2D	= total ( total ( e, 1), 1 )
	;f2D	= total ( total ( f, 1), 1 )

	old_dev = !D.name
	set_plot, 'ps'
	outfname	= 'output/b2D.eps'
	device, filename=outfname, preview=0, /color, bits_per_pixel=8, $
		xsize=14, ysize = 8,xoffset=.1, yoffset=.1, /encapsul
	!p.multi = 0
	loadct, 13, file='davect.tbl'
	levels	= (fIndGen(51)-25)/25 * 5d11
	colors	= bytScl ( fIndGen(51), top = 253 ) + 1
	plotI	= n_elements(b[*,0,0,0])/2
	plotJ	= n_elements(b[0,*,0,0])/2

	contour, abs(transpose (reform(b[plotI,plotJ,*,*]))), vPar*1d-6, vPer*1d-6, $
			lev = levels, $
			c_colo = colors, $
			/fill, $
			color = 255, $
			xTitle = 'vPar [m/s] x10!U6!N', $
			yTitle = 'vPer [m/s] x10!U6!N', $
			title = string(R[plotI])+string(z[plotJ])

	device, /close

	;contour, transpose (c2D ), vPar, vPer, $
	;		lev = levels, $
	;		c_colo = colors, $
	;		/fill, $
	;		color = 255, $
	;		title = 'b x10^5'
	;contour, transpose (e2D), vPar, vPer, $
	;		lev = levels, $
	;		c_colo = colors, $
	;		/fill, $
	;		color = 255, $
	;		title = 'b x10^5'
	;contour, transpose (f2D ), vPar, vPer, $
	;		lev = levels, $
	;		c_colo = colors, $
	;		/fill, $
	;		color = 255, $
	;		title = 'b x10^5'
	;
	;contour, transpose (b_cyl2D), vPar, vPer, $
	;		lev = levels, $
	;		c_colo = colors, $
	;		/fill, $
	;		color = 255 
	;contour, transpose (c_cyl2D ), vPar, vPer, $
	;		lev = levels, $
	;		c_colo = colors, $
	;		/fill, $
	;		color = 255
	;contour, transpose (e_cyl2D ), vPar, vPer, $
	;		lev = levels, $
	;		c_colo = colors, $
	;		/fill, $
	;		color = 255
	;contour, transpose (f_cyl2D ), vPar, vPer, $
	;		lev = levels, $
	;		c_colo = colors, $
	;		/fill, $
	;		color = 255

	;contour, transpose ( (c2D<max(levels))>min(levels) ), vPar, vPer, $
	;		lev = levels, $
	;		c_colo = colors, $
	;		/fill, $
	;		color = 0
	;contour, transpose ( (e2D<max(levels))>min(levels) ), vPar, vPer, $
	;		lev = levels, $
	;		c_colo = colors, $
	;		/fill, $
	;		color = 0
	;contour, transpose ( (f2D<max(levels))>min(levels) ), vPar, vPer, $
	;		lev = levels, $
	;		c_colo = colors, $
	;		/fill, $
	;		color = 0

stop
end
