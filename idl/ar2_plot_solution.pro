pro ar2_plot_solution, full = full, $
		scale = scale, $
		sav = sav, NoRunData=NoRunData, $
		nPhi = _nPhi, RHS=_RHS, $ 
        sumSpecies = sumSpecies

    if keyword_set(scale) then ScaleFac = scale else ScaleFac = 1
	if keyword_set(_RHS) then ThisRHS = _RHS else ThisRHS = 1

	ar2Input = ar2_read_namelist()
	ar2InputFile = ar2Input['AR2InputFileName']

	ThisGridNo = 1
	GridNoStr = string(ThisGridNo,format='(i3.3)')

	if keyword_set(_nPhi) then ThisNPhi = _nPhi else ThisNPhi = ar2Input['nPhi']

	nPhiStr = string(ThisNPhi,format='(i+4.3)')
	rhsStr = string(ThisRHS,format='(i6.6)')

	SolutionFile = 'output/solution_'+GridNoStr+'_'+nPhiStr+'_'+rhsStr+'.nc'
	RunDataFile = 'output/runData_'+GridNoStr+'_'+nPhiStr+'_'+rhsStr+'.nc'

	print, 'ar2InputFile: ', ar2InputFile
	print, 'SolutionFile: ', SolutionFile
	print, 'RunDataFile: ', RunDataFile

	ar2_read_ar2input, ar2InputFile, $
			rLim=rLim,zLim=zLim,LimMask=LimMask,ar2=ar2,$
            rlcfs=rlcfs, zlcfs=zlcfs

	@constants

	xAll = 0.0
	yAll = 0.0
	nAll = 0
	eAlphaAll = complex (0,0)
	eBetaAll = complex (0,0)
	eBAll = complex (0,0)
	
	cdfId = ncdf_open ( RunDataFile, /noWrite ) 
		nCdf_varGet, cdfId, 'nPhi', nPhi 
		nCdf_varGet, cdfId, 'freq', freq 
		nCdf_varGet, cdfId, 'capR', x 
		nCdf_varGet, cdfId, 'y', y 
		nCdf_varGet, cdfId, 'jr_re', jr_re 
		nCdf_varGet, cdfId, 'jr_im', jr_im 
		nCdf_varGet, cdfId, 'jt_re', jt_re 
		nCdf_varGet, cdfId, 'jt_im', jt_im 
		nCdf_varGet, cdfId, 'jz_re', jz_re 
		nCdf_varGet, cdfId, 'jz_im', jz_im 
		nCdf_varGet, cdfId, 'bmod', bmod
		nCdf_varGet, cdfId, 'brU', brU
		nCdf_varGet, cdfId, 'btU', btU
		nCdf_varGet, cdfId, 'bzU', bzU
		nCdf_varGet, cdfId, 'densitySpec', densitySpec
		nCdf_varGet, cdfId, 'tempSpec', tempSpec
		nCdf_varGet, cdfId, 'xx_re', xx_re 
		nCdf_varGet, cdfId, 'xx_im', xx_im
		nCdf_varGet, cdfId, 'yy_re', yy_re 
		nCdf_varGet, cdfId, 'yy_im', yy_im
		nCdf_varGet, cdfId, 'drBfn_re', dRbFn_bFn_re
		nCdf_varGet, cdfId, 'drBfn_im', dRbFn_bFn_im
		nCdf_varGet, cdfId, 'dzBfn_re', dzbFn_bFn_re
		nCdf_varGet, cdfId, 'dzBfn_im', dzbFn_bFn_im
		nCdf_varGet, cdfId, 'nuOmg', nuOmg
		nCdf_varGet, cdfId, 'LimMask', LimMask
        nCdf_varGet, cdfId, 'nZ_1D', nz_1D
	ncdf_close, cdfId

	xx	= complex ( xx_re, xx_im )
	yy	= complex ( yy_re, yy_im )

	jA_r = complex ( jr_re[*,*], jr_im[*,*] )
	jA_t = complex ( jt_re[*,*], jt_im[*,*] )
	jA_z = complex ( jz_re[*,*], jz_im[*,*] )

	dRbFn_bFn	= complex ( dRbFn_bFn_re, dRbFn_bFn_im )
	dZbFn_bFn	= complex ( dZbFn_bFn_re, dZbFn_bFn_im )

	nN	= n_elements ( xx[*,0] )
	nM	= n_elements ( yy[*,0] )

	cdfId = ncdf_open ( SolutionFile, /noWrite ) 

		nCdf_varGet, cdfId, 'r', r
		nCdf_varGet, cdfId, 'z', z
		nCdf_varGet, cdfId, 'nPhi', nPhi 
		nCdf_varGet, cdfId, 'freq', freq 
	
		nCdf_varGet, cdfId, 'ealpha_re', ealpha_re 
		ealpha_re = temporary(ealpha_re[*,*])
		nCdf_varGet, cdfId, 'ebeta_re', ebeta_re 
		ebeta_re = temporary(ebeta_re[*,*])
		nCdf_varGet, cdfId, 'eB_re', eB_re 
		eb_re = temporary(eb_re[*,*])

		nCdf_varGet, cdfId, 'ealpha_im', ealpha_im 
		ealpha_im = temporary(ealpha_im[*,*])
		nCdf_varGet, cdfId, 'ebeta_im', ebeta_im 
		ebeta_im = temporary(ebeta_im[*,*])
		nCdf_varGet, cdfId, 'eB_im', eB_im 
		eb_im = temporary(eb_im[*,*])

		nCdf_varGet, cdfId, 'ealphak_re', ealphak_re 
		ealphak_re = temporary(ealphak_re[*,*])
		nCdf_varGet, cdfId, 'ebetak_re', ebetak_re 
		ebetak_re = temporary(ebetak_re[*,*])
		nCdf_varGet, cdfId, 'eBk_re', eBk_re 
		ebk_re = temporary(ebk_re[*,*])

		nCdf_varGet, cdfId, 'ealphak_im', ealphak_im 
		ealphak_im = temporary(ealphak_im[*,*])
		nCdf_varGet, cdfId, 'ebetak_im', ebetak_im 
		ebetak_im = temporary(ebetak_im[*,*])
		nCdf_varGet, cdfId, 'eBk_im', eBk_im 
		ebk_im = temporary(ebk_im[*,*])

		nCdf_varGet, cdfId, 'er_re', er_re 
		er_re = temporary(er_re[*,*])
		nCdf_varGet, cdfId, 'et_re', et_re 
		et_re = temporary(et_re[*,*])
		nCdf_varGet, cdfId, 'ez_re', ez_re 
		ez_re = temporary(ez_re[*,*])

		nCdf_varGet, cdfId, 'er_im', er_im 
		er_im = temporary(er_im[*,*])
		nCdf_varGet, cdfId, 'et_im', et_im 
		et_im = temporary(et_im[*,*])
		nCdf_varGet, cdfId, 'ez_im', ez_im 
		ez_im = temporary(ez_im[*,*])

		nCdf_varGet, cdfId, 'jalpha_re', jalpha_re 
		nCdf_varGet, cdfId, 'jbeta_re', jbeta_re 
		nCdf_varGet, cdfId, 'jB_re', jB_re 

		nCdf_varGet, cdfId, 'jalpha_im', jalpha_im 
		nCdf_varGet, cdfId, 'jbeta_im', jbeta_im 
		nCdf_varGet, cdfId, 'jB_im', jB_im 

		nCdf_varGet, cdfId, 'jP_r_re', jPr_re 
		nCdf_varGet, cdfId, 'jP_t_re', jPt_re 
		nCdf_varGet, cdfId, 'jP_z_re', jPz_re 
	
		nCdf_varGet, cdfId, 'jP_r_im', jPr_im 
		nCdf_varGet, cdfId, 'jP_t_im', jPt_im 
		nCdf_varGet, cdfId, 'jP_z_im', jPz_im 
	
		nCdf_varGet, cdfId, 'jouleHeating', jouleHeating 

		ealpha	= complex ( ealpha_re, ealpha_im )
		ebeta	= complex ( ebeta_re, ebeta_im )
		eb	= complex ( eb_re, eb_im )

		ealphak	= complex ( ealphak_re, ealphak_im )
		ebetak	= complex ( ebetak_re, ebetak_im )
		ebk	= complex ( ebk_re, ebk_im )

		e_r	= complex ( er_re, er_im )
		e_t	= complex ( et_re, et_im )
		e_z	= complex ( ez_re, ez_im )

		jP_r_s	= complex ( jPr_re, jPr_im )
		jP_t_s	= complex ( jPt_re, jPt_im )
		jP_z_s	= complex ( jPz_re, jPz_im )

        jP_r = total(jP_r_s,3)
        jP_t = total(jP_t_s,3)
        jP_z = total(jP_z_s,3)

		jP_a_s = complex ( jAlpha_re, jAlpha_im )
		jP_b_s = complex ( jBeta_re, jBeta_im )
		jP_p_s = complex ( jB_re, jB_im )

        jP_a = total(jP_a_s,3)
        jP_b = total(jP_b_s,3)
        jP_p = total(jP_p_s,3)
        
	ncdf_close, cdfId


    w = 2*!pi*freq
    wRFc = complex(w+nuOmg*0,w*nuOmg)
	x = r
	y = z
    dr = r[1]-r[0]
	nx = n_elements(r)
	ny = n_elements(z)
	x2D	= rebin ( r, nX, nY )
	if nY gt 1 then begin
		y2D = transpose(rebin ( z, nY, nX ))
	endif else begin
		y2D = fltArr(1)+z
	endelse

	nSpec	= n_elements ( jp_r_s[0,0,*] )

	; jDotE
	; -----

	this_jDotE	= jp_r_s[*,*,0] * conj(e_r) $
				+ jp_t_s[*,*,0] * conj(e_t) $
				+ jp_z_s[*,*,0] * conj(e_z)

    jDotE_s = jP_r_s*0
    jDotE_s[*,*,0] = this_jDotE

	for s=1,nSpec-1 do begin

		this_jDotE	= jp_r_s[*,*,s] * conj(e_r) $
				+ jp_t_s[*,*,s] * conj(e_t) $
				+ jp_z_s[*,*,s] * conj(e_z)

        jDotE_s[*,*,s] = this_jDotE

	endfor

    for s=0,nSpec-1 do begin
        print, 'amu: ',round(ar2.amu[s]),' Z: ',round(ar2.atomicZ[s]),' Total(jDotE)['+string(s,format='(i1.1)')+']: ',$
            string(total(real_part(jDotE_s[*,*,s]))/total(real_part(jDotE_s))*100,format='(f4.1)')+'%'
    endfor


	; 1D catch
	; --------

	@load_colors

	; Plot all grid solutions
	; -----------------------

	if size(y,/dim) eq 0 then begin

		iix = n_elements(x)/2
		wrf = freq * 2 * !pi

		eRange = max(abs([e_r]))
		p_r = plot ( x, e_r, layout=[1,3,1],$
				title='Er',$
                yRange=[-eRange,eRange],$
                ytitle='Er [V/m]',$
                name='Re',$
                window_title='aorsa')
		p_i = plot ( x, imaginary(e_r), color='red',/over,name='Im')

		eRange = max(abs([e_t]))
		p_r = plot ( x, e_t, layout=[1,3,2],/current,$
				title='Et',$
                yRange=[-eRange,eRange],$
                ytitle='Et [V/m]',$
                name='Re')
		p_i = plot ( x, imaginary(e_t), color='red',/over,name='Im')

    	eRange = max(abs([e_z]))
		p_r = plot ( x, e_z, layout=[1,3,3],/current,$
				title='Ez',$
                yRange=[-eRange,eRange],$
                ytitle='Ez [V/m]',$
                name='Re')
		p_i = plot ( x, imaginary(e_z), color='red',/over,name='Im')

		p = plot (x,(total(JouleHeating,3))[*],color='b',thick=1,$
				title='J dot E',name='jP . E (Total)',font_size=10,$
				layout=[1,3,1],window_title='aorsa')
		colorArr = ['b','g','r','c','m','y','k']
        for s=0,n_elements(JouleHeating[0,0,*])-1 do begin
            p = plot(x,JouleHeating[*,0,s],color=colorArr[s],/over,thick=3,linestyle='--',transparency=50)
        endfor
		p_array = !NULL
		p_array = [p_array,p]

		for s=1,nSpec-1 do begin
			p = plot ( x, jDotE_s[*,*,0], /over,title='jDotE_0' )
			p_array = [p_array,p]
		endfor

        p = plot(x,jA_r,thick=2,/current,layout=[1,3,2],title='jA components')
        p = plot(x,jA_t,thick=2,/over)
        p = plot(x,jA_z,thick=2,/over)

        p = plot(x,imaginary(jA_r),thick=4,/over, lineStyle="--")
        p = plot(x,imaginary(jA_t),thick=4,/over, lineStyle="--")
        p = plot(x,imaginary(jA_z),thick=4,/over, lineStyle="--")


        nuOmgRange = [0,max(nuOmg)+0.02]
        p = plot(x,nuOmg,/current,title='nuOmg',layout=[1,3,3],yRange=nuOmgRange)

		; jP
		; --

    	yRange = max(abs(jp_r))

        plotJp_rtz = 1
    	if plotJp_rtz then begin
    
            nS = n_elements(jP_r_s[0,0,*])
   
            p=plot([0,0],[0,0],/noData, layout=[nS,3,1],window_title='aorsa_1d - jP_rtz',dimensions=[1200,800])
    
            for s=0,nS-1 do begin
                This_amu_str = ', amu: '+string(ar2.amu[s],format='(i1.1)')
                This_Z_str = ', Z: '+string(ar2.atomicZ[s],format='(i+2.1)')

    		    p_re = plot (r,jP_r_s[*,0,s],thick=2,$
    		    		title='jP_r'+This_amu_str+This_Z_str,name='Jp_re',font_size=10,$
    		    		layout=[nS,3,1+s],yRange=yRange,transparency=50,/current)
    		    p_im = plot (r,imaginary(jP_r_s[*,0,s]),color='r',thick=2,transparency=50,$
    		    		name='Ja_re',font_size=10,/over)
    
    		    p_re = plot (r,jP_t_s[*,0,s],thick=2,$
    		    		title='jP_t',name='Jp_re',font_size=10,$
    		    		layout=[nS,3,1+1*nS+s],/current,yRange=yRange,transparency=50)
    		    p_im = plot (r,imaginary(jP_t_s[*,0,s]),color='r',thick=2,transparency=50,$
    		    		name='Jp_re',font_size=10,/over)
    
    		    p_re = plot (r,jP_z_s[*,0,s],thick=2,$
    		    		title='jP_z',name='Jp_re',font_size=10,$
    		    		layout=[nS,3,1+2*nS+s],/current,yRange=yRange,transparency=50)
    		    p_im = plot (r,imaginary(jP_z_s[*,0,s]),color='r',thick=2,transparency=50,$
    		    		name='Jp_re',font_size=10,/over)
            endfor
    
            p=plot(r,jP_a,layout=[1,3,1], title='jP_rtz (summed over species)')
            p=plot(r,imaginary(jP_r),color='r',/over)
    
            p=plot(r,jP_t,layout=[1,3,2], /current)
            p=plot(r,imaginary(jP_t),color='r',/over)
    
            p=plot(r,jP_z,layout=[1,3,3], /current)
            p=plot(r,imaginary(jP_z),color='r',/over)
    
    	endif


        plotJp_abp = 1
    	if plotJp_abp then begin
    
            nS = n_elements(jP_r_s[0,0,*])
   
            p=plot([0,0],[0,0],/noData, layout=[nS,3,1],window_title='aorsa_1d - jP_abp',dimensions=[1200,800])
    
            for s=0,nS-1 do begin
                This_amu_str = ', amu: '+string(ar2.amu[s],format='(i1.1)')
                This_Z_str = ', Z: '+string(ar2.atomicZ[s],format='(i+2.1)')

    		    p_re = plot (r,jP_a_s[*,0,s],thick=2,$
    		    		title='jP_a'+This_amu_str+This_Z_str,name='Jp_re',font_size=10,$
    		    		layout=[nS,3,1+s],yRange=yRange,transparency=50,/current)
    		    p_im = plot (r,imaginary(jP_a_s[*,0,s]),color='r',thick=2,transparency=50,$
    		    		name='Ja_re',font_size=10,/over)
    
    		    p_re = plot (r,jP_b_s[*,0,s],thick=2,$
    		    		title='jP_b',name='Jp_re',font_size=10,$
    		    		layout=[nS,3,1+1*nS+s],/current,yRange=yRange,transparency=50)
    		    p_im = plot (r,imaginary(jP_b_s[*,0,s]),color='r',thick=2,transparency=50,$
    		    		name='Jp_re',font_size=10,/over)
    
    		    p_re = plot (r,jP_p_s[*,0,s],thick=2,$
    		    		title='jP_p',name='Jp_re',font_size=10,$
    		    		layout=[nS,3,1+2*nS+s],/current,yRange=yRange,transparency=50)
    		    p_im = plot (r,imaginary(jP_p_s[*,0,s]),color='r',thick=2,transparency=50,$
    		    		name='Jp_re',font_size=10,/over)
            endfor
    
            p=plot(r,jP_a,layout=[1,3,1], title='jP_abp (summed over species)')
            p=plot(r,imaginary(jP_a),color='r',/over)
    
            p=plot(r,jP_t,layout=[1,3,2], /current)
            p=plot(r,imaginary(jP_b),color='r',/over)
    
            p=plot(r,jP_z,layout=[1,3,3], /current)
            p=plot(r,imaginary(jP_p),color='r',/over)
    
    	endif



        ; Calculate the H vector & the Poynting vector

        II = complex(0,1)
       
        h_r = complexArr(nX)
        h_t = complexArr(nX)
        h_z = complexArr(nX)

        k_z = nZ_1D*0 ; this is NOT right, fix for non-zero nZ

        for i=2,nX-3 do begin

            dEz_dr = (1.0/12.0*e_z[i-2] - 2.0/3.0*e_z[i-1]$
                    +2.0/3.0*e_z[i+1] - 1.0/12.0*e_z[i+2])/dr
            drEt_dr = (1.0/12.0*r[i-2]*e_t[i-2] - 2.0/3.0*r[i-1]*e_t[i-1]$
                    +2.0/3.0*r[i+1]*e_t[i+1] - 1.0/12.0*r[i+2]*e_t[i+2])/dr

            h_r[i] = -II*k_z*e_t[i] + II*nPhi*e_z[i]/r[i]
            h_t[i] = II*k_z*e_r[i] - dEz_dr 
            h_z[i] = (-II*nPhi*e_r[i] + drEt_dr )/r[i]

        endfor

        h_r = h_r / (II*wRFc*u0)
        h_t = h_t / (II*wRFc*u0)
        h_z = h_z / (II*wRFc*u0)

        p=plot(r,h_r,layout=[1,3,1],title='h_r',window_title='aorsa')
        p=plot(r,imaginary(h_r),/over,color='r')
        p=plot(r,h_t,layout=[1,3,2],/current,title='h_t')
        p=plot(r,imaginary(h_t),/over,color='r')
        p=plot(r,h_z,layout=[1,3,3],/current,title='h_z')
        p=plot(r,imaginary(h_z),/over,color='r')

        ; Time average Poynting vector
        S_r = e_t*conj(h_z)-e_z*conj(h_t)
        S_t = -(e_r*conj(h_z)-e_z*conj(h_r))
        S_z = e_r*conj(h_t)-e_t*conj(h_r)

        S_r = 0.5*real_part(S_r)
        S_t = 0.5*real_part(S_t)
        S_z = 0.5*real_part(S_z)

        SRange = max(abs([S_r,S_t,S_z]))
        p=plot(r,S_r,layout=[1,3,1],title='Time average Poynting vector S_r',$
                window_title='aorsa',yRange=[-SRange,SRange])
        p=plot(r,S_t,layout=[1,3,2],/current,title='S_t',yRange=[-SRange,SRange])
        p=plot(r,S_z,layout=[1,3,3],/current,title='S_z',yRange=[-SRange,SRange])

	endif else begin ; Now 2D plotting

        !x.margin = !x.margin / 3
        !y.margin = !y.margin / 2
        dimensions = get_screen_size()*0.8
        plotpos = 1
        Layout = [3,2]
        aspect = 0
        xRange = [min(r),max(r)]
        yRange = [min(z),max(z)]

        LimColor = 'black'
        AntColor = 'black'
        AntCLevel = 0.1
        lcfsColor = 'black'
        cutoffColor = 'black'

        nPhiString = ' (nPhi: '+string(ThisNPHI,format='(i+4.3)')+')'

		nLevs = 10

		thisField = e_r[*,*]
        scale = max(abs(thisField))*scaleFac
        title = 'E_r'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
        aspect = 1.0
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
            aspect_ratio=aspect, dimensions=dimensions, title=title, xRange=xRange, yRange=yRange, layout=[[Layout,PlotPos]] )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
        p = plot ( rlim, zlim, /over, thick = 3, color=LimColor)
        p = plot ( rlcfs, zlcfs, /over, thick = 4, color=lcfsColor )
        ;c_zero_set = contour(ar2.kPerSq_F,ar2.r,ar2.z,/over,c_value=0.001,color=cutoffColor,C_LABEL_SHOW=0,c_thick=3)
        ;c_jant = contour(jA_z, x, y, /over, c_value = AntCLevel, color=AntColor, C_LABEL_SHOW=0, c_thick=5)
        ++PlotPos
        p.Save, "Er.png", border=0, height=600

		thisField = abs(e_r[*,*])
        scale = max(abs(thisField))*scaleFac
        title = 'abs(E_r)'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
            aspect_ratio=aspect, title=title, xRange=xRange, yRange=yRange, $
            layout=[[Layout,PlotPos]],/current)
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
        p = plot ( rlim, zlim, /over, thick = 3, color=LimColor)
        p = plot ( rlcfs, zlcfs, /over, thick = 4, color=lcfsColor )
        c_zero_set = contour(ar2.kPerSq_F,ar2.r,ar2.z,/over,c_value=0.001,color=cutoffColor,C_LABEL_SHOW=0,c_thick=3)
        c_jant = contour(jA_r, x, y, /over, c_value = AntCLevel, color=AntColor, C_LABEL_SHOW=0, c_thick=5)
        ++PlotPos
        p.Save, "absEr.png",  border=0, height=600

		nLevs = 50

		thisField = real_part(jDotE_s[*,*,0])
        title = 'jDotE (e)'
        scale = max(abs(thisField))*ScaleFac
        levels = fIndGen(nLevs)/(nLevs-1)*scale
        colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
            aspect_ratio=aspect, title=title, xRange=xRange, yRange=yRange, $
            font_size=16, Layout=[[Layout,PlotPos]],/current)
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
        p = plot ( rlim, zlim, /over, thick = 3, color=LimColor)
        p = plot ( rlcfs, zlcfs, /over, thick = 4, color=lcfsColor )
        c_zero_set = contour(ar2.kPerSq_F,ar2.r,ar2.z,/over,c_value=0.001,color=cutoffColor,C_LABEL_SHOW=0,c_thick=3)
        c_jant = contour(jA_r, x, y, /over, c_value = AntCLevel, color=AntColor, C_LABEL_SHOW=0, c_thick=5)
        ++PlotPos
        p.Save, "jDotE_e.png",  border=0, height=600

		thisField = nuOmg[*,*,0]
        title = 'nuOmg (e)'
        scale = max(abs(thisField))*ScaleFac
        levels = fIndGen(nLevs)/(nLevs-1)*scale
        colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
            aspect_ratio=aspect, title=title, xRange=xRange, yRange=yRange, $
            font_size=16, Layout=[[Layout,PlotPos]],/current)
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
        p = plot ( rlim, zlim, /over, thick = 3, color=LimColor)
        p = plot ( rlcfs, zlcfs, /over, thick = 4, color=lcfsColor )
        c_zero_set = contour(ar2.kPerSq_F,ar2.r,ar2.z,/over,c_value=0.001,color=cutoffColor,C_LABEL_SHOW=0,c_thick=3)
        c_jant = contour(jA_r, x, y, /over, c_value = AntCLevel, color=AntColor, C_LABEL_SHOW=0, c_thick=5)
        ++PlotPos
        p.Save, "nuOmg_e.png",  border=0, height=600


	endelse

	if keyword_set(sav) then begin

		xorig = x
		eAlphaOrig = eAlpha
		eBetaOrig = eBeta
		eBOrig = eB

		save, xorig, eAlphaOrig, eBetaOrig, eBOrig, $
				file = 'soln.sav'
	endif

	stop
end
