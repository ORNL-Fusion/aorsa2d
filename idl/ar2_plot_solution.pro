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

	nPhiStr = string(ThisNPhi,format='(i+7.6)')
	rhsStr = string(ThisRHS,format='(i6.6)')

	SolutionFile = 'output/solution_'+GridNoStr+'_'+nPhiStr+'_'+rhsStr+'.nc'
	RunDataFile = 'output/runData_'+GridNoStr+'_'+nPhiStr+'_'+rhsStr+'.nc'

	print, 'ar2InputFile: ', ar2InputFile
	print, 'SolutionFile: ', SolutionFile
	print, 'RunDataFile: ', RunDataFile

	ar2 = ar2_read_ar2input( './', $
			rLim=rLim,zLim=zLim,LimMask=LimMask,$
            rlcfs=rlcfs, zlcfs=zlcfs)


	@dlg_constants

	xAll = 0.0
	yAll = 0.0
	nAll = 0
	E_alpAll = complex (0,0)
	E_betAll = complex (0,0)
	E_prlAll = complex (0,0)

    arR = dlg_read_netcdf( RunDataFile )

	nN	= n_elements ( (arR['xx'])[*,0] )
	nM	= n_elements ( (arR['yy'])[*,0] )

    arS = dlg_read_netcdf( SolutionFile )

    arS['jP_r_sum'] = total(arS['jP_r'],3)
    arS['jP_t_sum'] = total(arS['jP_t'],3)
    arS['jP_z_sum'] = total(arS['jP_z'],3)

    arS['jP_alp_sum'] = total(arS['jP_alp'],3)
    arS['jP_bet_sum'] = total(arS['jP_bet'],3)
    arS['jP_prl_sum'] = total(arS['jP_prl'],3)

    nuOmg = arR['nuOmg']
    freq = arS['freq']

    w = 2*!pi*freq
    wRFc = complex(w+nuOmg*0,w*nuOmg)
	x = arS['r']
	y = arS['z']
    r = x
    z = y
    dr = r[1]-r[0]
	nx = n_elements(r)
	ny = n_elements(z)
	x2D	= rebin ( r, nX, nY )
	if nY gt 1 then begin
		y2D = transpose(rebin ( z, nY, nX ))
	endif else begin
		y2D = fltArr(1)+z
	endelse

	nSpec	= n_elements ( (arS['jP_r'])[0,0,*] )

	; jDotE
	; -----

    jP_r = arS['jP_r']
    jP_t = arS['jP_t']
    jP_z = arS['jP_z']

    jP_alp = arS['jP_alp']
    jP_bet = arS['jP_bet']
    jP_prl = arS['jP_prl']

    jP_r_sum = arS['jP_r_sum']
    jP_t_sum = arS['jP_t_sum']
    jP_z_sum = arS['jP_z_sum']

    jP_alp_sum = arS['jP_alp_sum']
    jP_bet_sum = arS['jP_bet_sum']
    jP_prl_sum = arS['jP_prl_sum']

    jA_r = arR['jA_r']
    jA_t = arR['jA_t']
    jA_z = arR['jA_z']
    
    E_r = arS['E_r']
    E_t = arS['E_t']
    E_z = arS['E_z']

    E_alp = arS['E_alp']
    E_bet = arS['E_bet']
    E_prl = arS['E_prl']


    jouleHeating = arS['jouleHeating']

	this_jDotE	= jP_r[*,*,0] * conj(E_r) $
				+ jP_t[*,*,0] * conj(E_t) $
				+ jP_z[*,*,0] * conj(E_z)

    jDotE_s = arS['jP_r']*0
    jDotE_s[*,*,0] = this_jDotE

	for s=1,nSpec-1 do begin

		this_jDotE	= jp_r[*,*,s] * conj(e_r) $
				+ jp_t[*,*,s] * conj(e_t) $
				+ jp_z[*,*,s] * conj(e_z)

        jDotE_s[*,*,s] = this_jDotE

	endfor

    for s=0,nSpec-1 do begin
        print, 'amu: ',round(ar2.amu[s]),' Z: ',round(ar2.atomicZ[s]),' Total(jDotE)['+string(s,format='(i1.1)')+']: ',$
            string( total(real_part(jDotE_s[*,*,s])) / total(real_part(jDotE_s)) *100,format='(f15.4)')+'%'
    endfor


	; 1D catch
	; --------

	@load_colors

	; Plot all grid solutions
	; -----------------------

    

	if size(y,/dim) eq 0 then begin

        margin = [0.2,0.2,0.1,0.2]
        thick = 2

		iix = n_elements(x)/2
		wrf = freq * 2 * !pi

		eRange = max(abs([e_r]))
		p_r = plot ( x, e_r, layout=[1,3,1],$
				title='Er',$
                yRange=[-eRange,eRange],$
                ytitle='Er [V/m]',$
                name='Re',$
                window_title='aorsa',margin=margin,thick=thick)
		p_i = plot ( x, imaginary(e_r), color='red',/over,name='Im',thick=thick)

		eRange = max(abs([e_t]))
		p_r = plot ( x, e_t, layout=[1,3,2],/current,$
				title='Et',$
                yRange=[-eRange,eRange],$
                ytitle='Et [V/m]',$
                name='Re',margin=margin,thick=thick)
		p_i = plot ( x, imaginary(e_t), color='red',/over,name='Im',thick=thick)

    	eRange = max(abs([e_z]))
		p_r = plot ( x, e_z, layout=[1,3,3],/current,$
				title='Ez',$
                yRange=[-eRange,eRange],$
                ytitle='Ez [V/m]',$
                name='Re',margin=margin,thick=thick,xtitle="r [m]")
		p_i = plot ( x, imaginary(e_z), color='red',/over,name='Im',thick=thick)

		p_r = plot ( x, E_alp, layout=[1,3,1],$
				title='E_alp',$
                yRange=[-eRange,eRange],$
                ytitle='E_alp [V/m]',$
                name='Re',$
                window_title='aorsa',margin=margin,thick=thick)
		p_i = plot ( x, imaginary(E_alp), color='red',/over,name='Im',thick=thick)

		eRange = max(abs([E_bet]))
		p_r = plot ( x, E_bet, layout=[1,3,2],/current,$
				title='E_bet',$
                yRange=[-eRange,eRange],$
                ytitle='E_bet [V/m]',$
                name='Re',margin=margin,thick=thick)
		p_i = plot ( x, imaginary(E_bet), color='red',/over,name='Im',thick=thick)

    	eRange = max(abs([E_prl]))
		p_r = plot ( x, E_prl, layout=[1,3,3],/current,$
				title='E_prl',$
                yRange=[-eRange,eRange],$
                ytitle='E_prl [V/m]',$
                name='Re',margin=margin,thick=thick,xtitle="r [m]")
		p_i = plot ( x, imaginary(E_prl), color='red',/over,name='Im',thick=thick)

        p_r.save, 'plot1.png', res=300


		p = plot (x,(total(JouleHeating,3))[*],color='b',thick=1,$
				title='J dot E',name='jP . E (Total)',font_size=10,$
				layout=[1,3,1],window_title='aorsa',xtitle="r [m]")
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
        current = 0
    	if plotJp_rtz then begin
    
            nS = n_elements(jP_r[0,0,*])
  

            for s=0,nS-1 do begin
                This_amu_str = ', amu: '+string(ar2.amu[s],format='(i1.1)')
                This_Z_str = ', Z: '+string(ar2.atomicZ[s],format='(i+2.1)')

    		    p_re = plot (r,jP_r[*,0,s],thick=2,window_title='aorsa_1d - jP_rtz',dimensions=[1200,800],$
    		    		title='jP_r'+This_amu_str+This_Z_str,name='Jp_re',font_size=10,$
    		    		layout=[nS,3,1+s],yRange=yRange,transparency=50,current=current<1)
    		    p_im = plot (r,imaginary(jP_r[*,0,s]),color='r',thick=2,transparency=50,$
    		    		name='Ja_re',font_size=10,/over)

                current++
    		    p_re = plot (r,jP_t[*,0,s],thick=2,$
    		    		title='jP_t',name='Jp_re',font_size=10,$
    		    		layout=[nS,3,1+1*nS+s],current=current<1,yRange=yRange,transparency=50)
    		    p_im = plot (r,imaginary(jP_t[*,0,s]),color='r',thick=2,transparency=50,$
    		    		name='Jp_re',font_size=10,/over)
    
    		    p_re = plot (r,jP_z[*,0,s],thick=2,$
    		    		title='jP_z',name='Jp_re',font_size=10,$
    		    		layout=[nS,3,1+2*nS+s],current=current<1,yRange=yRange,transparency=50)
    		    p_im = plot (r,imaginary(jP_z[*,0,s]),color='r',thick=2,transparency=50,$
    		    		name='Jp_re',font_size=10,/over)

            endfor
    
            p=plot(r,jP_r_sum,layout=[1,3,1], title='jP_rtz (summed over species)')
            p=plot(r,imaginary(jP_r_sum),color='r',/over)
    
            p=plot(r,jP_t_sum,layout=[1,3,2], /current)
            p=plot(r,imaginary(jP_t_sum),color='r',/over)
    
            p=plot(r,jP_z_sum,layout=[1,3,3], /current)
            p=plot(r,imaginary(jP_z_sum),color='r',/over)
    
    	endif


        plotJp_abp = 1
    	if plotJp_abp then begin
    
            nS = n_elements(jP_r[0,0,*])
   
            p=plot([0,0],[0,0],/noData, layout=[nS,3,1],window_title='aorsa_1d - jP_abp',dimensions=[1200,800])
    
            for s=0,nS-1 do begin
                This_amu_str = ', amu: '+string(ar2.amu[s],format='(i1.1)')
                This_Z_str = ', Z: '+string(ar2.atomicZ[s],format='(i+2.1)')

    		    p_re = plot (r,jP_alp[*,0,s],thick=2,$
    		    		title='jP_alp'+This_amu_str+This_Z_str,name='Jp_re',font_size=10,$
    		    		layout=[nS,3,1+s],yRange=yRange,transparency=50,/current)
    		    p_im = plot (r,imaginary(jP_alp[*,0,s]),color='r',thick=2,transparency=50,$
    		    		name='Ja_re',font_size=10,/over)
    
    		    p_re = plot (r,jP_bet[*,0,s],thick=2,$
    		    		title='jP_bet',name='Jp_re',font_size=10,$
    		    		layout=[nS,3,1+1*nS+s],/current,yRange=yRange,transparency=50)
    		    p_im = plot (r,imaginary(jP_bet[*,0,s]),color='r',thick=2,transparency=50,$
    		    		name='Jp_re',font_size=10,/over)
    
    		    p_re = plot (r,jP_prl[*,0,s],thick=2,$
    		    		title='jP_prl',name='Jp_re',font_size=10,$
    		    		layout=[nS,3,1+2*nS+s],/current,yRange=yRange,transparency=50)
    		    p_im = plot (r,imaginary(jP_prl[*,0,s]),color='r',thick=2,transparency=50,$
    		    		name='Jp_re',font_size=10,/over)
            endfor
    
            p=plot(r,jP_alp_sum,layout=[1,3,1], title='jP_abp (summed over species)')
            p=plot(r,imaginary(jP_alp_sum),color='r',/over)
    
            p=plot(r,jP_bet_sum,layout=[1,3,2], /current)
            p=plot(r,imaginary(jP_bet_sum),color='r',/over)
    
            p=plot(r,jP_prl_sum,layout=[1,3,3], /current)
            p=plot(r,imaginary(jP_prl_sum),color='r',/over)
    
    	endif



        ; Calculate the H vector & the Poynting vector

        II = complex(0,1)
       
        h_r = complexArr(nX)
        h_t = complexArr(nX)
        h_z = complexArr(nX)

        kz_1D = arR['kz_1d']
        nPhi = arS['nPhi']
        k_z = kz_1D*0 ; this is NOT right, fix for non-zero nZ

        for i=2,nX-3 do begin

            dEz_dr = (1.0/12.0*e_z[i-2] - 2.0/3.0*e_z[i-1]$
                    +2.0/3.0*e_z[i+1] - 1.0/12.0*e_z[i+2])/dr
            drEt_dr = (1.0/12.0*r[i-2]*e_t[i-2] - 2.0/3.0*r[i-1]*e_t[i-1]$
                    +2.0/3.0*r[i+1]*e_t[i+1] - 1.0/12.0*r[i+2]*e_t[i+2])/dr

            h_r[i] = -II*k_z*e_t[i] + II*nPhi*e_z[i]/r[i]
            h_t[i] = II*k_z*e_r[i] - dEz_dr 
            h_z[i] = (-II*nPhi*e_r[i] + drEt_dr )/r[i]

        endfor

        h_r = h_r / (_II*wRFc*_u0)
        h_t = h_t / (_II*wRFc*_u0)
        h_z = h_z / (_II*wRFc*_u0)

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
        Layout = [3,3]
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

        scale = 0.1

		thisField = e_r[*,*]
        scale = max(abs(real_part(thisField)))*scaleFac
        print, 'E_r scale: ', scale
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
        c_zero_set = contour(ar2.kPerSq_F,ar2.r,ar2.z,/over,c_value=0.001,color=cutoffColor,C_LABEL_SHOW=0,c_thick=3)
        c_jant = contour(jA_z, x, y, /over, c_value = AntCLevel, color=AntColor, C_LABEL_SHOW=0, c_thick=5)
        ++PlotPos
        p.Save, "Er.png", border=0, height=600

		thisField = e_t[*,*]
        scale = max(abs(real_part(thisField)))*scaleFac
        title = 'E_t'
        print, 'E_t scale: ', scale
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
        aspect = 1.0
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
            aspect_ratio=aspect, dimensions=dimensions, title=title, xRange=xRange, yRange=yRange, layout=[[Layout,PlotPos]], /current )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
        p = plot ( rlim, zlim, /over, thick = 3, color=LimColor)
        p = plot ( rlcfs, zlcfs, /over, thick = 4, color=lcfsColor )
        c_zero_set = contour(ar2.kPerSq_F,ar2.r,ar2.z,/over,c_value=0.001,color=cutoffColor,C_LABEL_SHOW=0,c_thick=3)
        c_jant = contour(jA_z, x, y, /over, c_value = AntCLevel, color=AntColor, C_LABEL_SHOW=0, c_thick=5)
        ++PlotPos
        p.Save, "Et.png", border=0, height=600

		thisField = e_z[*,*]
        scale = max(abs(real_part(thisField)))*scaleFac
        print, 'E_z scale: ', scale
        title = 'E_z'
		levels = fIndGen(nLevs)/(nLevs-1)*scale
		colors = reverse(bytScl(levels, top=253)+1)
		PlotField = (real_part(thisField)<max(levels))>min(-levels)
        aspect = 1.0
		c = contour ( PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=3, /fill, $
            aspect_ratio=aspect, dimensions=dimensions, title=title, xRange=xRange, yRange=yRange, layout=[[Layout,PlotPos]], /current )
		c = contour ( -PlotField, x, y, c_value=levels, rgb_indices=colors, rgb_table=1, /fill,/over )
        p = plot ( rlim, zlim, /over, thick = 3, color=LimColor)
        p = plot ( rlcfs, zlcfs, /over, thick = 4, color=lcfsColor )
        c_zero_set = contour(ar2.kPerSq_F,ar2.r,ar2.z,/over,c_value=0.001,color=cutoffColor,C_LABEL_SHOW=0,c_thick=3)
        c_jant = contour(jA_z, x, y, /over, c_value = AntCLevel, color=AntColor, C_LABEL_SHOW=0, c_thick=5)
        ++PlotPos
        p.Save, "Ez.png", border=0, height=600

		thisField = abs(e_r[*,*])
        ;scale = max(abs(thisField))*scaleFac
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

		thisField = abs(e_t[*,*])
        ;scale = max(abs(thisField))*scaleFac
        title = 'abs(E_t)'
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
        p.Save, "absEt.png",  border=0, height=600

		thisField = abs(e_z[*,*])
        ;scale = max(abs(thisField))*scaleFac
        title = 'abs(E_z)'
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
        p.Save, "absEz.png",  border=0, height=600

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
		E_alpOrig = E_alp
		E_betOrig = E_bet
		E_prlOrig = E_prl 

		save, xorig, E_alpOrig, E_betOrig, E_prlOrig, $
				file = 'soln.sav'
	endif

	stop
end
