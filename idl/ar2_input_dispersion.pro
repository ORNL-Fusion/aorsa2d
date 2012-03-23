pro ar2_input_dispersion

	stixp = 1d0
	stixl = 1d0
	stixr = 1d0

	resonances = fltarr(nx,ny,nspec)

	for s=0,nspec-1 do begin

		for i=0,nx-1 do begin
			for j=0,ny-1 do begin

				if mask_bbbs[i,j] eq 1 then begin

					xmap_n[i,j,s] = prof1(sqrt(xmap_psinorm[i,j]),nn[*,s])
					xmap_t[i,j,s] = prof1(sqrt(xmap_psinorm[i,j]),tt[*,s])

				endif else begin

					eta = xmap_dlcfs[i,j]*lambda
					xmap_n[i,j,s] = nn[0,s] * exp (-eta)
					xmap_t[i,j,s] = tt[0,s] * exp (-eta)

				endelse

			endfor
		endfor

		;xmap_n[*,*,s]=xmap_n[*,*,s]>nn[0,s]
		;xmap_t[*,*,s]=xmap_t[*,*,s]>tt[0,s]

		print, 'min ne after smoothing: ', min (xmap_n[*,*,s])

		xmap_r_a = ( xmap_r - eqdsk.rmaxis ) / (max(eqdsk.rbbbs) - eqdsk.rmaxis)

		;	calculate cutoffs etc

		;	these are in si not cgs units (i think)
		
		wp	= sqrt ( xmap_n[*,*,s] / abs(z[s]) * (z[s]*e)^2 / ( amu[s]*mi * epsilon0 ) )
		wc	=  z[s]*e * xmap_b / ( amu[s]*mi )

		for i=0,nx-1 do begin
			for j=0,ny-1 do begin
				resonances[i,j,s] = 1.0 / ( wc[i,j] mod w )
			endfor
		endfor

		stixp	= stixp - ( wp^2/w^2 )
		stixl	= stixl - ( wp^2 / ( w * ( w - wc ) ) )
		stixr	= stixr - ( wp^2 / ( w * ( w + wc ) ) )

	endfor

	stixs	= 0.5d0 * ( stixr + stixl )
	stixd	= 0.5d0 * ( stixr - stixl )
	stixa	= stixs

	if ~keyword_set(nphi) then nphi = -27
	print, 'nphi = ', nphi

	;	without upshift

	kpar2d	= nphi / xmap_r2d
	npar2d	= kpar2d * c / w

	stixb	= -1.0 * ( stixr * stixl + stixp * stixs - npar2d^2 * ( stixp + stixs ) )
	stixc	= stixp * ( npar2d^2 - stixr ) * ( npar2d^2 - stixl )
	b24ac_	= stixb^2 - 4.0 * stixa * stixc
	kperp2d_1	= sqrt ( ( -stixb + sqrt ( complex( b24ac_, b24ac_*0) ) ) / ( 2d0 * stixa ) ) * w / c
	kperp2d_2	= sqrt ( ( -stixb - sqrt ( complex( b24ac_, b24ac_*0) ) ) / ( 2d0 * stixa ) ) * w / c

   	s_ne = surface(xmap_n[*,*,0], layout=[3,1,1])
   	s_te = surface(xmap_t[*,*,0], layout=[3,1,2], /current)

	p = plot(xmap_r_a, kperp2d_2[*,ny/2], yrange = [0,200], color='black', thick=3)
	p = plot(xmap_r_a, kperp2d_1[*,ny/2], /over, color='blue', thick=3)
	p = plot(xmap_r_a, imaginary(kperp2d_2[*,ny/2]), /over,color='black',thick=3,linestyle='--')
	p = plot(xmap_r_a, imaginary(kperp2d_1[*,ny/2]), /over,color='blue',thick=3,linestyle='--')

	p = plot(xmap_r_a, xmap_n[*,ny/2,0], /ylog, yrange=[1e10,1e20])
	p = plot(xmap_r_a, xmap_n[*,ny/2,1], /over, thick = 3, color = 'red')
	p = plot(xmap_r_a, xmap_n[*,ny/2,2], /over, thick = 3, color = 'blue')

	ii=where(xmap_r gt 7.5) 
	p=plot(xmap_r[ii],xmap_n[ii,ny/2,0],/ylog)
	p=plot(xmap_dlcfs[ii,ny/2]/l_sol,xmap_n[ii,ny/2,0],/ylog)

	p = plot(xmap_dlcfs[ii,ny/2], kperp2d_2[ii,ny/2], yrange = [0,200], color='black', thick=3)
	p = plot(xmap_dlcfs[ii,ny/2], kperp2d_1[ii,ny/2], /over, color='blue', thick=3)
	p = plot(xmap_dlcfs[ii,ny/2], imaginary(kperp2d_2[ii,ny/2]), /over,color='black',thick=3,linestyle='--')
	p = plot(xmap_dlcfs[ii,ny/2], imaginary(kperp2d_1[ii,ny/2]), /over,color='blue',thick=3,linestyle='--')


end
