pro trace_ant_path

nX	= 20
nY	= 20

jx_int	= fltArr ( nX, nY )
jy_int	= fltArr ( nX, nY )

jx	= fltArr ( nX-1, nY )
jy	= fltArr ( nX, nY-1 )

pathX	= [10,10,10,10,10,09,09,09]
pathY	= [05,06,07,08,09,09,10,11]

xGrid	= fIndGen ( nX )
yGrid	= fIndGen ( nY )

xGrid_jx	= xGrid[1:*]-0.5
yGrid_jx	= yGrid

xGrid_jy	= xGrid
yGrid_jy	= yGrid[1:*]-0.5

x2D	= rebin ( xGrid, nX, nY )
y2D	= transpose ( rebin ( yGrid, nY, nX ) )

x2D_jx	= rebin ( xGrid_jx, nX-1, nY )
y2D_jx	= transpose ( rebin ( yGrid_jx, nY, nX-1 ) )

x2D_jy	= rebin ( xGrid_jy, nX, nY-1 )
y2D_jy	= transpose ( rebin ( yGrid_jy, nY-1, nX ) )

for i=0,n_elements ( pathX )-2 do begin

	jySgn	= (pathY[i+1]-pathY[i]) / abs (pathY[i+1]-pathY[i])
	jxSgn	= (pathX[i+1]-pathX[i]) / abs (pathX[i+1]-pathX[i])

	if jxSgn ne jxSgn then jxSgn = 0
	if jySgn ne jySgn then jySgn = 0

	jx[pathX[i]-1,pathY[i]]	+= 0.5 * jxSgn
	jx[pathX[i],pathY[i]]	+= 0.5 * jxSgn

	jy[pathX[i],pathY[i]-1]	+= 0.5 * jySgn
	jy[pathX[i],pathY[i]]	+= 0.5 * jySgn

endfor

;	calculate the divergence

divJ	= fltArr ( nX, nY )

for i=1,nX-2 do begin
	for j=1,nY-2 do begin

		divJ[i,j]	= ( jx[i,j]-jx[i-1,j] ) + ( jy[i,j-1] - jy[i,j] )

	endfor
endfor

window, 0, xSize = 800, ySize = 800
device, decomposed = 0
!p.background = 255
plot, x2D[*], y2D[*], psym = 4, color = 0

loadct, 3, /sil
for i=1,nX-2 do begin
	for j=1,nY-2 do begin

		polyFill,[xGrid[i]-1/2.0, $
			xGrid[i]+1/2.0, $
			xGrid[i]+1/2.0, $
			xGrid[i]-1/2.0, $
			xGrid[i]-1/2.0], $
			[yGrid[j]-1/2.0, $
			yGrid[j]-1/2.0, $
			yGrid[j]+1/2.0, $
			yGrid[j]+1/2.0, $
			yGrid[j]-1/2.0], $
			color = 255-divJ[i,j]*20, $
			/data


	endfor
endfor

plots, x2D[*], y2D[*], psym = 4, color = 0
plots, x2d_jx[*], y2d_jx[*], psym = 3, color = 0
plots, x2d_jy[*], y2d_jy[*], psym = 3, color = 0

loadct, 12, /sil
veloVect, jx, jx*0, xGrid_jx, yGrid_jx, /over, color = 12*16-1, length= 0.5
veloVect, jy*0, jy, xGrid_jy, yGrid_jy, /over, color = 12*16-1, length= 0.5

plots, pathX, pathY, psym = 4, color = 8*16-1

stop
end
