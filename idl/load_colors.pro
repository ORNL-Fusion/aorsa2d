common dlg_colors, ct12, $
	red, $
	blue, $
	green, $
	purple


loadct, 12, /sil, rgb_table = ct12

red	= transpose ( ct12[12*16-1,*] )
blue	= transpose ( ct12[8*16-1,*] )
green	= transpose ( ct12[2*16-1,*] )
purple	= transpose ( ct12[9*16-1,*] )
black	= red * 0


