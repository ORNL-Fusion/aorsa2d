	VorpalBox_rOffSet = 1.8
	VorpalBox_r = [0,0.8,0.8,0.0,0.0]+VorpalBox_rOffset
	VorpalBox_zOffSet = -0.36/2.0
	VorpalBox_z = [0,0,0.36,0.36,0.0]+VorpalBox_zOffset

	LeftSide_rlim = [ $
			4.066, $
			4.325, $
			5.000, $
			5.779, $
			6.951, $
			7.578, $
			7.993, $
			8.3, $
			8.3, $
			8.3, $
			7.924, $
			7.318, $
			6.315, $
			5.796, $
			4.516, $
			4.083, $
			4.083, $
			4.066 ]*ShrinkFac

	LeftSide_zlim = ([ $
			3.581, $
			4.256, $
			4.671, $
			4.498, $
			3.599, $	
			2.976, $
			2.215, $
			1.626, $
			0.588, $
			-0.519, $
			-1.401, $
			-2.318, $
			-3.253, $
			-3.374, $
			-3.235, $
			-2.561, $
			0.000, $
			3.581 ]-0.5)*ShrinkFac

    RightSide_rLim = [$
            +1.7,$
            +2.0,$
            max(LeftSide_rLim),$
            max(LeftSide_rLim),$
            +2.60,+$
            2.60,$
            max(LeftSide_rLim),$
            max(LeftSide_rLim),$
            +2.0,$
            +1.7,$
            +1.7] 

    RightSide_zLim = [-0.25,-0.25,-0.25,-0.18,-0.18,+0.18,+0.18,+0.25,+0.25,+0.25,-0.25] 


