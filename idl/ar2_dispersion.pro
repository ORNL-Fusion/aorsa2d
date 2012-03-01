pro ar2_dispersion

@constants

gFileName = 'Scen4_bn2.57_129x129'

freq 	= 53e6
ni		= 3.204e19 
Zi		= 1d0
B		= 5.3
nPhi	= -27
R		= 7.0
kPar	= nPhi / R

g = readgeqdsk(gFileName,/noTor)

B = g.bMag[*,g.nH/2]

wrf = freq * 2 * !pi;
wpi = sqrt ( ni * (Zi*e)^2 / (mi * e0) );
wci	= Zi*e*B/mi

kPerSq	= wrf^2/c^2*wpi^2/(wci*(wrf+wci))-kPar^2/2.0
kPer 	= sqrt(kPerSq)
lambda	= 2*!pi/kPer

p = plot(g.R,lambda)
stop
end
