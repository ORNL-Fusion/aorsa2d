pro zfunction

n = 101
rrange = 6d0
irange = 0.3d0
rp=(dIndGen(n)-n/2)/(n/2)*rrange
zp=(dIndGen(n)-n/2)/(n/2)*0+1d-5
rp2=rebin(rp,n,n)
zp2=transpose(rebin(zp,n,n))
zi=complex(0,1)
arg = complex(rp2,zp2)
;w=erfcx(arg)
w=exp(-arg^2)*erfc(-zi*arg)
z=zi*sqrt(!Pi)*w
s=surface(real_part(z),rp,zp)
s=surface(imaginary(z),rp,zp)
s=surface(arg*z,rp,zp)
stop
end
