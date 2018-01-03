function [ar2] = ar2_read_solution(folder)

files = dir(strcat(folder,'output/solution*.nc'));

[M,N] = size(files);

if (M > 1)
    print, 'Too many solution files, please cleanup';
    exit
end

file = strcat('output/',files(1).name);

ar2 = containers.Map();

ar2('r') = ncread(file,'r');

jr_re = ncread(file,'jP_r_re');
jr_im = ncread(file,'jP_r_im');
jt_re = ncread(file,'jP_t_re');
jt_im = ncread(file,'jP_t_im');
jz_re = ncread(file,'jP_z_re');
jz_im = ncread(file,'jP_z_im');

ar2('jr') = complex(jr_re,jr_im);
ar2('jt') = complex(jt_re,jt_im);
ar2('jz') = complex(jz_re,jz_im);

Er_re = ncread(file,'er_re');
Er_im = ncread(file,'er_im');
Et_re = ncread(file,'et_re');
Et_im = ncread(file,'et_im');
Ez_re = ncread(file,'ez_re');
Ez_im = ncread(file,'ez_im');

ar2('Er') = complex(Er_re,Er_im);
ar2('Et') = complex(Et_re,Et_im);
ar2('Ez') = complex(Ez_re,Ez_im);

end