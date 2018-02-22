function [arR] = ar2_read_rundata(folder)

files = dir(strcat(folder,'/output/runData*.nc'));

[M,N] = size(files);

if (M > 1)
    print, 'Too many runData files, please cleanup';
    exit
end

file = strcat(folder,'/output/',files(1).name);

arR = dlg_read_netcdf(file);

end