function [ar2] = ar2_read_solution(folder)

files = dir(strcat(folder,'/output/solution*.nc'));

[M,N] = size(files);

if (M > 1)
    print, 'Too many solution files, please cleanup';
    exit
end

file = strcat(folder,'/output/',files(1).name);


ar2 = dlg_read_netcdf(file);

end