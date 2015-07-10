pro ar2_plot_vorpal

    v = ar2_read_vorpal()

    layout = [3,1]
    plotPos = 1
    p=plot(v.r,v.E_r,layout=[[layout],[plotPos]]) 
    p=plot(v.r,imaginary(v.E_r),/over,color='r') 
    ++plotPos

    stop
end
