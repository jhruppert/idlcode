;
; Function for regressing a time-varying horizontal field y(x,y,t) onto an index x(t)
;
;   x = x(t)' is assumed to be standardized, such that variance(x) = 1.
;   y = y(i,j,t)
;   dims = [ni,nj,nt]
;
; James Ruppert
; 7/15/21

function regress_3d, x, y, dims

    dims2=dims[[2,0,1]] ; Need to place time first for rebinning
    x_rb = transpose(rebin(x,dims2),[1,2,0])

    ym = mean(y,dimension=3,/nan,/double)
    yp = y - rebin(ym,dims)

    nt=dims[2]
    return, total( (x_rb * yp) ,3,/nan,/double)/nt ; / variance (assumed to be =1)

end
