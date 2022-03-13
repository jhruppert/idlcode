; 
; Min, Mean, and Max
; 
pro stats, var

print,'Min, Mean, Max'
print,min(var,/nan),mean(var,/nan,/double),max(var,/nan)

end