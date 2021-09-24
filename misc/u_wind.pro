function u_wind,wspd,wdir
; 
; Calculate the u-wind from wind speed and wind direction (meteorological degrees).
; 
; Units returned are the same as those provided.
; 
; James Ruppert
; ruppert@atmos.colostate.edu
; 2/13/13
; 

u = (-1.) * wspd * sin( wdir * !dtor)

return,u

end