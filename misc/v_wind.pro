function v_wind,wspd,wdir
; 
; Calculate the v-wind from wind speed and wind direction (meteorological degrees).
; 
; Units returned are the same as those provided.
; 
; James Ruppert
; ruppert@atmos.colostate.edu
; 2/13/13
; 

v = (-1.) * wspd * cos( wdir * !dtor)

return,v

end