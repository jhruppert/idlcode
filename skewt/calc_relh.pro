; 
; Function to handle relative humidity calculation for both freezing
; and sub-freezing temperatures, given input variables in SI units.
; 
; qv  = mixing ratio (kg/kg)
; tmp = temperature (K)
; pr  = pressure (Pa)
; 

function calc_relh, qv, tmpk, pr, noice=noice

  relh = tmpk
  relh[*]=0.
  loc_ice = where(tmpk le 273.15,complement=loc_ni)

  if keyword_set(noice) then $
    relh=mixr2rh(qv*1.e3,pr*1e-2,tmpk) $ ; ignore ice RH
  else begin
    relh[loc_ice]=mixr2rh(qv[loc_ice]*1.e3,pr[loc_ice]*1e-2,tmpk[loc_ice],/ice)
    relh[loc_ni]=mixr2rh(qv[loc_ni]*1.e3,pr[loc_ni]*1e-2,tmpk[loc_ni])
  endelse


  return,relh

end
