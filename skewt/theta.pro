function theta,tmpc,pres,reverse=reverse

; FUNCTION TO CALCULATE POTENTIAL TEMPERATURE
; 
; INPUTS:
; 
;   TMPC: TEMPERATURE (C)
;   PRES: PRESSURE (hPa)
; 
; RETURNS:
; 
;   POTENTIAL TEMPERATURE (K)
; 
; OPTIONS:
; 
;   /REVERSE: if this is set, tmpc is assumed to be theta (in K), and temperature is instead returned (in C).
; 
; James Ruppert, ruppert@atmos.colostate.edu
; 8/4/14
; 


if ~keyword_set(reverse) then begin
  theta = (tmpc+273.15) * (1000./pres)^(2./7)
  return,theta
endif else begin
  theta = tmpc
  tmpc2 = theta * (1000./pres)^(-2./7)
  tmpc2 -= 273.15
  return,tmpc2
endelse


end