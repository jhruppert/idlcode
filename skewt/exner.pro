function exner,pres,reverse=reverse

; FUNCTION TO CALCULATE EXNER FUNCTION
; 
; INPUTS:
; 
;   PRES: PRESSURE (Pa or hPa)
; 
; RETURNS:
; 
;   EXNER FUNCTION (UNITLESS)
; 
; OPTIONS:
; 
;   REVERSE: When this is set, pressure is returned (Pa) assuming exner is input.
; 
; James Ruppert, ruppert@atmos.colostate.edu
; 2/20/15


; Check for Pa vs. hPa
p0=1e5
if ~keyword_set(reverse) then $
  if max(pres) le 2000 then p0=1e3

if keyword_set(reverse) then $
  out=(pres^(7./2))*p0 $
else $
  out=(pres/p0)^(2./7)

return,out


end