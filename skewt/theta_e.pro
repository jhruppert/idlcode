function theta_e,tmpk,pres,rv,rtot
; FUNCTION TO CALCULATE EQUIVALENT POTENTIAL TEMPERATURE
;
; BASED ON (34) OF BRYAN & FRITSCH 2002, OR (2.31) OF MARKOWSKI AND RICHARDSON,
; which is the "wet equivalent potential temperature" (BF02).
;
; INPUTS:
; 
;   TMPK: TEMPERATURE (K)
;   PRES: PRESSURE (Pa)
;   RV:   WATER VAPOR MIXING RATIO (KG/KG)
;   RTOT: TOTAL WATER (VAPOR+HYDROMETEOR) MIXING RATIO (KG/KG)
; 
; RETURNS:
; 
;   EQUIVALENT OTENTIAL TEMPERATURE (K)
; 
; OPTIONS:
; 
;   /REVERSE: if this is set, tmpc is assumed to be theta (in K), and temperature is instead returned (in C).
; 
; James Ruppert, ruppert@atmos.colostate.edu
; 8/4/14
; 

  ;CONSTANTS
    R=287.    ; J/K/kg
    lv0=2.5e6 ; J/kg
    cp=1004.  ; J/K/kg
    cpl=4186. ; J/k/kg
    cpv=1885. ; J/K/kg
    eps=18.0160/28.9660 ; Mw / Md (source: Brunner scripts)

  ;LATENT HEAT OF VAPORIZATION
    lv = lv0 - (cpl-cpv)*(tmpk-273.15)

  ;DRY AIR PRESSURE
    e = pres / ((eps/rv) + 1.)
    p_d = pres-e

  ;CALCULATE THETA-E
  th_e = tmpk * (1e5/p_d)^(R/(cp+cpl*rtot)) * exp( lv*rv / ((cp+cpl*rtot)*tmpk) )

  return,th_e

end
