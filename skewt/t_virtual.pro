;+
; NAME:
;	t_virtual
;
; PURPOSE:
;	Calculates the virtual temperature (K) for given
;       absolute temperature (K) and water vapor pressure e (hPa)
;   
; CATEGORY:
;	atmospheric physics
;
; CALLING SEQUENCE:
;	res=t_virtual(T,p,e)
;
; EXAMPLE:
;	T=293.   ; absolute temperature
;       U=[10.,20,30,40,50,60,70,80,90,100]    ; relative humidities
;       p=850.   ; ambient pressure
;       print,t_virtual(T,p,U/100.*esat(T))
;
; INPUTS:
;	T       temperature (K)
;       p       pressure (hPa)
;       e       water vapour pressure (hPa)
;
; OPTIONAL INPUT PARAMETERS:
;	
; KEYWORD INPUT PARAMETERS:
;
; OUTPUTS:
;	return virtual temperature in K
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;	
; MODIFICATION HISTORY:
;
; Dominik Brunner, March 2000 (brunner@atmos.umnw.ethz.ch)
; 
;-

FUNCTION t_virtual,T,p,e
  
ON_ERROR,2
;check requested parameters:
  IF N_PARAMS() NE 3 THEN BEGIN
    message,'Usage: tvirt=t_virtual(T,e)'
  ENDIF
  
return,T*(1+0.378*e/p)

end