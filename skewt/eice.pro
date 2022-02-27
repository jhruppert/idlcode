;+
; NAME:
;	eice
;
; PURPOSE:
;	compute saturation vapor pressure over ice given temperature
;	in K or C. The accuracy is stated by Marti & Mauersberger
;       as 2% in the temperature range of 170K to 250K.
;   
; CATEGORY:
;	atmospheric physics
;
; CALLING SEQUENCE:
;	result=eice(t)
;
; EXAMPLE:
;	print,eice(-40)
;
; INPUTS:
;	T	SCALAR OR VECTOR OF TEMPERATURES IN CELSIUS OR K
;
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD INPUT PARAMETERS:
;
; OUTPUTS:
;	returns the saturation vapor pressure over ice in hPa
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
;  Dominik Brunner (brunner@atmos.umnw.ethz.ch), March 2000
;       reference: Marti and Mauersberger, GRL 20, 363-366, 1993.
;-

FUNCTION  EICE, T
  
   ON_ERROR,2
  
   IF MIN(T, /Nan) LT 105. THEN T0=273.15 ELSE T0=0. ; degC or K?
   
   ; Define constants
   A=-2663.5D
   B=12.537
   logp=A/(T+T0)+B
   
   RETURN,10^logp/100.       ; conversion to hPa
    
END

