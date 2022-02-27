;+
; NAME:
;	esat
;
; PURPOSE:
;	compute saturation vapor pressure given temperature in K or C
;   
; CATEGORY:
;	atmospheric physics
;
; CALLING SEQUENCE:
;	result=esat(t)
;
; EXAMPLE:
;	print,esat(15)
;               prints    17.0523
;
; INPUTS:
;	T	SCALAR OR VECTOR OF TEMPERATURES IN CELSIUS OR K
;
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD INPUT PARAMETERS:
;
; OUTPUTS:
;	returns the saturation vapor pressure in hPa
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
;  Dominik Brunner (brunner@atmos.umnw.ethz.ch), Feb 2000
;       A good reference is Gibbins, C.J., Ann. Geophys., 8, 859-886, 1990
;-

FUNCTION  ESAT, T
  
  ON_ERROR,2
  IF T[0] LT 105. THEN T0=273.16 ELSE T0=0.
    
; Formula with T = temperature in K
;    esat = exp( -6763.6/(T+T0) - 4.9283*alog((T+T0)) + 54.2190 )
    
; Formula close to that of Magnus, 1844 with temperature TC in Celsius
;    ESAT = 6.1078 * EXP( 17.2693882 * TC / (TC + 237.3) ) ; TC in Celsius

; or Emanuel's formula (also approximation in form of Magnus' formula,
; 1844), which was taken from Bolton, Mon. Wea. Rev. 108, 1046-1053, 1980.
; This formula is very close to Goff and Gratch with differences of
; less than 0.25% between -50 and 0 deg C (and only 0.4% at -60degC)    
;    esat=6.112*EXP(17.67*TC/(243.5+TC))
    
; WMO reference formula is that of Goff and Gratch (1946), slightly
; modified by Goff in 1965:
    
    e1=1013.250
    TK=273.16
    esat=e1*10^(10.79586*(1-TK/(T+T0))-5.02808*alog10((T+T0)/TK)+$
                1.50474*1e-4*(1-10^(-8.29692*((T+T0)/TK-1)))+$
                0.42873*1e-3*(10^(4.76955*(1-TK/(T+T0)))-1)-2.2195983)    
    
    RETURN, ESAT  
    
 END
