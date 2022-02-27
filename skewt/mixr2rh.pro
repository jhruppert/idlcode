;+
; NAME:
;	mixr2rh.pro
;
; PURPOSE:
;	Convert mixing ratio (g H2O per kg of dry air) at given
;	temperature and pressure into relative humidity (%)
;   
; CATEGORY:
;	atmospheric physics
;
; CALLING SEQUENCE:
;	result=mixr2rh(mixr,p,T[,/ice])
;
; EXAMPLE:
;       mixr=findgen(5)+1
;	print,mixr2rh(mixr,1013.,273.15)
;         prints     0.750788      1.50339      2.25781      3.01407      3.77215
;;
; INPUTS:
;	MIXR: Float or FltArr(n) H2O mixing ratios in g H2O per kg dry air
;       p   : Float or FltArr(n) ambient pressure in hPa
;       T   : Float or FltArr(n) ambient Temperature in C or K
;
; KEYWORD INPUT PARAMETERS:
;     /ice  : Set keyword to return RH over ice
;
; OUTPUTS:
;	returns the relative humidity over liquid water or over ice
;	(if keyword /ice is set)
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
;  Dominik Brunner (brunner@atmos.umnw.ethz.ch), August 2001
;
;  Derivation:
;                                      Mw*e              e
;  W (mixing ratio) = m_h2o/m_dry = -------- = Mw/Md * ---
;                                    Md*(p-e)           p-e
;
;  RH (rel. hum.)    = e/esat(T)*100.
;
;-
FUNCTION  mixr2rh,MIXR,p,T,ice=ice
   ; check for input
   ON_ERROR,2
   IF n_params() NE 3 THEN message,'Usage: mixr2rh(RH,p,T)'
   ; note east accepts temperature in deg C or in Kelvin
   IF keyword_set(ice) THEN es=eice(T) ELSE es=esat(T)
   Mw=18.0160 ; molecular weight of water
   Md=28.9660 ; molecular weight of dry air
   fact=MIXR/1000.*Md/Mw
   return,p/es*fact/(1+fact)*100.
END
