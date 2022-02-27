;+
; 
; CREATED USING THE ATMOSPHERIC PHYSICS SCRIPTS OF DOMINIK BRUNNER
; BY JAMES RUPPERT
; 2/23/2015
; 
; 
; NAME:
;	rh2mixr.pro
;
; PURPOSE:
;	Convert relative humidity (%) at given
;	temperature and pressure into mixing ratio (kg H2O per kg of dry air)
;   
; CATEGORY:
;	atmospheric physics
;
; INPUTS:
;	RH: Float or FltArr(n) %
;       p   : Float or FltArr(n) ambient pressure in hPa
;       T   : Float or FltArr(n) ambient Temperature in C or K
;
; KEYWORD INPUT PARAMETERS:
;     /ice  : Set keyword to assess e_s over ice
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
FUNCTION  rh2mixr,RH,p,T,ice=ice
   ; check for input
   ON_ERROR,2
   IF n_params() NE 3 THEN message,'Usage: rh2mixr(RH,p,T)'
   ; note east accepts temperature in deg C or in Kelvin
   IF keyword_set(ice) THEN es=eice(T) ELSE es=esat(T)
   Mw=18.0160 ; molecular weight of water
   Md=28.9660 ; molecular weight of dry air
;   fact=MIXR/1000.*Md/Mw
;   return,p/es*fact/(1+fact)*100.
   fact=RH*es/100.
   return,(Mw/Md)*fact/(p-fact)
   
END
