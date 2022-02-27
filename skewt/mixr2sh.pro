;+
; NAME:
;	mixr2sh.pro
;
; PURPOSE:
;	Convert mixing ratio (g H2O per kg of dry air) into
;       specific humidity (g H2O per kg of humid air)
;   
; CATEGORY:
;	atmospheric physics
;
; CALLING SEQUENCE:
;	result=mixr2sh(q)
;
; EXAMPLE:
;	print,mixr2sh([5,10,20])
;
;       prints       4.97512      9.90099      19.6078
;
; INPUTS:
;	W: H2O mixing ratio (g H2O / kg dry air)
;
; OPTIONAL INPUT PARAMETERS:
;
; KEYWORD INPUT PARAMETERS:
;
; OUTPUTS:
;	returns the mixing ratio as g H2O per kg of dry air
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
;
;                                      Mw*e                e
;  W (mixing ration) = m_h2o/m_dry = -------- = Mw/Md * -------
;                                    Md*(p-e)           (p-e)/e

;                                                   Mw*e          W
;  Q (spec. hum.)    = m_h2o/(m_dry + m_h2o) =  ------------  = -----
;                                               Mw*e+Md*(p-e)   1 + W
;
;-
FUNCTION  mixr2sh,W
  return,W/(1.+W/1000.)
END
