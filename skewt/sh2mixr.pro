; Created by James Ruppert, 8/19/2014
; 
; INUPUT: q (specific humidity) in kg / kg
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
;  w = q / ( 1 - q )
;  
;-
FUNCTION  sh2mixr,q
  w = q / ( 1. - q )
  return,w
END
