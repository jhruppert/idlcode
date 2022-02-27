; 
; Function to calculate frozen moist static energy. See Wing and Emanuel 2014, JAMES
; 
; h = cp*T + g*z + Lv*qv - Lf*qi
; 
; qv   = vapor mixing ratio (kg/kg)
; qi   = (SET AS OPTION) TOTAL ice mixing ratio (kg/kg) (i.e., including cloud ice, snow, graupel, hail, ...)
; tmpk = temperature (K)
; z    = height (m)
; 

function calc_mse_frozen, qv, tmpk, z, qi=qi

  ;CONSTANTS
    g=9.81 ; m/s2
    cp=1004. ; J/K/kg
    rgas=287. ; J/kg
    xlv=2.5e6   ; J/kg
    xls=2.834e6 ; J/kg
    ;BELOW TAKEN FROM CM1, BRYAN & FRITSCH 2002, MWR
      cpv = 1870.0
      cpl = 4190.0
      cpi = 2106.0
      xlv0 = 2501000.0
      xls0 = 2834000.0
      lv1 = xlv0+(cpl-cpv)*273.15
      lv2 = cpl-cpv
      ls1 = xls0+(cpi-cpv)*273.15
      ls2 = cpi-cpv
      ;use as: lv = lv1 - lv2*T

  ;Latent heats
    lv = lv1 - lv2*tmpk
    ls = ls1 - ls2*tmpk
    lf = ls - lv

  if ~keyword_set(qi) then qi=0

  h = cp*tmpk + g*z + lv*qv - lf*qi

  return,h

end
