;*******************************************************************************
; eswat:
; procedure used in calculating the mixing ratio, called by
; rrmixv. 
; See Thermos (Flatau, Lee) 5.12 and 5.17.
; This procedure is the IDL version of the described functions.
; temp = temperature, K
; res  = result, returned by this procedure
;
; history:
;               1-9-93          created, Tony Darugar

pro _eswat, temp, res

;Clear the math status value:
math_status = check_math(1,1)

rr = 373.16/temp

es = -7.90298 * (rr - 1) + 5.02808 * alog10(rr) - $
        1.3816e-7 * (10 ^ (11.344 * (1-(1/rr)))-1) + $
        8.1328e-3 * (10 ^ (-3.49149 * (rr - 1)) -1 )

res = 1013.246 * 10^es

math_status = check_math(0,0)
;if (math_status ne 0) then begin $
;        print, 'Math Errors were detected in the procedure <eswat>'
;        print, 'It is likely that there are missing data points'
;        print, 'in the temperature variable passed in.'
;endif

end

;*******************************************************************************
; rrmixv:
; Procedure for calculating mixing ratio, given the pressure, 
; temperature, and relative humidity.
; See Thermos (Flatau, Lee) 5.12 and 5.17.
; Parameters:
;       pres:   pressure
;       temp:   temperature, Kelvin
;               keyword
;       rh:     relative humidity, 0.0-1.0
;       res:    result, mixing ratio, g/g
;       celsius: keyword, if set to 1, temperature if taken to
;               be in celsius, if not present, or set to 0,
;               temp is assumed to be Kelvin.
;
; history:
;               1-9-93          created, Tony Darugar

pro _rrmixv, pres, temp, rh, res

ktemp = temp

res = ktemp
_eswat, ktemp, esat

for i=long(0),n_elements(ktemp) - 1 do $
  if (esat(i) ge 0.616*pres(i)) THEN $
        res(i) = 0 $
  else $
        res(i) = rh(i)*0.62198*esat(i)/(pres(i) - rh(i) * esat(i))

end


;*******************************************************************************
pro _esice, temp, res            ; returns e_sat over ice, mb

;Clear the math status value:
math_status = check_math(1,1)

rr = temp/273.16

es = (.876793 - 9.09718/rr) * (1.-rr) + 3.56654*alog10(rr)

res = 6.1071 * 10^es

math_status = check_math(0,0)
if (math_status ne 0) then begin $
        print, 'Math Errors were detected in the procedure <eswat>'
        print, 'It is likely that there are missing data points'
        print, 'in the temperature variable passed in.'
endif

end


;************************************************************************
pro _rrmixi, pres, temp, rh, res    ; returns mixing ratio over ice, g/g

ktemp = temp

res = ktemp
_esice, ktemp, esat

for i=0,n_elements(ktemp) - 1 do $
  if (esat(i) ge 0.616*pres(i)) THEN $
        res(i) = 0 $
  else $
        res(i) = rh(i)*0.62198*esat(i)/(pres(i) - rh(i) * esat(i))

end

pro _rrmix, pres, temp, rh, res    ; returns mixing ratio over ice/liquid.

res = temp*0
wh = where(temp lt 273.15, nice)
if nice gt 0 then begin
    _rrmixi, pres(wh), temp(wh), rh(wh), resice
    res(wh) = resice
endif
wh = where(temp ge 273.15, nliq)
if nliq gt 0 then begin
    _rrmixv, pres(wh), temp(wh), rh(wh), resliq
    res(wh) = resliq
endif
end

;***********************************************************************
function tsat, mixrat  ;  Saturation temp, mixing ratio is in g/kg
tsat =  255.33+12.46*alog(mixrat*(1.0+.0265*mixrat));
return, tsat
end
;*************

;**********************************************************************
; program pott-- calculates O and Oe-- mixing ratio is in g/kg
pro _pott, pres, temp, mixrat, pott, epott

pott = temp*0+999.
epott = pott
wh = where(temp lt 900 and temp gt 100)
pott(wh) = temp(wh)*(1000./pres(wh))^.286
wh = where(mixrat lt 90 and mixrat gt 0)
;mixrat2 = mixrat(wh)
;epott(wh) = pott(wh) * $
;  exp((3.376/Tsat(mixrat2)-.00254)*mixrat2*(1.0+.00081*mixrat2))
epott(wh) = pott(wh) * $
  exp((3.376/Tsat(mixrat[wh])-.00254)*mixrat[wh]*(1.0+.00081*mixrat[wh]))
end
;***********


