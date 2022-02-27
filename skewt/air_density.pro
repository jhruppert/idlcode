;+
; NAME:
;    air_density (function)
;
; PURPOSE:
;
;        Compute air mass density for given pressure(s) and 
;        temperature(s) in SI units, i.e. kg/m3
;
; CATEGORY:
;
;        Atmospheric physics, thermodynamics
;
; CALLING SEQUENCE:
;
;        air_density(p,T)
;
; INPUTS:
;
;        p     : Float or FltArr(n): pressures in hPa (mbar)
;        T     : Float or FltArr(m): temperatures in Kelvin
;      
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;
;        The dry air mass density in kg/m3. The result type will 
;        match the type and array dimensions of p unless p is a scalar
;        and T an array.
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;
; EXAMPLE:
;
;    print,air_density(1013.25,273.15)
;    ; prints 2.69000e+19
;
;    p = findgen(5)*100+500
;    print,air_density(p,273.15)
;    ; prints 1.32741e+19  1.59289e+19  1.85838e+19  2.12386e+19  2.38934e+19
;
; MODIFICATION HISTORY:
;
;   (c) Dominik Brunner, Aug 2001: first implementation
;
;-
function air_density,p,T
   
   ; some constants
   Rd=287. ; Gas constant for dry air
   ; create result depending on dimensions of p,T
   result = make_array(size=size([p]),value=-9999.)
   if (n_elements(p) eq 1 AND n_elements(T) gt 1) then  $
      result = make_array(size=size(T),value=-9999.)
   
   ; compute dens only for valid values of p and T 
   okind = where( T gt 0. AND p gt 0.)
   
   result[okind] = p*1e2/(Rd*T[okind])
   
   ; convert vector to scalar if only on element 
   if (n_elements(result) eq 1) then result = result[0]
 
   return,result
 
end
