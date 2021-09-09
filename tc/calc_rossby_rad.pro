; 
; Estimate Rossby radius of deformation.
; 
; James Ruppert
; 1/10/19
; 
pro calc_rossby_rad

; Coriolis Parameter: f = 2 * Omega * sin(phi), Omega = 7.292e-5 /s

; Rossby Radius Lr = N*H / f*pi

;INPUTS

  phi=10.0  ; Latitude
  bvs=0.01 ; Brunt-Vaisala frequency (/s)
  h=10.    ; Scale height (km)

;CALCULATION

  f = 2. * 7.292e-5 * sin(phi*!dtor)

  lr = bvs * h / ( f )

;PRINT

  print,'Lr = ',lr

end
