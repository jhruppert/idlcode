; 
; Function to calculate a profile of vertical motion in Weak Temperature
;   Gradient (WTG) balance according to the spectral-WTG method of
;   Herman and Raymond (2014, JAMES).
; 
; If deriv(hght) is non-constant (i.e., for a stretched vertical input),
;   grid, interpolates input onto a regular temporary grid. w_wtg is
;   returned on original input grid.
; 
; INPUT:
; 
;   hght = hght(z): height [m]
;   rho  = rho(z) : density [kg/m**3]
;   th_p = th_p(z): perturbation potential* temperature [K]
;   th_b = th_b(z): mean potential* temperature [K]
; 
;   * Should be equivalent potential temperature.
; 
; James Ruppert
; james.ruppert@mpimet.mpg.de
; 25.07.17
; 
function spectral_wtg, hght, rho, th_p, th_b

  ;SETTINGS (See HR16 for nomenclature)

    ;Horizontal relaxation length scale L (m)
    ;L_adj = 500e3
    L_adj = 250e3

    ;BVS frequency**2 thresholds
    bvs_sq_trop = 3e-4 ; >= this value in stratosphere
    bvs_sq_pbl  = 1e-4 ; <= this value in PBL


  ;VERTICAL GRID

    nz = n_elements(hght)

    ;check vertical consistency
    dz = deriv(hght)
    dz = dz[sort(dz)]
    if dz[0] lt 0 then message,'Inputs must have surface at first index!'
    if min(dz) ne max(dz) then is_stretched=1 else begin
      is_stretched=0
      dz = dz[0]
    endelse

    if is_stretched then begin

      ;save input grid
      hght_out = hght

      ;new grid spacing
      dz = 100.

      ;overwrite for new grid
      nz = round( (max(hght) - min(hght)) / dz )
      hght = findgen(nz) * dz + min(hght)
      rho  = interpol(rho,hght_out,hght)
      th_p = interpol(th_p,hght_out,hght)
      th_b = interpol(th_b,hght_out,hght)

    endif


  ;REQUIRED DIAGNOSTICS

    dthdz = deriv(hght,th_b) ; [ K/m ]

    d_term = th_p / dthdz

    ;BVS frequency
    bvs_sq = 9.81 * dthdz / th_b

    ;tropopause
    iz_trop = (where(bvs_sq ge bvs_sq_trop))[0]
    h_trop = hght[iz_trop]
    if h_trop lt 12e3 then stop

    ;PBL
    loc_pbl = max(where(bvs_sq[0:50] le bvs_sq_pbl))
    nz_valid = iz_trop-loc_pbl - 1
    loc_valid = indgen(nz_valid)+loc_pbl+1
    tmp=intarr(nz)
    tmp[loc_valid]=1
    loc_invalid=where(tmp eq 0)

    ;troposphere-mean BVS frequency
    bvs_sq_mean = total(bvs_sq[loc_valid] * rho[loc_valid],/double) $
                / total(rho[loc_valid],/double)
    bvs_mean = sqrt(bvs_sq_mean)


  ;POWER SPECTRUM OF D_TERM

  nmodes = round(nz_valid/2.)-1 ; Maximum resolvable signal
;  nmodes = 40 ; this yields tau <= 113 hours

  w_wtg_tmp = fltarr(nz)
;  w_wtg_tmp[*] = !values.f_nan
;  big_th = w_wtg_tmp

;  m = fltarr(nmodes)
;  tau = fltarr(nmodes)

  for ij=0,nmodes-1 do begin

;    m[ij] = ij * !pi / h_trop
    m = (ij+1) * !pi / h_trop
    tau = m * L_adj / bvs_mean

    big_th=0.
    for iz = loc_valid[0] , max(loc_valid) do $
      big_th += d_term[iz] * sin(m*hght[iz]) * dz
;    big_th[*,iz] = total( d_term[iz] * sin(m*hght[iz]) ,/double )

    big_th *= 2./h_trop

    for iz = loc_valid[0] , max(loc_valid) do $
      w_wtg_tmp[iz] += big_th * sin(m*hght[iz]) / tau

  endfor

  w_wtg_tmp[loc_invalid] = !values.f_nan


  ;PUT BACK ONTO ORIGINAL GRID

  w_wtg = interpol(w_wtg_tmp,hght,hght_out)


return, w_wtg

end