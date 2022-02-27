;---------------------------------------------------------------------------------------------------------------
;  Subroutine to locate the tropopause from a sounding.
;  
;  FINDTROP, PRES, HGHT, TEMPC, MAX, TPLEV, TTROP, TMIN, PMIN, KS, LTRP
;  
;  PROVIDE 1-D SOUNDING: PRES, HGHT, TEMPC
;  
;  RETURNED:
;  
;  TPLEV:     Tropopause pressure
;  TTROP:     Tropopause temperature
;  LTRP:      Vertical index of tropopause location
;  
;  OPTIONS:
;  THRESH:    Threshold (minimum) pressure value required for determining tropopause height (returns NANs if
;             not met)
;  
;  James Ruppert (ruppert@atmos.colostate.edu)
;  4/25/2012
;---------------------------------------------------------------------------------------------------------------
pro findtrop,pres,hght,tempc,tplev,ttrop,ztrop,thresh=thresh

  on_error,0
  
  ;IGNORE NANS
  loc_valid=where(finite(pres))
  pres=pres[loc_valid]
  hght=hght[loc_valid]
  tempc=tempc[loc_valid]
  
  loc_min=(where(tempc eq min(tempc)))[0]
  
  dims=size(pres)
  nz=dims[1]
  top=nz-1
  
  ;CHECK THRESHOLD
  if keyword_set(thresh) then $
    if pres[top] gt thresh then begin
      nan=!values.f_nan
      ltrp=nan
      tplev=nan
      ttrop=nan
      return
    endif
  
  dtdz=deriv(hght,tempc)*1000.
  idxstable=where(dtdz ge -2)
  loc=min(idxstable)
  
  izz=loc+1
  hght_check=hght[izz]-hght[loc]
  
  while hght[top]-hght[loc] ge 2000 and hght_check lt 2000 and izz+1 le top do begin
    
    if dtdz[loc] ge -2 and finite(hght[loc]) then begin
      
      izz=loc+1
      while ~finite(hght[izz]) do izz+=1
      ave_dtdz = (tempc[izz]-tempc[loc]) / (hght[izz]-hght[loc]) * 1000.
      
      while ave_dtdz ge -2 and hght_check lt 2000 and izz+1 le top do begin
        
        izz+=1
        while ~finite(hght[izz]) do izz+=1
        ave_dtdz = (tempc[izz]-tempc[loc]) / (hght[izz]-hght[loc]) * 1000.
        hght_check = hght[izz]-hght[loc]
        
      endwhile
      
      if ave_dtdz lt -2 then loc+=1
      
      if hght_check ge 2000 and hght[izz]-hght[izz-1] gt 500 then begin
        loc+=1
        hght_check=0
      endif
      
      if hght[loc_min]-hght[loc] gt 1000 then begin
        loc+=1
        hght_check=0
      endif
      
    endif else loc+=1
    
    while ~finite(hght[loc]) do loc+=1
    
  endwhile
  
  ltrp=loc
  tplev=pres[ltrp]
  ttrop=tempc[ltrp]
  ztrop=hght[ltrp]

end