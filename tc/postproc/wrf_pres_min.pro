; 
; Identify center of pressure minimum based on directory containing WRF output.
;
; Set /local to indicate local reference frame (i.e., of lat/lon from time zero).
;   In that case, use as lon[pmin.lon[it]] and lat[pmin.lat[it]].
;
; James Ruppert
; 3/12/19
; 
function wrf_pres_min, output_dir, time, hurdat, dims, local=local


;READ AND SMOOTH SLP

  slp=reform(read_nc_var(output_dir+'post/SLP.nc',' ',varid='0'))

  specs=size(slp,/dimensions)
  nx=specs[0]
  ny=specs[1]
  nt=specs[2]

  ismooth=[3,3,0]
  for i=0,3 do $
    slp=smooth(temporary(slp),ismooth,/edge_truncate,/nan)

; CHECKED HOW THIS LOOKS; DOES AN ACCURATE JOB, AND EFFECTIVELY REMOVES NOISE AT SMALLEST SCALES

;NEXT MASK SLP WITH NANS
  slp=wrf_maria_mask(temporary(slp),time,hurdat,dims)

;LOCATE P-MIN

  pres=fltarr(nt)
  minloc=fltarr(2,nt)

  for it=0,nt-1 do begin
    pres[it]=min(reform(slp[*,*,it]),loc,/nan)
    minloc[*,it]=array_indices([nx,ny],loc,/dimensions)
  endfor


if keyword_set(local) then begin

;INNER DOMAIN REFERENCE FRAME

  pmin=create_struct('pres',pres,'lon',reform(minloc[0,*]),'lat',reform(minloc[1,*]))

endif else begin

;FIND ACTUAL LAT/LON

  spawn,'ls '+output_dir+'wrfout_d03*',rawfils

  lon=fltarr(nt)
  lat=lon

  count=[1,1,1]
  for i=0,nt-1 do begin
    offset=[minloc[0,i],minloc[1,i],0]
    lon[i]=reform(read_nc_var(rawfils[i],'XLONG',count=count,offset=offset))
    lat[i]=reform(read_nc_var(rawfils[i],'XLAT',count=count,offset=offset))
  endfor

  pmin=create_struct('pres',pres,'lon',lon,'lat',lat)

endelse

return,pmin

end
