;
; Convert u/v wind to azimuthal u (radial) and v (tangential) wind.
;
; Assumes u,v are fltarr( nradius , nazimuth [,...] ) and azimuth = fltarr( nazimuth )
;
; James Ruppert
; 12 May 2019
; Bayshore, NY
function azim_wind_conv, u, v, azimuth

d2r = !pi/180.

specs=size(u)
ndims=specs[0]
nrad=specs[1]
nazim=specs[2]

;AZIMUTH EXPANDED
  if ndims eq 2 then begin
    azim=fltarr(nrad,nazim)
    for iaz=0,nazim-1 do azim[*,iaz]=azimuth[iaz]
  endif
  if ndims eq 3 then begin
    azim=fltarr(nrad,nazim,specs[3])
    for iaz=0,nazim-1 do azim[*,iaz,*]=azimuth[iaz]
  endif
  if ndims eq 4 then begin
    azim=fltarr(nrad,nazim,specs[3],specs[4])
    for iaz=0,nazim-1 do azim[*,iaz,*,*]=azimuth[iaz]
  endif

;CONVERT WINDS
  wdir = atan( v / u ) / d2r
  loc_neg = where(u lt 0)
  wdir[loc_neg] += 180.
  wspd = sqrt( u^2 + v^2 )
  u_rad = wspd * cos( (wdir-azim) * d2r)
  v_tan = wspd * sin( (wdir-azim) * d2r)
;ind=90
;print,u[ind,ind,10,10],v[ind,ind,10,10]
;print,u_rad[ind,ind,10,10],v_tan[ind,ind,10,10],azim[ind,ind,10,10]

wnd = create_struct('u_rad',u_rad,'v_tan',v_tan)

return, wnd

end
