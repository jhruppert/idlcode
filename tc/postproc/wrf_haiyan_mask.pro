; 
; Mask out region beyond some Lat/Lon radius threshold from a desired location (for Maria simulations).
;
; James Ruppert
; 7/7/19
; 
function wrf_haiyan_mask, invar, time, hurdat, dims


outvar=invar

vardims=size(outvar)
ndim=vardims[0]
if ndim eq 3 then ntvar=vardims[3] else ntvar=vardims[4]

;SUBSET HURDAT TO MODEL TIME AND SPACE DOMAIN
  subset=where( (hurdat.jultim ge time[0]-0.01d and hurdat.jultim le max(time)+0.042d) and $
                (hurdat.lon ge dims.lon[0] and hurdat.lon le max(dims.lon)) and $
                (hurdat.lat ge dims.lat[0] and hurdat.lat le max(dims.lat)) )
  hdtim=hurdat.jultim[subset]
  hdlon=hurdat.lon[subset]
  hdlat=hurdat.lat[subset]

;INTERPOLATE HURDAT ONTO MODEL TIME
  t0=hdtim[0]
  tmax=max(hdtim)
  t_ind=where(time ge t0-0.01d and time le tmax+0.042d,nt_hd)

  ;Simply replace hdtim with time
  hdtim0=hdtim
  hdtim=time[t_ind]
  hdlon=interpol(temporary(hdlon),hdtim0,hdtim)
  hdlat=interpol(temporary(hdlat),hdtim0,hdtim)

if ndim ne 3 and ndim ne 4 then message,'Check variable format!'

;PRIOR TO STORM ID

  t_prior=where(time lt hdtim[0],ntprior)
;print,n_elements(where(finite(outvar[*,*,t_prior])))

  if ntprior gt 0 then begin

    ynan=where((dims.lat ge 10) or (dims.lat le 3),nnan,complement=ngood)
    if ndim eq 3 then outvar[*,ynan,t_prior]=!values.f_nan
    if ndim eq 4 then outvar[*,ynan,*,t_prior]=!values.f_nan
    xnan=where(dims.lon le 158,nnan,complement=ngood)
    if ndim eq 3 then  outvar[xnan,*,t_prior]=!values.f_nan
    if ndim eq 4 then outvar[xnan,*,*,t_prior]=!values.f_nan

  endif

;OTHER SIMPLE FILTERS

  ;CUT OFF INFLOW WALL

  it_xnan=where(time gt julday(11,1,2013,20,0,0),count)
;  ;indgen(ntvar-10)+10
  if count gt 0 then begin
    xnan=where(dims.lon gt 165,nnan,complement=ngood)
    if nnan gt 0 then begin
      if ndim eq 3 then  outvar[xnan,*,it_xnan]=!values.f_nan
      if ndim eq 4 then outvar[xnan,*,*,it_xnan]=!values.f_nan
    endif
  endif

;  it_xnan=where(time gt julday(9,15,2017,12,0,0),count)
;  ;indgen(ntvar-24)+24
;  if count gt 0 then begin
;    xnan=where(dims.lon gt -38,nnan,complement=ngood)
;    if nnan gt 0 then begin
;      if ndim eq 3 then  outvar[xnan,*,it_xnan]=!values.f_nan
;      if ndim eq 4 then outvar[xnan,*,*,it_xnan]=!values.f_nan
;    endif
;  endif

;  it_xnan=where(time gt julday(9,16,2017,1,0,0),count)
;  ;indgen(ntvar-37)+37
;  if count gt 0 then begin
;    xnan=where(dims.lon gt -41,nnan,complement=ngood)
;    if nnan gt 0 then begin
;      if ndim eq 3 then  outvar[xnan,*,it_xnan]=!values.f_nan
;      if ndim eq 4 then outvar[xnan,*,*,it_xnan]=!values.f_nan
;    endif
;  endif

;  endif

;HURDAT LOCATIONS

  buff=5. ; distance limit (in degrees) of x,y distance from Best-track center
;  t_match=where((time ge hdtim[0]) and (time le max(hdtim)),nt)
;  if max(t_match) gt ntvar-1 then nt-=max(t_match)-(ntvar-1)

  for it=0,nt_hd-1 do begin

    ynan=where(abs(dims.lat-hdlat[it]) gt buff)
    if ndim eq 3 then outvar[*,ynan,t_ind[it]]=!values.f_nan
    if ndim eq 4 then outvar[*,ynan,*,t_ind[it]]=!values.f_nan
    xnan=where(abs(dims.lon-hdlon[it]) gt buff)
    if ndim eq 3 then  outvar[xnan,*,t_ind[it]]=!values.f_nan
    if ndim eq 4 then outvar[xnan,*,*,t_ind[it]]=!values.f_nan

  endfor

;IN ADDITION TO HURDAT LOCATION, LOCALIZE TO MARIA'S CORNER OF DOMAIN

;    t_sel=indgen(ntvar-t_ind[0])+t_ind[0]
;
;    ynan=where(dims.lat le 11)
;    if ndim eq 3 then outvar[*,ynan,t_sel]=!values.f_nan
;    if ndim eq 4 then outvar[*,ynan,*,t_sel]=!values.f_nan
;    xnan=where(dims.lon ge -52)
;    if ndim eq 3 then  outvar[xnan,*,t_sel]=!values.f_nan
;    if ndim eq 4 then outvar[xnan,*,*,t_sel]=!values.f_nan

return,outvar

end
