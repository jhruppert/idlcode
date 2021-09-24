;
; Program to TRMM 3B42V7 data from post-processed file (from read_trmm_hdf.pro).
;
; Units of MM / H
;
; James Ruppert
; 6 Nov 2018
;
function read_trmm_post, t_sel=t_sel, lon=lon, lat=lat, sel_lon=sel_lon, sel_lat=sel_lat,$
  read_file=read_file, time=time

  read_dir='/work/06040/tg853394/stampede2/trmm/'

  if keyword_set(read_file) then begin
    spawn,'ls '+read_file,out,err
    if keyword_set(err) then message,err
  endif else $
    read_file=read_dir+'trmm_3b42v7_01011998H00:00:00-01012014H00:00:00.nc'

;  if keyword_set(filt) then begin
;    if ~keyword_set(wave) then message,'Select wave for TRMM read'
;    read_dir=read_dir+'filtered/'+wave+'/'
;    read_file=read_dir+'trmm_rain_rate_mmph.nc'
;  endif

;TIME SPAN

  d_start=julday(1,1,1998,0,0,0)
  d_end=julday(1,1,2014,0,0,0)
  times=timegen(start=d_start,final=d_end,step_size=3,units='H')
  nt=n_elements(times)

  if keyword_set(t_sel) then begin
    ;CHECK BOUNDS
    julbuf=0.5d/24 ; 30-min safeguard against tiny fractions
    if (t_sel[0]+julbuf lt times[0]) or (t_sel[1]-julbuf gt max(times)) then $
      message,'T_SEL outside the available time range!'
    it_sel=where((times ge t_sel[0]) and (times le t_sel[1]),nt)
  endif else $
    it_sel=indgen(nt)

;DATASET SPECS

  nx=1440
  ny=400
  lon=findgen(nx)*0.25-180;+0.25/2
;  lon=shift(lon,-nx/2)
  lat=findgen(ny)*0.25-50;+0.25/2

  if keyword_set(sel_lon) then begin
    ;CAREFUL - TWO WAYS TO SPECIFY LON FOR WESTERN HEMISPHERE
    for i=0,1 do if sel_lon[i] lt 0 then sel_lon[i]+=360
    xsel=where(lon ge sel_lon[0] and lon le sel_lon[1],nxsel)
    ysel=where(lat ge sel_lat[0] and lat le sel_lat[1],nysel)
  endif else begin
    nxsel=nx
    nysel=ny
    xsel=indgen(nx)
    ysel=indgen(ny)
  endelse

  lon=lon[xsel]
  lat=lat[ysel]

;READ

  count=[nxsel,nysel,nt] & offset=[xsel[0],ysel[0],it_sel[0]] ; x,y,t
  rain=read_nc_var(read_file,'trmm_rain_rate_mmph',count=count,offset=offset) ; MM / H

  nan=where(rain lt -1000,count)
  if count gt 0 then rain[nan]=!values.f_nan

;RETURN

time=times[it_sel]

return,reform(rain)

end
