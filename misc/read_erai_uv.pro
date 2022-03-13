;
; Program to ERAi U,V from daily-averaged file.
;
; Units of m/s
;
; James Ruppert
; 16 April 2019
;
function read_erai_uv, t_sel=t_sel, ilev=ilev, lon=lon, lat=lat, sel_lon=sel_lon, sel_lat=sel_lat, $
  time=time

  read_dir='/work/06040/tg853394/stampede2/era/'
  read_file=read_dir+'erai_uv.1998-01-01_2014-12-31_daily.nc'


;TIME SPAN

  d_start=julday(1,1,1998,0,0,0)
  d_end=julday(12,31,2014,0,0,0)
  times=timegen(start=d_start,final=d_end,step_size=1,units='D')
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

  ;0.75 DEGREE GRID
  nx=480
  ny=241
  lon=findgen(nx)*0.75
  lat=reverse(findgen(ny)*0.75-90)

  if keyword_set(sel_lon) then begin
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

  levs=[ 10000, 15000, 20000, 30000, 50000, 70000, 85000, 100000 ] * 1e-2 ; hPa
  if ~keyword_set(ilev) then message,'Pick a pressure level (set ilev): ',strtrim(levs,2)
  print,'Reading ERAi, pressure: ',strtrim(levs[ilev],2)


;READ

  count=[nxsel,nysel,1,nt] & offset=[xsel[0],ysel[0],ilev,it_sel[0]] ; x,y,t
  u=reform( read_nc_var(read_file,'var131',count=count,offset=offset) ) ; U: m/s
  v=reform( read_nc_var(read_file,'var132',count=count,offset=offset) ) ; V: m/s

  nan=where(abs(u) gt 1e4,count)
  if count gt 0 then u[nan]=!values.f_nan
  nan=where(abs(v) gt 1e4,count)
  if count gt 0 then v[nan]=!values.f_nan

;RETURN

time=times[it_sel]

era_wind=create_struct('u',u,'v',v)

return,era_wind

end
