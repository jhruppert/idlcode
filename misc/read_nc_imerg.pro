;
; Function to read and return IMERG data from full-length NetCDF files.
;
; Input: time_read - Julian time array of times desired
;        ncfil - filename
;
; Optional input: bounds - [W,S,E,N] specifications for the read-in
;
; Optional returns: lon, lat
;
; James Ruppert
; 20 July 2021
;
function read_nc_imerg, time_read, ncfile, bounds=bounds, lon=lon, lat=lat

  imonthly=0
  if strmatch(ncfile,'*monthly*') then imonthly=1
  idaily=0
  if strmatch(ncfile,'*daily*') then idaily=1

  varname='precipitationCal'
  if imonthly then varname='precipitation'

  ;IMERG LAT/LON
    lon=read_nc_var(ncfile,'lon')
    lat=read_nc_var(ncfile,'lat')
    nx=n_elements(lon)
    ny=n_elements(lat)

    ;TIME SUBSET
      time=reform(read_nc_var(ncfile,'time')) ; Seconds since 1970-01-01 00:00:00Z
      if imonthly then time=julday(1,1,1970,0,0,time) $
      else if idaily then time=julday(1,1,1970,0,0,0)+time $
      else time=julday(1,1,1970,0,0,time);time+=julday(1,1,1970,0,0,0)
      diftime=abs(time-time_read[0])
      it0=where(diftime eq min(diftime))
      it0=it0[0]
      diftime=abs(time-max(time_read))
      it1=where(diftime eq min(diftime))
      it1=it1[0]

      nt_read=it1-it0+1
      it=indgen(nt_read)+it0

    ;SPACE SUBSET
    ix0=0 & iy0=0
    if keyword_set(bounds) then begin
      ix=where((lon ge bounds[0]) and (lon le bounds[2]),nx)
      iy=where((lat ge bounds[1]) and (lat le bounds[3]),ny)
      ix0=ix[0] & iy0=iy[0]
      lon=lon[ix] & lat=lat[iy]
    endif

    count=[ny,nx,nt_read] & offset=[iy0,ix0,it[0]] ; x,y,t
    rain=reform(read_nc_var(ncfile,varname,count=count,offset=offset))
    rain=transpose(temporary(rain),[1,0,2])

  return,rain

end
