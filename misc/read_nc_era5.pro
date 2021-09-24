;
; Function to read and return DAILY ERA5 data from full-length NetCDF files.
;
; Inputs:
;   time_read - Julian time array of times desired
;   ncfil - filename
;   varname - string containing the variable name in netcdf file
;
; Optional inputs:
;   bounds - [W,S,E,N] specifications for the read-in
;   plev - pressure (hPa) level to read, if a 3d file
;
; Optional returns: lon, lat
;
; James Ruppert
; 22 July 2021
;
function read_nc_era5, time_read, ncfile, varname, $
  bounds=bounds, plev=plev, $
  lon=lon, lat=lat

  ;DIMS
    lon=read_nc_var(ncfile,'lon')
    lat=read_nc_var(ncfile,'lat')

    nx=n_elements(lon)
    ny=n_elements(lat)

  ;CHECK FOR VERTICAL
    cdfid=ncdf_open(ncfile)
    varid=ncdf_varid(cdfid,varname)
    varinq=ncdf_varinq(cdfid,varid)
    ndims=varinq.ndims
    if ndims eq 4 then i3d=1 else i3d=0
    ncdf_close,cdfid

    ;TIME SUBSET
      time=reform(read_nc_var(ncfile,'time')) ; Hours since 2000-01-01 00:00:00Z
      nt=n_elements(time)
      if time[0] gt 1e7 then $
        time=timegen(n_elements(time),start=julday(1,1,2000,0,0,0),step_size=1,units="Days") $
      else $
        time=julday(1,1,2000,time,0,0)
      diftime=abs(time-time_read[0])
      it0=where(diftime eq min(diftime))
      it0=it0[0]
      diftime=abs(time-max(time_read))
      it1=where(diftime eq min(diftime))
      it1=it1[0]
      ntread=it1-it0+1
      it=indgen(ntread)+it0

    ;SPACE SUBSET
    ix0=0 & iy0=0
    if keyword_set(bounds) then begin
      ix=where((lon ge bounds[0]) and (lon le bounds[2]),nx)
      iy=where((lat ge bounds[1]) and (lat le bounds[3]),ny)
      ix0=ix[0] & iy0=iy[0]
      lon=lon[ix] & lat=lat[iy]
    endif

    ;VERTICAL
    if i3d then begin

      iz=0 ; set a default

      ;PLEVS
      p_era=reform(read_nc_var(ncfile,'plev'))*1d-2 ; Pa --> hPa
      iz=(where(p_era eq plev,count))[0]
      if count eq 0 then message,'Pressure level not found!'

      count=[nx,ny,1,ntread]
      offset=[ix0,iy0,iz,it0] ; x, y, p, t

    endif else begin

      count=[nx,ny,ntread]
      offset=[ix0,iy0,it0] ; x, y, t

    endelse

    var=reform(read_nc_var(ncfile,varname,count=count,offset=offset))

    ;NEED TO REVERSE LATITUDE
    var=reverse(var,2,/overwrite)
    lat=reverse(lat,/overwrite)

  return,var

end
