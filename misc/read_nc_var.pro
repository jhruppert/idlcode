; 
; Read an individual variable from a NetCDF file.
; 
; fname: complete NetCDF file path
; 
; var_str: string of variable name as it appears in NetCDF file.
; 
; count,
; offset,
; stride: see info page for ncdf_varget
; 
; units_str: set this to return string of the variable units.
; 
; James Ruppert
; 24.08.16
; 
function read_nc_var , fname, var_str, $
  count=count, offset=offset, stride=stride, $
  units_str=units_str, varid=varid

  if keyword_set(count) then $
    s_count=count
  if keyword_set(offset) then $
    s_offset=offset
  if keyword_set(stride) then $
    s_stride=stride

  fid = ncdf_open(fname,/nowrite)
;  dinfo=NCDF_INQUIRE(fid)
  if ~keyword_set(varid) then $
    varid = ncdf_varid(fid,var_str) $
  else varid=fix(varid)

  ncdf_varget, fid, varid, var, count=s_count, offset=s_offset, stride=s_stride

  ;CHECK FOR SCALING AND OFFSETS
  ;CODE TAKEN AND MINIMALLY MODIFIED FROM "READ_NETCDF.PRO" AUTHORED
  ;BY Gunnar Spreen, PROVIDED BY STEFAN KERN (ICDC, HAMBURG) ON 20.04.2017

    varinf = ncdf_varinq(fid,varid)

    stest=0
    otest=0
    utest=0
    filltest=0
    misstest=0
    FOR ii=0,varinf.natts-1 DO BEGIN
      attname = NCDF_ATTNAME(fid,varid,ii)
      CASE attname OF
        'scale_factor': stest=1
        'add_offset': otest=1
        'units': utest=1
        '_FillValue': filltest=1
        'missing_value': misstest=1
        ELSE:
      ENDCASE
    ENDFOR

    count_fill=0
    IF filltest THEN BEGIN
      NCDF_ATTGET,fid,varid,'_FillValue',fill
      ifill=where(var eq fill,count_fill)
    ENDIF
    count_miss=0
    IF misstest THEN BEGIN
      NCDF_ATTGET,fid,varid,'missing_value',miss
      imiss=where(var eq miss,count_miss)
    ENDIF
    IF stest THEN BEGIN
      NCDF_ATTGET,fid,varid,'scale_factor',sfact
      var=sfact*TEMPORARY(var)
    ENDIF
    IF otest THEN BEGIN
      NCDF_ATTGET,fid,varid,'add_offset',add_offset
      var=TEMPORARY(var)+add_offset
    ENDIF
    IF utest THEN BEGIN
      NCDF_ATTGET,fid,varid,'units',unit
      units_str=STRING(unit)
    ENDIF

    IF filltest AND (count_fill gt 0) THEN $
      var[ifill]=!values.f_nan
    IF misstest AND (count_miss gt 0) THEN $
      var[imiss]=!values.f_nan

  ncdf_close,fid

  return,var

end
