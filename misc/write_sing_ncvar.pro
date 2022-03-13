; 
; Write an individual variable to a NetCDF file.
; 
; fname: complete NetCDF file path
; 
; var: variable (**WILL ASSUME FLOAT**)
; var_str: desired NETCDF variable name (string).
; 
; var1, var2, var3, ... var5: to add variables to file
;   with an identical shape to var
; var1_str, ...
; 
; James Ruppert
; 6.3.18
; 

pro write_nc_end,fid

  ncdf_control,fid,/endef
  ncdf_close,fid

  message,'Check that input variables match!'

end

pro write_sing_ncvar, fname, var, var_str, $
  var1=var1, var2=var2, var3=var3, var4=var4, var5=var5, $
  str_var1=str_var1, str_var2=str_var2, str_var3=str_var3, $
  str_var4=str_var4, str_var5=str_var5, $
  dim1=dim1, dimtag1=dimtag1, $
  dim2=dim2, dimtag2=dimtag2, $
  dim3=dim3, dimtag3=dimtag3, $
  dim4=dim4, dimtag4=dimtag4, $
  dim5=dim5, dimtag5=dimtag5

  ;VAR SPECS
  specs=size(var)
  ndims=specs[0]

  filename=fname
  rm_str='rm '+filename
  spawn,rm_str,rm_out,rm_err
  fid=ncdf_create(filename,/clobber,/netcdf4_format)

  ;DEFINE DIMENSIONS

  dimid=lonarr(ndims)
  for idim=0,ndims-1 do begin

    dimtag='dim'+strtrim(idim+1,2)

    if idim eq 0 and keyword_set(dimtag1) then dimtag=dimtag1
    if idim eq 1 and keyword_set(dimtag2) then dimtag=dimtag2
    if idim eq 2 and keyword_set(dimtag3) then dimtag=dimtag3
    if idim eq 3 and keyword_set(dimtag4) then dimtag=dimtag4
    if idim eq 4 and keyword_set(dimtag5) then dimtag=dimtag5

    dimid[idim]=ncdf_dimdef(fid,dimtag,specs[idim+1])

  endfor

  ;TO WRITE DIMENSIONS
    if keyword_set(dim1) then $
      varid_dim1=ncdf_vardef(fid,dimtag1,dimid[0],/float)
    if keyword_set(dim2) then $
      varid_dim2=ncdf_vardef(fid,dimtag2,dimid[1],/float)
    if keyword_set(dim3) then $
      varid_dim3=ncdf_vardef(fid,dimtag3,dimid[2],/float)
    if keyword_set(dim4) then $
      varid_dim4=ncdf_vardef(fid,dimtag4,dimid[3],/float)
    if keyword_set(dim5) then $
      varid_dim5=ncdf_vardef(fid,dimtag5,dimid[4],/float)

  ;DEFINE VARIABLES

  varid=ncdf_vardef(fid,var_str,dimid,/float)
  if keyword_set(var1) then begin
    if (size(var1))[0] ne ndims then write_nc_end,fid
    varid1=ncdf_vardef(fid,str_var1,dimid,/float)
  endif
  if keyword_set(var2) then begin
    if (size(var2))[0] ne ndims then write_nc_end,fid
    varid2=ncdf_vardef(fid,str_var2,dimid,/float)
  endif
  if keyword_set(var3) then begin
    if (size(var3))[0] ne ndims then write_nc_end,fid
    varid3=ncdf_vardef(fid,str_var3,dimid,/float)
  endif
  if keyword_set(var4) then begin
    if (size(var4))[0] ne ndims then write_nc_end,fid
    varid4=ncdf_vardef(fid,str_var4,dimid,/float)
  endif
  if keyword_set(var5) then begin
    if (size(var5))[0] ne ndims then write_nc_end,fid
    varid5=ncdf_vardef(fid,str_var5,dimid,/float)
  endif

  ;ADD ATTRIBUTES

;  ncdf_attput,fid,varid_pres,'def','surf pressure'
;  ncdf_attput,fid,varid_pres,'units','Pa'

  ;END DEFINE MODE

  ncdf_control,fid,/endef

  ;INSERT VARIABLES

  ncdf_varput,fid,varid,float(var)

  ;EXTRA VARS
  if keyword_set(var1) then ncdf_varput,fid,varid1,float(var1)
  if keyword_set(var2) then ncdf_varput,fid,varid2,float(var2)
  if keyword_set(var3) then ncdf_varput,fid,varid3,float(var3)
  if keyword_set(var4) then ncdf_varput,fid,varid4,float(var4)
  if keyword_set(var5) then ncdf_varput,fid,varid5,float(var5)

  ;DIMENSIONS
  if keyword_set(dim1) then ncdf_varput,fid,varid_dim1,float(dim1)
  if keyword_set(dim2) then ncdf_varput,fid,varid_dim2,float(dim2)
  if keyword_set(dim3) then ncdf_varput,fid,varid_dim3,float(dim3)
  if keyword_set(dim4) then ncdf_varput,fid,varid_dim4,float(dim4)
  if keyword_set(dim5) then ncdf_varput,fid,varid_dim5,float(dim5)

  ;CLOSE FILE

  ncdf_close,fid


  print,'Done writing out!'


end
