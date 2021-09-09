;
; Simple program to plot the mean field for a variable, used frquently
;  for testing.
;
; James Ruppert
; 7/22/2021
;
pro plot_mean_map_myanmar, dirs, var, lon, lat, var_plot, figname, $
  u=u, v=v, eralon=eralon, eralat=eralat, cvar=cvar, cross=cross

if var_plot eq 'rain' then begin
  var_str='RAINNC'
  setmax=25 & setmin='0.'
  cbform='(i2)'
endif else if var_plot eq 'pw' then begin
  var_str=var_plot
  setmax=70 & setmin=30
  cbform='(i2)'
endif

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=2
  figspecs=create_struct(figspecs,'figname',figname)
  figspecs.cbar_format=cbform
  figspecs.ndivs-=1
  figspecs.title='Mean map test'

  if keyword_set(u) then begin
    wind=create_struct('u',u,'v',v,'x',eralon,'y',eralat)
    wspd=sqrt(u^2+v^2)
  endif

  if ~keyword_set(cvar) and keyword_set(u) then $
    cvar=create_struct('cvar',wspd,'x',eralon,'y',eralat)

  wrf_myanmar_map_plot, dirs, var, lon, lat, figspecs, wind=wind, cvar=cvar, cross=cross

end
