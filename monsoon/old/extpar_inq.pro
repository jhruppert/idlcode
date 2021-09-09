; 
; Get specs from EXTPAR file for running ICON.
; 
; James Ruppert
; 15.12.17
; 
pro extpar_inq

;SOME SETTINGS

;  extparfil='/home/mpim/m300462/mistwork/extpar_lem/exp/'+$
;    'icon_extpar_0039_R02B07_G_20170202_tiles.nc'
  homedir='/home/mpim/m300462/'
  extparfil=homedir+'Downloads/icon_extpar_0039_R02B07_G_20170202_tiles.nc'
  figdir=homedir+'misthome/IDLWorkspace85/Default/figures/monsoon/'
  extpar_oceanland_fil=figdir+'extpar_ocean_land_fromglobal.txt'

  ncdf_list,extparfil,/variables,vname=all_vars
  nvars=n_elements(all_vars)

  ;disregard TIME and MLEV
  all_vars=all_vars[2:nvars-1]
  nvars-=2

;OPEN FILE
  fid=ncdf_open(extparfil)

  ;LAT and LON
    vid=ncdf_varid(fid,'lon')
    ncdf_varget,fid,vid,lon
    vid=ncdf_varid(fid,'lat')
    ncdf_varget,fid,vid,lat

  ;TOPO
    vid=ncdf_varid(fid,'topography_c')
    ncdf_varget,fid,vid,topog

  ;OCEAN AND LAND POINT
    ocean_point=[77.5,0] ; [x,y] Just E of Maldives
    land_point=[75,15]
    xtmp=abs(lon-ocean_point[0])
    ytmp=abs(lat-ocean_point[1])
    ind_ocean=where((xtmp le 0.1) and (ytmp le 0.1))
    ind_ocean=ind_ocean[0]
    xtmp=abs(lon-land_point[0])
    ytmp=abs(lat-land_point[1])
    ind_land=where((xtmp le 0.1) and (ytmp le 0.1))
    ind_land=ind_land[0]

  ;WRITE DATA TO ASCII FILE
    openw,1,extpar_oceanland_fil

  for ivar=0,nvars-1 do begin
    varstr=all_vars[ivar]
    vid=ncdf_varid(fid,varstr)
    ncdf_varget,fid,vid,var
    if varstr eq 'LU_CLASS_FRACTION' then $
      for i=0,22 do $
        printf,1,varstr,' class',i,'  ocean:',var[ind_ocean,i],'  land: ',var[ind_land,i] $
    else $
      printf,1,varstr,'  ocean:',var[ind_ocean],'  land: ',var[ind_land]
  endfor

  close,1

stop

;CLOSE FILE
  ncdf_close,fid

;stop
;PLOT MAP

;  lon+=180
  x=lon
  y=lat
  area = [min(y),min(x),max(y),max(x)]
  area = [-10,40,30,120]

  xloc = where( (lon ge area[1] and lon le area[3]) and $
                (lat ge area[0] and lat le area[2]) )

;  set_plot,'x'
  set_plot,'ps'
  figname=figdir+'map'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=5,ysize=3,/inches,$
    /helvetica

  position=[0.1,0.1,0.9,0.9]
  map_set,0,0,/mercator,limit=area,xmargin=0,ymargin=0,position=position;,/isotropic

  loadct,8,/silent
;  polyfill,position[[0,0,2,2]],position[[1,3,3,1]],color=200,/normal
;  polyfill,area[[0,0,2,2]],area[[1,3,3,1]],color=200,/data

  loadct,0,/silent
  map_continents,/coasts,color=0,mlinethick=0.4,/hires,/countries;,fill_continents=1

  ;TOPOGRAPHY
  itopo=0
  if itopo then begin
    print,'Contouring topo'
    maxv=max(topog)
    nlevs=20;255
    levels=findgen(nlevs)/(nlevs-1)*maxv
    colors=indgen(nlevs)
    colors=reverse(colors)
    loadct,0,/silent
;    contour,topog[xloc],x[xloc],y[xloc],/cell_fill,/overplot,/irregular,$
;      levels=levels,c_colors=colors
    cint=100
    clevs=indgen(50)*cint+cint
    contour,topog[xloc],x[xloc],y[xloc],/follow,/overplot,/irregular,$
      levels=clevs,c_colors=0
  endif

  plots,lon[ind_ocean],lat[ind_ocean],psym=1,thick=2.0,/data
;  plots,land_point[0]+180,land_point[1],psym=2,thick=2.0,/data
  plots,lon[ind_land],lat[ind_land],psym=2,thick=2.0,/data

;  lats=indgen(5)*5-5 & latnames=['5!9%!XS','0!9%!X','5!9%!XN','10!9%!XN','15!9%!XN']
;  lons=indgen(5)*15-60 & lonnames=['60!9%!XW','45!9%!XW','30!9%!XW','15!9%!XW','0!9%!X']
  map_grid,/box_axes,londel=5,latdel=5,glinestyle=0;,lats=lats,lons=lons,charsize=csize*0.8,color=0,/no_grid
;  map_grid,/box_axes,londel=20,latdel=20,glinestyle=0

  device,/close

  convert_png,figname,res=200,/remove_eps

print,'Done!'

end
