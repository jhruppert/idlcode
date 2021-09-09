; 
; Plot maps from DYNAMO WRF simulations.
; 
; James Ruppert
; 20.10.16
; 
pro gfs_map_plot, var, lon, lat, figspecs, $
  cvar=cvar, wind=wind, tcloc=tcloc

iwrf_dom=1 ; Plot WRF domain?
wrf_dom=[-63,2,-35,18] ; L,B,R,T

;CREATE FIGURE

  ;PLOT SPECS
    csize=0.6;1.0
    position=[0.05,0.01,0.84,0.95]
    xsize=4.2 & ysize=2.9
    xtitle='Longitude'
    ytitle='Latitude'

  ;AXES
    ix=reform(lon[*,0])
    iy=reform(lat[0,*])
    x=ix;lon
    y=iy;lat
;    if ~keyword_set(xrange) then $
;      xrange=[min(x),max(x)]
;    if ~keyword_set(xrange) then $
;      yrange=[min(y),max(y)]

  set_plot,'ps'
  epsname=figspecs.figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

;  area = [min(y),min(x),max(y),max(x)]
area=[0,-80,29,-25]

  map_set,0,0,/mercator,/isotropic,limit=area,xmargin=0,ymargin=0,position=position;,color=254;,$
;    charsize=csize,title=figspecs.title

  loadct,figspecs.col_table,/silent

  ;FILL SHADING
    for i=0,1 do $
      contour,var*figspecs.scale,x,y,/cell_fill,/overplot,$
        levels=figspecs.levels,c_colors=figspecs.colors

  loadct,0,/silent

  ;CONTOURS
    if keyword_set(cvar) then begin
      contour,cvar,x,y,/follow,/overplot,$
        levels=figspecs.clevs,c_colors=0,c_thick=1,c_charsize=csize*0.7;,c_labels=replicate(0,30)
    endif ;else if icont then $
;      contour,cvar,x,y,/follow,/overplot,$
;        levels=clevels,c_charsize=csize*0.3;0.7

  ;OVERLAY TITLE
    xyouts,0.5,0.95,figspecs.title,align=0.5,charsize=csize*1.,/normal

  ;LAND
    map_continents,/coasts,limit=area,color=figspecs.landcol,mlinethick=1.0,/hires
;    map_continents,color=0,/coasts,/hires,mlinethick=0.1

  ;WIND
    if keyword_set(wind) then begin
      xskip=5
      yskip=5
      nx=n_elements(ix) & ny=n_elements(iy)
      xplt=lindgen(500)*xskip & xplt=xplt[where(xplt lt nx)]
      yplt=lindgen(500)*yskip & yplt=yplt[where(yplt lt ny)]
      u=wind.u[xplt,*] & u=u[*,yplt]
      v=wind.v[xplt,*] & v=v[*,yplt]
      velovect,u,v,ix[xplt],iy[yplt],/overplot,length=1.3
    endif

  ;TC SYMBOL
  if keyword_set(tcloc) then begin
     psav=!P.Font
     !P.Font=-1
     xyouts,tcloc[0],tcloc[1],'!20<!X',charsize=0.5,charthick=2.5,align=0.5,color=figspecs.landcol,/data
     !P.Font=psav
;    tx=text(tcloc[0],tcloc[1],"<",font_name="Hershey 20",/overplot)
  endif

;WRF MODEL DOMAIN
  if iwrf_dom then begin
    plots,wrf_dom[[0,2]],replicate(wrf_dom[1],2),linestyle=0,thick=2,color=figspecs.landcol,/data
    plots,wrf_dom[[0,2]],replicate(wrf_dom[3],2),linestyle=0,thick=2,color=figspecs.landcol,/data
    plots,replicate(wrf_dom[0],2),wrf_dom[[1,3]],linestyle=0,thick=2,color=figspecs.landcol,/data
    plots,replicate(wrf_dom[2],2),wrf_dom[[1,3]],linestyle=0,thick=2,color=figspecs.landcol,/data
  endif

  ;LATLON GRID
;;    lats=indgen(5)*10-20 & latnames=['10!9%!XS','Eq.','10!9%!XN']
;;    lons=indgen(8)*10+50 & lonnames=strtrim(lons,2)+'!9%!XE'
;;    ;for i=1,14,2 do lonnames[i]=''
;;    map_grid,/box_axes,lats=lats,lons=lons,charsize=csize*0.8,color=0
    map_grid,/box_axes,latdel=5,londel=5,color=0,charsize=csize*0.7,glinethick=1.5,glinestyle=1,/no_grid

  ;COLOR BAR
  if figspecs.icbar then begin
    cpos= [ position[2]+0.045 ,$
            position[1]+0.3 ,$
            position[2]+0.06 ,$
            position[3]-0.3 ]
    loadct,figspecs.col_table,/silent
    colorbar2, colors=figspecs.colors, range=[min(figspecs.levels),max(figspecs.levels)],divisions=figspecs.ndivs,$
      charsize=csize*0.7, position=cpos, /right, /vertical, title=figspecs.cbar_tag, $
      annotatecolor='black',format=figspecs.cbar_format
    loadct,0,/silent
  endif

  device,/close

  convert_png,figspecs.figname,/remove_eps,res=400

end
