; 
; Plot maps from DYNAMO WRF simulations.
; 
; James Ruppert
; 20.10.16
; 
pro wrf_monthly_myanmar_plot, dirs, rain_all, rain_sum, u_all, v_all, u_sum, v_sum, $
  lon, lat, eralon, eralat, figspecs, $
  cvar=cvar;, wind=wind, box=box, bounds=bounds, $
;  noscalewind=noscalewind, cross=cross, zoomin=zoomin

;CREATE FIGURE

  ;PLOT SPECS
    csize=0.5
    position=[0.03,0.00,0.92,0.95]
    xsize=8. & ysize=2.5

  ;MULTI-PANEL
    position=[0.02,0.06,0.935,0.94]
    dx=position[2]-position[0]
    dy=position[3]-position[1]
    ygap=0.08;0.11
    xmid=9./16
    xgap=0.025

    position1 = [position[0], position[1]+dy/2+ygap/2, position[0]+dx*xmid, position[3]]
    position2 = [position[0], position[1],             position[0]+dx*xmid, position[1]+dy/2-ygap/2]

    position3 = [position[0]+dx*xmid+xgap, position[1], position[2], position[3]]

    position=fltarr(3,4)
    position[0,*]=position1
    position[1,*]=position2
    position[2,*]=position3

  ;AXES
    x=lon
    y=lat
    if ~keyword_set(xrange) then $
      xrange=[min(x),max(x)]
    if ~keyword_set(xrange) then $
      yrange=[min(y),max(y)]

  !p.multi=[3,3,1]
;  !p.multi=[2,2,1]

  set_plot,'ps'
  epsname=figspecs.figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

for iplot=0,2 do begin

  if iplot eq 0 then begin
    rain=rain_all
    u=u_all
    v=v_all
  endif else if iplot ge 1 then begin
    rain=rain_sum
    u=u_sum
    v=v_sum
  endif

  loadct,0,/silent

  if iplot lt 2 then begin
    area = [min(lat),min(lon),max(lat),max(lon)]
    area[[0,2]]=[-38,38]
    y0=0
    x0=0;80
  endif else begin
    area = [0.,60,30,109]
    y0=0
    x0=mean(area[[1,3]])
  endelse

  map_set,y0,x0,/mercator,limit=area,xmargin=0,ymargin=0,position=reform(position[iplot,*]),/noerase;,/isotropic

  loadct,figspecs.col_table,/silent,file=dirs.ctfil

  ;FILL SHADING
    for i=0,1 do $
      contour,rain,x,y,/cell_fill,/overplot,$
        levels=figspecs.levels,c_colors=figspecs.colors

  loadct,0,/silent

  ;OVERLAY TITLE
;    xyouts,mean(area[[1,3]]),area[2]+8,figspecs.title,align=0.5,charsize=csize,/data

  ;LAND
    landcol=0
    landthick=1.3
;    map_continents,/coasts,limit=area,color=landcol,mlinethick=landthick;,/hires
    if iplot eq 2 then begin
      landthick=1.1
      map_continents,/coasts,limit=area,color=landcol,fill_continents=0,mlinethick=landthick,/hires
;      map_continents,/continents,limit=area,color=255,fill_continents=1,mlinethick=landthick,/hires
;      map_continents,/continents,limit=area,color=landcol,fill_continents=0,mlinethick=landthick,/hires
    endif else $
      map_continents,/continents,limit=area,color=landcol,mlinethick=landthick;,/hires

  ;CONTOURS
    if iplot eq 2 and keyword_set(cvar) then begin
ccolor=254;0
loadct,11,/silent
ccolor=220
      contour,cvar.cvar,cvar.x,cvar.y,/follow,/overplot,$
        levels=figspecs.clevs,c_colors=ccolor,c_thick=2,c_charsize=csize;,c_labels=replicate(0,30)
;      contour,cvar.cvar,cvar.x,cvar.y,/follow,/overplot,$
;        levels=-1*reverse(figspecs.clevs),c_colors=0,c_thick=2,c_linestyle=2,c_charsize=csize;,c_labels=replicate(0,30)
loadct,0,/silent
    endif

  ;PLOT BOX
;    if iplot eq 2 then begin
;      thick=2
;      bounds=[81.5,9.5,97,23.5]
;      for i=0,2,2 do plots,bounds[[i,i]],bounds[[1,3]],thick=thick,/data
;      for i=1,3,2 do plots,bounds[[0,2]],bounds[[i,i]],thick=thick,/data
;    endif

  ;WIND
  iwind=1
    if iwind then begin
      ix=where((eralon ge area[1] and eralon le area[3]),nx)
      iy=where((eralat ge area[0] and eralat le area[2]),ny)
      nskip=round(1.*nx/50.)
  if iplot eq 2 then nskip=9
      xplt=lindgen(500)*nskip & xplt=xplt[where(xplt lt nx)]
      yplt=lindgen(500)*nskip & yplt=yplt[where(yplt lt ny)]
      ilon=eralon[ix[xplt]] & ilat=eralat[iy[yplt]]
      u=u[ix[xplt],*] & u=u[*,iy[yplt]]
      v=v[ix[xplt],*] & v=v[*,iy[yplt]]
      ;SCALE VECTORS BY A CONSTANT
        maxw=10. ; Speed for vector of unit size
        ws=max(sqrt(u^2+v^2))
      scale=ws/maxw
      if keyword_set(noscalewind) then scale=1.
      length=1.0*scale
      thick=1.2
      velovect,u,v,ilon,ilat,/overplot,length=length,thick=thick
print,'Max wind:',max(sqrt(u^2+v^2))
    endif

      ;LATLON GRID
      latdel=30
      londel=60
      if iplot eq 2 then begin
        latdel=10
        londel=15
      endif
      map_grid,/box_axes,latdel=latdel,londel=londel,color=0,charsize=csize,glinethick=1.5,glinestyle=1,/no_grid

endfor ; iplot

      ;COLOR BAR
      cpos= [ position[2,2]+0.022 ,$
      position[2,1]+0.20 ,$
      position[2,2]+0.03 ,$
      position[2,3]-0.20 ]
    loadct,figspecs.col_table,/silent,file=dirs.ctfil
    colorbar2, colors=figspecs.colors, range=[min(figspecs.levels),max(figspecs.levels)],divisions=figspecs.ndivs,$
      charsize=csize*2, position=cpos, /right, /vertical, title=figspecs.cbar_tag,$
      annotatecolor='black',format=figspecs.cbar_format,$
      setlevels=setlevs
    loadct,0,/silent

  device,/close

  convert_png,figspecs.figname,res=300,/remove_eps

end
