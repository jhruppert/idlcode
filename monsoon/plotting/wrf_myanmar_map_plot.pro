; 
; Plot maps from DYNAMO WRF simulations.
; 
; James Ruppert
; 20.10.16
; 
pro wrf_myanmar_map_plot, dirs, var, lon, lat, figspecs, $
  cvar=cvar, wind=wind, box=box, bounds=bounds, $
  noscalewind=noscalewind, cross=cross

;CREATE FIGURE

  ;PLOT SPECS
    csize=0.75;6
    position=[0.05,0.00,0.84,0.95]
    xsize=4.2 & ysize=3.8
    xtitle='Longitude'
    ytitle='Latitude'

  ;AXES
    x=lon
    y=lat
    if ~keyword_set(xrange) then $
      xrange=[min(x),max(x)]
    if ~keyword_set(xrange) then $
      yrange=[min(y),max(y)]

  set_plot,'ps'
  epsname=figspecs.figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  area = [min(y),min(x),max(y),max(x)]
;  area=[5,74,26,103]
;  area=[4,67,30,104]
;  area=[0,64,30,108]
;  area=[6,78,29,104]

  map_set,0,0,/mercator,/isotropic,limit=area,xmargin=0,ymargin=0,position=position;,color=254;,$
;    charsize=csize,title=figspecs.title

  loadct,figspecs.col_table,/silent,file=dirs.ctfil

  ;FILL SHADING
    for i=0,1 do $
      contour,var*figspecs.scale,x,y,/cell_fill,/overplot,$
        levels=figspecs.levels,c_colors=figspecs.colors

  loadct,0,/silent

  ;CONTOURS
    if keyword_set(cvar) then begin
      contour,cvar.cvar,cvar.x,cvar.y,/follow,/overplot,$
        levels=figspecs.clevs,c_colors=0,c_thick=2,c_charsize=csize;,c_labels=replicate(0,30)
      contour,cvar.cvar,cvar.x,cvar.y,/follow,/overplot,$
        levels=-1*reverse(figspecs.clevs),c_colors=0,c_thick=2,c_linestyle=2,c_charsize=csize;,c_labels=replicate(0,30)
    endif ;else if icont then $
;      contour,cvar,x,y,/follow,/overplot,$
;        levels=clevels,c_charsize=csize

  ;OVERLAY TITLE
;    xyouts,0.5,0.93,figspecs.title,align=0.5,charsize=csize,/normal
    xyouts,mean(area[[1,3]]),area[2]+2.3,figspecs.title,align=0.5,charsize=csize,/data

  ;LAND
    landcol=0
;    landcol=255
    landthick=2.4
    hires=1;0
    if hires then landthick=0.8
    map_continents,/coasts,/countries,limit=area,color=landcol,mlinethick=landthick*0.35,hires=hires
    ;FILL LAND
    do_fill=0
;    do_fill=1
    if do_fill then begin
      map_continents,/coasts,limit=area,color=254,mlinethick=landthick,hires=hires,fill_continents=1
      map_continents,/coasts,limit=area,color=landcol,mlinethick=landthick,hires=hires
    endif else $
      map_continents,/coasts,limit=area,color=landcol,mlinethick=landthick,hires=hires

;    map_continents,color=0,/coasts,hires=hires,mlinethick=0.1

  ;CROSS SECTIONS
;  if keyword_set(cross) then $
;  xcross=[77.8,96.2]; -13. ; lon,lat
;  ycross=[12.,21.5] ;- 2.5  ; lon,lat
;  cross=[xcross,ycross]
;    plots,cross[0:1],cross[2:3],linestyle=2,thick=2,/data

  ;PLOT BOX
    if keyword_set(box) then begin
      thick=2
      if ~keyword_set(bounds) then $
      ;COASTAL INDICES BOX
      bounds=[81.5,9.5,97,23.5]
bounds=[90.,15.5,95.5,22.5] ; Northern Myanmar coastline
bounds=[69.5,8.,77.5,20.] ; Western Ghats
      ;POWER SPECTRUM BOX
  ;    bounds=[85.,13,96,23.5]
      if size(bounds,/n_dim) eq 2 then begin
        thick=3
        nb=(size(bounds,/dim))[1]
        for ib=0,nb-1 do begin
          ibounds=reform(bounds[*,ib])
          for i=0,2,2 do plots,ibounds[[i,i]],ibounds[[1,3]],thick=thick,linestyle=ib,/data
          for i=1,3,2 do plots,ibounds[[0,2]],ibounds[[i,i]],thick=thick,linestyle=ib,/data
          xyouts,ibounds[0]+0.15,ibounds[1]+0.2,'Box '+strtrim(ib+1,2),align=0,charsize=csize,/data
       endfor
      endif else begin
        for i=0,2,2 do plots,bounds[[i,i]],bounds[[1,3]],thick=thick,/data
        for i=1,3,2 do plots,bounds[[0,2]],bounds[[i,i]],thick=thick,/data
      endelse
    endif

  ;WIND
    if keyword_set(wind) then begin
      ix=where((wind.x ge area[1] and wind.x le area[3]),nx)
      iy=where((wind.y ge area[0] and wind.y le area[2]),ny)
      nskip=round(1.*nx/30.)
      xplt=lindgen(500)*nskip & xplt=xplt[where(xplt lt nx)]
      yplt=lindgen(500)*nskip & yplt=yplt[where(yplt lt ny)]
      ilon=wind.x[ix[xplt]] & ilat=wind.y[iy[yplt]]
      u=wind.u[ix[xplt],*] & u=u[*,iy[yplt]]
      v=wind.v[ix[xplt],*] & v=v[*,iy[yplt]]
      ;SCALE VECTORS BY A CONSTANT
        maxw=10. ; Speed for vector of unit size
        ws=max(sqrt(u^2+v^2))
      scale=ws/maxw
      if keyword_set(noscalewind) then scale=1.
      length=1.0*scale
      thick=1.2
      velovect,u,v,ilon,ilat,/overplot,length=length,thick=thick
    endif

      ;LATLON GRID
;;    lats=indgen(5)*10-20 & latnames=['10!9%!XS','Eq.','10!9%!XN']
;;    lons=indgen(8)*10+50 & lonnames=strtrim(lons,2)+'!9%!XE'
;;    ;for i=1,14,2 do lonnames[i]=''
;;    map_grid,/box_axes,lats=lats,lons=lons,charsize=csize,color=0
      map_grid,/box_axes,latdel=10,londel=10,color=0,charsize=csize,glinethick=1.5,glinestyle=1,/no_grid

      ;COLOR BAR
    if figspecs.icbar then begin
      dy=0.3
      ;dy=0.35
      cpos= [ position[2]+0.048 ,$
      position[1]+dy ,$
      position[2]+0.064 ,$
      position[3]-dy ]
    loadct,figspecs.col_table,/silent,file=dirs.ctfil
;    if strmatch(figspecs.figname,'*wspd*') then setlevs=strtrim(figspecs.levels,2) else setlevs=''
;    if strmatch(figspecs.figname,'*rainrate*') then setlevs=strtrim(figspecs.levels,2) else setlevs=''
    colorbar2, colors=figspecs.colors, range=[min(figspecs.levels),max(figspecs.levels)],divisions=figspecs.ndivs,$
      charsize=csize*0.7, position=cpos, /right, /vertical, title=figspecs.cbar_tag,$
      annotatecolor='black',format=figspecs.cbar_format,$
      setlevels=setlevs
    loadct,0,/silent
    endif

  device,/close

  convert_png,figspecs.figname,res=300,/remove_eps

end
