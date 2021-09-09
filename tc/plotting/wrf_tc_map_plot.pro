; 
; Plot maps from DYNAMO WRF simulations.
; 
; James Ruppert
; 20.10.16
; 
pro wrf_tc_map_plot, dirs, var, lon, lat, figspecs, $
  cvar=cvar, wind=wind, loc_pmin=loc_pmin, itcloc=itcloc, $
  shr=shr

imaria=1 ; zoom in around Maria or Haiyan?

;CREATE FIGURE

  ;PLOT SPECS
    csize=1.3;0.8
;    if imaria then csize=1.1
    position=[0.05,0.00,0.84,0.95]
    xsize=4.2 ;& ysize=3.9
    ysize=3.0
    xtitle='Longitude'
    ytitle='Latitude'

  ;AXES
    x=lon
    y=lat
    if ~keyword_set(xrange) then $
      xrange=[min(x),max(x)]
    if ~keyword_set(xrange) then $
      yrange=[min(y),max(y)]

  area = [min(y),min(x),max(y),max(x)]

  if imaria then $
;;    area=[10,-69,19.8,-53]
;    area=[11,-68,21,-52] ; WSPD HISTORY ZOOM
;    area=[12.5,-64.5,19.5,-55] ; WSPD HISTORY ZOOM (SMALL)
    area=[12.5,-68,20,-55.5] ; PNAS PAPER
    ;FOR HAIYAN
    if strmatch(figspecs.figname,'*haiy*') then begin
      area[[0,2]]-=5.1       ; PNAS PAPER - TURN IMARIA ON
      area[[1,3]]+=127.+63.8
      ;EXPAND FOR HAIYAN
;        dx=area[3]-area[1]
;        dy=area[2]-area[0]
;        fac=1.4
;        area[0]=6 & area[2]=area[0]+dy*fac
;        area[1]=min(x)+0.2 & area[3]=area[1]+dx*fac
    endif
;    area=[7,-73.2,24,-49]
;    area=[10,-70,22,-52] ; OLR ZOOM

  ;AREA BASED ON TC LOC (COMMENT THIS OUT FOR WSPD HISTORY)
;    if imaria and keyword_set(loc_pmin) then begin
;;      dy=8.5
;      dy=4.5; OLR SMALL ANIM
;;      dy=5.0
;;      dy=7. ; OLR LARGE ANIM
;      dx=dy/0.68
;      area=[loc_pmin[1]-dy,loc_pmin[0]-dx,loc_pmin[1]+dy,loc_pmin[0]+dx]
;    endif

  set_plot,'ps'
  epsname=figspecs.figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

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
      contour,cvar,x,y,/follow,/overplot,$
        levels=figspecs.clevs,c_colors=0,c_thick=2;,c_labels=replicate(0,30)
    endif ;else if icont then $
;      contour,cvar,x,y,/follow,/overplot,$
;        levels=clevels,c_charsize=csize*0.3;0.7

  ;BOX AROUND PLOT
;This is problematic due to the inexact plotting bounding box.
;    for i=0,1 do begin
;      plots,xrange,replicate(yrange[i],2),linestyle=0,color=0,thick=1,/data
;      plots,replicate(xrange[i],2),yrange,linestyle=0,color=0,thick=1,/data
;    endfor

  ;OVERLAY TITLE
;    xyouts,0.5,0.93,figspecs.title,align=0.5,charsize=csize*1.,/normal
    xyouts,mean(area[[1,3]]),area[2]+1.3,figspecs.title,align=0.5,charsize=csize*0.7,/data

  ;LAND
    landcol=0
;    landcol=255
    landthick=1.5;1.8
    map_continents,/coasts,limit=area,color=landcol,mlinethick=landthick,/hires
;    map_continents,color=0,/coasts,/hires,mlinethick=0.1

  ;WIND
    if keyword_set(wind) then begin
      xskip=30;10
      yskip=30;10
      nx=n_elements(x) & ny=n_elements(y)
      xplt=lindgen(500)*xskip & xplt=xplt[where(xplt lt nx)]
      yplt=lindgen(500)*yskip & yplt=yplt[where(yplt lt ny)]
      u=wind.u[xplt,*] & u=u[*,yplt]
      v=wind.v[xplt,*] & v=v[*,yplt]
      velovect,u,v,x[xplt],y[yplt],/overplot,length=1.3
    endif

  ;SHEAR VECTOR
    if keyword_set(shr) then begin
      print,'Shear magnitude: ',sqrt(shr[0]^2+shr[1]^2)
      mag=sqrt(shr[0]^2+shr[1]^2)
      magset=2.2
      fac=magset/mag
      arrow,loc_pmin[0],loc_pmin[1],loc_pmin[0]+shr[0]*fac,loc_pmin[1]+shr[1]*fac,thick=2.5,/overplot,/data;,hsize=!D.X_SIZE/60.
    endif


  ;LOCATION OF SLP MINIMIUM
  if itcloc and keyword_set(loc_pmin) then begin
    loadct,11,/silent
;    plots,lon[loc_pmin[0]],lat[loc_pmin[1]],psym=7,symsize=0.9,thick=3.5,color=250,/data
    plots,loc_pmin[0],loc_pmin[1],psym=7,symsize=0.9,thick=3.5,color=250,/data
    loadct,0,/silent
  endif

  ;LATLON GRID
;;;    lats=indgen(5)*10-20 & latnames=['10!9%!XS','Eq.','10!9%!XN']
;;;    lons=indgen(8)*10+50 & lonnames=strtrim(lons,2)+'!9%!XE'
;;;    ;for i=1,14,2 do lonnames[i]=''
;;;    map_grid,/box_axes,lats=lats,lons=lons,charsize=csize*0.8,color=0
;    map_grid,/box_axes,latdel=5,londel=5,color=0,charsize=csize*0.60,glinethick=1.5,glinestyle=1,/no_grid

  ;COLOR BAR
  if figspecs.icbar then begin
    cpos= [ position[2]+0.048 ,$
            position[1]+0.30 ,$
            position[2]+0.064 ,$
            position[3]-0.30 ]
    loadct,figspecs.col_table,/silent,file=dirs.ctfil
;    if strmatch(figspecs.figname,'*wspd*') then setlevs=strtrim(figspecs.levels,2) else setlevs=''
    if strmatch(figspecs.figname,'*rainrate*') then setlevs=strtrim(figspecs.levels,2) else setlevs=''
    colorbar2, colors=figspecs.colors, range=[min(figspecs.levels),max(figspecs.levels)],divisions=figspecs.ndivs,$
      charsize=csize*0.6, position=cpos, /right, /vertical, title=figspecs.cbar_tag,$
      annotatecolor='black',format=figspecs.cbar_format,$
      setlevels=setlevs
    loadct,0,/silent
  endif

  device,/close

  convert_png,figspecs.figname,res=400,/remove_eps
;  convert_png,figspecs.figname,res=300,/remove_eps

end
