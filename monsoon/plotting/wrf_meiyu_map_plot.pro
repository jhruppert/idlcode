; 
; Plot maps from DYNAMO WRF simulations.
; 
; James Ruppert
; 20.10.16
; 
pro wrf_meiyu_map_plot, dirs, var, lon, lat, figspecs, $
  cvar=cvar, wind=wind 

;CREATE FIGURE

  ;PLOT SPECS
    csize=0.75
    position=[0.05,0.00,0.84,0.95]
    xsize=4.2 & ysize=4.1
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

;    area=[10,-69,19.8,-53]

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
;        levels=clevels,c_charsize=csize

  ;OVERLAY TITLE
;    xyouts,0.5,0.93,figspecs.title,align=0.5,charsize=csize,/normal
    xyouts,mean(area[[1,3]]),area[2]+2.3,figspecs.title,align=0.5,charsize=csize,/data

  ;LAND
    landcol=0
;    landcol=255
    landthick=1.5
    map_continents,/coasts,/countries,limit=area,color=landcol,mlinethick=landthick,/hires
;    map_continents,color=0,/coasts,/hires,mlinethick=0.1

  ;WIND
    if keyword_set(wind) then begin
      nx=n_elements(x) & ny=n_elements(y)
      nskip=round(1.*nx/30.)
      xplt=lindgen(500)*nskip & xplt=xplt[where(xplt lt nx)]
      yplt=lindgen(500)*nskip & yplt=yplt[where(yplt lt ny)]
      ix=x[xplt] & iy=y[yplt]
      u=wind.u[xplt,*] & u=u[*,yplt]
      v=wind.v[xplt,*] & v=v[*,yplt]
      ;SCALE VECTORS BY A CONSTANT
        maxw=10 ; Speed for vector of unit size
        ws=sqrt(u^2+v^2)
        u*=ws/maxw
        v*=ws/maxw
      length=1.3
      thick=1.2
      velovect,u,v,ix,iy,/overplot,length=length,thick=thick
    endif

      ;LATLON GRID
      ;;    lats=indgen(5)*10-20 & latnames=['10!9%!XS','Eq.','10!9%!XN']
      ;;    lons=indgen(8)*10+50 & lonnames=strtrim(lons,2)+'!9%!XE'
      ;;    ;for i=1,14,2 do lonnames[i]=''
      ;;    map_grid,/box_axes,lats=lats,lons=lons,charsize=csize,color=0
      map_grid,/box_axes,latdel=5,londel=5,color=0,charsize=csize,glinethick=1.5,glinestyle=1,/no_grid

      ;COLOR BAR
      if figspecs.icbar then begin
      cpos= [ position[2]+0.048 ,$
      position[1]+0.30 ,$
      position[2]+0.064 ,$
      position[3]-0.30 ]
    loadct,figspecs.col_table,/silent,file=dirs.ctfil
;    if strmatch(figspecs.figname,'*wspd*') then setlevs=strtrim(figspecs.levels,2) else setlevs=''
;    if strmatch(figspecs.figname,'*rainrate*') then setlevs=strtrim(figspecs.levels,2) else setlevs=''
    colorbar2, colors=figspecs.colors, range=[min(figspecs.levels),max(figspecs.levels)],divisions=figspecs.ndivs,$
      charsize=csize, position=cpos, /right, /vertical, title=figspecs.cbar_tag,$
      annotatecolor='black',format=figspecs.cbar_format,$
      setlevels=setlevs
    loadct,0,/silent
  endif

  device,/close

  convert_png,figspecs.figname,res=200,/remove_eps

end
