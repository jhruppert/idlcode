; 
; Plot maps from DYNAMO WRF simulations with both Myanmar and IMERG side-by-side.
; 
; James Ruppert
; 20.10.16
; 
pro plot_myanmar_obcomp_maps, dirs, figname, s_wrf, s_obs, figspecs

iwind=0
if keyword_set(s_obs.u) then iwind=1

;CREATE FIGURE

  ;PLOT SPECS
    csize=0.75
    xsize=8.2 & ysize=3.8
    xtitle='Longitude'
    ytitle='Latitude'

    gap=0.05
    edge=0.05
    position1=[edge,0.0,0.5-gap/2,0.95]
    position2=[0.5+gap/2,0.0,1.-edge,0.95]

  ;AXES
    x=s_wrf.lon
    y=s_wrf.lat
    if ~keyword_set(xrange) then $
      xrange=[min(x),max(x)]
    if ~keyword_set(xrange) then $
      yrange=[min(y),max(y)]

  area = [min(y),min(x),max(y),max(x)]
;    area=[10,-69,19.8,-53]

  !p.multi=[2,2,2]

  set_plot,'ps'
  epsname=figspecs.figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica


for iplot=0,1 do begin

  if iplot eq 0 then position=position1 else position=position2

  loadct,0,/silent

  map_set,0,0,/mercator,/isotropic,limit=area,xmargin=0,ymargin=0,position=position,/NOERASE;,color=254;,$
;    charsize=csize,title=figspecs.title
;  plot,x,y,/nodata,xrange=area[[1,3]],xmargin=0,ymargin=0,position=position;,color=254;,$

  loadct,figspecs.col_table,/silent,file=dirs.ctfil

  if iplot eq 0 then begin
    ititle='JJAS13-17';'WRF'
    var=s_wrf.var*figspecs.scale
    if iwind then begin
      xskip=26
      u=s_wrf.u
      v=s_wrf.v
      wlon=s_wrf.lon
      wlat=s_wrf.lat
    endif
  endif else begin
    ititle='COMP(800km,12mm/d)';'IMERG/ERA5'
    x=s_obs.lon
    y=s_obs.lat
    var=s_obs.var*figspecs.scale
    if iwind then begin
      xskip=3
      u=s_obs.u
      v=s_obs.v
      wlon=s_obs.lonwind
      wlat=s_obs.latwind
    endif
  endelse

  ;FILL SHADING
    for i=0,1 do $
      contour,var,x,y,/cell_fill,/overplot,$
        levels=figspecs.levels,c_colors=figspecs.colors

  loadct,0,/silent

  ;OVERLAY TITLE
;    xyouts,0.5,0.93,figspecs.title,align=0.5,charsize=csize,/normal
    xyouts,mean(area[[1,3]]),area[2]+1.9,ititle,align=0.5,charsize=csize,/data

  ;LAND
    landcol=0
;    landcol=255
    landthick=1.5
    map_continents,/coasts,/countries,limit=area,color=landcol,mlinethick=landthick,/hires
;    map_continents,color=0,/coasts,/hires,mlinethick=0.1

  ;WIND
    if iwind then begin
      nx=n_elements(wlon) & ny=n_elements(wlat)
      nskip=round(1.*nx/30.)
      xplt=lindgen(500)*nskip & xplt=xplt[where(xplt lt nx)]
      yplt=lindgen(500)*nskip & yplt=yplt[where(yplt lt ny)]
      iu=u[xplt,*] & iu=iu[*,yplt]
      iv=v[xplt,*] & iv=iv[*,yplt]
      maxw=10
      ws=sqrt(iu^2+iv^2)
      iu*=ws/maxw
      iv*=ws/maxw
      length=1.3
      thick=1.2
      velovect,iu,iv,wlon[xplt],wlat[yplt],/overplot,length=length,thick=thick
    endif

      ;LATLON GRID
      ;;    lats=indgen(5)*10-20 & latnames=['10!9%!XS','Eq.','10!9%!XN']
      ;;    lons=indgen(8)*10+50 & lonnames=strtrim(lons,2)+'!9%!XE'
      ;;    ;for i=1,14,2 do lonnames[i]=''
      ;;    map_grid,/box_axes,lats=lats,lons=lons,charsize=csize,color=0
      map_grid,/box_axes,latdel=5,londel=5,color=0,charsize=csize,glinethick=1.5,glinestyle=1,/no_grid

endfor ; iplot

  ;OVERLAY TITLE
    xyouts,0.5,position1[3],figspecs.title,align=0.5,charsize=csize,/normal

      ;COLOR BAR
      if figspecs.icbar then begin
      cpos= [ position1[0]+0.02 ,$
      position2[1]+0.11 ,$
      position1[0]+0.032 ,$
      position2[1]+0.35 ]
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

  !p.multi=0

  convert_png,figspecs.figname,res=200,/remove_eps

end
