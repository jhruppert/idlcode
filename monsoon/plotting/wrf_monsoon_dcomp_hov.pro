; 
; Plot diurnal composite Hovmoller.
;
; James Ruppert
; 10/13/20
; 
pro wrf_monsoon_dcomp_hov, dirs, figspecs, var, time, xsav, wind=wind, cvar=cvar

nt=n_elements(time)

  ;PLOT SPECS
    csize=0.6
    position=[0.16,0.13,0.84,0.93]
    xsize=2.7 & ysize=3.2
    ytitle='Hour [ LT ]'
    xtitle='Distance from coast [ km ]'

  ;AXES
    npd=(size(var,/dim))[1]
    y=findgen(nt)*24/npd
    x=xsav
    nx=n_elements(x)
    xrange=[x[0],x[nx-1]]
    if ~keyword_set(xrange) then $
      xrange=[min(x),max(x)]
    if ~keyword_set(yrange) then $
      yrange=[min(y),max(y)]

  ;REPEAT DIURNAL CYCLE
    nt=npd*2+1
    y=findgen(nt)*48/(nt-1)
    yrange=[min(y),max(y)]
    var=[[var],[var],[var[*,0]]]
    yticks=16
    ytickv=indgen(yticks+1)*3
    ytickname=strtrim([indgen(8),indgen(9)]*3,2)

  set_plot,'ps'
  epsname=figspecs.figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  plot,x,y,/nodata,position=position,$
    xstyle=9,ystyle=9,$
    yticklen=-0.022,xticklen=-0.018,$
    yticks=yticks,ytickv=ytickv,ytickname=ytickname,$
    xrange=xrange,yrange=yrange,yminor=2,$
    xtitle=xtitle,ytitle=ytitle,$
    charsize=csize;,$
    ;title=figspecs.title

  loadct,figspecs.col_table,/silent,file=dirs.ctfil

  ;FILL SHADING
    for i=0,1 do $
      contour,var*figspecs.scale,x,y,/cell_fill,/overplot,$
        levels=figspecs.levels,c_colors=figspecs.colors

  loadct,0,/silent

  ;CONTOURS
    if keyword_set(cvar) then begin
      message,'still need to contour settings'
      contour,cvar,x,y,/follow,/overplot,$
        levels=clevs,c_colors=0,c_thick=1.5,c_charsize=0.9*csize;,c_labels=replicate(1,200)
      contour,cvar,x,y,/follow,/overplot,c_linestyle=1,$
        levels=-1.*reverse(clevs),c_colors=0,c_thick=1.5,c_charsize=0.9*csize;,c_labels=replicate(1,200)
    endif

  ;COASTLINE
  plots,[0,0],!y.crange,/data

  ;BOX AROUND PLOT
    for i=0,1 do begin
      plots,xrange,replicate(yrange[i],2),linestyle=0,color=0,thick=1,/data
      plots,replicate(xrange[i],2),yrange,linestyle=0,color=0,thick=1,/data
    endfor

  ;OVERLAY TITLE
    xyouts,mean(position[[1,3]]),position[3]+0.03,figspecs.title,align=0.5,charsize=csize,/normal

  ;COLOR BAR
  if figspecs.icbar then begin
    ybuff=0.21
    cpos= [ position[2]+0.020 ,$
            position[1]+ybuff ,$
            position[2]+0.039 ,$
            position[3]-ybuff ]
    loadct,figspecs.col_table,/silent,file=dirs.ctfil
    colorbar2, colors=figspecs.colors, range=[min(figspecs.levels),max(figspecs.levels)],divisions=figspecs.ndivs,$
      charsize=csize*0.9, position=cpos, /right, /vertical, title=figspecs.cbar_tag,$
      annotatecolor='black',format=figspecs.cbar_format;,$
      ;setlevels=setlevs
    loadct,0,/silent
  endif

  device,/close

  convert_png,figspecs.figname,res=300,/remove_eps


end
