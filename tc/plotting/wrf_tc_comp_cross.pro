; 
; Plot cross sections from DYNAMO WRF simulations.
; 
; James Ruppert
; 20.10.16
; 
pro wrf_tc_comp_cross, var, x, pres, figspecs, $
  var2=var2, cvar=cvar, xtitle=xtitle, clevs=clevs, $
  olr=olr, csol=csol, irad=irad

ishallow=0 ; zoomed into PBL?

;CREATE FIGURE

  ;PLOT SPECS
    ylog=1
    csize=0.85
    position=[0.14,0.12,0.85,0.97]
    xsize=4.8 & ysize=3.7
    ytitle='Pressure [ hPa ]'

  ;DUAL PLOT
    dy_lower=0.14
    buff=0.04
    position2=position
    position[1]+=dy_lower+buff
    position2[3]=position[1]-buff

  ;AXES
    y=pres
    ;xrange=[min(x),max(x)]
    nx=n_elements(x)
    fin=where(finite(x))
    if ~keyword_set(xrange) then $
      xrange=[x[fin[0]],x[max(fin)]]
    if ~keyword_set(yrange) then $
      yrange=[max(y),min(y)]

  ytickv=[1000,700,500,400,300,200,100];,50]
  if ishallow then ytickv=[1000,925,850,800];,700]
;  ytickv=[1000,700,500,400,300,200];,100];,50]
  yrange=[max(y),min(ytickv)]
  yticks=n_elements(ytickv)-1

  set_plot,'ps'
  if ishallow then figspecs.figname+='_shall'
  epsname=figspecs.figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  !p.multi=[2,0,2]

  loadct,0,/silent

  plot,x,y,/nodata,position=position,ylog=ylog,$
    xstyle=9,ystyle=9,xminor=4,yticks=yticks,ytickv=ytickv,$
    yticklen=-0.011,xticklen=-0.022,$
    xrange=xrange,yrange=yrange,$
    ;xtitle=xtitle,
    xtickname=replicate(' ',10),$
    ytickname=ytickname,yminor=0,$
    ytitle=ytitle,charsize=csize;,title=ititle

;for i=0,1 do $
;var=smooth(var,3,/edge_truncate,/nan)

  ;UNDERLAY EXTREMES FOR RAD
    if irad then begin
      ;6-16
;        nlevs_add=11;6
;        levs_add=findgen(nlevs_add)-16
      ;4-10
        nlevs_add=7
        levs_add=findgen(nlevs_add)-10
      ;5-15
;        nlevs_add=11
;        levs_add=findgen(nlevs_add)-15
      cols_add=findgen(nlevs_add)/(nlevs_add-1)*255
      ;cols_add=reverse(cols_add)
      ctadd=0
      loadct,ctadd,/silent
vartmp=var*figspecs.scale
vartmp[where(vartmp lt min(levs_add))]=min(levs_add)
      for i=0,1 do $
        contour,vartmp,x,y,/cell_fill,/overplot,$
          levels=levs_add,c_colors=cols_add
    endif

  loadct,figspecs.col_table,/silent

  ;FILL SHADING
    for i=0,1 do $
      contour,var*figspecs.scale,x,y,/cell_fill,/overplot,$
        levels=figspecs.levels,c_colors=figspecs.colors

  loadct,0,/silent

  ;CONTOURS
    if keyword_set(cvar) then begin
      if ~keyword_set(clevs) then clevs=figspecs.clevs
      contour,cvar,x,y,/follow,/overplot,$
        levels=clevs,c_colors=0,c_thick=1.5,c_labels=replicate(1,200),c_charsize=0.9*csize
      contour,cvar,x,y,/follow,/overplot,c_linestyle=1,$
        levels=-1.*reverse(clevs),c_colors=0,c_thick=1.5,c_labels=replicate(1,200),c_charsize=0.9*csize
    endif ;else if icont then $
;      contour,cvar,x,y,/follow,/overplot,$
;        levels=clevels,c_charsize=csize*0.3;0.7

  ;OLR
;    if keyword_set(olr) then begin
;      axis,yaxis=1,ystyle=1,yminor=1,charsize=csize,yrange=[-200,350],ylog=0,/save
;      oplot,x,olr,linestyle=0,thick=4,color=0
;    endif

  ;BOX AROUND PLOT
    for i=0,1 do begin
      plots,xrange,replicate(yrange[i],2),linestyle=0,color=0,thick=1,/data
      plots,replicate(xrange[i],2),yrange,linestyle=0,color=0,thick=1,/data
    endfor

  ;OLR

    iolr=0
    if iolr then begin
      ytitle='OLR!C[ W m!U-2!N ]'
      yrange=[60,310]
      ytickv=[100,200,300]
      yticks=2
      ylog=0
    endif else begin
      ytitle='N'
      yrange=[1e1,1e6]
      ytickv=[1e1,1e3,1e5]
      yticks=2
      ylog=1
    endelse

    plot,x,y,/nodata,position=position2,ylog=ylog,$
      xstyle=9,ystyle=9,xminor=4,$
      yticklen=0.015,xticklen=0.05,$
      xrange=xrange,yticks=yticks,ytickv=ytickv,yrange=yrange,$
      xtitle=xtitle,yminor=4,$
      ytitle=ytitle,charsize=csize;,title=ititle

    oplot,x,olr,linestyle=0,thick=3,color=0
    if keyword_set(csol) then oplot,x,csol,linestyle=1,thick=3,color=0

  ;COLOR BAR
  if figspecs.icbar then begin

    loadct,figspecs.col_table,/silent
    cpos= [ position[2]+0.018 ,$
            position[1]+0.1 ,$
            position[2]+0.031 ,$
            position[3]-0.1 ]

    if irad then $
      cpos[1]+=0.15

;    if iw then begin
;      dy=position[3]-position[1]
;      icpos= [ cpos[0] , position[1]+0.51*dy , cpos[2] , cpos[3] ]
;      colorbar2, colors=colors_up, range=[min(levels_up),max(levels_up)],$
;        charsize=csize*1., position=icpos, /right, /vertical, title=figspecs.cbar_tag, $
;        DIVISIONS=ndivs_up,annotatecolor='black',format='(i2)'
;      icpos= [ cpos[0] , cpos[1]             , cpos[2] , position[3]-0.51*dy ]
;      colorbar2, colors=colors_dn, range=[min(levels_dn),max(levels_dn)],$
;        charsize=csize*1., position=icpos, /right, /vertical, title=figspecs.cbar_tag, $
;        DIVISIONS=ndivs_dn,annotatecolor='black',format='(i2)'
;    endif else begin
    colorbar2, colors=figspecs.colors, range=[min(figspecs.levels),max(figspecs.levels)],divisions=figspecs.ndivs,$
      charsize=csize*0.9, position=cpos, /right, /vertical, title=figspecs.cbar_tag, $
      annotatecolor='black',format=figspecs.cbar_format
;    endelse

    if irad then begin
      loadct,ctadd,/silent
      cpos[1]=cpos[1]-0.15
      cpos[3]=cpos[1]+0.13
      colorbar2, colors=cols_add, range=[min(levs_add),max(levs_add)],divisions=2,$
        charsize=csize*0.9, position=cpos, /right, /vertical, title=' ', $
        annotatecolor='black',format='(i3)';figspecs.cbar_format
    endif

    loadct,0,/silent

  endif

  device,/close

  !p.multi=0

  convert_png,figspecs.figname,/remove_eps,res=400

end
