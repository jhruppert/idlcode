pro cm1_island_flux_tser,figname,var,x,y,$
  exp=exp

  col_table=70;13

  maxv=60.
  levels=findgen(255)/254*2*maxv-maxv
  colors=indgen(255)
  colors=reverse(colors)

;  cint=0.5 & nlevs=30
;  clevels=(findgen(nlevs)+1)*cint

  ;PLOT SETTINGS
  icbar=1 ; color bar?
  ileg=0 ; legend?
  icont=0 ; line contours?
  linethick=3
  csize=0.7
  position=[0.12,0.17,0.96,0.91]
  xticklen=-0.03 & yticklen=-0.012
  xsize=4.5 & ysize=2.5
  xtitle='x [ km ]'
  ytitle='Time [ d ]'

  ;AXES
  if ~keyword_set(xrange) then $
    xrange=[min(x),max(x)]
;    xrange=[299.0,301.5]
  if ~keyword_set(yrange) then $
    yrange=[min(y),max(y)]
;    yrange=[0,10]
;    yrange=[0,14]

  ;TIME SERIES

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  plot,x,y,/nodata,position=position,$
    xstyle=9,ystyle=9,color=0,$
    xtick_get=xtickv,ytick_get=ytickv,$
    xticklen=xticklen,yticklen=yticklen,yminor=4,$
    xrange=xrange,yrange=yrange,$;yticks=4,ytickname=strtrim([0,12,0,12,24],2),
    ;    xrange=[x03,x03+4],yrange=[0,60],$
    xtitle=xtitle,ytitle=ytitle,$
    charsize=csize,title='Sensible Heat Flux ('+strupcase(exp)+')'

  loadct,col_table,/silent

  ;FILL SHADING
  ;  for i=0,1 do $
  contour,var,x,y,/cell_fill,/overplot,$
    levels=levels,c_colors=colors

  ;  oplot,x,y,linestyle=2,thick=2.0,color=0
  ;  oplot,x2,y,linestyle=2,thick=2.0,color=220

  loadct,0,/silent

  ;CONTOURS
  if icont then $
    contour,pdf,x,y,/follow,/overplot,$
    levels=clevels,c_charsize=csize*0.7

  ;OVERPLOT AXIS LINES
  for i=0,1 do $
    plots,!x.crange[[i,i]],!y.crange
  plots,!x.crange,!y.crange[[0,0]]

  ;PRECIP RATE AND AREA FRAC
  plot_agg_size=0
  if plot_agg_size then begin

    ;    axis,yaxis=1,ystyle=1,charsize=csize,yticklen=yticklen,$
    ;      yrange=[0,10],ytitle='Rainfall [ mm d!U-1!N ]',/save,/data
    ;    oplot,times_ave,precip,linestyle=2,thick=2.0
    ;
    ;    ;DIURNAL AMPLITUDE
    ;    if diurn_pcp_amp then $
    ;      oplot,time_days,d_amp,color=0,linestyle=3,thick=thick_2

    xloc=103
    ;    axis,xloc,yaxis=1,ystyle=1,charsize=csize,yticklen=yticklen,$
    ;      yrange=[0,100],ytitle='Ascent Area Fraction [ % ]',/save,/data
    axis,yaxis=1,ystyle=1,charsize=csize,yticklen=yticklen,$
      yrange=[0,100],ytitle='Ascent Area Fraction [ % ]',/save,/data
    oplot,times_ave,frac_up,linestyle=2,thick=2.0

  endif

  ;COLOR BAR
  if icbar then begin
    loadct,col_table,/silent
    ;    cpos= [ position[2]+0.03 ,$
    ;    cpos= [ 0.91 ,$
    ;      position[1] ,$
    ;;      position[2]+0.04 ,$
    ;      0.922 ,$
    ;      position[3] ]
    ;    if plot_agg_size then $
    ;      cpos[[0,2]]+=0.02
    cpos= [ position[0]+0.03 , 0.54 , position[0]+0.03+0.01 , position[3]-0.04 ]
    colorbar2, colors=colors, range=[min(levels),max(levels)], divisions=5, $
      charsize=csize*0.7, position=cpos, /right, /vertical, title=cbar_tag, $
      annotatecolor='black',format='(f5.1)'
    loadct,0,/silent
  endif

  ;LEGEND
  if ileg then begin
    ;  if ileg then begin
    csize_fac=0.7
    margin=0.1
    pspacing=2.0 ; length of lines
    spacing=0.8 ; between lines
    if plot_agg_size then begin
      ;leg_str=['Mean','Standard Dev.','Area Frac.']
      leg_str=['Mean','Precip','Area Frac.']
      leg_str=['Mean','Area Frac.']
      leg_style=[0,2]
      leg_thick=[2,2]
      leg_color=[0,0]
    endif else begin
      leg_str=['Mean','Standard Dev.']
      leg_style=[0,2]
      leg_thick=[2,2]
      leg_color=[0,0]
    endelse
    legend2,leg_str,linestyle=leg_style,thick=leg_thick,COLORS=leg_color,$
      charsize=csize*csize_fac,/top_legend,/right_legend,/normal,$
      margin=margin,pspacing=pspacing,spacing=spacing,box=0
  endif

  device,/close
  convert_png,figname,res=200,/remove_eps


print,'Done plotting!'


end