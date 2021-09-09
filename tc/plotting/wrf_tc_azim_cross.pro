; 
; Plot cross sections from DYNAMO WRF simulations.
; 
; James Ruppert
; 20.10.16
; 
pro wrf_tc_azim_cross, var, x, y, figspecs, $
  var2=var2, cvar=cvar, clevs=clevs, irad=irad, $
  var1d=var1d, var2_1d=var2_1d, $
  cv2=cv2, cl2=cl2, $
  pr_seq=pr_seq, rad_seq=rad_seq, var_seq=var_seq, cint_seq=cint_seq, $
  rmw=rmw


if keyword_set(var_seq) then doseq=1 else doseq=0

ishallow=0 ; zoomed into PBL?

;CHECK FOR SUB-PANEL VARS
  do1d=0
  nnan=where(~finite(var1d),ncomplement=ngood)
  if ngood gt 10 then do1d=1
  nnan=where(~finite(var2_1d),ncomplement=ngood2)
  if ngood2 gt 10 then do1dv2=1

;CREATE FIGURE

  ;PLOT SPECS
    csize=0.8
    xtitle='Radius [ km ]'
    ixtitle=xtitle

;  xsize=3.4 & ysize=3.0
  xsize=3.0 & ysize=3.2
  position=[0.21,0.15,0.82,0.97]

  if do1d then begin
  ;INCLUDE SUB-PANEL
    ixtitle=' '
;    xsize=3.4 & ysize=3.7
    dy_lower=0.14
    buff=0.04
;    position=[0.14,0.12,0.85,0.97]
;    position=[0.19,0.15,0.86,0.97]
    position2=position
    position[1]+=dy_lower+buff
    position2[3]=position[1]-buff
    xtickname=replicate(' ',10)
  endif; else begin
;  ;SINGLE PANEL
;    position=[0.19,0.15,0.86,0.97]
;  endelse

  ;AXES
    xrange=[min(x),max(x)]
    xrange=[xrange[0],600];500];800];400]
    nx=n_elements(x)
    fin=where(finite(x))
    if ~keyword_set(xrange) then $
      xrange=[x[fin[0]],x[max(fin)]]
    if ~keyword_set(yrange) then $
      yrange=[max(y),min(y)]

  if y[0] le 100 then ihght=1 else ihght=0

  if ihght then begin
    ylog=0
    ytitle='Height [ km ]'
    yrange=[0,max(y)]
  endif else begin
    ylog=1
    ytitle='Pressure [ hPa ]'
    ytickv=[1000,700,500,400,300,200,150];,100];,50]
    if ishallow then ytickv=[1000,925,850,800];,700]
  ;  ytickv=[1000,700,500,400,300,200];,100];,50]
    yrange=[max(y),min(ytickv)]
    yticks=n_elements(ytickv)-1
  endelse

  set_plot,'ps'
  if ishallow then figspecs.figname+='_shall'
  epsname=figspecs.figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  !p.multi=0
  if do1d then !p.multi=[2,0,2]

  loadct,0,/silent

  plot,x,y,/nodata,position=position,ylog=ylog,$
    xstyle=9,ystyle=9,xminor=2,yticks=yticks,ytickv=ytickv,$
    yticklen=-0.022,xticklen=-0.018,$
    xrange=xrange,yrange=yrange,$
    xtitle=ixtitle,$
    xtickname=xtickname,$
    ytickname=ytickname,yminor=0,$
    ytitle=ytitle,charsize=csize;,title=ititle

;for i=0,1 do $
;var=smooth(var,3,/edge_truncate,/nan)

  ;UNDERLAY EXTREMES FOR RAD
    if irad then begin
      ;5-15
;        nlevs_add=11
;        levs_add=findgen(nlevs_add)-15
      ;6-16
;        nlevs_add=11;6
;        levs_add=findgen(nlevs_add)-16
      ;8-22
;        nlevs_add=8;15
;        levs_add=findgen(nlevs_add)*2-22
      ;4-10
        nlevs_add=9;7
        levs_add=findgen(nlevs_add)-12;10
      ;----
      cols_add=findgen(nlevs_add)/(nlevs_add-1)*255
cols_add=findgen(nlevs_add)/(nlevs_add)*255*0.95
      ;cols_add=reverse(cols_add)
      ctadd=0
      loadct,ctadd,/silent
      for i=0,1 do $
        contour,var*figspecs.scale,x,y,/cell_fill,/overplot,$
          levels=levs_add,c_colors=cols_add
    endif

  loadct,figspecs.col_table,/silent

  ;FILL SHADING
    for i=0,1 do $
      contour,var*figspecs.scale,x,y,/cell_fill,/overplot,$
        levels=figspecs.levels,c_colors=figspecs.colors

  loadct,0,/silent

  ;CONTOURS
    if keyword_set(cvar) and ~doseq then begin
      if ~keyword_set(clevs) then clevs=figspecs.clevs
      contour,cvar,x,y,/follow,/overplot,$
        levels=clevs,c_colors=0,c_thick=1.5,c_charsize=0.9*csize,c_labels=replicate(1,200)
;        levels=clevs,c_colors=0,c_thick=2,c_charsize=1.0*csize,c_labels=replicate(1,200)
      contour,cvar,x,y,/follow,/overplot,c_linestyle=1,$
        levels=-1.*reverse(clevs),c_colors=0,c_thick=1.5,c_charsize=0.9*csize,c_labels=replicate(1,200)
    endif ;else if icont then $
    if keyword_set(cv2) then begin
;      loadct,7,/silent
;      if ~keyword_set(clevs2) then clevs2=figspecs.clevs
      contour,cv2,x,y,/follow,/overplot,$
        levels=cl2,c_colors=0,c_thick=1.7,c_linestyle=1,c_labels=replicate(0,200),c_charsize=0.9*csize
;      contour,cvar2,x,y,/follow,/overplot,c_linestyle=1,$
;        levels=-1.*reverse(clevs2),c_colors=200,c_thick=1.5,c_labels=replicate(1,200),c_charsize=0.9*csize
    endif
;      contour,cvar,x,y,/follow,/overplot,$
;        levels=clevels,c_charsize=csize*0.3;0.7

  ;OVERLAY RMW
    if keyword_set(rmw) then begin
      oplot,rmw,y,linestyle=0,thick=4,color=0
    endif

  ;OVERLAY SEQ VARIABLE
    if doseq then begin
      clevs_seq=findgen(60)*cint_seq+cint_seq
      contour,var_seq,rad_seq,pr_seq,/follow,/overplot,$
        levels=clevs_seq,c_colors=0,c_thick=1.5,c_charsize=0.9*csize,c_labels=replicate(1,200)
      contour,var_seq,rad_seq,pr_seq,/follow,/overplot,c_linestyle=2,$
        levels=-1.*reverse(clevs_seq),c_colors=0,c_thick=1.5,c_charsize=0.9*csize,c_labels=replicate(1,200)
    endif

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

  ;1D VAR

  if do1d then begin

    ixtitle=xtitle

    if strmatch(figspecs.ytitle1d,'OLR*') then begin
      ytickv=[100,200,300]
      yticks=2
    endif else begin
      yticks=0
      ytickv=0
    endelse

    plot,x,y,/nodata,position=position2,ylog=0,$
      xstyle=9,ystyle=9,xminor=4,$
      yticklen=0.015,xticklen=0.05,$
      xrange=xrange,yticks=yticks,ytickv=ytickv,yrange=figspecs.yrange1d,$
      xtitle=ixtitle,yminor=4,$
      ytitle=figspecs.ytitle1d,charsize=csize;,title=ititle

    oplot,x,var1d,linestyle=0,thick=2,color=0
    if do1dv2 then oplot,x,var2_1d,linestyle=1,thick=2,color=0

  endif

  ;COLOR BAR
  if figspecs.icbar then begin

    loadct,figspecs.col_table,/silent

    ybuff=0.17
    if do1d then ybuff=0.08

    cpos= [ position[2]+0.020 ,$
            position[1]+ybuff ,$
            position[2]+0.039 ,$
            position[3]-ybuff ]

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
      cpos[1]=cpos[1]-0.16
      cpos[3]=cpos[1]+0.14
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
