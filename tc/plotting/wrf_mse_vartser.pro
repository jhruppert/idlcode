; 
; Time series of TC output

; James Ruppert
; 3/26/20
; 
pro wrf_mse_vartser, tcname, hvar, halfddt, hpsefp, hplwp, hpswp, nc, ntime, figdir, cases, idaily


for ivar=0,5 do begin
;for ivar=0,0 do begin
;ivar=2 ; 0-var budget, 1-ratio, 2-MSEVAR, 3-SEF, 4-LW, 5-SW


;----CREATE PLOT--------------------

if ivar le 1 then begin
  if ivar eq 0 then begin
    figname=figdir+'tser_varbud'
    ytitle1='Tendency [ day!U-1!N ]'
    yrange=[-0.75,0.75]
;yrange=[-0.2,0.5]
    iv2=0
    ytitle2='Variance [ J!U2!N m!U-4!N ]'
    yrange2=[1e13,1e15]
    ylog=0
    ylog2=1
    yminor=2
  endif else begin
    figname=figdir+'tser_varbud_ratio'
    ytitle1='WISHE [ J!U2!N m!U-4!N s!U-1!N ]'
    yrange=[-0.5,1.5]
    iv2=1
    ytitle2='Ratio LW:SEF'
    yrange2=[-0.5,1.5]
    ylog=0
    ylog2=0
    yminor=2
    var2=abs(reform(hplwp[0,*]/hpsefp[0,*]))
  endelse
endif else begin
  iv2=0
  ylog=0
  yminor=2;''
  if ivar eq 2 then begin
    figname=figdir+'tser_msevar'
    scale=1e-14
    var1=hvar*scale
    ytitle1='Variance [ 10!U14!N J!U2!N m!U-4!N ]'
;    yrange=[1e13,1e15]
    yrange=[0,10];4]
    ylog=0;1
  endif else if ivar eq 3 then begin
    figname=figdir+'tser_hpsefp'
    scale=1e-8
    var1=hpsefp*scale
    ytitle1='Tendency [ 10!U8!N J!U2!N m!U-4!N s!U-1!N ]'
    yrange=[-5,25]
  endif else if ivar eq 4 then begin
    figname=figdir+'tser_hplwp'
    scale=1e-8
    var1=hplwp*scale
    ytitle1='Tendency [ 10!U8!N J!U2!N m!U-4!N s!U-1!N ]'
    yrange=[-5,25]
  endif else if ivar eq 5 then begin
    figname=figdir+'tser_hpswp'
    scale=1e-8
    var1=hpswp*scale
    ytitle1='Tendency [ 10!U8!N J!U2!N m!U-4!N s!U-1!N ]'
    yrange=[-5,25]
  endif
endelse

  ;PLOT SPECS
    csize=0.8
    position=[0.14,0.18,0.89,0.89]
    xsize=4.2 & ysize=2
    title=''
    icplot=indgen(nc)

  ncplot=n_elements(icplot)

  ;AXES
    x=findgen(ntime)/24
    y=indgen(10)
    if ~keyword_set(yrange) then $
      yrange=[min(y),max(y)]

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  if tcname eq 'maria' then begin
    xticks=5
    xtickv=indgen(xticks+1)+15
    xrange=[14.5,20.5]
;xticks=4
;xtickv=indgen(xticks+1)+16
;xrange=[15.5,20.5]
    x+=14.5
    ;xtickname=strtrim(indgen(xticks+1)+15,2)
    xtitle='Day in Sept [ UTC ]'
  endif else if tcname eq 'haiyan' then begin
    xticks=7
    xtickv=indgen(xticks+1)+1
    x+=1
    xrange=[1,8]
    ;xtickname=strtrim(indgen(xticks+1)+1,2)
    xtitle='Day in Nov [ UTC ]'
;xticks=6
;xtickv=indgen(xticks+1)+2
;xrange=[2,8]
  endif

  if idaily then $
    xplot=findgen(ntime/24)+x[0]+0.5 $
  else $
    xplot=x

  if ~keyword_set(xrange) then $
    xrange=[min(x),max(x)]

  plot,x,y,/nodata,position=position,ylog=ylog,$
    xstyle=9,ystyle=9,$
    xrange=xrange,yrange=yrange,yminor=yminor,$
    xticks=xticks,xtickv=xtickv,xticklen=0.034,xminor=2,$;xtickname=xtickname,$
    xtitle=xtitle,ytitle=ytitle1,$
    charsize=csize,$
    title=title

;DAY LINES
;  for id=1,nd-1 do $
;    plots,replicate(24*id,2),!y.crange,linestyle=0,thick=1,color=col

if ivar le 1 then begin
;MSE VARIANCE BUDGET

  loadct,6,/silent

  colmse=0
  coladv=0
  colsh=125 ; green
  collh=125
  colsw=80 ; red
  collw=80

  stmse=0
  stadv=1
  stsh=1
  stlh=0
  stsw=1
  stlw=0

  lthick=2.5

  if ivar eq 0 then begin
  
    plots,!x.crange,[0,0],linestyle=0,thick=1.5,color=0

    adv = halfddt - (hpsefp + hplwp + hpswp)

    oplot,xplot,reform(halfddt[0,*]/hvar[0,*]*3600*24),linestyle=0,thick=lthick,color=0
    oplot,xplot,reform(adv[0,*]/hvar[0,*]*3600*24),linestyle=1,thick=lthick,color=0
  
    oplot,xplot,reform(hpsefp[0,*]/hvar[0,*]*3600*24),linestyle=stlh,thick=lthick,color=collh
    oplot,xplot,reform(hplwp[0,*]/hvar[0,*]*3600*24),linestyle=stlw,thick=lthick,color=collw
    oplot,xplot,reform(hpswp[0,*]/hvar[0,*]*3600*24),linestyle=stsw,thick=lthick,color=colsw
  
  endif else begin
  
    plots,!x.crange,[0,0],linestyle=0,thick=1.5,color=0
  
    oplot,xplot,reform(hpsefp[0,*]/hvar[0,*]*3600*24),linestyle=stlh,thick=lthick,color=collh
;    oplot,xplot,abs(reform(hplwp[0,*]/hpsefp[0,*])),linestyle=stlw,thick=lthick,color=collw
  
  endelse

endif else begin

;VARIABLES

  cols=(findgen(nc)+1)/nc*220
  cols[0]=[0]
  lstyle=replicate(0,nc);indgen(nc)
;  lstyle[[4,5]]=[1,2]
  lthick=reverse((findgen(nc)+3)*4.5/nc)
lthick*=0.6
;  lthick_ctl=2.4
;  lthick[0]=lthick_ctl

  loadct,3,/silent

  plots,!x.crange,[0,0],linestyle=0,thick=1.5,color=0

  ;VAR1
  for ic=0,nc-1 do $
;  for ic=0,3 do $
    oplot,xplot,var1[ic,*],linestyle=lstyle[ic],thick=lthick[ic],color=cols[ic]

  ;ADD HASHES FOR START TIME
;  if ivar ne 1 then $
  for ic=1,nc-1 do begin
    iv=reform(var1[ic,*])
    x0=min(where(finite(iv)))
    iv0=reform(iv[x0])
;    wd=3
;    plots,replicate(x0,2),[iv0-wd,iv0+wd],linestyle=lstyle[ic],thick=lthick[ic],color=cols[ic]
    wd=(yrange[1]-yrange[0])*0.034
;    datc=convert_coord([0.2,position[1]+wd],/normal,/to_data)
    plots,replicate(xplot[x0],2),[iv0-wd,iv0+wd],linestyle=lstyle[ic],thick=lthick[ic],color=cols[ic]
  endfor

endelse


if iv2 then begin

  ;VAR2
  axis,yaxis=1,ystyle=9,ytitle=ytitle2,charsize=csize,yrange=yrange2,yminor=yminor2,ylog=ylog2,/save

if ivar le 1 then begin
  ;MSE VARIANCE
  oplot,xplot,var2,linestyle=1,thick=lthick,color=0
endif else begin
  loadct,3,/silent
  for ic=0,nc-1 do $
    oplot,xplot,var2[ic,*],linestyle=lstyle[ic],thick=lthick[ic],color=cols[ic]
endelse

endif; else begin


;endelse

  loadct,0,/silent

  ;LEGEND
  ileg=1
  if ileg and ivar gt 1 then begin
  loadct,3,/silent
;  if ileg then begin
    csize_fac=0.7
    margin=0.1
    pspacing=2.0 ; length of lines
    spacing=0.8 ; between lines
    leg_str=strupcase(cases[icplot])
    leg_style=lstyle
    leg_thick=lthick
    leg_color=cols;replicate(0,ncplot);cols
    legend2,leg_str,linestyle=leg_style,thick=leg_thick,COLORS=leg_color,$
      charsize=csize*csize_fac,/top_legend,/left_legend,/normal,$
      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.16,0.75]
;      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.56,0.35]
;  loadct,0,/silent
  endif

  device,/close

  convert_png,figname,res=200;,/remove_eps


endfor

print,'DONE!!'
end
