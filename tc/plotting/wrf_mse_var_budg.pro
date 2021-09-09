; 
; Plot MSE variance budget from WRF TC output.
; 
; James Ruppert
; 3/9/2020
; 
pro wrf_mse_var_budg, figname, mseddt, radsw, radlw, lh, sh, adv

tcname='maria'

;CREATE FIGURE

  ;PLOT SPECS
    csize=0.8
    position=[0.14,0.18,0.89,0.89]
    xsize=4.2 & ysize=2
    ytitle='Tendency [ day!U-1!N ]'
    yrange=[-1,1]

  if tcname eq 'maria' then begin
    xticks=5
    xtickv=indgen(xticks+1)*24+12
    xtickname=strtrim(indgen(xticks+1)+15,2)
    xtitle='Date in Sept [ UTC ]'
  endif else if tcname eq 'haiyan' then begin
    xticks=6
    xtickv=indgen(xticks+1)*24
    xtickname=strtrim(indgen(xticks+1)+1,2)
    xtitle='Date in Nov [ UTC ]'
  endif

  ;AXES
    nt=n_elements(mseddt)
    x=indgen(nt)
    y=mseddt
    if ~keyword_set(yrange) then $
      yrange=[max(y),min(y)]

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

;  !p.multi=[2,0,2]

  loadct,0,/silent

  plot,x,y,/nodata,position=position,$
    xstyle=9,ystyle=9,$
    xrange=xrange,yrange=yrange,yminor=2,$
    xticks=xticks,xtickv=xtickv,xtickname=xtickname,xticklen=0.034,xminor=1,$
    xtitle=xtitle,ytitle=ytitle,$
    charsize=csize,$
    title=title

  oplot,x,mseddt*0,linestyle=0,thick=1.5,color=0

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

  oplot,x,mseddt,linestyle=stmse,thick=lthick,color=colmse
  oplot,x,adv,linestyle=stadv,thick=lthick,color=coladv
;  oplot,x,sh,linestyle=stsh,thick=lthick,color=colsh
  oplot,x,lh,linestyle=stlh,thick=lthick,color=collh
  oplot,x,radsw,linestyle=stsw,thick=lthick,color=colsw
  oplot,x,radlw,linestyle=stlw,thick=lthick,color=collw

;  loadct,0,/silent

  device,/close

  !p.multi=0

  convert_png,figname,/remove_eps,res=200

end
