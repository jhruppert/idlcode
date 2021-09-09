; 
; Plot cross sections from DYNAMO WRF simulations.
; 
; James Ruppert
; 20.10.16
; 
pro wrf_tc_comp_prof, prof1, prof2, pres, figspecs, xtitle=xtitle

;CREATE FIGURE

  ;PLOT SPECS
    ylog=1
    csize=0.85
    position=[0.14,0.12,0.85,0.97]
    xsize=4.8 & ysize=3.7
    ytitle='Pressure [ hPa ]'

  ;AXES
    y=pres
    x=indgen(10)
;    xrange=[-10,3]
    xrange=[-8,8]
    if ~keyword_set(xrange) then $
      xrange=[x[fin[0]],x[max(fin)]]

  ytickv=[1000,700,500,400,300,200,100];,50]
  yrange=[max(y),min(ytickv)]
  yticks=n_elements(ytickv)-1

  set_plot,'ps'
  epsname=figspecs.figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  plot,x,y,/nodata,position=position,ylog=ylog,$
    xstyle=9,ystyle=9,yticks=yticks,ytickv=ytickv,$
;    yticklen=-0.011,xticklen=-0.022,$
    xrange=xrange,yrange=yrange,$
    xminor=4,$
    ;xtitle=xtitle,
;    xtickname=replicate(' ',10),$
    ytickname=ytickname,yminor=0,$
    ytitle=ytitle,charsize=csize;,title=ititle

  loadct,0,/silent

  lthick=2.

  oplot,prof1,y,linestyle=0,thick=lthick,color=0
  oplot,prof2,y,linestyle=2,thick=lthick,color=0

  device,/close

  convert_png,figspecs.figname;,/remove_eps,res=200

end
