; 
; Plot time series of spatially averaged stddev of IMERG rainfall.
;
; Called by run_imerg_coast_dc_regress.pro
; 
; James Ruppert
; 1/6/2021
; 
pro plot_stddev_tser, dirs, figname, var

;CREATE FIGURE

  ;PLOT SPECS
    csize=0.75
    xsize=8.2 & ysize=2.8

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  nt=n_elements(var)
  x=findgen(nt*2)/2
  y=[var,var]

  position=[0.07,0.13,0.9,0.93]

  plot,x,y,/nodata,position=position,$
    xstyle=9,ystyle=9,$
    xrange=[0,48],yrange=[0,80],$
    xtitle='Hour',ytitle=' ',$
    charsize=csize;,$

  oplot,x,y,linestyle=0,thick=2,color=0

  plots,[24,24],!y.crange

  device,/close

  convert_png,figname,res=200,/remove_eps

end
