; 
; Plot rainfall histogram from IMERG rainfall.
;
; Called by run_imerg_spectrum.pro
; 
; James Ruppert
; 5/26/2021
; 
pro plot_imerg_hist, figdir, rain

;CREATE FIGURE

  figname=figdir+'imerg_histogram'

  ;PLOT SPECS
    csize=0.75
    xsize=7 & ysize=2.8
    position=[0.07,0.12,0.97,0.93]

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  ;GENERATE HISTOGRAM
  hist=histogram(rain,locations=x);nbins=100)
help,hist
stats,hist
  ;GAUSSIAN
  s=stddev(rain)

;  x=findgen(nt*2)/2
  y=hist

  plot,x,y,/nodata,position=position,$
    xstyle=9,ystyle=9,$
    xrange=[0,48],yrange=[0,80],$
    xtitle='Rain',ytitle='N',$
    charsize=csize;,$

  oplot,x,y,linestyle=0,thick=2,color=0

  device,/close

  convert_png,figname,res=200,/remove_eps

end
