; 
; Plot IMERG power spectrum as either 1-dimenional or map.
;
; Called by run_imerg__spectrum.pro
; 
; plot_imerg_pspec_map
; plot_imerg_pspec_series
; plot_imerg_pspec
;
; James Ruppert
; 5/26/2021
; 
; ------------------------
; 
; Map of spectral power
; 
pro plot_imerg_pspec_map, dirs, figname, lon, lat, band

  var_str='RAINNC'
  setmax=4;8 & setmin='0.'
  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figspecs=create_struct(figspecs,'figname',figname)
  figspecs.title=' '
  figspecs.cbar_format='(f3.1)';'(i1)'
  figspecs.cbar_tag='Power [ % ]'

  ;NEW COLOR TABLE
    figspecs.col_table=70;71
    ncols=n_elements(figspecs.colors)
    colors=findgen(ncols)/(ncols-1)*255.;/2;+255/2
    irev=1
    if irev then colors=reverse(colors)
    figspecs.colors=colors
;    figspecs.ndivs-=1

  wrf_myanmar_map_plot, dirs, band, lon, lat, figspecs, bounds=bounds

  print,'DONE PLOTTING!'

end
; ------------------------
; 
; 1-dimensional power spectrum
; 
pro plot_imerg_pspec_spectrum, figname, spec, freq, $
  mspec=mspec, redspec=redspec, sigspec=sigspec

  ;PLOT SPECS
    csize=0.65
;    xsize=7 & ysize=2.4
;    position=[0.09,0.15,0.97,0.93]
  xsize=2.2 & ysize=1.8
  position = [0.17,0.19,0.96,0.87]

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  x=freq
  nx=n_elements(x)
  y=spec
  xtitle='Frequency [ cycles day!U-1!N ]'
  ytitle='Power [ % ]'

;  xrange=[min(x),.142];0.5];max(x)]
  xrange=[1./365,1./2];6.8];1./2]
xrange[0]=1./128
xrange[1]=[1./7]
;xrange=[2*365,2]

  yrange=[0,20];max(y)]

;  ;REPLACE AXIS ANNOTATIONS WITH PERIOD
;  ;PLOT IN FREQUENCY AND REPLACE VALUES WITH PERIOD
;;    persel=[50,25,15,10,7,5,4,3,2]
;    persel=[round(2.*nx/3/10)*10,90,60,30,20,15,10,7];,5];,4,3,2]
;persel[0]=180;[300,150,100,50,25,15,10,7]
;;    xtickv=[0,1./float(persel)]
;;    tickstr=['Inf',strtrim(persel,2)]
;    xtickv=1./float(persel)
;    tickstr=strtrim(persel,2)
;    xticks=n_elements(xtickv)-1
    xtitle2='Period [ days ]'
;;  axis,0,0,xaxis=0,xstyle=1,charsize=csize,xrange=xrange,xtitle='Period [ days ]',$
;;    xticks=xticks,xtickv=xtickv,xtickname=tickstr

period=1./x
        period2 = FIX(ALOG(period)/ALOG(2))   ; integer powers of 2 in period
        xtickv = 2.^(period2[UNIQ(period2)])  ; unique powers of 2
;Overwrite
xtickv=[120,60,30,20,15,7,5]
        xtickv = xtickv[where(xtickv le 1./xrange[0] and xtickv ge 1./xrange[1])]
        xticks=n_elements(xtickv)-1
    xtickname=strtrim(fix(xtickv),2)
    xtickv=1./xtickv

  plot,x,y,/xlog,/nodata,position=position,$
    xstyle=9,ystyle=9,xminor=1,$
    xrange=xrange,yrange=yrange,$
    xtitle=xtitle2,ytitle=ytitle,yticklen=0.01,$
    xticks=xticks,xtickv=xtickv,xtickname=xtickname,$
    charsize=csize

  if keyword_set(mspec) then begin
;    nspec=(size(mspec,/dimen))[1]
;    for isp=0,nspec-1 do $
;      oplot,x,reform(mspec[*,isp]),linestyle=0,thick=1.0,color=200
    mins=min(mspec,dimension=2,/nan)
    maxs=max(mspec,dimension=2,/nan)
    polyfill,[x,reverse(x)],[mins,reverse(maxs)],color=220,/data,$
      clip=[xrange[0],yrange[0],xrange[1],yrange[1]],noclip=0
  endif

  thick=2.5

  oplot,x,y,linestyle=0,thick=thick,color=0

  if keyword_set(redspec) then $
;print,'Doing red'
;print,redspec
    oplot,x,redspec,linestyle=1,thick=thick,color=0

  if keyword_set(sigspec) then $
    oplot,x,sigspec,linestyle=2,thick=thick,color=0

;TESTS FOR BUTTERWORTH
;nf=2562
;
;t_sel=14
;co=50
;nw=7
;
;t_sel=31
;co=35
;nw=7
;
;f3=findgen(nf)/nf
;period3=1./f3
;t2=abs(t_sel-period3)
;loc=(where(t2 eq min(t2)))[0]
;
;butter_temp=butterworth(nf,cutoff=co,order=nw,/origin)
;butter_temp=shift(butter_temp, -1*((nf-1)/2-loc) )
;
;print,period3[where(butter_temp ge 0.25)]
;
;oplot,f3,butter_temp*38,thick=3.5

  device,/close

  convert_png,figname,res=200;,/remove_eps

  print,'DONE PLOTTING!'

end

; ------------------------
; 
; 1-dimensional filtered time series
; 
pro plot_imerg_pspec_series, figname, series1, $
  series2=series2, series3=series3

  ;PLOT SPECS
    csize=0.65
;    xsize=7 & ysize=2.4
;    position=[0.09,0.15,0.97,0.93]
  xsize=3.8 & ysize=1.8
  position = [0.17,0.19,0.96,0.87]

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  nx=(size(series1,/dimen))[0]
  x=indgen(nx)
  xtitle='Day'
  ytitle=' '

yrange=[min(series1),max(series1)]
;yrange=[0,35]
yrange=[-12,16]

  plot,x,series1,/nodata,position=position,$
    xstyle=9,ystyle=9,xminor=1,$
    xrange=xrange,yrange=yrange,$
    xtitle=xtitle,ytitle=ytitle,yticklen=0.01,$
    charsize=csize

  thick=2.5

  oplot,x,series1,linestyle=0,thick=thick,color=0
  if keyword_set(series2) then $
    oplot,x,series2,linestyle=1,thick=thick,color=0
  if keyword_set(series3) then $
    oplot,x,series3,linestyle=2,thick=thick,color=0

  device,/close

  convert_png,figname,res=200,/remove_eps

  print,'DONE PLOTTING!'

end


;
; SETUP FOR POWER SPECTRUM PLOTS
;
pro plot_imerg_pspec, figdir

end

