; ------------------------------------------------------------------------------------------
; 
; Plotting routines for wavelet analysis.
; 
; plot_imerg_wavelet.pro:
; 
; plot_imerg_wave_scaleavg
; plot_imerg_wave_global
; plot_imerg_wavelet
;
; James Ruppert
; 6/12/21
; 
; ------------------------------------------------------------------------------------------
; 
; Plot time series of scale-averaged wavelet spectrum, based on wavelet software of:
; "http://paos.colorado.edu/research/wavelets/"
; Written January 1998 by C. Torrence
;
; James Ruppert
; 6/13/21
; 
pro plot_imerg_wave_scaleavg, figname, time, iavg, iavg2, avg_tser, avg_tser2, avg_signif, avg_signif2

  nb=(size(avg_tser,/dim))[1]

  !P.CHARSIZE = 0.6;75
  !P.MULTI = 0
  !X.STYLE = 1
  !Y.STYLE = 1

  LOADCT,0,/silent

;--- Plot scale-average time series

  labels=label_date(date_format=['%y'])
  labeldays=timegen(start=time[0],final=max(time),step_size=1,units='years')

  xrange = [min(time),max(time)] ; plotting range
  yrange = [0,2.9]

  dy=0.11
  pos1 = [0.1,0.75,0.7,0.95]
  pos2 = [pos1[0],0.39,pos1[2],pos1[1]-dy]
  pos3 = [pos1[0],0.08,pos1[2],pos2[1]-dy]

for i=0,1 do begin

  if i eq 0 then begin
    ifigname=figname+'_avg1'
    i_iavg=iavg
    i_avg_tser=avg_tser
    i_avg_signif=avg_signif
  endif else begin
    ifigname=figname+'_avg2'
    i_iavg=iavg2
    i_avg_tser=avg_tser2
    i_avg_signif=avg_signif2
  endelse

  xex=0.02

  set_plot,'ps'
  epsname=ifigname+'.eps'
  !p.font=0
  xsize=6 & ysize=4
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

    title_scaleavg=strtrim(i_iavg[0],2)+'-'+strtrim(i_iavg[1],2)+' d Scale-average Time Series'

        PLOT,time,reform(i_avg_tser[*,0]),/NOERASE,POSITION=pos1, $
                xtickv=labeldays,xtickunits='time',xtickformat='label_date',$
                XRANGE=xrange,YRANGE=yrange, $
                XTITLE='Time (year)',YTITLE='Avg variance', $
                TITLE=title_scaleavg;'d) 8-16 d Scale-average Time Series'
        OPLOT,xrange,i_avg_signif[0]+[0,0],LINES=1
        oplot,time,reform(i_avg_tser[*,0]),linestyle=ib,thick=2,color=0
xyouts,pos1[0]+xex,pos1[3]-xex*2,string(mean(i_avg_tser[*,0]),format='(f6.3)'),/normal

        PLOT,time,reform(i_avg_tser[*,1]),/NOERASE,POSITION=pos2, $
                xtickv=labeldays,xtickunits='time',xtickformat='label_date',$
                XRANGE=xrange,YRANGE=yrange, $
                XTITLE='Time (year)',YTITLE='Avg variance', $
                TITLE=title_scaleavg;'d) 8-16 d Scale-average Time Series'
        OPLOT,xrange,i_avg_signif[1]+[0,0],LINES=1
        oplot,time,reform(i_avg_tser[*,1]),linestyle=ib,thick=2,color=0
xyouts,pos2[0]+xex,pos2[3]-xex*2,string(mean(i_avg_tser[*,1]),format='(f6.3)'),/normal

        PLOT,time,reform(i_avg_tser[*,2]),/NOERASE,POSITION=pos3, $
                xtickv=labeldays,xtickunits='time',xtickformat='label_date',$
                XRANGE=xrange,YRANGE=yrange, $
                XTITLE='Time (year)',YTITLE='Avg variance', $
                TITLE=title_scaleavg;'d) 8-16 d Scale-average Time Series'
        OPLOT,xrange,i_avg_signif[2]+[0,0],LINES=1
        oplot,time,reform(i_avg_tser[*,2]),linestyle=ib,thick=2,color=0
xyouts,pos3[0]+xex,pos3[3]-xex*2,string(mean(i_avg_tser[*,2]),format='(f6.3)'),/normal

        device,/close
        convert_png,ifigname,res=200,/remove_eps

endfor

end

; ------------------------------------------------------------------------------------------
; 
; Plot global averaged wavelet spectrum, based on wavelet software of:
; "http://paos.colorado.edu/research/wavelets/"
; Written January 1998 by C. Torrence
;
; James Ruppert
; 6/13/21
; 
pro plot_imerg_wave_global, figname, period, avg_var

  nb=(size(avg_var,/dim))[1]

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  xsize=2.2 & ysize=1.8
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  !P.CHARSIZE = 0.6;75
  !P.MULTI = 0
  !X.STYLE = 1
  !Y.STYLE = 1

  LOADCT,0,/silent

;--- Plot global spectra

  xrange = [365,2] ; plotting range
xrange[0]=360*2

  yrange = [0,25]

  position = [0.17,0.19,0.96,0.87]

        period2 = FIX(ALOG(period)/ALOG(2))   ; integer powers of 2 in period
        xtickv = 2.^(period2[UNIQ(period2)])  ; unique powers of 2
        xtickv = xtickv[where(xtickv le xrange[0])]
        xticks=n_elements(xtickv)-1

        PLOT,period,reform(avg_var[*,0]),XRANGE=xrange,yrange=yrange,/nodata,position=position, $
                xtickv=xtickv,xticks=xticks,$
                xstyle=9,ystyle=9,/xtype,$
                XTITLE='Period (days)',YTITLE='Variance', $
                TITLE='Global variance'

        for ib=0,nb-1 do $
          oplot,period,reform(avg_var[*,ib]),linestyle=ib,thick=2,color=0

  ;PLOT LINES FOR 7-25D BAND
    ib=7
    plots,[ib,ib],yrange,linestyle=0,thick=0.5,color=0
    ib=25
    plots,[ib,ib],yrange,linestyle=0,thick=0.5,color=0

        device,/close
        convert_png,figname,res=200,/remove_eps

end

; ------------------------------------------------------------------------------------------
; 
; Plot all fields in combination using framework of WAVETEST.PRO from wavelet software of:
; "http://paos.colorado.edu/research/wavelets/"
; Written January 1998 by C. Torrence
;
; James Ruppert
; 6/12/21
; 
pro plot_imerg_wavelet, figname, time, rain, iavg, wavelet

  period=wavelet.period
  recon=wavelet.recon_var
  power =wavelet.power
  signif=wavelet.signif
  coi   =wavelet.coi
  global_ws      =wavelet.global_ws
  global_signif  =wavelet.global_signif
  scale_avg      =wavelet.scale_avg
  scaleavg_signif=wavelet.scaleavg_signif

  labels=label_date(date_format=['%y'])
  labeldays=timegen(start=time[0],final=max(time),step_size=1,units='years')

  xrange = [min(time),max(time)] ; plotting range

  dy=0.11
  pos1 = [0.1,0.75,0.7,0.95]
  pos2 = [pos1[0],0.39,pos1[2],pos1[1]-dy]
  pos3 = [0.74,pos2[1],0.95,pos2[3]]
  pos4 = [pos1[0],0.08,pos1[2],pos2[1]-dy]

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  xsize=6 & ysize=4
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  !P.CHARSIZE = 0.6;75
  !P.MULTI = 0
  !X.STYLE = 1
  !Y.STYLE = 1

  LOADCT,0,/silent

;--- Plot time series
        PLOT,time,rain,XRANGE=xrange, $
                xtickv=labeldays,xtickunits='time',xtickformat='label_date',$
                xstyle=1,ystyle=0,$
                XTITLE='Time (year)',YTITLE='Rainfall (sigma)', $
                TITLE='a) IMERG Rainfall (BP filter: 5-90 d)', $
                POSITION=pos1
;        IF (N_ELEMENTS(recon) GT 1) THEN OPLOT,time,recon,COLOR=144,thick=0.5

;--- Contour plot wavelet power spectrum
        yrange = [365,2] ; period range
yrange[0]=2*365
        levels = [0,2.^(indgen(9)-2)];1,2,4,8,16,32]
;        colors = [64,128,208,254]
        nlev=n_elements(levels)
        colors=findgen(nlev)/(nlev-1)*254
;stats,power
        period2 = FIX(ALOG(period)/ALOG(2))   ; integer powers of 2 in period
        ytickv = 2.^(period2[UNIQ(period2)])  ; unique powers of 2
        ytickv = ytickv[where(ytickv le yrange[0])]

        CONTOUR,power,time,period,/NOERASE,POSITION=pos2,/nodata, $
                xtickv=labeldays,xtickunits='time',xtickformat='label_date',$
                xstyle=1,ystyle=0,$
                XRANGE=xrange,YRANGE=yrange,/YTYPE, $
                YTICKS=N_ELEMENTS(ytickv)-1,YTICKV=ytickv, $
                LEVELS=levels,C_COLORS=colors,/FILL, $
                XTITLE='Time (year)',YTITLE='Period (d)', $
                TITLE='b) Wavelet Power Spectrum'

  LOADCT,64,/silent
        CONTOUR,power,time,period,/overplot, $
                LEVELS=levels,C_COLORS=colors,/FILL
  LOADCT,0,/silent

        CONTOUR,power,time,period,/NOERASE,POSITION=pos2,/nodata, $
                xtickv=labeldays,xtickunits='time',xtickformat='label_date',$
                xstyle=1,ystyle=0,$
                XRANGE=xrange,YRANGE=yrange,/YTYPE, $
                YTICKS=N_ELEMENTS(ytickv)-1,YTICKV=ytickv
        CONTOUR,signif,time,period,/OVERPLOT,LEVEL=1,THICK=2, $
                C_LABEL=0,C_CHARSIZE=1;,C_ANNOT='90%'
; cone-of-influence, anything "below" is dubious
        x = [time[0],time,MAX(time)]
        y = [MAX(period),coi,MAX(period)]
        color = 4
        POLYFILL,x,y,ORIEN=+45,SPACING=0.5,COLOR=color,NOCLIP=0,THICK=1
        POLYFILL,x,y,ORIEN=-45,SPACING=0.5,COLOR=color,NOCLIP=0,THICK=1
        PLOTS,time,coi,COLOR=color,NOCLIP=0,THICK=1

;--- Plot global wavelet spectrum
        blank = REPLICATE(' ',29)
        PLOT,global_ws,period,/NOERASE,POSITION=pos3, $
                THICK=2,XSTYLE=9,YSTYLE=9, $
                YRANGE=yrange,/YTYPE,YTICKLEN=-0.02, $
                ;xrange=[0,max(global_ws)],XTICKS=2,XMINOR=2, $
                YTICKS=N_ELEMENTS(ytickv)-1,YTICKV=ytickv,YTICKNAME=blank, $
                xrange=[0,25],XTITLE='Power',TITLE='c) Global'
        OPLOT,global_signif,period,LINES=1
        XYOUTS,17,60,'95%'

;--- Plot scale-average time series
    title_scaleavg='d) '+strtrim(iavg[0],2)+'-'+strtrim(iavg[1],2)+' d Scale-average Time Series'

        PLOT,time,scale_avg,/NOERASE,POSITION=pos4, $
                xtickv=labeldays,xtickunits='time',xtickformat='label_date',$
                XRANGE=xrange,YRANGE=[0,MAX(scale_avg)*1.25],THICK=2, $
                XTITLE='Time (year)',YTITLE='Avg variance', $
                TITLE=title_scaleavg;'d) 8-16 d Scale-average Time Series'

        ;IF PRESENT, OVERLAY 2ND SCALE AVG
        if keyword_set(wavelet.scale_avg2) then $
          OPLOT,time,wavelet.scale_avg2,thick=3,color=120

        OPLOT,xrange,scaleavg_signif+[0,0],LINES=1

        device,/close
        convert_png,figname,res=200,/remove_eps

;print,'Done plotting!!'

end
