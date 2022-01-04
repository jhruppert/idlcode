;
; Spectral plot of the butterworth filter employed in filtering.
;
pro plot_butter, freq, butter_band, figname, s0=s0, s1=s1

  set_plot,'ps'

  epsname=figname+'.eps'
  !p.font=0

  device,filename=epsname,/encapsulated,/color,bits=8,xsize=3.5,ysize=2,/inches,/helvetica
  loadct,0,/silent

  period = 1./freq

  xrange=[1./365,1./2]

    period2 = FIX(ALOG(period)/ALOG(2))   ; integer powers of 2 in period
    xtickv = 2.^(period2[UNIQ(period2)])  ; unique powers of 2
    xtickv = xtickv[where(xtickv le 1./xrange[0])]
    xticks=n_elements(xtickv)-1
    xtickname=strtrim(fix(xtickv),2)
    xtickv=1./xtickv

  plot,freq,butter_band,/xlog,position=[0.07,0.14,0.95,0.95],$
    xstyle=9,ystyle=9,xminor=1,$
    xrange=xrange,yrange=[0,1.2],$
    xtitle='Period',xticks=xticks,xtickv=xtickv,xtickname=xtickname,$
    charsize=0.7

  if keyword_set(s0) then plots,replicate(1./s0,2),!y.crange,linestyle=1,/data
  if keyword_set(s1) then plots,replicate(1./s1,2),!y.crange,linestyle=1,/data

  device,/close

  convert_png,figname,res=200,/remove_eps

end
;============================================================================
; 
; Program to run filtering procedures for a given input variable for BSISO
; and Biweekly bands (monsoon study).
;
; Assumes either a 1d [t] or 3d [x,y,t] input variable.
;
; TYPE = 'rmean' or 'fft' for using running mean or FFT filter approach
;
; Will also generate spectral plots of Butterworth filters applied in run
; directory.
;
; Returns:
;   var_bw, var_intra: filtered variables
;
; James Ruppert
; 7/22/21
; 
pro filter_monsoon, var_sav, type, var_bw=var_bw, var_intra=var_intra

;if type eq 'fft' then $
;  print,'NOTE: Butterworth filter is valid for a daily time'+$
;    'series of length 2000-2020!'

;Don't overwrite input variable.
var=var_sav

dims=size(var)
ndims=dims[0]

i3d=0
if ndims eq 3 then begin
  i3d=1
  nx=dims[1]
  ny=dims[2]
  nt=dims[3]
endif else nt=dims[1]


;----FILTER PT 1--------------------

;BAND-PASS FILTER USING SIMPLE RUNNING MEAN
;THIS WILL DETREND AND REMOVE MOST VARIABILITY OUTSIDE THESE BOUNDS

;  max_ts=121
;  min_ts=5
;
;  ;SUBTRACT RUNNING MEAN
;    ism=max_ts
;    if i3d then ism=[0,0,max_ts]
;    var -= smooth(var,ism,/edge_truncate)
;
;  ;KEEP BAND-PASS
;    min_ts = 5 ; n-days of running average
;    ism=min_ts
;    if i3d then ism=[0,0,min_ts]
;    var = smooth(temporary(var),ism,/edge_truncate)
;
;;----REMOVE MEAN--------------------

if i3d then begin
  varm = mean(var,dimension=3,/nan,/double)
  varm = rebin(varm,nx,ny,nt)
  var -= varm
endif else $
  var -= mean(var,/nan,/double)

;----FILTER PT 2--------------------

  ;See e.g. Lee et al. 2013, Clim Dym; Fujinami et al. 2014, Clim Dyn
  ;  for time scales

for iband=1,2 do begin

  ;SETTINGS

  if iband eq 1 then begin

  ;BIWEEKLY: 7-25 days
    tag='biweekly'

    ;RUNNING AVERAGE
      s0=7
      s1=25
;s0=7
;s1=60
s0=10
s1=20
;    ;BUTTERWORTH FILTER
;;      t_sel=14 ; selected midpoint in period
;;      cutoff=50 ; scales the width of the filter
;;      nwts=7 ; scales the steepness at cutoff
;      t_sel=11 ; selected midpoint in period
;      cutoff=390 ; scales the width of the filter
;      nwts=9 ; scales the steepness at cutoff

    ;WILL SIMPLY ZERO OUT OUTSIDE OF BOUNDS
      iper=[10,20]

  endif else if iband eq 2 then begin

  ;BSISO: 30-60 days
    tag='intra'

    ;RUNNING AVERAGE
      s0=30
      s1=60
;    ;FOR FFT BUTTERWORTH
;;      t_sel=31 ; selected midpoint in period
;;      cutoff=35 ; scales the width of the filter
;;      nwts=7 ; scales the steepness at cutoff
;      t_sel=39 ; selected midpoint in period
;      cutoff=70 ; scales the width of the filter
;      nwts=9 ; scales the steepness at cutoff

    ;WILL SIMPLY ZERO OUT OUTSIDE OF BOUNDS
      iper=[30,60]

  endif


  ;RUN FILTERING PROCEDURES

  if type eq 'rmean' then begin

    ;RUNNING AVERAGE

    ism=s0
    if i3d then ism=[0,0,s0]
    var_filt = smooth(var,ism,/edge_truncate)

    ism=s1
    if i3d then ism=[0,0,s1]
    var_filt -= smooth(var_filt,ism,/edge_truncate)

  endif else if type eq 'fft' then begin

    freq = calc_freq(nt,1) ; N, DeltaT
    per  = 1./freq
    per  = abs(per)

    ;GENERATE BUTTERWORTH FILTER
;      t2=abs(t_sel-per)
;      loc=where(t2 eq min(t2))
;
;      butter_temp=butterworth(nt,cutoff=cutoff,order=nwts,/origin)
;      butter_band=shift(butter_temp, -1*((nt-1)/2-loc) )

    ;SIMPLY ZERO OUT FREQUENCIES OUTSIDE DESIRED RANGE
      ifreq=where(per ge iper[0] and per le iper[1])

      band=fltarr(nt)
      band[ifreq]=1.

    figname='butterworth_'+tag
;    plot_butter,freq,butter_band,figname,s0=s0,s1=s1

    ;PRINT PERIODS RETAINED BY FILTER
;      print,per[where(butter_band ge 0.25)]
;      print

    if i3d then band=transpose(rebin(band,nt,nx,ny),[1,2,0])

    if iband eq 1 then begin
      if i3d then temp_fft1=fft(var,-1,dimension=3) $; DO NOT NEED TO RECALCULATE THIS A 2ND TIME
      else        temp_fft1=fft(var,-1)
    endif

    if i3d then var_filt = real_part ( fft(temp_fft1*band,1,dimension=1) ) $
    else        var_filt = real_part ( fft(temp_fft1*band,1) )

  endif else message,'Bad TYPE specified'

  ;REMOVE MEAN AGAIN
  if i3d then begin
    varm = mean(var_filt,dimension=3,/nan,/double)
    varm = rebin(varm,nx,ny,nt)
    var_filt -= varm
  endif else $
    var_filt -= mean(var_filt,/nan,/double)

  if iband eq 1 then var_bw = var_filt $
  else if iband eq 2 then var_intra = var_filt

endfor ; iband


end
