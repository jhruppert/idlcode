; 
; Calculate the power spectrum (one dimensional) for coast-normal wind index.
;
; Coastal distance downloaded from https://oceancolor.gsfc.nasa.gov/docs/distfromcoast/dist2coast.txt.bz2
;   See run_imerg_coast_index.pro for plot of coastal distance.
; 
; James Ruppert
; 5/25/21
; 
pro run_coast_spectrum_series

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

;BOUNDS FOR READ-IN
;  bounds=[85.,13,96,23.5]
;  tag='box1'
;  bounds=[85.,10,92,15]
;  tag='box2'
  bounds=[78,6,104,29]
  tag='box3'

;SELECT DATE RANGE
  yy_plot=[2000,2020]
  mm_plot=[1,12]
  dd_plot=[1,31] ; inclusive

  ;DATE STRING
    form2='(i2.2)'
    form4='(i4)'
    dat_str=string(mm_plot[0],format=form2)+string(dd_plot[0],format=form2)+strmid(strtrim(yy_plot[0],2),2,2)+'-'+$
            string(mm_plot[1],format=form2)+string(dd_plot[1],format=form2)+strmid(strtrim(yy_plot[1],2),2,2)

;----OB DIRECTORIES--------------------

  era_dir=dirs.wkdir+'era5/'
  era_fil=era_dir+'ERA5-20000101-20201231-pl_dayavg.nc'

  ifigdir=dirs.figdir+'myanmar/imerg/power_spec/'

;----TIME ARRAY FOR DAILY AVERAGE TIME SERIES--------------------

  ;SELECTED TIME ARRAY
    time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
      final=julday(mm_plot[1],dd_plot[1],yy_plot[1],23,59,59),step_size=1,units='Days');30,units='Minutes')
    nt=n_elements(time)
    nd=nt
    nyr=yy_plot[1]-yy_plot[0]+1

  caldat,time,mm,dd,yy
  jjas=where(mm ge 6 and mm le 9)


;=====BEGIN READING=========================================================

;----READ ERA5--------------------

  ;USE COAST-NORMAL WIND AS INDEX INSTEAD

  ;850-HPA WIND
  psel_era=850 ; hPa
  u=read_nc_era5(time,era_fil,'var131',plev=psel_era,lon=eralon,lat=eralat,bounds=bounds)
  v=read_nc_era5(time,era_fil,'var132',plev=psel_era,lon=eralon,lat=eralat,bounds=bounds)

  icross=3
  coast_thresh=200 ; km
  coastnormal, icross, coast_thresh, u, v, eralon, eralat, uindex=uindex

  ;STANDARDIZE
  sttu=stddev(uindex,/nan,/double)
  uindex = (uindex-mean(uindex,/nan,/double))/sttu

;----POWER SPECTRA--------------------

  ;BOUNDS FOR AVERAGING RAINFALL

  ;STEPS FOR GENERATING POWER SPECTRUM
  ;OVERVIEW: GENERATE POWER SPECTRA FOR MULTIPLE OVERLAPPING SEGMENTS
  ;THEN AVERAGE THE POWER ACROSS THEM (WHEELER AND HENDON 2004, P.3)

  ;STEP 2 - BAND-PASS FILTER
    ;REMOVE PERIODS >= MIN_TS DAYS (THIS ALSO ACTS TO DETREND)
      min_ts = 121 ; n-days of running average
      lp = smooth(uindex,min_ts,/edge_truncate)
    ;SUBTRACT LP
      uindex -= lp
    ;KEEP BAND-PASS
      min_ts = 7 ; n-days of running average
      uindex = smooth(temporary(uindex),min_ts,/edge_truncate)

  ;SUBSET TO JJAS
    do_jjas=1
    if do_jjas then begin
      uindex=uindex[jjas]
      nd=n_elements(jjas)
    endif

  ;STEP 2 - CREATE (OVERLAPPING?) SEGMENTS OF LENGTH = NSEG
    npseg=nd/nyr ; (days)
    extra=(nd mod nyr)
    nband=npseg/2
    ngap=npseg;10
    seg=uindex[0:npseg-1]
    count=0
    while count+ngap+npseg le (nd-extra) do begin
      count+=ngap
      seg=[[seg],[uindex[count:count+npseg-1]]]
    endwhile
   nseg=(size(seg,/dimen))[1]
   print,''
   print,'NP = ',strtrim(npseg,2),', Yields ',strtrim(nseg,2),' segments'
   print,''

  ;STEP 7 - CALCULATE SPECTRUM
    spec_seg = fft(seg,-1,dimension=1,/double)
    spec_seg = abs(temporary(spec_seg))^2

  ;STEP 7.1 - 3-POINT RUNNING MEAN TO SMOOTH
;    spec_seg = smooth(spec_seg,[3,0],/edge_truncate)

  ;STEP 7.2 - REMOVE 1ST 3 HARMONICS TO REMOVE ANNUAL CYCLE
;    spec_seg[indgen(4),*]=0
;    spec_seg[indgen(2),*]=0

  ;STEP 7.3 - AVERAGE ACROSS SEGMENTS
    spec_seg = spec_seg[0:nband-1,*]
    spec = mean(spec_seg,dimension=2,/double,/nan)

  ;STEP 8 - NORMALIZE
    freq = calc_freq(npseg,1)
    freq = freq[0:nband-1]
    t    = 1./freq
    area = total(spec,/double);*(deriv(freq))[0]
    spec = 1e2*temporary(spec)/area

;----SPECIAL CASES--------------------

  ;SPECTRA OF INDIVIDUAL SEGMENTS
    area = total(spec_seg,1,/double);*(deriv(freq))[0]
    area = transpose(rebin(temporary(area),[nseg,nband]))
    spec_seg = 1e2*temporary(spec_seg)/area

;  ;SPECTRUM OF ENTIRE SEQUENTIAL DATASET
;    ;STEP 4 - DETREND
;      itim  = findgen(nd) ; simply an array of ascending values
;      itim  = (itim - mean(itim,/nan,/doubl))/stddev(itim) ; standardize array
;      trend = total((itim * uindex),/nan,/double)/nd
;      trend_line = itim * trend
;      u_dt = uindex; - trend_line
;    ;STEP 6 - TAPER EDGES BY MIN_TS
;;      u_tap=u_dt
;;      ntap=7;min_ts
;;      for it=0,ntap-1 do begin
;;        scale=1.*it/float(ntap)
;;        rain_tap[it] *= scale
;;        rain_tap[npseg-1-it] *= scale
;;      endfor
;    ;REPLACE WITH JOINED TAPERED SEGMENTS
;;      rain_tap=reform(seg_tap,nd)
;    ;TEST WITH SINUSOID
;;      period=100.
;;      rain_tap=sin( findgen(nd)/nd*2*!pi*(1.*nd/period) )
;    ;STEP 7 - CALCULATE SPECTRUM
;      nband_sing=nd/2
;      spec_sing = abs(fft(u_dt,-1,/double))^2
;      spec_sing = spec_sing[0:nband_sing-1]
;    ;STEP 8 - NORMALIZE
;      freq_sing = (findgen(nd)/nd)[0:nband_sing-1]
;      t    = 1./freq_sing
;      area = total(spec_sing,/double);*(deriv(freq_sing))[0]
;      spec_sing = 1e2*temporary(spec_sing)/area


;----RED NOISE AND SIGNIF--------------------

  ;LAG-1 AUTOCORRELATION
    a=fltarr(nseg)
    for iseg=0,nseg-1 do $
      a[iseg]=a_correlate(reform(seg[*,iseg]),1)
;      a[iseg]=0.5*(a_correlate(reform(seg[*,iseg]),1) + sqrt(a_correlate(reform(seg[*,iseg]),2)))
    efold=-1./alog(a) ; T = -Delta_t / ln(a)
    print,''
    print,'e-folding timescale (days):'
    print,'  ',efold
    efold=-1./alog(mean(a,/nan,/double))

  ;RED SPEC
    freq2p = (2*!pi)*(findgen(npseg)/npseg)[1:nband-1]
    redspec=2.*efold/(1. + efold^2 * freq2p^2 )
    ;NORMALIZE
      area = total(redspec,/double);*deriv(freq),/double)
      redspec = 1e2*temporary(redspec)/area
;print,redspec
  ;SIGSPEC
    nsamp = nseg;(2*efold);1.*nd/npseg
print,'DOF:',nsamp
    dof=nsamp
    siglvl=0.95
;    chisq=chisqr_cvf(1.-siglvl,dof)
;    f_stat=1.
    t_stat=t_cvf(1.-siglvl,dof)
    sigspec=redspec*t_stat;chisq

;----PLOTS--------------------

  plot_imerg_pspec, ifigdir
  figname=ifigdir+'uindex_pspec_series'
  plot_imerg_pspec_spectrum, figname, spec, freq, mspec=spec_seg, redspec=redspec, sigspec=sigspec
;  figname=ifigdir+'uindex_pspec_series_sing'
;  plot_imerg_pspec_spectrum, figname, spec_sing, freq;, mspec=spec_seg, redspec=redspec, sigspec=sigspec


print,'DONE!!'
end
