; 
; Calculate the power spectrum (one dimensional) for IMERG rainfall.
;
; Coastal distance downloaded from https://oceancolor.gsfc.nasa.gov/docs/distfromcoast/dist2coast.txt.bz2
;   See run_imerg_coast_index.pro for plot of coastal distance.
; 
; James Ruppert
; 5/25/21
; 
pro run_imerg_spectrum_series

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

;BOUNDS FOR READ-IN
;  bounds=[85.,13,96,23.5]
;  tag='box1'
;  bounds=[85.,10,92,15]
;  tag='box2'
  bounds=[78,6,104,29]
  tag='box3'

do_coast_subset=0 ; subset rainfall data by coast/land-sea mask?
iplot_mn_rain=0   ; plot time-mean rainfall?

;SELECT DATE RANGE
  yy_plot=[2000,2020]
  mm_plot=[6,12]
  dd_plot=[1,31] ; inclusive
yy_plot[0]=2001
mm_plot[0]=1
;    yy_plot=[2013,2017]
;    mm_plot=[1,12]
;    dd_plot=[1,31] ; inclusive

  ;DATE STRING
    form2='(i2.2)'
    form4='(i4)'
    dat_str=string(mm_plot[0],format=form2)+string(dd_plot[0],format=form2)+strmid(strtrim(yy_plot[0],2),2,2)+'-'+$
            string(mm_plot[1],format=form2)+string(dd_plot[1],format=form2)+strmid(strtrim(yy_plot[1],2),2,2)

;----OB DIRECTORIES--------------------

  im_fil=dirs.wkdir+'imerg/imerg/data/jjas_2013-2017_daymean_3B-HHR.MS.MRG.3IMERG.V06B.nc4'
;    im_fil=dirs.scdir+'imerg/imerg_3B-HHR.MS.MRG.3IMERG.V06_daily_2000-2021.nc4'
    im_fil=dirs.wkdir+'imerg/imerg/data/imerg_3B-HHR.MS.MRG.3IMERG.V06_daily_2000-2021_latlonsubset.nc4'
;im_fil=dirs.wkdir+'imerg/imerg/data/imerg_3B-HHR.MS.MRG.3IMERG.V06_daily_2000-2021_latlonsubset_85-96-13-23.5mean.nc4'
    ;FILE IS INDEED IDENTICAL TO ORIG JJAS-2013-17 FILE OVER THAT DATE RANGE
;    im_fil=dirs.wkdir+'imerg/imerg/data/imerg_3B-HHR.MS.MRG.3IMERG.V06_daily_2000-2021_latlonsubset_jjas2013-2017.nc4' ; test match to other file

  npd_imerg=48
  era_dir=dirs.wkdir+'wrf_wps/wrfv4/output_myanmar/JJAS_2013-2017/era5/'
  coast_dist_fil=dirs.home+'idl/code/misc/dist2coast.nc'

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

;----READ RAIN--------------------

  rain_sav=read_nc_imerg(time,im_fil,lon=lonim,lat=latim,bounds=bounds) ; already in mm/d
  nx=n_elements(lonim)
  ny=n_elements(latim)

;----READ ERA5--------------------

;  ;500-HPA WIND
;  varstr='var131' ; u
;  varstr='var132' ; v
;  psel_era=500 ; hPa
;  u=read_nc_era5(time,era_fil,varstr,plev=psel_era,lon=eralon,lat=eralat,bounds=bounds)
;
;;OVERWRITE
;rain_sav=u
;
;  ;SST OR PW
;;  varstr='var34' ; SST
;;  ;varstr='var137' ; PW
;;  cvar=read_nc_era5(time,era_fil_s,varstr,lon=eralon,lat=eralat,bounds=bounds)
;
;lonim=eralon
;latim=eralat
;  nx=n_elements(lonim)
;  ny=n_elements(latim)

;----PLOT TIME-MEAN RAINFALL--------------------

if iplot_mn_rain then begin

  var_str='RAINNC'
  setmax=30 & setmin='0.'
  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figdir=dirs.figdir+'myanmar/imerg/'
  figspecs=create_struct(figspecs,'figname',figdir+'imerg_mean_'+dat_str)
  figspecs.cbar_format='(i2)'
  
  rain_mn=mean(rain_dm,dimension=3,/nan,/double)
  wrf_myanmar_map_plot, dirs, rain_mn, lonim, latim, figspecs, bounds=bounds

endif

;----READ LAND/SEA MASK--------------------

if do_coast_subset then begin

  spawn,'ls '+era_dir+'met_em.d01*',era_fil
  fil=era_fil[0]
;  topo=reform(read_nc_var(fil,'HGT_M'))*1e3 ; m --> km
  lsmask=reform(read_nc_var(fil,'LANDMASK'))
  dimsw=size(lsmask,/dimen)
  nxw=dimsw[0] & nyw=dimsw[1]

  lonw=read_nc_var(fil,'XLONG_M') ; deg
  lonw=reform(lonw[*,0])
  latw=read_nc_var(fil,'XLAT_M') ; deg
  latw=reform(latw[0,*])

  ;INTERPOLATE ONTO IMERG GRID
    lsmaskx=fltarr(nx,nyw)
    for iy=0,nyw-1 do lsmaskx[*,iy]=interpol(reform(lsmask[*,iy]),lonw,lonim)
    lsm2=fltarr(nx,ny)
    for ix=0,nx-1 do lsm2[ix,*]=interpol(reform(lsmaskx[ix,*]),latw,latim)
    lsmask=lsm2

;endif

;----READ COASTAL DISTANCE--------------------

;if do_coast_subset then begin

  dist_fil=read_nc_var(coast_dist_fil,'coast_dist')
  londist=reform(dist_fil[0,*])
  latdist=reform(dist_fil[1,*])
  coast_dist=reform(dist_fil[2,*])

  ;KEEP ONLY IMERG DOMAIN
  ixkeep=where((londist ge min(lonim)) and (londist le max(lonim)) and (latdist ge min(latim)) and (latdist le max(latim)))
  londist=londist[ixkeep]
  latdist=latdist[ixkeep]
  coast_dist=coast_dist[ixkeep]

  ;INTERPOLATE ONTO IMERG GRID
    triangulate,londist,latdist,tri
    coast_dist_im=griddata(londist,latdist,coast_dist,/linear,triangles=tri,xout=lonim,yout=latim,/grid)
    coast_dist=coast_dist_im

endif

;----RAINFALL POWER SPECTRA--------------------

  ;BOUNDS FOR AVERAGING RAINFALL
    bounds_all=[ $
;      [85.,13,96,23.5], $
      [85.,11,100,23.5], $  ; Large BoB
      [87.,15.5,96,23.5], $ ; Northern BoB
      [88.,23.5,93.5,27] ]  ; Bangladesh
    nb=(size(bounds_all,/dim))[1]

;  for ib=0,nb-1 do begin
  for ib=0,1 do begin

    ;SUBSET TO BOUNDING BOX
      ix=where((lonim ge bounds_all[0,ib]) and (lonim le bounds_all[2,ib]),countx)
      iy=where((latim ge bounds_all[1,ib]) and (latim le bounds_all[3,ib]),county)
      rain_dm=rain_sav[ix,*,*] & rain_dm=rain_dm[*,iy,*]

  ;STEPS FOR GENERATING POWER SPECTRUM
  ;OVERVIEW: GENERATE POWER SPECTRA FOR MULTIPLE OVERLAPPING SEGMENTS
  ;THEN AVERAGE THE POWER ACROSS THEM (WHEELER AND HENDON 2004, P.3)

  ;STEP 0 OPTION
  if do_coast_subset then begin
    coast_dist=coast_dist[ix,*] & coast_dist=coast_dist[*,iy]
    lsmask=lsmask[ix,*] & lsmask=lsmask[*,iy]
    ;ADD TIME INDICES TO THESE TO SIMPLIFY WHERE FUNCTIONS
      coast_dist2=rain_dm
      for id=0,nd-1 do coast_dist2[*,*,id]=coast_dist
      lsmask2=rain_dm
      for id=0,nd-1 do lsmask2[*,*,id]=lsmask
    ;USE COASTAL DISTANCE VARIABLE
      coast_thresh=100.
    isea=where((coast_dist2 ge coast_thresh) and (lsmask2 le 0.1),complement=inot)
    ;icoast=where((coast_dist2 le coast_thresh),complement=inot)
    rain_dm[inot]=!values.f_nan
    ;NEW RAINFALL MAP
      rain_mn=mean(rain_dm,dimension=3,/nan,/double)
      wrf_myanmar_map_plot, dirs, rain_mn, lonim[ix], latim[iy], figspecs, bounds=bounds
  endif

  ;STEP 1 - SPATIAL AVERAGE
;    rain=mean(mean(rain_dm,dimension=1,/nan,/double),dimension=1,/nan,/double)
    ;PLOT HISTOGRAM
      ;plot_imerg_hist, dirs.figdir+'myanmar/imerg/', rain;reform(segp,nd)

;REPLACE WITH COAST-NORMAL WIND
;restore,'uindex'
;help,uindex
;;  filter_monsoon, uindex, 'fft', var_bw=u_bw, var_intra=u_intra
;;uindex=u_intra
;;uindex=u_bw
;rain=uindex

  ;STEP 2 - BAND-PASS FILTER
    ;REMOVE PERIODS >= MIN_TS DAYS (THIS ALSO ACTS TO DETREND)
      min_ts = 121 ; n-days of running average
      rain_lp = smooth(rain,min_ts,/edge_truncate)
    ;SUBTRACT LP
      rain -= rain_lp
    ;KEEP BAND-PASS
      min_ts = 7 ; n-days of running average
      rain = smooth(temporary(rain),min_ts,/edge_truncate)

  ;SUBSET TO JJAS
    do_jjas=0
;do_jjas=1
    if do_jjas then begin
      rain=rain[jjas]
      rain_sav=rain_sav[*,*,jjas]
      nd=n_elements(jjas)
    endif

  ;STEP 2 - CREATE (OVERLAPPING?) SEGMENTS OF LENGTH = NSEG
    npseg=nd/nyr;100 ; (days)
    extra=(nd mod nyr)
    nband=npseg/2
    ngap=npseg;10
    seg=rain[0:npseg-1]
    count=0
    while count+ngap+npseg le (nd-extra) do begin
      count+=ngap
      seg=[[seg],[rain[count:count+npseg-1]]]
    endwhile
   nseg=(size(seg,/dimen))[1]
   print,''
   print,'NP = ',strtrim(npseg,2),', Yields ',strtrim(nseg,2),' segments'
   print,''

  ;STEP 3 - STANDARDIZE
segp=seg
help,seg
exit
;    segm = mean(seg,dimension=1,/nan,/double)
;    segs = stddev(seg,dimension=1,/nan,/double)
;    segm = rebin(temporary(segm),nseg,npseg)
;    segs = rebin(temporary(segs),nseg,npseg)
;    segm = transpose(temporary(segm))
;    segs = transpose(temporary(segs))
;    segp = (seg - segm);/segs ; Now each segment is r'
;    seg=0 & segm=0 & segs=0

  ;STEP 3.5 - REMOVE 1ST 3 HARMONICS TO REMOVE ANNUAL CYCLE
;    filter=fltarr(npseg,nseg) & filter[*]=1.
;    ind0=[indgen(4),npseg-indgen(3)]
;    filter[ind0,*]=0
;    spec=fft(segp,-1,dimension=1,/double)
;    spec*=filter
;    segp = real_part( fft(spec,1,dimension=1,/double) )
;    spec=0
;    filter=0

      ;PLOT FIRST 3 HARMONICS
;        seg1=reform(segp[*,4])
;        spec=fft(seg1,-1,dimension=1,/double)
;    
;        filt=fltarr(npseg)
;        filt[*]=0
;        filt[[1,npseg-1]]=1
;        filt[[npseg-1:npseg-3]]=[1,0,0]
;        h1=real_part(fft(spec*filt,1,/double))
;    
;        filt[1:3]=[0,1,0]
;        h2=real_part(fft(spec*filt,1,/double))
;    
;        filt[1:3]=[0,0,1]
;        filt[[npseg-1:npseg-3]]=[0,0,1]
;        h3=real_part(fft(spec*filt,1,/double))
;    
;        freq=calc_freq(npseg,1)
;        t=1./freq
;        t=abs(t)
;    
;        iper=[30,60]
;        iband=where(t ge iper[0] and t le iper[1])
;        filt[*]=0
;        filt[iband]=1
;        h3=real_part(fft(spec*filt,1,/double))
;    
;    ;    filt[*]=1
;    ;    filt[1:3]=0
;    ;    resid=real_part(fft(spec*filt,1,/double))
;      
;        plot_imerg_pspec, ifigdir
;        figname=ifigdir+'harmonics'
;        plot_imerg_pspec_series, figname, h1, $
;          series2=h2, series3=h3
;        figname=ifigdir+'harmonics_actual'
;        plot_imerg_pspec_series, figname, seg1, $
;          ;series2=h0 + h1+h2+h3
;          series2=h3
;    ;    figname=ifigdir+'harmonics_resid'
;    ;    plot_imerg_pspec_series, figname, seg1, $
;    ;      series2=resid

  ;STEP 4 - DETREND
;    itim  = findgen(npseg) ; simply an array of ascending values
;    itim  = (itim - mean(itim,/nan,/doubl))/stddev(itim) ; standardize array
;    itim  = rebin(itim,npseg,nseg)
;    trend = total((itim * segp),1,/nan,/double)/npseg
;    trend_line = itim * transpose(rebin(trend,nseg,npseg))
;    segdt = segp - trend_line
;    segp=0 & trend_line=0 & itim=0


  ;STEP 6 - TAPER EDGES BY MIN_TS
    seg_tap=segp;seg_sm; & seg_sm=0
    ntap=min_ts
;    for it=0,ntap-1 do begin
;      scale=1.*it/float(ntap)
;      seg_tap[it,*] *= scale
;      seg_tap[npseg-1-it,*] *= scale
;    endfor
    ;COSINE WINDOW
    edge=0.1 ; fraction of time series to taper at both sides
    ntap=round(edge*npseg)
    wind=sin(0.5*!pi*indgen(ntap)/(ntap-1))
    wind=[wind,fltarr(npseg-ntap*2)+1,reverse(wind)]
    wind=rebin(wind,npseg,nseg)
;    seg_tap*=wind

  ;STEP 7 - CALCULATE SPECTRUM
    spec_seg = fft(seg_tap,-1,dimension=1,/double)
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
;    spec = 1e2*temporary(spec)/area

;----SPECIAL CASES--------------------

  ;SPECTRA OF INDIVIDUAL SEGMENTS
    area = total(spec_seg,1,/double);*(deriv(freq))[0]
    area = transpose(rebin(temporary(area),[nseg,nband]))
;    spec_seg = 1e2*temporary(spec_seg)/area

;  ;SPECTRUM OF ENTIRE SEQUENTIAL DATASET
;    ;STEP 3 - STANDARDIZE
;      segm = mean(rain,/nan,/double)
;      segs = stddev(rain,/nan,/double)
;      rainp = (rain - segm)/segs
;      segm=0 & segs=0
;    ;STEP 4 - DETREND
;      itim  = findgen(nd) ; simply an array of ascending values
;      itim  = (itim - mean(itim,/nan,/doubl))/stddev(itim) ; standardize array
;      trend = total((itim * rainp),/nan,/double)/nd
;      trend_line = itim * trend
;      raindt = rainp - trend_line
;      rainp=0 & trend_line=0 & itim=0
;    ;STEP 6 - TAPER EDGES BY MIN_TS
;      rain_tap=raindt; & seg_sm=0
;      ntap=7;min_ts
;      for it=0,ntap-1 do begin
;        scale=1.*it/float(ntap)
;        rain_tap[it] *= scale
;        rain_tap[npseg-1-it] *= scale
;      endfor
;    ;REPLACE WITH JOINED TAPERED SEGMENTS
;;      rain_tap=reform(seg_tap,nd)
;    ;TEST WITH SINUSOID
;;      period=100.
;;      rain_tap=sin( findgen(nd)/nd*2*!pi*(1.*nd/period) )
;    ;STEP 7 - CALCULATE SPECTRUM
;      nband_sing=nd/2
;      spec_sing = abs(fft(rain_tap,-1,/double))^2
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
      a[iseg]=a_correlate(reform(seg_tap[*,iseg]),1)
;      a[iseg]=0.5*(a_correlate(reform(seg_tap[*,iseg]),1) + sqrt(a_correlate(reform(seg_tap[*,iseg]),2)))
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
  figname=ifigdir+'imerg_pspec_series_'+tag+'_b'+strtrim(ib+1,2)
;figname=ifigdir+'uindex_pspec_series_'+tag+'_b'+strtrim(ib+1,2)
  plot_imerg_pspec_spectrum, figname, spec, freq, mspec=spec_seg, redspec=redspec, sigspec=sigspec
;  plot_imerg_pspec_spextrum, figname, spec_sing, freq_sing;, mspec=spec_sav, sigspec=sigspec

endfor ; ib (bounding box)

print,'DONE!!'
end
