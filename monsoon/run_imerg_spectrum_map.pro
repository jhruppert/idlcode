; 
; Calculate maps of spectral power from IMERG rainfall.
;
; Coastal distance downloaded from https://oceancolor.gsfc.nasa.gov/docs/distfromcoast/dist2coast.txt.bz2
;   See run_imerg_coast_index.pro for plot of coastal distance.
; 
; James Ruppert
; 5/27/21
; 
pro run_imerg_spectrum_map

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

;BOUNDS FOR READ-IN
  bounds=[62,5,103,28]

do_coast_subset=0 ; subset rainfall data by coast/land-sea mask?
iplot_mn_rain=0   ; plot time-mean rainfall?

do_seg=1 ; Segment into individual years?
 ; Overrides to 1 if JJAS

do_jjas=1 ; subset to JJAS?

;SELECT DATE RANGE
  yy_plot=[2000,2020]
  mm_plot=[6,12]
  dd_plot=[1,31] ; inclusive
yy_plot[0]=2001
mm_plot[0]=1
;   yy_plot=[2013,2017]
;   mm_plot=[1,12]
;   dd_plot=[1,31] ; inclusive

  ;DATE STRING
    form2='(i2.2)'
    form4='(i4)'
    dat_str=string(mm_plot[0],format=form2)+string(dd_plot[0],format=form2)+strmid(strtrim(yy_plot[0],2),2,2)+'-'+$
            string(mm_plot[1],format=form2)+string(dd_plot[1],format=form2)+strmid(strtrim(yy_plot[1],2),2,2)

;----OB DIRECTORIES--------------------

;    im_fil=dirs.wkdir+'imerg/imerg/data/jjas_2013-2017_daymean_3B-HHR.MS.MRG.3IMERG.V06B.nc4'
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
;  era_fil_p=dirs.wkdir+'era5/ERA5-20000101-20201231-pl_monthly.nc4'
;    ; u,v [ m/s ]
;    ; RH [ % ]
;  era_fil_s=dirs.wkdir+'era5/ERA5-20000101-20201231-sl_monthly.nc4'
;    ; SST [ K ]
;    ; PW [ mm ]

  ifigdir=dirs.figdir+'myanmar/imerg/power_spec/'

;----TIME ARRAY FOR DAILY AVERAGE TIME SERIES--------------------

  ;SELECTED TIME ARRAY
    time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
      final=julday(mm_plot[1],dd_plot[1],yy_plot[1],23,59,59),step_size=1,units='Days');30,units='Minutes')
    nt=n_elements(time)
    nd=nt
    nyr=yy_plot[1]-yy_plot[0]+1

  caldat,time,mm,dd,yy
;  caldat,time[0:364],mm,dd,yy
  jjas=where(mm ge 6 and mm le 9)

;=====BEGIN READING=========================================================

;----READ RAIN--------------------

  rain=read_nc_imerg(time,im_fil,lon=lonim,lat=latim,bounds=bounds) ; already in mm/d
  rain_sav=rain

  nx=n_elements(lonim)
  ny=n_elements(latim)

;----READ ERA5--------------------

;  ;500-HPA WIND
;  varstr='var131' ; u
;;  varstr='var132' ; v
;  psel_era=500 ; hPa
;  cvar=read_nc_era5(time,era_fil,varstr,plev=psel_era,lon=eralon,lat=eralat,bounds=bounds)
;
;  ;SST OR PW
;;  varstr='var34' ; SST
;;  ;varstr='var137' ; PW
;;  cvar=read_nc_era5(time,era_fil_s,varstr,lon=eralon,lat=eralat,bounds=bounds)
;
;  nx2=n_elements(eralon)
;  ny2=n_elements(eralat)

;----PLOT TIME-MEAN RAINFALL--------------------

if iplot_mn_rain then begin

  var_str='RAINNC'
  setmax=24 & setmin='0.'
  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figspecs=create_struct(figspecs,'figname',ifigdir+'imerg_mean_'+dat_str)
  figspecs.cbar_format='(i2)'
  
  rain_mn=mean(rain_sav[*,*,jjas],dimension=3,/nan,/double)
  wrf_myanmar_map_plot, dirs, rain_mn, lonim, latim, figspecs, bounds=bounds
exit
endif

;----RAINFALL POWER SPECTRA--------------------

  ;STEPS FOR GENERATING POWER SPECTRUM
  ;OVERVIEW: GENERATE POWER SPECTRA FOR MULTIPLE OVERLAPPING SEGMENTS

  ;PLOT HISTOGRAM
    ;plot_imerg_hist, dirs.figdir+'myanmar/imerg/', rain;reform(segp,nd)

  ;STEP 1 - BAND-PASS FILTER
    ;REMOVE PERIODS >= MIN_TS DAYS (THIS ALSO ACTS TO DETREND)
;      min_ts = 121 ; n-days of running average
;      rain_lp = smooth(rain,[0,0,min_ts],/edge_truncate)
;    ;SUBTRACT LP
;      rain -= rain_lp
;    ;KEEP BAND-PASS
;      min_ts = 7 ; n-days of running average
;      rain = smooth(temporary(rain),min_ts,/edge_truncate)

  if do_jjas then begin
  ;STEP 1.5 - SUBSET TO JJAS
    rain=rain[*,*,jjas]
    nd=n_elements(jjas)
    if keyword_set(cvar) then cvar=cvar[*,*,jjas]
    do_seg=1
  endif

if do_seg then begin

  ;STEP 2 - CREATE (OVERLAPPING?) SEGMENTS OF LENGTH = NSEG
    npseg=nd/nyr;100 ; (days)
    extra=(nd mod nyr)
    nband=npseg/2
    ngap=npseg;10
    seg=reform(rain[*,*,0:npseg-1],[1,nx,ny,npseg])
;    seg2=reform(cvar[*,*,0:npseg-1],[1,nx2,ny2,npseg])
    count=0
    while count+ngap+npseg le (nd-extra) do begin
      count+=ngap
      iseg=reform(rain[*,*,count:count+npseg-1],[1,nx,ny,npseg])
      seg=[seg,iseg]
;      iseg=reform(cvar[*,*,count:count+npseg-1],[1,nx2,ny2,npseg])
;      seg2=[seg2,iseg]
    endwhile
   nseg=(size(seg,/dimen))[0]
   print,''
   print,'NP = ',strtrim(npseg,2),', Yields ',strtrim(nseg,2),' segments'
   print,''

  ;STEP 3 - STANDARDIZE
;    segm = mean(seg,dimension=4,/nan,/double)
;    segs = stddev(seg,dimension=4,/nan,/double)
;    segm = rebin(temporary(segm),nseg,nx,ny,npseg)
;    segs = rebin(temporary(segs),nseg,nx,ny,npseg)
;    segp = (seg - segm)/segs ; Now each segment is r'
;    seg=0 & segm=0 & segs=0
segp=seg
;segp2=seg2

  ;STEP 4 - DETREND
;    itim  = findgen(npseg) ; simply an array of ascending values
;    itim  = (itim - mean(itim))/stddev(itim) ; standardize array
;    itim  = rebin(itim,npseg,nseg,nx,ny)
;    itim = transpose(temporary(itim),[1,2,3,0])
;    trend = total((itim * segp),4,/nan,/double)/npseg
;    trend_line = itim * rebin(trend,nseg,nx,ny,npseg)
;    segdt = segp - trend_line
;    segp=0 & trend=0 & trend_line=0 & itim=0
segdt=segp
;segdt2=segp2

  ;STEP 5 - SELECT MIN TIME SCALE AND RUN RUNNING AVERAGE
    ;min_ts = 3 ; n-days of running average and tapering
;    seg_sm = segdt;smooth(segdt,[min_ts,0],/edge_truncate)
;    segdt=0

  ;STEP 6 - TAPER EDGES BY MIN_TS
;    seg_tap=seg_dt; & seg_sm=0
;    ntap=min_ts
;    for it=0,ntap-1 do begin
;      scale=1.*it/float(ntap)
;      segdt[*,*,*,it] *= scale
;      segdt[*,*,*,npseg-1-it] *= scale
;    endfor

  ;STEP 7 - CALCULATE SPECTRUM, AVERAGE ACROSS SEGMENTS
    segdt=transpose(temporary(segdt),[3,0,1,2]) ; [ day, year, x, y ]
;    segdt2=transpose(temporary(segdt2),[3,0,1,2]) ; [ day, year, x, y ]

    spec = fft(segdt,-1,dimension=1,/double)
    spec = abs(spec)^2
;    spec2 = fft(segdt2,-1,dimension=1,/double)
;    spec2 = abs(spec2)^2

    freq=calc_freq(npseg,1) ; N, DeltaT
    t=1./freq
    t=abs(t)

    iper=[27,64] ; Aiming for 30-60
;    print,'Band1:',iper
    iband=where(t ge iper[0] and t le iper[1])
print,'Band1:',min(t[iband]),'-',max(t[iband])
;    rain_band1=total(spec[iband,*,*,*],1,/nan,/double)
    rain_band1=mean(spec[iband,*,*,*],dimension=1,/nan,/double)
    rain_band1=mean(rain_band1,dimension=1,/nan,/double) / npseg ; --> mm^2 / d
stats,rain_band1
exit
;    var2_band1=mean(spec2[iband,*,*,*],dimension=1,/nan,/double)
;    var2_band1=mean(var2_band1,dimension=1,/nan,/double)

    ;iper=[7,25]
    iper=[10,21]
;    print,'Band2:',iper
    iband=where(t ge iper[0] and t le iper[1])
print,'Band2:',min(t[iband]),'-',max(t[iband])
    print,''
    rain_band2=total(spec[iband,*,*,*],1,/nan,/double)
    rain_band2=mean(rain_band2,dimension=1,/nan,/double) / npseg ; --> mm^2 / d
;    var2_band2=mean(spec2[iband,*,*,*],dimension=1,/nan,/double)
;    var2_band2=mean(var2_band2,dimension=1,/nan,/double)

endif else begin

;message,"Haven't processed CVAR for this"

  ;CALCULATE SPECTRUM ACROSS ENTIRE RAINFALL TIME SERIES

    rain=transpose(temporary(rain),[2,0,1]) ; [ day, x, y ]
print,'running fft'
    spec = fft(rain,-1,dimension=1,/double)
    spec = abs(spec)^2

    freq=calc_freq(nd,1) ; N, DeltaT
    t=1./freq
    t=abs(t)

    iper=[30,60]
    print,'Band1:',iper
    iband=where(t ge iper[0] and t le iper[1])
print,min(t[iband]),'-',max(t[iband])
    rain_band1=mean(spec[iband,*,*],dimension=1,/nan,/double); / nd ; --> mm^2 / d

    iper=[7,25]
    print,'Band2:',iper
    iband=where(t ge iper[0] and t le iper[1])
print,min(t[iband]),'-',max(t[iband])
    print,''
    rain_band2=mean(spec[iband,*,*],dimension=1,/nan,/double); / nd ; --> mm^2 / d

endelse

stats,rain_band1
stats,rain_band2
;stats,var2_band1
;stats,var2_band2

scale=1.
if ~do_seg then scale=1e2

rain_band1*=scale
rain_band2*=scale

  ;STEP 8 - NORMALIZE
;    area = total(spec,3,/double);*(deriv(freq))[0]
;    area=rebin(temporary(area),nx,ny,npseg)
;    spec = 1e2*temporary(spec)/area

;----PLOTS--------------------

;  plot_imerg_pspec, figdir
;  figname=ifigdir+'imerg_pspec_map_band1'
;  plot_imerg_pspec_map, dirs, figname, lonim, latim, band1
;  figname=ifigdir+'imerg_pspec_map_band2'
;  plot_imerg_pspec_map, dirs, figname, lonim, latim, band2

  figname=' ';ifigdir+'imerg_pspec_map_ratio'
;  plot_imerg_pspec_map, dirs, figname, lonim, latim, band2/band1-1

  if do_seg then begin
    setmax=5
    cbar_format='(i1)'
    cbar_tag='[ (mm/d)!U2!N ]'
    figtag='_yrseg'
    if do_jjas then begin
      setmax=36
      cbar_format='(i2)'
      figtag='_yrseg_jjas'
    endif
  endif else begin
    setmax=24
    cbar_format='(i2)'
    cbar_tag='[ 10!U-2!N (mm/d)!U2!N ]'
    figtag='_fullseries'
  endelse

      var_str='RAINNC'
      setmin=0.
      myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
      figspecs=create_struct(figspecs,'figname',figname)
      figspecs.title=' '
      figspecs.cbar_format=cbar_format
      figspecs.cbar_tag=cbar_tag

      ;NEW COLOR TABLE
;        figspecs.col_table=70;71
;        ncols=n_elements(figspecs.colors)
;        colors=findgen(ncols)/(ncols-1)*255.;/2;+255/2
;        irev=1
;        if irev then colors=reverse(colors)
;        figspecs.colors=colors
;    ;    figspecs.ndivs-=1

;  if keyword_set(cvar) then begin
;;    cvar=create_struct('cvar',mean(cvar,dimension=3,/nan,/double),'x',eralon,'y',eralat)
;    cvar=create_struct('cvar',var2_band1,'x',eralon,'y',eralat)
;    cint=1
;    figspecs.clevs=(findgen(50)+1)*cint + 290
;  endif

  figname=ifigdir+'imerg_pspec_map_band1'+figtag
  figspecs.figname=figname
      wrf_myanmar_map_plot, dirs, rain_band1, lonim, latim, figspecs, cvar=cvar;, bounds=bounds

;  if keyword_set(cvar) then begin
;;    cvar=create_struct('cvar',mean(cvar,dimension=3,/nan,/double),'x',eralon,'y',eralat)
;    cvar=create_struct('cvar',var2_band2,'x',eralon,'y',eralat)
;    cint=1
;    figspecs.clevs=(findgen(50)+1)*cint + 290
;  endif

  figname=ifigdir+'imerg_pspec_map_band2'+figtag
  figspecs.figname=figname
      wrf_myanmar_map_plot, dirs, rain_band2, lonim, latim, figspecs, cvar=cvar;, bounds=bounds


      var_str='RAINNC'
      setmax=24 & setmin=0.
      myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
      figspecs=create_struct(figspecs,'figname',figname)
      figspecs.title=' '
  figspecs.cbar_format='(i2)'
  figspecs.cbar_tag='[ mm/d ]'
  figname=ifigdir+'mean'
  figspecs.figname=figname
;rainm=mean(rain_sav[*,*,jjas],dimension=3,/nan,/double)
rainm=mean(mean(segdt,dimension=1,/nan,/double),dimension=1,/nan,/double)
      wrf_myanmar_map_plot, dirs, rainm, lonim, latim, figspecs;, bounds=bounds

;stop

print,'DONE!!'
end
