; 
; Calculate the either the FFT or wavelet spectra for IMERG rainfall.
;
; Now with coastal/offshore separation.
;
; Based on WAVETEST and utilizes the wavelet code package developed by:
; See "http://paos.colorado.edu/research/wavelets/"
; Written January 1998 by C. Torrence
;
; James Ruppert
; 7/2/21
; 
pro run_imerg_wavelet_coastal

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

;BOUNDS FOR READ-IN
  bounds=[78,6,104,29]
  tag='box3'

var_plot='rain';
var_plot='pw';'rain';

iera=1 ; read/plot ERA?
if var_plot ne 'rain' and ~iera then message,'Must include ERA if you want any other field'

;Standard FFT parsing by year or Wavelet?
option='fft';'wave';
;option='wave'

do_coast_subset=1 ; subset rainfall data by coast/land-sea mask?
iplot_mn_rain=0   ; plot time-mean rainfall?

;SELECT DATE RANGE
  yy_plot=[2000,2020]
;  yy_plot=[2015,2020]
  mm_plot=[6,12]
;mm_plot[0]=1
  dd_plot=[1,31] ; inclusive

  ;DATE STRING
    form2='(i2.2)'
    form4='(i4)'
    dat_str=string(mm_plot[0],format=form2)+string(dd_plot[0],format=form2)+strmid(strtrim(yy_plot[0],2),2,2)+'-'+$
            string(mm_plot[1],format=form2)+string(dd_plot[1],format=form2)+strmid(strtrim(yy_plot[1],2),2,2)

;----OB DIRECTORIES--------------------

;  im_fil=dirs.scdir+'imerg/imerg_3B-HHR.MS.MRG.3IMERG.V06_daily_2000-2021.nc4'
  im_fil=dirs.wkdir+'imerg/imerg/data/imerg_3B-HHR.MS.MRG.3IMERG.V06_daily_2000-2021_latlonsubset.nc4'
;im_fil=dirs.wkdir+'imerg/imerg/data/imerg_3B-HHR.MS.MRG.3IMERG.V06_daily_2000-2021_latlonsubset_85-96-13-23.5mean.nc4'
  time_imerg=timegen(start=julday(6,1,2000,0,0,0),final=julday(12,31,2020,23,59,59),step_size=1,units='Days')
  ;FILE IS INDEED IDENTICAL TO ORIG JJAS-2013-17 FILE OVER THAT DATE RANGE
;  im_fil=dirs.wkdir+'imerg/imerg/data/imerg_3B-HHR.MS.MRG.3IMERG.V06_daily_2000-2021_latlonsubset_jjas2013-2017.nc4' ; test match to other file
  npd_imerg=48
  era_dir=dirs.wkdir+'wrf_wps/wrfv4/output_myanmar/JJAS_2013-2017/era5/'
  coast_dist_fil=dirs.home+'idl/code/misc/dist2coast.nc'

  maindir=dirs.wkdir
  era_dir2=maindir+'era5/jjas_2013-2017/';ERAi-JJAS13-17-pl.nc4'
  era_fil=era_dir2+'ERAi-JJAS13-17-pl_dayavg.nc4'
  era_sfil=era_dir2+'ERAi-JJAS13-17-sl_dayavg.nc4'
;  npd_era=24
era_dir2=dirs.scdir+'era5/'
era_fil=era_dir2+'ERA5-20000101-20201231-pl_dayavg.nc'
era_pw=era_dir2+'ERA5-20000101-20201231-pw_dayavg.nc'
  time_era=timegen(start=julday(1,1,2000,0,0,0),final=julday(12,31,2020,23,59,59),step_size=1,units='Days')

  figdir=dirs.figdir+'myanmar/imerg/'

;----READ RAIN--------------------

  ;SELECTED TIME ARRAY
    time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
      final=julday(mm_plot[1],dd_plot[1],yy_plot[1],23,59,59),step_size=1,units='Days');30,units='Minutes')
    nt=n_elements(time)
    nd=nt;/npd_imerg
    nyr=yy_plot[1]-yy_plot[0]+1

  rain_sav=read_nc_imerg(time_imerg,im_fil,lon=lon,lat=lat,bounds=bounds) ; already in mm/d
  nx=n_elements(lon)
  ny=n_elements(lat)

  ;SAVE BASIC TIME-MEAN
  rain_mn_sav=mean(rain_sav,dimension=3,/nan,/double)

  ;JJAS INDICES
    caldat,time,mm,dd,yy
    jjas=where((mm ge 6) and (mm le 9))
    ntjjas=n_elements(jjas)

  if var_plot eq 'rain' then begin
    lon=lon
    lat=lat
  endif

;----READ ERA--------------------

if iera then begin

  ;DIMENSIONS
  eralon=read_nc_var(era_fil,'lon')
  eralat=read_nc_var(era_fil,'lat')
  nxera=n_elements(eralon)
  nyera=n_elements(eralat)

  ;PLEVS
;  p_era=reform(read_nc_var(era_fil,'plev'))*1d-2 ; Pa --> hPa
;  nzera=n_elements(p_era)

  ;SELECTED LEVEL
;  izlev_era=(where(p_era eq psel_era,count))[0]
;  if count eq 0 then stop

  ;DIFFERENT START TIME FROM IMERG
  tdiff=abs(time_era-time[0])
  it0=(where(tdiff eq min(tdiff)))[0]
  ntread=n_elements(time_era)-it0
  if ntread ne nt then stop

  ix=where((eralon ge bounds[0]) and (eralon le bounds[2]),nxera)
  iy=where((eralat ge bounds[1]) and (eralat le bounds[3]),nyera)
  eralon=eralon[ix]
  eralat=eralat[iy]

  ;WINDS
;    count=[nxera,nyera,1,ntread] & offset=[0,0,izlev_era,it0] ; x, y, p, t
;    u=reform(read_nc_var(era_fil,'var131',count=count,offset=offset))
;    v=reform(read_nc_var(era_fil,'var132',count=count,offset=offset))

  ;PW
  if var_plot eq 'pw' then begin;icalc_pw=1 else icalc_pw=0
    count=[nxera,nyera,ntread] & offset=[ix[0],iy[0],it0] ; x, y, t
    var=reform(read_nc_var(era_pw,'var137',count=count,offset=offset))
    rain_sav=var
    nx=nxera
    ny=nyera
    lon=eralon
    lat=eralat
  endif

endif ; iera

;----PLOT TIME-MEAN RAINFALL--------------------

if iplot_mn_rain then begin

  var_str='RAINNC'
  setmax=30 & setmin='0.'
setmax=10
  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figspecs=create_struct(figspecs,'figname',figdir+'imerg_mean_'+dat_str)
  figspecs.cbar_format='(i2)'

  wrf_myanmar_map_plot, dirs, rain_mn_sav, lon, lat, figspecs;, bounds=bounds

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
    for iy=0,nyw-1 do lsmaskx[*,iy]=interpol(reform(lsmask[*,iy]),lonw,lon)
    lsm2=fltarr(nx,ny)
    for ix=0,nx-1 do lsm2[ix,*]=interpol(reform(lsmaskx[ix,*]),latw,lat)
    lsmask=lsm2

;endif

;----READ COASTAL DISTANCE--------------------

;if do_coast_subset then begin

  dist_fil=read_nc_var(coast_dist_fil,'coast_dist')
  londist=reform(dist_fil[0,*])
  latdist=reform(dist_fil[1,*])
  coast_dist=reform(dist_fil[2,*])

  ;KEEP ONLY IMERG DOMAIN
  ixkeep=where((londist ge min(lon)) and (londist le max(lon)) and (latdist ge min(lat)) and (latdist le max(lat)))
  londist=londist[ixkeep]
  latdist=latdist[ixkeep]
  coast_dist=coast_dist[ixkeep]

  ;INTERPOLATE ONTO IMERG GRID
    triangulate,londist,latdist,tri
    coast_dist_im=griddata(londist,latdist,coast_dist,/linear,triangles=tri,xout=lon,yout=lat,/grid)
    coast_dist=coast_dist_im

  ;ADD TIME INDICES TO THESE TO SIMPLIFY WHERE FUNCTIONS
      coast_dist2=rebin(coast_dist,[nx,ny,nd])
      lsmask2=rebin(lsmask,[nx,ny,nd])

endif

;----PREP FOR WAVELET SPECTRA--------------------

  dt = 1 ; time interval (1 = 1 per time step)
;  pad = 1
  s0 = 2*dt  ; this says start at a scale of 2*dt
  dj = 0.25  ; this will do 4 sub-octaves per octave
  j1 = 9./dj ; this says do 9 powers-of-two with dj sub-octaves each
  mother = 'Morlet'
  ;mother = 'Dog'
  iavg=[7,25] ; averaging bounds (period; units of dt)
  iavg2=[25,80] ; averaging bounds (period; units of dt)

;----SUBSETS OF RAINFALL--------------------

if var_plot eq 'rain' then begin
    var_str='RAINNC'
;    setmax=10 & setmin='0.'
    setmax=25 & setmin='0.'
endif else if var_plot eq 'pw' then begin
    var_str=var_plot
    setmax=70 & setmin=20
endif

  ;PLOT ALL RAINFALL
    myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
    figspecs=create_struct(figspecs,'figname',figdir+'imerg_mean_all_'+dat_str)
    figspecs.cbar_format='(i2)'
    figspecs.ndivs-=1
    wrf_myanmar_map_plot, dirs, mean(rain_sav[*,*,jjas],dimension=3,/nan,/double), lon, lat, figspecs, bounds=bounds_all

  ;0 - BOB
  ;1 - BOB OFFSHORE
  ;2 - BOB OFFSHORE N
  ;3 - BOB OFFSHORE S
  ;4 - BOB COASTAL
  ;5 - BOB COASTAL N
  ;6 - BOB COASTAL S
  ;7 - BANGLADESH

  nb=8

  ;BOUNDS FOR AVERAGING RAINFALL
;    bounds_all=[ $
;;      [85.,11,100,23.5], $  ; Large BoB
;      [85.,9,100,23.5], $   ; Large BoB
;;      [87.,15.5,96,23.5], $ ; Northern BoB
;      [88.,23.5,93.5,27] ]  ; Bangladesh
;    nb=(size(bounds_all,/dim))[1]

  ;NORTH-SOUTH BOB CUTOFF
    iy_cutoff=17. ; latitude separation point for two rainfall maxima
  ;COASTAL THRESHOLD
    coast_thresh=100. ; km

  ;SAVE ARRAYS
  avg_tser=fltarr(nt,nb)
  avg_tser2=fltarr(nt,nb)
  avg_signif=fltarr(nb)
  avg_signif2=fltarr(nb)

  for ib=0,nb-1 do begin
;  for ib=nb-1,nb-1 do begin

    ;SUBSET REGIONS
    if ib le 6 then $
      bounds=[85.,7,100,23.5] $ ; Large BoB
    else $
      bounds=[88.,23,93.5,27] ; Bangladesh
    if ib eq 1 then $
      bounds=[85.,14,94,23.5] ; Different box for "OFFSHORE"

    ;SUBSET TO BOUNDING BOX
      ix=where((lon ge bounds[0]) and (lon le bounds[2]))
      iy=where((lat ge bounds[1]) and (lat le bounds[3]))
      rain=rain_sav
      rain=rain[ix,*,*] & rain=rain[*,iy,*]
      rainjj=rain[*,*,jjas]
      rainjm=mean(rainjj,dimension=3,/nan,/double)

    ;SUBSET OFFSHORE/COASTAL
      cd3=coast_dist2[ix,*,*] & cd3=cd3[*,iy,*]
      ls3=lsmask2[ix,*,*] & ls3=ls3[*,iy,*]
      ;OFFSHORE
      if ib ge 1 and ib le 3 then begin
        ioffshore=where((cd3 ge coast_thresh) and (ls3 le 0.1),complement=inan)
        rain[inan]=!values.f_nan
      ;COASTAL
      endif else if ib ge 4 and ib le 6 then begin
        icoastal=where(cd3 le coast_thresh,complement=inan)
        rain[inan]=!values.f_nan
      endif

    ;SUBSET N/S
      ;NORTH
      if ib eq 2 or ib eq 5 then begin
        inorth=where(lat[iy] ge iy_cutoff,complement=inan)
        rain[*,inan,*]=!values.f_nan
      ;SOUTH
      endif else if ib eq 3 or ib eq 6 then begin
        isouth=where(lat[iy] lt iy_cutoff,complement=inan)
        rain[*,inan,*]=!values.f_nan
      endif

      ;PLOT RAINFALL SUBSET
        figspecs.figname=figdir+'imerg_mean_b'+strtrim(ib,2)+'_'+dat_str
        rainm=mean(rain,dimension=3,/nan,/double)
        igood=where(finite(rainm),complement=inan)
        rainjm[inan]=!values.f_nan

        wrf_myanmar_map_plot, dirs, rainjm, lon[ix], lat[iy], figspecs, bounds=bounds_all

;continue ; if only desire plots

;----CALCULATE SPECTRA--------------------

    ;STEP 1 - SPATIAL AVERAGE
      rain=mean(mean(temporary(rain),dimension=1,/nan,/double),dimension=1,/nan,/double)

    ;STEP 2 - BAND-PASS FILTER

      ;REMOVE PERIODS >= MIN_TS DAYS (THIS ALSO ACTS TO DETREND)
        min_ts = 121 ; n-days of running average

      ;RUNNING MEAN
;        rain_lp = smooth(rain,min_ts,/edge_truncate)
        rain_lp = fltarr(nd)
        it=indgen(min_ts)-(min_ts-1)/2
        ;MAIN
        for id=-1.*it[0],nd-max(it) do $
          rain_lp[id]=mean(rain[id+it],/nan,/double)
        ;ENDS
        for id=0,(min_ts-1)/2-1 do begin
          iit=id+it
          iit[where(iit lt 0)]=0
          rain_lp[id]=mean(rain[iit],/nan,/double)
        endfor
        for id=nd-max(it),nd-1 do begin
          iit=id+it
          iit[where(iit gt nd-1)]=nd-1
          rain_lp[id]=mean(rain[iit],/nan,/double)
        endfor

      ;SUBTRACT LP
        rain -= rain_lp

      ;KEEP BAND-PASS
        min_ts = 5 ; n-days of running average
        rain = smooth(temporary(rain),min_ts,/edge_truncate)

    ;STEP 3 - STANDARDIZE
;      rain = (rain - mean(rain,/nan,/double)) / stddev(rain,/nan,/double)
;DO THIS AFTER JJAS

    ;STEP 4 - CALCULATE POWER SPECTRUM

    if option eq 'fft' then begin

      ;STEP 4.1 - PARSE INTO INDIVIDUAL YEARS
      ;JJAS ONLY
        rain=rain[jjas]

    ;STEP 3 - STANDARDIZE
      rain = (rain - mean(rain,/nan,/double)) / stddev(rain,/nan,/double)

        ntpy=ntjjas/nyr
        rain_yr=fltarr(ntpy,nyr)
        caldat,time[jjas],mm,dd,yy
        years=yy[uniq(yy)]
        for iyr=0,nyr-1 do $
          rain_yr[*,iyr]=rain[where(yy eq (years[iyr]))]

      ;FREQUENCY INFORMATION
        nseg=nyr
        npseg=ntpy
        nband=npseg/2
        freq = (findgen(npseg)/npseg)[1:nband-1]
        t    = 1./freq

      ;STEP 4.2 - TAPER EDGES BY MIN_TS
        seg_tap=rain_yr
        ntap=min_ts
        for it=0,ntap-1 do begin
          scale=1.*it/float(ntap)
          seg_tap[it,*] *= scale
          seg_tap[npseg-1-it,*] *= scale
        endfor

      ;STEP 4.3 - CALCULATE SPECTRUM XX, AVERAGE ACROSS SEGMENTS
        spec_seg = abs(fft(seg_tap,-1,dimension=1,/double))^2
        spec_seg = spec_seg[0:nband-1,*]
;        spec = mean(spec_seg,dimension=2,/double,/nan)
      ;STEP 4.4 - 1-2-1 FILTER
        coef=0.25*[1.,2,1]
        coef=rebin(coef,3,nseg)
        spec_filt=spec_seg
        for iband=1,nband-2 do $
          spec_filt[iband,*] = total(spec_seg[iband-1:iband+1,*]*coef,1,/double)
        spec_seg=spec_filt
      ;STEP 4.5 - NOW AVERAGE OVER SEGMENTS
        spec = mean(spec_seg,dimension=2,/double,/nan)

      ;STEP 4.6 - NORMALIZE
        freq = (findgen(npseg)/npseg)[1:nband-1]
        t    = 1./freq
        area = total(spec,/double);*(deriv(freq))[0]
        spec = 1e2*temporary(spec)/area
;spec = spec[1:nband-1]

        area = total(spec_seg,1,/double);*(deriv(freq))[0]
        area = transpose(rebin(temporary(area),[nseg,nband]))
        spec_seg = 1e2*temporary(spec_seg)/area

      ;STEP 4.7 - RED NOISE AND SIGNIF
      
        ;LAG-1 AUTOCORRELATION
          a=fltarr(nseg)
          for iseg=0,nseg-1 do $
            ;a[iseg]=a_correlate(reform(seg_tap[*,iseg]),1)
            a[iseg]=0.5*(a_correlate(reform(seg_tap[*,iseg]),1) + sqrt(a_correlate(reform(seg_tap[*,iseg]),2)))
          efold=-1./alog(a) ; T = -Delta_t / ln(a)
;          print,''
;          print,'e-folding timescale (days):'
;          print,'  ',efold
          efold=-1./alog(mean(a,/nan,/double))
      
        ;RED SPEC
          freq2p = (2*!pi)*(findgen(npseg)/npseg)[1:nband-1]
          redspec=2.*efold/(1. + efold^2 * freq2p^2 )
          ;NORMALIZE
            area = total(redspec,/double);*deriv(freq),/double)
            redspec = 1e2*temporary(redspec)/area

        ;SIGSPEC
          nsamp = round( nyr * (1.*npseg/min_ts) ) ;nd/min_ts;(2*efold);1.*nd/npseg
          dof=nsamp
          siglvl=0.95
          chisq=chisqr_cvf(1.-siglvl,dof)
          t_stat=t_cvf(1.-siglvl,dof)
          sigspec=redspec*t_stat;chisq

      ;SAVE GLOBAL SPECTRUM
        if ib eq 0 then spec_sav=spec $
        else spec_sav = [[spec_sav],[spec]]

      ;PLOT WAVELET MULTI-PANEL
        plot_imerg_pspec, figdir
        figname=figdir+'imerg_pspec_series_'+tag+'_b'+strtrim(ib,2)
        plot_imerg_pspec_spectrum, figname, spec, freq, mspec=spec_seg, redspec=redspec, sigspec=sigspec

    endif else if option eq 'wave' then begin

      ;CALCULATE WAVELET SPECTRUM
        wavelet=wavelet_structure(time,rain,dt,pad,s0,dj,j1,mother,iavg,iavg2=iavg2)
  
      ;SAVE GLOBAL SPECTRUM
        if ib eq 0 then avg_var=wavelet.global_ws $
        else avg_var = [[avg_var],[wavelet.global_ws]]
  
;      ;SAVE SCALE-AVERAGED SPECTRUM
;        avg_tser[*,ib]=wavelet.scale_avg
;        avg_tser2[*,ib]=wavelet.scale_avg2
;        avg_signif[ib]=wavelet.scaleavg_signif

      ;PLOT WAVELET MULTI-PANEL
        figname=figdir+'wavelet_b'+strtrim(ib,2)
        plot_imerg_wavelet, figname, time, rain, iavg, wavelet

    endif

;------------------------------------------------------ Plotting

  endfor ; ibounding box

;  figname=figdir+'wavelet_global'
;  period=wavelet.period
;  plot_imerg_wave_global, figname, period, avg_var
;
;  figname=figdir+'wavelet_scaleavg'
;  period=wavelet.period
;  plot_imerg_wave_scaleavg, figname, time, iavg, iavg2, avg_tser, avg_tser2, avg_signif, avg_signif2

print,'Done!!'

end
