; 
; Calculate the wavelet spectrum for IMERG rainfall.
;
; Based on WAVETEST and utilizes the wavelet code package developed by:
; See "http://paos.colorado.edu/research/wavelets/"
; Written January 1998 by C. Torrence
;
; James Ruppert
; 6/11/21
; 
pro run_imerg_wavelet

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
  time_fil=timegen(start=julday(6,1,2000,0,0,0),final=julday(12,31,2020,23,59,59),step_size=1,units='Days')
  ;FILE IS INDEED IDENTICAL TO ORIG JJAS-2013-17 FILE OVER THAT DATE RANGE
;  im_fil=dirs.wkdir+'imerg/imerg/data/imerg_3B-HHR.MS.MRG.3IMERG.V06_daily_2000-2021_latlonsubset_jjas2013-2017.nc4' ; test match to other file
  npd_imerg=48
  era_dir=dirs.wkdir+'wrf_wps/wrfv4/output_myanmar/JJAS_2013-2017/era5/'
  coast_dist_fil=dirs.home+'idl/code/misc/dist2coast.nc'

  figdir=dirs.figdir+'myanmar/imerg/'

;----READ RAIN--------------------

  ;SELECTED TIME ARRAY
    time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
      final=julday(mm_plot[1],dd_plot[1],yy_plot[1],23,59,59),step_size=1,units='Days');30,units='Minutes')
    nt=n_elements(time)
    nd=nt;/npd_imerg
    nyr=yy_plot[1]-yy_plot[0]+1

  rain_sav=read_nc_imerg(time_imerg,im_fil,lon=lonim,lat=latim,bounds=bounds) ; already in mm/d
  nx=n_elements(lonim)
  ny=n_elements(latim)

  ;SAVE BASIC TIME-MEAN
  rain_mn_sav=mean(rain_sav,dimension=3,/nan,/double)

;----PLOT TIME-MEAN RAINFALL--------------------

if iplot_mn_rain then begin

  var_str='RAINNC'
  setmax=30 & setmin='0.'
setmax=10
  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figspecs=create_struct(figspecs,'figname',figdir+'imerg_mean_'+dat_str)
  figspecs.cbar_format='(i2)'

  wrf_myanmar_map_plot, dirs, rain_mn_sav, lonim, latim, figspecs;, bounds=bounds

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

;----PREP RAINFALL ARRAY--------------------

  ;STEP 0 OPTION
  if do_coast_subset then begin
    coast_dist=coast_dist[ix,*] & coast_dist=coast_dist[*,iy]
    lsmask=lsmask[ix,*] & lsmask=lsmask[*,iy]
    ;ADD TIME INDICES TO THESE TO SIMPLIFY WHERE FUNCTIONS
      coast_dist2=rain_sav
      for id=0,nd-1 do coast_dist2[*,*,id]=coast_dist
      lsmask2=rain_sav
      for id=0,nd-1 do lsmask2[*,*,id]=lsmask
    ;USE COASTAL DISTANCE VARIABLE
      coast_thresh=100.
    isea=where((coast_dist2 ge coast_thresh) and (lsmask2 le 0.1),complement=inot)
    ;icoast=where((coast_dist2 le coast_thresh),complement=inot)
    rain[inot]=!values.f_nan
    ;NEW RAINFALL MAP
      rain_mn=mean(rain,dimension=3,/nan,/double)
      wrf_myanmar_map_plot, dirs, rain_mn, lonim[ix], latim[iy], figspecs, bounds=bounds
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

;----BOX-AVERAGED SPECTRA--------------------

  ;BOUNDS FOR AVERAGING RAINFALL
    bounds_all=[ $
;      [85.,13,96,23.5], $
      [85.,11,100,23.5], $  ; Large BoB
      [87.,15.5,96,23.5], $ ; Northern BoB
      [88.,23.5,93.5,27] ]  ; Bangladesh
    nb=(size(bounds_all,/dim))[1]

  ;PLOT BOUNDING BOXES
    var_str='RAINNC'
    setmax=10 & setmin='0.'
    myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
    figspecs=create_struct(figspecs,'figname',figdir+'imerg_mean_avgboxes_'+dat_str)
    figspecs.cbar_format='(i2)'
;    wrf_myanmar_map_plot, dirs, rain_mn_sav, lonim, latim, figspecs, bounds=bounds_all

  ;SAVE ARRAYS
  avg_tser=fltarr(nt,nb)
  avg_tser2=fltarr(nt,nb)
  avg_signif=fltarr(nb)
  avg_signif2=fltarr(nb)

  for ib=0,nb-1 do begin

    ;SUBSET TO BOUNDING BOX
      ix=where((lonim ge bounds_all[0,ib]) and (lonim le bounds_all[2,ib]))
      iy=where((latim ge bounds_all[1,ib]) and (latim le bounds_all[3,ib]))
      rain=rain_sav[ix,*,*] & rain=rain[*,iy,*]

    ;STEP 1 - SPATIAL AVERAGE
      rain=mean(mean(temporary(rain),dimension=1,/nan,/double),dimension=1,/nan,/double)

    ;STEP 2 - BAND-PASS FILTER
      ;REMOVE PERIODS >= MIN_TS DAYS (THIS ALSO ACTS TO DETREND)
        min_ts = 91;121 ; n-days of running average
        rain_lp = smooth(rain,min_ts,/edge_truncate)
      ;SUBTRACT LP
        rain -= rain_lp
      ;KEEP BAND-PASS
        min_ts = 5 ; n-days of running average
        rain = smooth(temporary(rain),min_ts,/edge_truncate)

    ;STEP 3 - STANDARDIZE
      rain = (rain - mean(rain,/nan,/double)) / stddev(rain,/nan,/double)

    ;STEP 4 - CALCULATE WAVELET SPECTRUM
      wavelet=wavelet_structure(time,rain,dt,pad,s0,dj,j1,mother,iavg,iavg2=iavg2)

    ;SAVE GLOBAL SPECTRUM
      if ib eq 0 then avg_var=wavelet.global_ws $
      else avg_var = [[avg_var],[wavelet.global_ws]]

    ;SAVE SCALE-AVERAGED SPECTRUM
      avg_tser[*,ib]=wavelet.scale_avg
      avg_signif[ib]=wavelet.scaleavg_signif

    ;SAVE SCALE-AVERAGED SPECTRUM FOR BAND2
;      wavelet2=wavelet_structure(time,rain,dt,pad,s0,dj,j1,mother,iavg2)
      avg_tser2[*,ib]=wavelet.scale_avg2
;      avg_signif2[ib]=wavelet2.scaleavg_signif

;------------------------------------------------------ Plotting

    figname=figdir+'wavelet_b'+strtrim(ib+1,2)
    plot_imerg_wavelet, figname, time, rain, iavg, wavelet

  endfor ; ibounding box

  figname=figdir+'wavelet_global'
  period=wavelet.period
  plot_imerg_wave_global, figname, period, avg_var

  figname=figdir+'wavelet_scaleavg'
  period=wavelet.period
  plot_imerg_wave_scaleavg, figname, time, iavg, iavg2, avg_tser, avg_tser2, avg_signif, avg_signif2

print,'Done!!'

end
