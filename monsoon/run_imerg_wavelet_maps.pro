; 
; Maps based on wavelet spectrum for IMERG rainfall.
;
; Based on WAVETEST and utilizes the wavelet code package developed by:
; See "http://paos.colorado.edu/research/wavelets/"
; Written January 1998 by C. Torrence
;
; James Ruppert
; 6/13/21
; 
pro run_imerg_wavelet_maps

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

;BOUNDS FOR READ-IN
  bounds=[78,6,104,29]
;bounds=[88.,23.5,93.5,27]
;bounds=[87.,15.5,96,23.5]
  tag='box3'

iread_old=1 ; read in from old instead of running transforms?
do_coast_subset=0 ; subset rainfall data by coast/land-sea mask?
iplot_mn_rain=0   ; plot time-mean rainfall?

;SELECT DATE RANGE
  yy_plot=[2000,2020]
  mm_plot=[6,12]
;  yy_plot=[2015,2020]
;  mm_plot=[1,12]
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
  time_fil=timegen(start=julday(6,1,2000,0,0,0),final=julday(12,31,2021,23,59,59),step_size=1,units='Days')
  ;FILE IS INDEED IDENTICAL TO ORIG JJAS-2013-17 FILE OVER THAT DATE RANGE
;  im_fil=dirs.wkdir+'imerg/imerg/data/imerg_3B-HHR.MS.MRG.3IMERG.V06_daily_2000-2021_latlonsubset_jjas2013-2017.nc4' ; test match to other file
  npd_imerg=48
  era_dir=dirs.wkdir+'wrf_wps/wrfv4/output_myanmar/JJAS_2013-2017/era5/'
  coast_dist_fil=dirs.home+'idl/code/misc/dist2coast.nc'

  figdir=dirs.figdir+'myanmar/imerg/'

;----PREP FOR WAVELET SPECTRA--------------------

  dt = 1 ; time interval (1 = 1 per time step)
;  pad = 1
  s0 = 2*dt  ; this says start at a scale of 2*dt
  dj = 0.25  ; this will do 4 sub-octaves per octave
  j1 = 9./dj ; this says do 9 powers-of-two with dj sub-octaves each
  mother = 'Morlet'
  ;mother = 'Dog'
  iavg1=[7,25] ; averaging bounds (period; units of dt)
  iavg2=[25,80] ; averaging bounds (period; units of dt)

;----READ RAIN--------------------

  ;SELECTED TIME ARRAY
    time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
      final=julday(mm_plot[1],dd_plot[1],yy_plot[1],23,59,59),step_size=1,units='Days');30,units='Minutes')
    nt=n_elements(time)
    nd=nt;/npd_imerg
    nyr=yy_plot[1]-yy_plot[0]+1

  ;IMERG LAT/LON
    lonim=read_nc_var(im_fil,'lon')
    latim=read_nc_var(im_fil,'lat')
    nx=n_elements(lonim)
    ny=n_elements(latim)

  ;USE DAILY-MEAN RAINFALL
    ;TIME SUBSET
      temp=reform(read_nc_var(im_fil,'time')) ; Seconds since 1970-01-01 00:00:00Z
      it=where(time_fil ge time[0] and time_fil le max(time))
      nt_read=n_elements(it)
      if nt_read ne nt then stop
    ;SPACE SUBSET
      ix=where((lonim ge bounds[0]) and (lonim le bounds[2]),nx)
      iy=where((latim ge bounds[1]) and (latim le bounds[3]),ny)
      lonim=lonim[ix] & latim=latim[iy]

  if iread_old then begin

    outfile=dirs.wkdir+'imerg/idl_iavg1_'+dat_str+'.nc'
outfile=dirs.wkdir+'imerg/idl_iavg1_'+dat_str+'_unnorm.nc'
    avg1=reform(read_nc_var(outfile,'avg1'))
    outfile=dirs.wkdir+'imerg/idl_iavg2_'+dat_str+'.nc'
outfile=dirs.wkdir+'imerg/idl_iavg2_'+dat_str+'_unnorm.nc'
    avg2=reform(read_nc_var(outfile,'avg2'))

  endif else begin

    count=[ny,nx,nt] & offset=[iy[0],ix[0],it[0]] ; x,y,t
    rain_sav=reform(read_nc_var(im_fil,'precipitationCal',count=count,offset=offset)) ; already in mm/d
    rain_sav=transpose(temporary(rain_sav),[1,0,2])

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

;----CALCULATE WAVELET SPECTRA--------------------

  ;SAVE ARRAYS
  avg1=fltarr(nx,ny)
  avg2=fltarr(nx,ny)

  for ix=0,nx-1 do begin
  for iy=0,ny-1 do begin

    ;STEP 1 - SPATIAL AVERAGE
      rain=reform(rain_sav[ix,iy,*])

    ;STEP 2 - BAND-PASS FILTER
      ;REMOVE PERIODS >= MIN_TS DAYS (THIS ALSO ACTS TO DETREND)
        min_ts = 91 ; n-days of running average
        rain_lp = smooth(rain,min_ts,/edge_truncate)
      ;SUBTRACT LP
        rain -= rain_lp
      ;KEEP BAND-PASS
        min_ts = 5 ; n-days of running average
        rain = smooth(temporary(rain),min_ts,/edge_truncate)

    ;STEP 3 - STANDARDIZE
      rain = (rain - mean(rain,/nan,/double)) / stddev(rain,/nan,/double)

    ;STEP 4 - CALCULATE WAVELET SPECTRUM
      wavelet=wavelet_structure(time,rain,dt,pad,s0,dj,j1,mother,iavg1,iavg2=iavg2)
      ;SAVE SCALE-AVERAGED SPECTRUM
        avg1[ix,iy]=mean(wavelet.scale_avg,/nan,/double)
        avg2[ix,iy]=mean(wavelet.scale_avg2,/nan,/double)

;------------------------------------------------------ Plotting

  endfor ; iy
  endfor ; ix

  ;WRITE OUT AVERAGES TO AVOID RECALCULATING
  outfile=dirs.wkdir+'imerg/idl_iavg1_'+dat_str+'_unnorm.nc'
  write_sing_ncvar,outfile,avg1,'avg1'
  outfile=dirs.wkdir+'imerg/idl_iavg2_'+dat_str+'_unnorm.nc'
  write_sing_ncvar,outfile,avg2,'avg2'

  endelse
stats,avg1
stats,avg2
  ;BOUNDS FOR AVERAGING RAINFALL (just to include in plots)
    bounds_all=[ $
;      [85.,13,96,23.5], $
      [85.,11,100,23.5], $  ; Large BoB
      [87.,15.5,96,23.5], $ ; Northern BoB
      [88.,23.5,93.5,27] ]  ; Bangladesh
    nb=(size(bounds_all,/dim))[1]

  var_str='RAINNC'
;  setmax=0.5 & setmin='0.'
;  setmax=0.6 & setmin=0.4
setmax=10 & setmin=5
  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figspecs=create_struct(figspecs,'figname',figdir+'wavelet_varmap_'+dat_str)
  ifigname=figspecs.figname
  figspecs.cbar_format='(f4.2)'
  figspecs.cbar_tag='Variance'

  ;NEW COLOR TABLE
;    figspecs.col_table=66;70;71
;    ncols=n_elements(figspecs.colors)
;    colors=findgen(ncols)/(ncols-1)*255./2+255/2
;    irev=0;1
;    if irev then colors=reverse(colors)
;    figspecs.colors=colors
    figspecs.cbar_format='(i2)'
;;    figspecs.ndivs-=1

  avgtag=strtrim(iavg1[0],2)+'-'+strtrim(iavg1[1],2)+'d'
  figspecs.title=avgtag
  figspecs.figname=ifigname+'_'+avgtag
  wrf_myanmar_map_plot, dirs, avg1, lonim, latim, figspecs, bounds=bounds_all

  var_str='RAINNC'
;  setmax=0.5 & setmin='0.'
;  setmax=0.4 & setmin=0.25
setmax=20 & setmin=10
  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figspecs=create_struct(figspecs,'figname',figdir+'wavelet_varmap_'+dat_str)
  ifigname=figspecs.figname
  figspecs.cbar_format='(f4.2)'
  figspecs.cbar_tag='Variance'

  ;NEW COLOR TABLE
;    figspecs.col_table=66;70;71
;    ncols=n_elements(figspecs.colors)
;    colors=findgen(ncols)/(ncols-1)*255./2+255/2
;    irev=0;1
;    if irev then colors=reverse(colors)
;    figspecs.colors=colors
    figspecs.cbar_format='(i2)'
;;    figspecs.ndivs-=1

  avgtag=strtrim(iavg2[0],2)+'-'+strtrim(iavg2[1],2)+'d'
  figspecs.title=avgtag
  figspecs.figname=ifigname+'_'+avgtag
  wrf_myanmar_map_plot, dirs, avg2, lonim, latim, figspecs, bounds=bounds_all
;stop
print,'Done!!'

end
