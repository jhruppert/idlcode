; 
; Create lag-regression plots to show patterns associated with heavy rainfall
; in target areas, using IMERG and ERA5
;
; James Ruppert
; 7/20/21
; 
pro run_imerg_lagreg

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

var_maps=1 ; Plot maps of variance?
  ifft=0 ; do FFT & Butterworth OR simple running means

var_plot='rain';
var_plot='pw';'rain';

iera=1 ; read/plot ERA?
;ERAi SETTINGS
;LEVEL SELECTION
  psel_era=850;500;700;925

if var_plot ne 'rain' and ~iera then message,'Must include ERA if you want any other field'

;BOUNDS FOR READ-IN
  bounds=[78,6,104,29]

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

  maindir=dirs.wkdir
;  im_fil=maindir+'imerg/imerg/data/jjas_2013-2017_daymean_3B-HHR.MS.MRG.3IMERG.V06B.nc4'
  im_fil=dirs.wkdir+'imerg/imerg/data/imerg_3B-HHR.MS.MRG.3IMERG.V06_daily_2000-2021_latlonsubset.nc4'
  time_imerg=timegen(start=julday(6,1,2000,0,0,0),final=julday(12,31,2020,23,59,59),step_size=1,units='Days')
;  npd_imerg=48

  era_dir=maindir+'era5/jjas_2013-2017/';ERAi-JJAS13-17-pl.nc4'
  era_fil=era_dir+'ERAi-JJAS13-17-pl_dayavg.nc4'
  era_sfil=era_dir+'ERAi-JJAS13-17-sl_dayavg.nc4'
;  npd_era=24
era_dir=dirs.scdir+'era5/'
era_fil=era_dir+'ERA5-20000101-20201231-pl_dayavg.nc'
era_pw=era_dir+'ERA5-20000101-20201231-pw_dayavg.nc'
  time_era=timegen(start=julday(1,1,2000,0,0,0),final=julday(12,31,2020,23,59,59),step_size=1,units='Days')

  irain_dir=dirs.figdir+'myanmar/imerg/'

;Local SOLAR time conversion
;local=6;round(mean(dims.lon)/360.*24.) ; deg lon --> hours
;print,'Adding +'+strtrim(local,2)+' for LT'

;----ONE TIME ARRAY FOR DAILY AVERAGE TIME SERIES--------------------

  ;SELECTED TIME ARRAY
    time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
      final=julday(mm_plot[1],dd_plot[1],yy_plot[1],23,59,59),step_size=1,units='Days');30,units='Minutes')
    nt=n_elements(time)
    nd=nt;/npd_imerg
    nyr=yy_plot[1]-yy_plot[0]+1

  ;SAVE JJAS INDICES
    caldat,time,mm,dd,yy
    jjas=where((mm ge 6) and (mm le 9))

;----READ RAIN--------------------

  rain_sav=read_nc_imerg(time_imerg,im_fil,lon=lonim,lat=latim,bounds=bounds) ; already in mm/d
  nx=n_elements(lonim)
  ny=n_elements(latim)

  ;SAVE BASIC TIME-MEAN
;  rain_mn_sav=mean(rain_sav,dimension=3,/nan,/double)

  if var_plot eq 'rain' then begin
    lon=lonim
    lat=latim
    var=rain_sav
  endif

;----READ ERA--------------------

if iera then begin

  ;DIMENSIONS
  eralon=read_nc_var(era_fil,'lon')
  eralat=read_nc_var(era_fil,'lat')
  nxera=n_elements(eralon)
  nyera=n_elements(eralat)

  ;PLEVS
  p_era=reform(read_nc_var(era_fil,'plev'))*1d-2 ; Pa --> hPa
  nzera=n_elements(p_era)

  ;SELECTED LEVEL
  izlev_era=(where(p_era eq psel_era,count))[0]
  if count eq 0 then stop

  ;DIFFERENT START TIME FROM IMERG
  tdiff=abs(time_era-time[0])
  it0=(where(tdiff eq min(tdiff)))[0]
  ntread=n_elements(time_era)-it0
  if ntread ne nt then stop

  ;WINDS
    count=[nxera,nyera,1,ntread] & offset=[0,0,izlev_era,it0] ; x, y, p, t
    u=reform(read_nc_var(era_fil,'var131',count=count,offset=offset))
    v=reform(read_nc_var(era_fil,'var132',count=count,offset=offset))

  ;PW
  if var_plot eq 'pw' then begin;icalc_pw=1 else icalc_pw=0
    count=[nxera,nyera,ntread] & offset=[0,0,it0] ; x, y, t
    var=reform(read_nc_var(era_pw,'var137',count=count,offset=offset))
    var_sav=var
    nx=nxera
    ny=nyera
    lon=eralon
    lat=eralat
  endif

endif ; iera

;----READ LAND/SEA MASK--------------------

;  spawn,'ls '+era_dir+'met_em.d01*',era_fil
;  fil=era_fil[0]
;;  topo=reform(read_nc_var(fil,'HGT_M'))*1e3 ; m --> km
;  lsmask=reform(read_nc_var(fil,'LANDMASK'))
;  dimsw=size(lsmask,/dimen)
;  nxw=dimsw[0] & nyw=dimsw[1]
;
;  lonw=read_nc_var(fil,'XLONG_M') ; deg
;  lonw=reform(lonw[*,0])
;  latw=read_nc_var(fil,'XLAT_M') ; deg
;  latw=reform(latw[0,*])
;
;  ;INTERPOLATE ONTO IMERG GRID
;    lsmaskx=fltarr(nx,nyw)
;    for iy=0,nyw-1 do lsmaskx[*,iy]=interpol(reform(lsmask[*,iy]),lonw,lonim)
;    lsm2=fltarr(nx,ny)
;    for ix=0,nx-1 do lsm2[ix,*]=interpol(reform(lsmaskx[ix,*]),latw,latim)
;    lsmask=lsm2
;
;;----READ COASTAL DISTANCE--------------------
;
;  dist_fil=read_nc_var(coast_dist_fil,'coast_dist')
;  londist=reform(dist_fil[0,*])
;  latdist=reform(dist_fil[1,*])
;  coast_dist=reform(dist_fil[2,*])
;
;  ;KEEP ONLY IMERG DOMAIN
;  ixkeep=where((londist ge min(lonim)) and (londist le max(lonim)) and (latdist ge min(latim)) and (latdist le max(latim)))
;  londist=londist[ixkeep]
;  latdist=latdist[ixkeep]
;  coast_dist=coast_dist[ixkeep]
;
;  ;INTERPOLATE ONTO IMERG GRID
;    triangulate,londist,latdist,tri
;    coast_dist_im=griddata(londist,latdist,coast_dist,/linear,triangles=tri,xout=lonim,yout=latim,/grid)
;    coast_dist=coast_dist_im

;----READ IN RAINFALL INDICES--------------------

  irain_all=fltarr(nd)
  irain_coast=irain_all
  irain_offshore=irain_all
;  irain_onshore=irain_all
  irain_bang=irain_all
;  irain_smd=irain_all
  openr,1,irain_dir+'rain_all_'+dat_str+'.txt'
  readf,1,irain_all & close,1
  openr,1,irain_dir+'rain_coast_'+dat_str+'.txt'
  readf,1,irain_coast & close,1
  openr,1,irain_dir+'rain_offshore_'+dat_str+'.txt'
  readf,1,irain_offshore & close,1
;  openr,1,irain_dir+'rain_onshore_'+dat_str+'.txt'
;  readf,1,irain_onshore & close,1
  openr,1,irain_dir+'rain_bangladesh_'+dat_str+'.txt'
  readf,1,irain_bang & close,1
;  openr,1,irain_dir+'rain_smd_'+dat_str+'.txt'
;  readf,1,irain_smd & close,1

  ;OFFSHORE - ONSHORE; COAST
;  irain_cdiff = irain_offshore - irain_onshore

;----FILTER--------------------

;BAND-PASS FILTER

  ;REMOVE PERIODS >= MIN_TS DAYS (THIS ALSO ACTS TO DETREND)
    min_ts = 121 ; n-days of running average
  ;SUBTRACT RUNNING MEAN
    var -= smooth(var,[0,0,min_ts],/edge_truncate)
    if iera then u -= smooth(u,[0,0,min_ts],/edge_truncate)
    if iera then v -= smooth(v,[0,0,min_ts],/edge_truncate)

    irain_all -= smooth(irain_all,min_ts,/edge_truncate)
    irain_coast -= smooth(irain_coast,min_ts,/edge_truncate)
    irain_offshore -= smooth(irain_offshore,min_ts,/edge_truncate)
    irain_bang -= smooth(irain_bang,min_ts,/edge_truncate)

  ;KEEP BAND-PASS
    min_ts = 5 ; n-days of running average
    var = smooth(temporary(var),[0,0,min_ts],/edge_truncate)
    if iera then u = smooth(temporary(u),[0,0,min_ts],/edge_truncate)
    if iera then v = smooth(temporary(v),[0,0,min_ts],/edge_truncate)

    irain_all = smooth(temporary(irain_all),min_ts,/edge_truncate)
    irain_coast = smooth(temporary(irain_coast),min_ts,/edge_truncate)
    irain_offshore = smooth(temporary(irain_offshore),min_ts,/edge_truncate)
    irain_bang = smooth(temporary(irain_bang),min_ts,/edge_truncate)

  ;REMOVE BASIC MEAN
    varm = mean(var,dimension=3,/nan,/double)
    varm = rebin(varm,nx,ny,nd)
    var -= varm


;----KEEP JJAS ONLY--------------------

  nd=n_elements(jjas)

  var=var[*,*,jjas]
  if iera then begin
    u=u[*,*,jjas]
    v=v[*,*,jjas]
  endif

  irain_all=irain_all[jjas]
  irain_coast=irain_coast[jjas]
  irain_offshore=irain_offshore[jjas]
  irain_bang=irain_bang[jjas]

;----VARIANCE MAPS--------------------

if var_maps then begin

if var_plot eq 'rain' then begin
  var_str='RAINNC'
  setmax=200. & setmin='0.'
  cbform='(i2)'
;  if ifft then begin
    setmax=1 & setmin='0.'
;setmax=0.6 & setmin=0.2
    cbform='(f4.2)'
;  endif
endif else if var_plot eq 'pw' then begin
  var_str=var_plot
  setmax=30 & setmin='0.'
  cbform='(i2)'
setmax=0.4 & setmin='0.'
cbform='(f4.2)'
endif

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=3
  figdir=dirs.figdir+'myanmar/imerg/'
  figspecs=create_struct(figspecs,'figname',figdir)
  figspecs.cbar_format=cbform
  figspecs.ndivs-=1
  if ifft then figspecs.cbar_tag=''
;figspecs.cbar_tag='%'

  ;SAVE FOR NORMALIZATION
  total_var = variance(var,dimension=3,/nan,/double)

for iband=1,2 do begin

if iband eq 1 then begin
  figtag='biweekly'
  ;BIWEEKLY
    s0=7
    s1=25
  ;FOR FFT BUTTERWORTH
    t_sel=14
    cutoff=50
    nwts=7
endif else if iband eq 2 then begin
  figtag='intra'
  ;INTRASEASONAL
    s0=30
    s1=60
  ;FOR FFT BUTTERWORTH
    t_sel=31
    cutoff=35
    nwts=7
endif

    var_plt=var

    if ~ifft then begin
      var_plt = smooth(var_plt,[0,0,s0],/edge_truncate)
      var_plt -= smooth(var_plt,[0,0,s1],/edge_truncate)
    endif

;figspecs.figname=figdir+'imerg_var_biweekly_'+dat_str
;figname=figspecs.figname
;var2=mean(mean(var_sav,dimension=1,/nan,/double),dimension=1,/nan,/double)
;nband=nd/2
;spec = abs(fft(var2,-1,dimension=1,/double))^2
;spec = spec[0:nband-1,*]
;stop
;freq=(findgen(nd)/nd)[1:nband-1]
;plot_imerg_pspec, figdir
;plot_imerg_pspec_series, figname, spec, freq

  if ifft then begin
  
    ;nband=nd/2
    freq = (findgen(nd)/nd);[1:nband-1]
    per    = 1./freq
    
    t2=abs(t_sel-per)
    loc=where(t2 eq min(t2))
    
    butter_temp=butterworth(nd,cutoff=cutoff,order=nwts,/origin)
    butter_band=shift(butter_temp, -1*((nd-1)/2-loc) )
    ;PRINT PERIODS RETAINED BY FILTER
      print,per[where(butter_band ge 0.25)]
    butter_band=transpose(rebin(butter_band,nd,nx,ny),[1,2,0])
    
    if iband eq 1 then $
      temp_fft1=fft(var_plt,-1,dimension=3) ; DO NOT NEED TO RECALCULATE THIS A 2ND TIME
    
    var_plt=fft(temp_fft1*butter_band,1,dimension=3)
    var_plt=abs(var_plt)
  
  endif

    ;CALCULATE VARIANCE
    vari=variance(var_plt,dimension=3,/nan,/double)
    var_plt=vari & vari=0

;if ifft then 
;var_plt/=max(var_plt)
var_plt /= total_var

    stats,var_plt

    figspecs.figname=figdir+'imerg_var_'+var_str+'_'+figtag+'_'+dat_str
if ifft then figspecs.figname+='_fft'
    figspecs.title='Variance'

    wrf_myanmar_map_plot, dirs, var_plt, lon, lat, figspecs

endfor

exit
endif ; var_maps

;----STANDARDIZE--------------------

  std_var = stddev(var,dimension=3,/nan,/double)
  std_var = rebin(temporary(std_var),[nx,ny,nd])
  if iera then begin
    std_u = stddev(u,dimension=3,/nan,/double)
    std_u = rebin(temporary(std_u),[nxera,nyera,nd])
    std_v = stddev(v,dimension=3,/nan,/double)
    std_v = rebin(temporary(std_v),[nxera,nyera,nd])
  endif

  std_all=stddev(irain_all,/nan,/double)
  std_coast=stddev(irain_coast,/nan,/double)
  std_offshore=stddev(irain_offshore,/nan,/double)
;  std_onshore=stddev(irain_onshore,/nan,/double)
  std_bang=stddev(irain_bang,/nan,/double)
;  std_smd=stddev(irain_smd,/nan,/double)
;  std_cdiff=stddev(irain_cdiff,/nan,/double)

  print,'Standard Deviations:'
  print,'All:',std_all
  print,'Coast:',std_coast
  print,'Offshore:',std_offshore
;  print,'Onshore:',std_onshore
  print,'Bangladesh:',std_bang
;  print,'SMD:',std_smd
;  print,'Offshore-Onshore:',std_cdiff

  varm = mean(var,dimension=3,/nan,/double)
  varm = rebin(varm,nx,ny,nd)
  var = (var - varm)/std_var
  if iera then begin
    um = mean(u,dimension=3,/nan,/double)
    um = rebin(um,nxera,nyera,nd)
    u = (u - um)/std_u
    vm = mean(v,dimension=3,/nan,/double)
    vm = rebin(vm,nxera,nyera,nd)
    v = (v - vm)/std_v
  endif

  irain_all = (irain_all - mean(irain_all,/nan,/double)) / std_all
  irain_coast = (irain_coast - mean(irain_coast,/nan,/double)) / std_coast
  irain_offshore = (irain_offshore - mean(irain_offshore,/nan,/double)) / std_offshore
;  irain_onshore = (irain_onshore - mean(irain_onshore,/nan,/double)) / std_onshore
  irain_bang = (irain_bang - mean(irain_bang,/nan,/double)) / std_bang
;  irain_smd = (irain_smd - mean(irain_smd,/nan,/double)) / std_smd
;  irain_cdiff = (irain_cdiff - mean(irain_cdiff,/nan,/double)) / std_cdiff

  ;N-CASES where |irain| > 1 sigma
  iall = where(abs(irain_all) ge 1,nall)
  icoast = where(abs(irain_coast) ge 1,ncoast)
  ioffshore = where(abs(irain_offshore) ge 1,noffshore)
;  ionshore = where(abs(irain_onshore) ge 1,nonshore)
  ibang = where(abs(irain_bang) ge 1,nbang)
;  ismd = where(abs(irain_smd) ge 1,nsmd)
;  icdiff = where(abs(irain_cdiff) ge 1,ncdiff)

  print,'N where |index| >= 1:'
  print,'All:',nall
  print,'Coast:',ncoast
  print,'Offshore:',noffshore
;  print,'Onshore:',nonshore
  print,'Bangladesh:',nbang
;  print,'SMD:',nsmd
;  print,'Offshore-Onshore:',ncdiff

;----CALCULATE REGRESSIONS--------------------



end
