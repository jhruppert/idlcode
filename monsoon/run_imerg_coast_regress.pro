; 
; Regression analysis of daily-mean rainfall (Myanmar/Meiyu IMERG data) and ERA-i
; using area-averaged rainfall time series calculated in run_imerg_coast_index.pro
;
; IMERG data: 2013-2017 JJAS
;
; James Ruppert
; 12/07/20
; 
pro run_imerg_coast_regress

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

;BOUNDS FOR READ-IN
;  bounds=[78,6,104,29]
  bounds=[70,2,108,32]

iplot_mn_rain=0

jjas_only=1 ; Subset time to JJAS?

dolag=0 ; Do lagged regressions?

var_plot='rain';
;var_plot='pw';'rain';

iera=1 ; read/plot ERA?
;ERAi SETTINGS
;LEVEL SELECTION
  psel_era=850;500;700;925

if var_plot ne 'rain' and ~iera then message,'Must include ERA if you want any other field'

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

;USE WRF MET_EM DATA FOR LAND/SEA MASK
;expname='myanmar'
;cases=['ctl']
;idir='2dom'
;if idir eq '2dom' then $
;  domtag='d01'
;casedir=dirs.scdir+expname+'/WRF/experiment/'+idir+'/'
;dirs.figdir+=expname+'/'
;config_wrfexp, casedir=casedir,cases=cases,dirs=dirs,$
;  dims=dims, vars=vars, nfils_set=nfils_set, domtag=domtag;, /verbose


;----OB DIRECTORIES--------------------

  maindir=dirs.wkdir
;  im_fil=maindir+'imerg/imerg/data/jjas_2013-2017_daymean_3B-HHR.MS.MRG.3IMERG.V06B.nc4'
  im_fil=dirs.wkdir+'imerg/imerg/data/imerg_3B-HHR.MS.MRG.3IMERG.V06_daily_2000-2021_latlonsubset.nc4'
;  npd_imerg=48

  era_dir=maindir+'era5/'
  era_fil=era_dir+'ERA5-20000101-20201231-pl_dayavg.nc'
  era_pw=era_dir+'ERA5-20000101-20201231-pw_dayavg.nc'

  irain_dir=dirs.figdir+'myanmar/imerg/rainfall_indices/'

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

  rain=read_nc_imerg(time,im_fil,lon=lonim,lat=latim,bounds=bounds) ; already in mm/d
  nx=n_elements(lonim)
  ny=n_elements(latim)

  ;SAVE BASIC TIME-MEAN
;  rain_mn_sav=mean(rain_sav,dimension=3,/nan,/double)

;;READ TOPOGRAPHY AND LAND MASK FROM ERA5 MET_EM FILES
;  spawn,'ls '+era_dir+'met_em.'+domtag+'*',era_fil
;  fil=era_fil[0]
;  topo=reform(read_nc_var(fil,'HGT_M'))*1e3 ; m --> km
;  lsmask=reform(read_nc_var(fil,'LANDMASK'))

  if var_plot eq 'rain' then begin
    lon=lonim
    lat=latim
    var=rain
  endif

;----READ ERA--------------------

if iera then begin

  u=read_nc_era5(time,era_fil,'var131',plev=psel_era,lon=eralon,lat=eralat,bounds=bounds)
  v=read_nc_era5(time,era_fil,'var132',plev=psel_era,lon=eralon,lat=eralat,bounds=bounds)
  nxera=n_elements(eralon)
  nyera=n_elements(eralat)

  ;PW
  if var_plot eq 'pw' then begin
    var=read_nc_era5(time,era_pw,'var137',lon=lon,lat=lat,bounds=bounds)
    nx=nxera
    ny=nyera
  endif

endif ; iera


;----MEAN MAPS--------------------

if iplot_mn_rain then begin

  var_plt=var[*,*,jjas]

if var_plot eq 'rain' then begin
  var_str='RAINNC'
  setmax=25 & setmin='0.'
  cbform='(i2)'
endif else if var_plot eq 'pw' then begin
  var_str=var_plot
  setmax=70 & setmin=30
  cbform='(i2)'
endif

  var_mn=mean(var_plt,dimension=3,/nan,/double)

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=3
  figdir=dirs.figdir+'myanmar/imerg/'
  figspecs=create_struct(figspecs,'figname',figdir+'imerg_mean_'+var_plot+'_'+dat_str)
  figspecs.cbar_format=cbform
  figspecs.ndivs-=1
  figspecs.title='IMERG Rainfall (JJAS, 2000-2020)'

  if iera then begin
    u_mn=mean(u[*,*,jjas],dimension=3,/nan,/double)
    v_mn=mean(v[*,*,jjas],dimension=3,/nan,/double)
    wind=create_struct('u',u_mn,'v',v_mn,'x',eralon,'y',eralat)
    wspd=sqrt(u_mn^2+v_mn^2)
    cvar=create_struct('cvar',wspd,'x',eralon,'y',eralat)
  endif

  wrf_myanmar_map_plot, dirs, var_mn, lon, lat, figspecs, wind=wind, cvar=cvar
exit
endif

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

  filter_monsoon, irain_all, 'rmean', var_bw=irain_all_bw, var_intra=irain_all_intra
  filter_monsoon, irain_coast, 'rmean', var_bw=irain_coast_bw, var_intra=irain_coast_intra
  filter_monsoon, irain_offshore, 'rmean', var_bw=irain_offshore_bw, var_intra=irain_offshore_intra
  filter_monsoon, irain_bang, 'rmean', var_bw=irain_bang_bw, var_intra=irain_bang_intra

;  filter_monsoon, var, 'rmean', var_bw=var_bw, var_intra=var_intra
;  if iera then filter_monsoon, u, 'rmean', var_bw=u_bw, var_intra=u_intra
;  if iera then filter_monsoon, v, 'rmean', var_bw=v_bw, var_intra=v_intra

  iband=1
  if iband eq 1 then begin
    band_tag='bw'
    irain_all=irain_all_bw
    irain_coast=irain_coast_bw
    irain_offshore=irain_offshore_bw
    irain_bang=irain_bang_bw
;    var=var_bw
;    if iera then u=u_bw
;    if iera then v=v_bw
  endif else if iband eq 2 then begin
    band_tag='intra'
    irain_all=irain_all_intra
    irain_coast=irain_coast_intra
    irain_offshore=irain_offshore_intra
    irain_bang=irain_bang_intra
;    var=var_intra
;    if iera then u=u_intra
;    if iera then v=v_intra
  endif

;BAND-PASS FILTER

  ;REMOVE PERIODS >= MIN_TS DAYS (THIS ALSO ACTS TO DETREND)
    min_ts = 121 ; n-days of running average
  ;SUBTRACT RUNNING MEAN
    var -= smooth(var,[0,0,min_ts],/edge_truncate)
    if iera then u -= smooth(u,[0,0,min_ts],/edge_truncate)
    if iera then v -= smooth(v,[0,0,min_ts],/edge_truncate)

;    irain_all -= smooth(irain_all,min_ts,/edge_truncate)
;    irain_coast -= smooth(irain_coast,min_ts,/edge_truncate)
;    irain_offshore -= smooth(irain_offshore,min_ts,/edge_truncate)
;    irain_bang -= smooth(irain_bang,min_ts,/edge_truncate)

  ;KEEP BAND-PASS
    min_ts = 5 ; n-days of running average
    var = smooth(temporary(var),[0,0,min_ts],/edge_truncate)
    if iera then u = smooth(temporary(u),[0,0,min_ts],/edge_truncate)
    if iera then v = smooth(temporary(v),[0,0,min_ts],/edge_truncate)

;    irain_all = smooth(temporary(irain_all),min_ts,/edge_truncate)
;    irain_coast = smooth(temporary(irain_coast),min_ts,/edge_truncate)
;    irain_offshore = smooth(temporary(irain_offshore),min_ts,/edge_truncate)
;    irain_bang = smooth(temporary(irain_bang),min_ts,/edge_truncate)

;----KEEP JJAS ONLY--------------------

if jjas_only then begin

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

endif

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

ido_reg=1
if ido_reg then begin

;----LAG RAINFALL INDICES FOR LAG-REGRESSIONS

  lagall=0

  if dolag then begin
    lagall=indgen(4)*5+5
    lagall=[reverse(-1*lagall),0,lagall]
  endif

  nlag=n_elements(lagall)

  irain_all_sav=irain_all
  irain_coast_sav=irain_coast
  irain_offshore_sav=irain_offshore
  irain_bang_sav=irain_bang

for ilag=0,nlag-1 do begin

  lag=lagall[ilag]

  irain_all = shift(irain_all_sav,-1*lag)
  irain_coast = shift(irain_coast_sav,-1*lag)
  irain_offshore = shift(irain_offshore_sav,-1*lag)
  irain_bang = shift(irain_bang_sav,-1*lag)

  ;REGRESS MAIN VARIABLE

    dims=[nx,ny,nd]

    var_all      = regress_3d(irain_all,var,dims)
    var_coast    = regress_3d(irain_coast,var,dims)
    var_offshore = regress_3d(irain_offshore,var,dims)
    var_bang     = regress_3d(irain_bang,var,dims)

  ;REGRESS ERA

  if iera then begin

    dims=[nxera,nyera,nd]

    u_all      = regress_3d(irain_all,u,dims)
    u_coast    = regress_3d(irain_coast,u,dims)
    u_offshore = regress_3d(irain_offshore,u,dims)
    u_bang     = regress_3d(irain_bang,u,dims)

    v_all      = regress_3d(irain_all,v,dims)
    v_coast    = regress_3d(irain_coast,v,dims)
    v_offshore = regress_3d(irain_offshore,v,dims)
    v_bang     = regress_3d(irain_bang,v,dims)

  endif

  ;BACK INTO UNITS OF Y

;  var_all*=std_all
;
;  var_coast*=std_coast
;  var_offshore*=std_offshore
;;  var_onshore*=std_onshore
;  var_bang*=std_bang
;;  var_smd*=std_smd
;;  var_cdiff*=std_cdiff
;
;  if iera then begin
;    u_all*=std_all
;    v_all*=std_all
;    u_coast*=std_coast
;    v_coast*=std_coast
;    u_offshore*=std_offshore
;    v_offshore*=std_offshore
;;    u_onshore*=std_onshore
;;    v_onshore*=std_onshore
;    u_bang*=std_bang
;    v_bang*=std_bang
;;    u_smd*=std_smd
;;    v_smd*=std_smd
;;    u_cdiff*=std_cdiff
;;    v_cdiff*=std_cdiff
;  endif else begin
;    u_all=0 & v_all=0
;    u_coast=0 & v_coast=0
;    u_offshore=0 & v_offshore=0
;;    u_onshore=0 & u_offshore=0
;    u_bang=0 & v_bang=0
;;    u_smd=0 & v_smd=0
;;    u_cdiff=0 & v_cdoff=0
;  endelse


;----REGRESSION PLOTS--------------------

;iplot_reg=0
;if iplot_reg then begin

  if var_plot eq 'rain' then begin
    var_str='RAINNC'
;    setmax=100 & setmin='0.'
;;    setmax=25 & setmin='0.'
;    cbform='(i3)'
;    cbtag='Rain [ mm d!U-1!N ]'
  endif else if var_plot eq 'pw' then begin
    var_str=var_plot
;    setmax=15 & setmin=-15
;    cbform='(i3)'
;    cbtag='PW [ mm ]'
  endif

;  setmax=max(var_all)
  setmax=0.5
  setmin=-1.*setmax
  cbform='(f4.1)'

  myan_figspecs, 'regress', figspecs, setmin=setmin, setmax=setmax, setndivs=4, set_cint=3
  figdir=dirs.figdir+'myanmar/imerg/regression/'
  figspecs=create_struct(figspecs,'figname',' ')
  figspecs.cbar_format=cbform
;  figspecs.cbar_tag=cbtag
;  if var_plot eq 'rain' then figspecs.ndivs+=1

  figspecs.cbar_tag=' '
  if var_plot eq 'rain' then begin
    figspecs.col_table=71
    figspecs.colors=reverse(figspecs.colors)
  endif

for itest=0,3 do begin

  if itest eq 0 then begin
    var_plt=var_all
    u_plt=u_all
    v_plt=v_all
    itag='All'
  endif else if itest eq 1 then begin
    var_plt=var_coast
    u_plt=u_coast
    v_plt=v_coast
    itag='Coast'
  endif else if itest eq 2 then begin
    var_plt=var_offshore
    u_plt=u_offshore
    v_plt=v_offshore
    itag='Offshore'
;  endif else if itest eq 3 then begin
;    if var_plot eq 'rain' then $
;      var_plt=rain_onshore $
;    else if var_plot eq 'pw' then $
;      var_plt=pw_onshore
;    u_plt=u_onshore
;    v_plt=v_onshore
;    itag='Onshore'
  endif else if itest eq 3 then begin
    var_plt=var_bang
    u_plt=u_bang
    v_plt=v_bang
    itag='Bangladesh'
;  endif else if itest eq 4 then begin
;    if var_plot eq 'rain' then $
;      var_plt=rain_smd $
;    else if var_plot eq 'pw' then $
;      var_plt=pw_smd
;    u_plt=u_smd
;    v_plt=v_smd
;    itag='SMD'
;  endif else if itest eq 5 then begin
;    if var_plot eq 'rain' then $
;      var_plt=rain_cdiff $
;    else if var_plot eq 'pw' then $
;      var_plt=pw_cdiff
;    u_plt=u_cdiff
;    v_plt=v_cdiff
;    itag='cdiff'
  endif

  print,'ITEST: ',itag

  if iera then begin
    wind=create_struct('u',u_plt,'v',v_plt,'x',eralon,'y',eralat)
    wspd=sqrt(u_plt^2+v_plt^2);*10.
    cvar=create_struct('cvar',wspd,'x',eralon,'y',eralat)
  endif

  figspecs.title='Regression: "'+itag+'"'
  figspecs.figname=figdir+'imerg_reg_'+var_plot+'_'+strlowcase(itag)+'_'+strtrim(psel_era,2)

  if dolag then begin
    lagtag=strtrim(lag,2)
    if lag lt 0 then lagtag='m'+strtrim(-1*lag,2)
    figspecs.figname+='_lag'+lagtag
  endif

  wrf_myanmar_map_plot, dirs, var_plt, lon, lat, figspecs, wind=wind, cvar=cvar, /noscalewind

endfor ; plots
;endif

endfor ; ilag

endif ; All regression analysis (ido_reg)


;----THRESHOLD AVERAGES--------------------

ido_thresh=0
if ido_thresh then begin

  thresh=1.5

  it_all=where(irain_all ge thresh,count)
  print,'N-all: ',count
  it_coast=where(irain_coast ge thresh,count)
  print,'N-coast: ',count
  it_offshore=where(irain_offshore ge thresh,count)
  print,'N-offshore: ',count
  it_onshore=where(irain_onshore ge thresh,count)
  print,'N-onshore: ',count

  it_all_m=where(irain_all le -1.*thresh,count)
  print,'N-all-m: ',count
  it_coast_m=where(irain_coast le -1.*thresh,count)
  print,'N-coast-m: ',count
  it_offshore_m=where(irain_offshore le -1.*thresh,count)
  print,'N-offshore-m: ',count
  it_onshore_m=where(irain_onshore le -1.*thresh,count)
  print,'N-onshore-m: ',count

  if var_plot eq 'rain' then begin
    var_str='RAINNC'
    pvar=rain
    setmax=60 & setmin='0.'
    cbform='(i3)'
    lon=lonim
    lat=latim
  endif else if var_plot eq 'pw' then begin
    var_str=var_plot
    pvar=pw
    setmax=70 & setmin=30
    cbform='(i2)'
    lon=eralon
    lat=eralat
  endif

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, setndivs=6, set_cint=3
  figdir=dirs.figdir+'myanmar/imerg/'
  figspecs=create_struct(figspecs,'figname',' ')
  figspecs.cbar_format=cbform

for itest=0,7 do begin

  if itest eq 0 then begin
    it_sel=it_all
    itag='All'
  endif else if itest eq 1 then begin
    it_sel=it_coast
    itag='Coast'
  endif else if itest eq 2 then begin
    it_sel=it_offshore
    itag='Offshore'
  endif else if itest eq 3 then begin
    it_sel=it_onshore
    itag='Onshore'
  endif else if itest eq 4 then begin
    it_sel=it_all_m
    itag='All-m'
  endif else if itest eq 5 then begin
    it_sel=it_coast_m
    itag='Coast-m'
  endif else if itest eq 6 then begin
    it_sel=it_offshore_m
    itag='Offshore-m'
  endif else if itest eq 7 then begin
    it_sel=it_onshore_m
    itag='Onshore-m'
  endif

  if n_elements(it_sel) le 1 then continue

  var_plt=mean(pvar[*,*,it_sel],dimension=3,/nan,/double)
  u_plt=mean(u[*,*,it_sel],dimension=3,/nan,/double)
  v_plt=mean(v[*,*,it_sel],dimension=3,/nan,/double)

  wind=create_struct('u',u_plt,'v',v_plt,'x',eralon,'y',eralat)
  wspd=sqrt(u_plt^2+v_plt^2)
  cvar=create_struct('cvar',wspd,'x',eralon,'y',eralat)

  figspecs.title='Average: "'+itag+'"'
  figspecs.figname=figdir+'imerg_thresh_'+var_plot+'_'+strlowcase(itag)+'_'+strtrim(psel_era,2)

  wrf_myanmar_map_plot, dirs, var_plt, lon, lat, figspecs, wind=wind, cvar=cvar

endfor
  

endif ; ido_thresh

print,'DONE!!'
end
