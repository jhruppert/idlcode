; 
; Regression analysis of diurnally varying rainfall (Myanmar/Meiyu TRMM data) and ERA5
; using area-averaged rainfall time series calculated in run_trmm_coast_index.pro
;
; James Ruppert
; 2/5/21
; 
pro run_trmm_coast_dc_regress

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

var_plot='rain';
;var_plot='pw';'rain';

do_wind=0 ; include ERA5 wind?

;ERAi SETTINGS
;LEVEL SELECTION
  psel_era=925;500;700;925

;USE FULL JJAS 2013-2017
  yy_plot=[2013,2017]
  mm_plot=[6,9]
  dd_plot=[1,30] ; inclusive

;GIVEN THE ABOVE, THIS IS ALL WE NEED FOR TIME
  npd_imerg=8
  npd_era=24
  nd=610
  nt=nd*npd_imerg
  nt_era=nd*npd_era

  ;DATE STRING
    form2='(i2.2)'
    form4='(i4)'
    dat_str=string(mm_plot[0],format=form2)+string(dd_plot[0],format=form2)+strmid(strtrim(yy_plot[0],2),2,2)+'-'+$
            string(mm_plot[1],format=form2)+string(dd_plot[1],format=form2)+strmid(strtrim(yy_plot[1],2),2,2)

;----WRF STUFF FOR LAND/SEA MASK------

  expname='myanmar'
  cases=['ctl']
  idir='2dom'
  if idir eq '2dom' then $
    domtag='d01'
  casedir=dirs.scdir+expname+'/WRF/experiment/'+idir+'/'
  dirs.figdir+=expname+'/'
  config_wrfexp, casedir=casedir,cases=cases,dirs=dirs,$
    dims=dims, vars=vars, nfils_set=nfils_set, domtag=domtag;, /verbose
  
  ;READ TOPOGRAPHY AND LAND MASK FROM ERA5 MET_EM FILES
    met_dir='/work/06040/tg853394/stampede2/wrf_wps/wrfv4/output_myanmar/JJAS_2013-2017/era5/'
    spawn,'ls '+met_dir+'met_em.'+domtag+'*',met_fil
    met_fil=met_fil[0]
    topo=reform(read_nc_var(met_fil,'HGT_M'))*1e3 ; m --> km
    lsmask=reform(read_nc_var(met_fil,'LANDMASK'))


;----DATA DIRECTORIES--------------------

  maindir=dirs.wkdir+'trmm/'
  im_fil_dm=maindir+'TRMM_3B42v7.2013-2017_JJAS_daymean.nc4'

  era_dir=maindir+'era5/jjas_2013-2017/';ERAi-JJAS13-17-pl.nc4'
;  era_fil=era_dir+'ERAi-JJAS13-17-pl_dayavg.nc4'
;  era_sfil=era_dir+'ERAi-JJAS13-17-sl_dayavg.nc4'

  irain_dir=dirs.figdir+'imerg/'

;Local SOLAR time conversion
  local=6;round(mean(dims.lon)/360.*24.) ; deg lon --> hours
  print,'Adding +'+strtrim(local,2)+' for LT'
  ltim_imerg=findgen(npd_imerg)*24/npd_imerg+local
  for it=0,npd_imerg-1 do ltim_imerg[it]-=24*(ltim_imerg[it] ge 24)
  ltim_era=findgen(npd_era)*24/npd_era+local
  for it=0,npd_era-1 do ltim_era[it]-=24*(ltim_era[it] ge 24)


;----TIME ARRAYS--------------------

  ;IMERG TIME
;    if yy_plot[0] eq yy_plot[1] then $
;      time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
;        final=julday(mm_plot[1],dd_plot[1],yy_plot[1],23,59,59),step_size=30,units='minutes') $
;    else begin
;      time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
;        final=julday(mm_plot[1],dd_plot[1],yy_plot[0],23,59,59),step_size=30,units='minutes')
;      for iyy=yy_plot[0]+1,yy_plot[1] do begin
;        itime=timegen(start=julday(mm_plot[0],dd_plot[0],iyy,0,0,0),$
;          final=julday(mm_plot[1],dd_plot[1],iyy,23,59,59),step_size=30,units='minutes')
;        time=[time,itime]
;      endfor
;    endelse
;    nt=n_elements(time)
;    nd=nt/npd_imerg
;    nyy=yy_plot[1]-yy_plot[0]+1

;    rain_dm=reform(read_nc_var(im_fil_dm,'precipitationCal'))*24. ; change to mm/d
;    rain_dm=transpose(temporary(rain_dm),[1,0,2])


;----READ IMERG RAINFALL--------------------

  print,'Reading TRMM...'

  spawn,'ls '+maindir+'TRMM_3B42v7.2013-2017_JJAS.nc4',trmmfil

  ;LAT/LON
    lonim=read_nc_var(trmmfil,'nlon')
    latim=read_nc_var(trmmfil,'nlat')
    nx=n_elements(lonim)
    ny=n_elements(latim)

    ;READ RAIN AND PUT INTO F(X,Y,D,HOUR)
    tmp=reform(read_nc_var(trmmfil,'precipitation'))
    tmp=transpose(temporary(tmp),[1,0,2])
    rain=fltarr(nx,ny,npd_imerg,nd)
    idarr=indgen(npd_imerg)
    for id=0,nd-1 do $
      rain[*,*,*,id]=reform(tmp[*,*,idarr+id*npd_imerg])
    rain*=24.  ; mm/hr --> mm/d
    rain_dc=mean(rain,/nan,/double,dimension=4) ; all-time diurnal composite


;----READ ERA--------------------

if do_wind then begin

  print,'Beginning read of ERA...'

  spawn,'ls '+era_dir+'ERA5-201[0-9]*-pl.nc',era_fils_pl
  spawn,'ls '+era_dir+'ERA5-201[0-9]*-sl.nc',era_fils_sl
  nfil=n_elements(era_fils_pl)

  ;DIMENSIONS
    era_fil=era_fils_sl[0]
    eralon=read_nc_var(era_fil,'lon')
    eralat=read_nc_var(era_fil,'lat')
    nxera=n_elements(eralon)
    nyera=n_elements(eralat)

  ;PLEVS
    era_fil=era_fils_pl[0]
    p_era=reform(read_nc_var(era_fil,'plev'))*1d-2 ; Pa --> hPa
    nzera=n_elements(p_era)
    ;SELECTED LEVEL
    izlev_era=(where(p_era eq psel_era,count))[0]
    if count eq 0 then stop

  ;READ DATA FROM YEARLY FILES
    u=fltarr(nxera,nyera,nt_era)
    v=u
    for iyy=0,nfil-1 do begin
      print,'  Reading year ',string(iyy+yy_plot[0],format='(i4)')
      it_read=where((yyy_era eq yy_plot[0]+iyy), nt_read)
      ;READ ENTIRE SET OF DAYS
        count=[nxera,nyera,1,nt_read] & offset=[0,0,izlev_era,0] ; x, y, p, t
        u[*,*,it_read]=reform(read_nc_var(era_fils_pl[iyy],'var131',count=count,offset=offset)) ; m/s
        v[*,*,it_read]=reform(read_nc_var(era_fils_pl[iyy],'var132',count=count,offset=offset)) ; m/s
    endfor

endif


;----HOVMOLLER PREP------------------------

  ;REORDER LOCAL TIME TO PUT 0 LT AT FIRST INDEX
  ltim_imerg_shift=shift(ltim_imerg,local*npd_imerg/24)

  if expname eq 'myanmar' then begin
    ;CROSS SECTION
      xcross=[77.8,97.4] ; lon,lat
      ycross=[12.,20.]   ; lon,lat
      width=3.5 ; width of cross section (degrees)
  endif else if expname eq 'meiyu' then begin
  endif

  cross_topo=cross_diag(topo,dims.lon,dims.lat,width,x_bounds=xcross,y_bounds=ycross,lonout=xlon,latout=xlat)
  topo=mean(cross_topo,dimension=2,/nan,/double) & cross_topo=0
  cross_lsmask=cross_diag(lsmask,dims.lon,dims.lat,width,x_bounds=xcross,y_bounds=ycross)
  lsmask=mean(cross_lsmask,dimension=2,/nan,/double) & cross_lsmask=0

  ;INTERPOLATE ONTO HOVMOLLER CROSS
  tmp=cross_diag(reform(rain[*,*,0,0]),lonim,latim,width,x_bounds=xcross,y_bounds=ycross,lonout=xlon_im,latout=xlat_im)

  ;DISTANCE RELATIVE TO COASTLINE (KM)
  ;DISTANCE FROM LAT/LON
    nxhov=n_elements(xlon)
    xhov_era=fltarr(nxhov)
    for ix=0,nxhov-1 do $
      xhov_era[ix]=111.*sqrt( (xlon[ix]*cos(xlat[ix]*!dtor))^2 + xlat[ix]^2 ) ; Converts to km
  ;ADJUST BY COASTLINE LOCATION
    ixcheck=indgen(nxhov/2)+nxhov/2
    ixcoast=(where(lsmask[ixcheck] gt 0))[0]
    xhov_era=reverse(xhov_era)
    xcoast=xhov_era[ixcheck[ixcoast]]
    xhov_era-=xcoast

  ;IMERG: DISTANCE FROM LAT/LON
    nxhov=n_elements(xlon_im)
    xhov_im=fltarr(nxhov)
    for ix=0,nxhov-1 do $
      xhov_im[ix]=111.*sqrt( (xlon_im[ix]*cos(xlat_im[ix]*!dtor))^2 + xlat_im[ix]^2 ) ; Converts to km
    xhov_im=reverse(xhov_im)
    xhov_im-=xcoast


;----MEAN MAPS AND COMPOSITES AS A SANITY CHECK--------------------

;Map

iplot_mn_rain=0
if iplot_mn_rain then begin

if var_plot eq 'rain' then begin
  var_str='RAINNC'
  setmax=30 & setmin='0.'
  cbform='(i2)'
  var_mn=mean(rain_dc,dimension=3,/nan,/double)
  lon=lonim
  lat=latim
endif else if var_plot eq 'pw' then begin
  var_str=var_plot
  setmax=70 & setmin=30
  cbform='(i2)'
  var_mn=mean(pw,dimension=3,/nan,/double)
  lon=eralon
  lat=eralat
endif

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=3
  figdir=dirs.figdir+'trmm/'
  figname=figdir+'trmm_mean_'+var_plot+'_'+dat_str
  figspecs=create_struct(figspecs,'figname',figname)
  figspecs.cbar_format=cbform

;  u_mn=mean(u,dimension=3,/nan,/double)
;  v_mn=mean(v,dimension=3,/nan,/double)
;  wind=create_struct('u',u_mn,'v',v_mn,'x',eralon,'y',eralat)
;  wspd=sqrt(u_mn^2+v_mn^2)
;  cvar=create_struct('cvar',wspd,'x',eralon,'y',eralat)

  wrf_myanmar_map_plot, dirs, var_mn, lon, lat, figspecs, wind=wind, cvar=cvar

endif

;Diurnal composite hovmoller

iplot_dhov=0
if iplot_dhov then begin

  if var_plot eq 'rain' then begin
    var_str='RAINNC'
    setmax=1.5 & setmin='0.'
    cbform='(f3.1)'
    cbtag='Rain [ mm h!U-1!N ]'
    lon=lonim
    lat=latim
  endif else if var_plot eq 'pw' then begin
    var_str=var_plot
    setmax=15 & setmin=-15
    cbform='(i3)'
    cbtag='PW [ mm ]'
    lon=eralon
    lat=eralat
  endif

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, setndivs=3, set_cint=3
  figdir=dirs.figdir+'trmm/'
  figspecs=create_struct(figspecs,'figname',' ')
  figspecs.cbar_format=cbform
  figspecs.cbar_tag=cbtag

  ;INTERPOLATE ONTO HOVMOLLER CROSS
  var_plt=rain_dc
  cross=cross_diag(var_plt,lonim,latim,width,x_bounds=xcross,y_bounds=ycross,lonout=xlon_im,latout=xlat_im)
  imerg=mean(cross,dimension=2,/nan,/double) & cross=0
  imerg*=1./24 ; mm/d --> mm/h

  figspecs.figname=figdir+'../hov/trmm_dcomptest_'+var_plot
  figspecs.title='Diurnal Composite (JJAS 2013-17)'

  ;REORDER TO PUT 0 LT AT FIRST INDEX
  imerg=shift(temporary(imerg),[0,1.*local*npd_imerg/24])

;  if iwind then begin
;    ltim_era=findgen(npd_era)*24/npd_era+local
;    for it=0,npd_era-1 do ltim_era[it]-=24*(ltim_era[it] ge 24)
;    ltim_era=shift(temporary(ltim_era),local*npd_era/24)
;    u=shift(temporary(u),[0,local*npd_era/24])
;    v=shift(temporary(v),[0,local*npd_era/24])
;  endif

;  if iwind then wind=create_struct('u',iu,'v',iv)

  stats,imerg

  wrf_monsoon_dcomp_hov, dirs, figspecs, imerg, ltim_imerg_shift, xhov_im, wind=wind, cvar=cvar

endif


;----RAINFALL INDICES--------------------

  irain_all=fltarr(nd) ; units of mm/d
  irain_coast=irain_all
  irain_offshore=irain_all
  irain_onshore=irain_all
  openr,1,irain_dir+'rain_all_'+dat_str+'.txt'
  readf,1,irain_all & close,1
  openr,1,irain_dir+'rain_coast_'+dat_str+'.txt'
  readf,1,irain_coast & close,1
  openr,1,irain_dir+'rain_offshore_'+dat_str+'.txt'
  readf,1,irain_offshore & close,1
  openr,1,irain_dir+'rain_onshore_'+dat_str+'.txt'
  readf,1,irain_onshore & close,1

  ;STANDARDIZE

  ;remove means
  irain_all -= mean(irain_all,/double,/nan)
  irain_coast -= mean(irain_coast,/nan,/double)
  irain_offshore -= mean(irain_offshore,/nan,/double)
  irain_onshore -= mean(irain_onshore,/nan,/double)

  ;standard deviations
  std_all = sqrt( total(irain_all^2,/nan,/double)/(nd-1) )
  std_coast=sqrt( total(irain_coast^2,/nan,/double)/(nd-1) )
  std_offshore=sqrt( total(irain_offshore^2,/nan,/double)/(nd-1) )
  std_onshore=sqrt( total(irain_onshore^2,/nan,/double)/(nd-1) )

  print,'  Standard deviations:'
  print,'    All:',std_all
  print,'    Coast:',std_coast
  print,'    Offshore:',std_offshore
  print,'    Onshore:',std_onshore

  ;divide by stddev
  irain_all /= std_all
  irain_coast /= std_coast
  irain_offshore /= std_offshore
  irain_onshore /= std_onshore


;----CALCULATE REGRESSIONS--------------------

ido_reg=1
if ido_reg then begin

  ;REGRESS RAIN

  print,'Beginning regression of rain...'
stop
tic
  calc_dc_regress, rain, $
    irain_all, irain_coast, irain_offshore, irain_onshore, $ ; input indices
    a1_all=a1_all, a1_coast=a1_coast, a1_offshore=a1_offshore, a1_onshore=a1_onshore
toc

;PLOT HOVMOLLER OF STDDEV
;stdev = stddev(rain,dimension=4,/nan,/double) ; [x,y,hour]
;;stdevmn=mean(mean(stdev,dimension=1,/nan,/double),dimension=1,/nan,/double)
;;  stdevmn2=shift(stdevmn,local*npd_imerg/24)
;;for it=0,npd_imerg-1 do print,stdevmn2[it]
;;SPATIALLY AVERAGE STDDEV OVER HOVMOLLER AREA
;  xloc=where((lonim ge xcross[0]) and (lonim le xcross[1]))
;  yloc=where((latim ge ycross[0]) and (latim le ycross[1]))
;  stdevmn=mean((mean(stdev[xloc,*,*],dimension=1,/nan,/double))[yloc,*],dimension=1,/nan,/double)
;figdir=dirs.figdir+'imerg/'
;figname=figdir+'stddev_tser'
;plot_stddev_tser, dirs, figname, stdevmn

  ;Stddev of rainfal index
  a1_all*=std_all
  a1_coast*=std_coast
  a1_offshore*=std_offshore
  a1_onshore*=std_onshore

  ;Stddev of point-rainfall as a function of hour
;for it=0,npd_imerg-1 do a1_all[*,*,it]*=stdevmn[it]
;for it=0,npd_imerg-1 do a1_coast[*,*,it]*=stdevmn[it]
;for it=0,npd_imerg-1 do a1_offshore[*,*,it]*=stdevmn[it]
;for it=0,npd_imerg-1 do a1_onshore[*,*,it]*=stdevmn[it]

  stats,a1_all
  stats,a1_coast
  stats,a1_offshore
  stats,a1_onshore

  ;REGRESS WIND

  if do_wind then begin

    print,'Beginning regression of ERA wind...'
  
    calc_dc_regress, u, npd_era, $
      irain_all, irain_coast, irain_offshore, irain_onshore, $ ; input indices
      a1_all=u_all, a1_coast=u_coast, a1_offshore=u_offshore, a1_onshore=u_onshore

    calc_dc_regress, v, npd_era, $
      irain_all, irain_coast, irain_offshore, irain_onshore, $ ; input indices
      a1_all=v_all, a1_coast=v_coast, a1_offshore=v_offshore, a1_onshore=v_onshore

    u_all*=std_all
    v_all*=std_all
    u_coast*=std_coast
    v_coast*=std_coast
    u_offshore*=std_offshore
    v_offshore*=std_offshore
    u_onshore*=std_onshore
    v_onshore*=std_onshore

  endif


;----PLOTTING-------------------------------

;stop

rtag='trmm'

dc_regress_plots, dirs, rtag, dims, var_plot, npd_imerg, local, ltim_imerg, psel_era, $
  xcross, ycross, width, xhov_im, $
  a1_all, a1_coast, a1_offshore, a1_onshore, $
  do_wind, u_all=u_all, u_coast=u_coast, u_offshore=u_offshore, u_onshore=u_onshore, $
  lonim=lonim, latim=latim, $
  lonera=lonera, latera=latera


endif ; All regression analysis (ido_reg)


print,'DONE!!'
end
