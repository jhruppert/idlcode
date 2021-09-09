; 
; Create histograms of rainfall as a function of filtered coast-normal wind speed.
;
; James Ruppert
; 7/22/21
; 
pro run_imerg_coast_windhist

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

coast_thresh=150 ; km

min_ts=5 ; n-days of running average

;ERAi SETTINGS
;LEVEL SELECTION
  psel_era=850;500;700;925

;SELECT DATE RANGE
  yy_plot=[2000,2020]
  mm_plot=[6,12]
  dd_plot=[1,31] ; inclusive
;yy_plot[0]=2001
;mm_plot[0]=1
dat_str='2000-2020'
;   yy_plot=[2013,2017]
;   mm_plot=[1,12]
;   dd_plot=[1,31] ; inclusive
;dat_str='2013-2017'

;----DIRECTORIES--------------------

  figdir=dirs.figdir+'myanmar/imerg/wind_histogram/'
  ifigdir=figdir+'ind_unorm/'+dat_str+'/'

  maindir=dirs.wkdir
;  im_fil=dirs.wkdir+'imerg/imerg/data/imerg_3B-HHR.MS.MRG.3IMERG.V06_daily_2000-2021_latlonsubset.nc4'
  im_fil=dirs.wkdir+'imerg/imerg/easthem/3B-HHR.MS.MRG.3IMERG.2000-2020_JJAS_dayavg.V06B.nc4'
;  npd_imerg=48

  era_dir=maindir+'era5/'
  era_fil=era_dir+'ERA5-20000101-20201231-pl_dayavg.nc'
  era_pw=era_dir+'ERA5-20000101-20201231-pw_dayavg.nc'
;  npd_era=24

;Local SOLAR time conversion
;local=6;round(mean(dims.lon)/360.*24.) ; deg lon --> hours
;print,'Adding +'+strtrim(local,2)+' for LT'

;----TIME ARRAY FOR DAILY AVERAGE TIME SERIES--------------------

  ;SELECTED TIME ARRAY
    time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
      final=julday(mm_plot[1],dd_plot[1],yy_plot[1],23,59,59),step_size=1,units='Days');30,units='Minutes')
    nt=n_elements(time)
    nd=nt;/npd_imerg
    nyr=yy_plot[1]-yy_plot[0]+1

  ;SAVE JJAS INDICES
    caldat,time,mm,dd,yy
    jjas=where((mm ge 6) and (mm le 9))

;=====BEGIN READING=========================================================

;WESTERN GHATS

  bounds=[69.5,8.,77.5,20.] ; Western Ghats

;----READ RAIN--------------------

  rain_wg=read_nc_imerg(time,im_fil,lon=lonwg,lat=latwg,bounds=bounds)*24. ; mm/h --> mm/d

  coast_dist_wg = read_coastal_dist(lonwg, latwg)
  lsmask_wg = read_lsmask(lonwg, latwg)

;----READ ERA--------------------

  u_wg=read_nc_era5(time,era_fil,'var131',plev=psel_era,lon=eralonwg,lat=eralatwg,bounds=bounds)
  v_wg=read_nc_era5(time,era_fil,'var132',plev=psel_era,bounds=bounds)

;=====BOB================

  bounds=[90.,15.5,95.5,22.5] ; Northern Myanmar coastline

;----READ RAIN--------------------

  rain_bob=read_nc_imerg(time,im_fil,lon=lonbob,lat=latbob,bounds=bounds)*24. ; mm/h --> mm/d

  coast_dist_bob = read_coastal_dist(lonbob, latbob)
  lsmask_bob = read_lsmask(lonbob, latbob)

;----READ ERA--------------------

  u_bob=read_nc_era5(time,era_fil,'var131',plev=psel_era,lon=eralonbob,lat=eralatbob,bounds=bounds)
  v_bob=read_nc_era5(time,era_fil,'var132',plev=psel_era,bounds=bounds)


;----KEEP JJAS ONLY--------------------

  nd=n_elements(jjas)

  ;RAIN IS ALREADY JJAS
;  rain_wg=rain_wg[*,*,jjas]
;  rain_bob=rain_bob[*,*,jjas]

  u_wg=u_wg[*,*,jjas]
  v_wg=v_wg[*,*,jjas]
  u_bob=u_bob[*,*,jjas]
  v_bob=v_bob[*,*,jjas]


;=====AVERAGE RAINFALL AROUND COASTLINE=========================================================

  ;REBIN VARIABLES FOR WHERE STATEMENTS
    nx=n_elements(lonwg)
    ny=n_elements(latwg)
    coast_dist_wg=rebin(coast_dist_wg,nx,ny,nd)
    lsmask_wg=rebin(lsmask_wg,nx,ny,nd)

    nx=n_elements(lonbob)
    ny=n_elements(latbob)
    coast_dist_bob=rebin(coast_dist_bob,nx,ny,nd)
    lsmask_bob=rebin(lsmask_bob,nx,ny,nd)

  ;IMPOSE RULES
;    isea=where((coast_dist_wg le coast_thresh) and (lsmask_wg le 0.1),complement=inot)
    isea=where((coast_dist_wg le coast_thresh),complement=inot)
    rain_wg[inot]=!values.f_nan
;    isea=where((coast_dist_bob le coast_thresh) and (lsmask_bob le 0.1),complement=inot)
    isea=where((coast_dist_bob le coast_thresh),complement=inot)
    rain_bob[inot]=!values.f_nan

  rain_wg = mean(mean(rain_wg,dimension=1,/nan,/double),dimension=1,/nan,/double)
  rain_bob = mean(mean(rain_bob,dimension=1,/nan,/double),dimension=1,/nan,/double)


;=====COAST-NORMAL WIND INDEX=========================================================

  icross=2
  coastnormal, icross, coast_thresh, u_wg, v_wg, eralonwg, eralatwg, uindex=uindex_wg
  icross=1
  coastnormal, icross, coast_thresh, u_bob, v_bob, eralonbob, eralatbob, uindex=uindex_bob

;----FILTER--------------------

  ;SMOOTH TO REMOVE HF VARIANCE
    uindex_wg = smooth(temporary(uindex_wg),min_ts,/edge_truncate)
    uindex_bob = smooth(temporary(uindex_bob),min_ts,/edge_truncate)
    rain_wg = smooth(temporary(rain_wg),min_ts,/edge_truncate)
    rain_bob = smooth(temporary(rain_bob),min_ts,/edge_truncate)

  ;UNORM INDEX
;  filter_monsoon, uindex, 'fft', var_bw=u_bw, var_intra=u_intra
;  filter_monsoon, uindex, 'rmean', var_bw=u_bw, var_intra=u_intra


;=====GENERATE HISTOGRAM=========================================================

  bins = [0,2,4,6,8,10]
  nbin=n_elements(bins)

for idom=0,1 do begin

  if idom eq 0 then begin
    figname=ifigdir+'rain_histogram_wg'
    rain=rain_wg
    uindex=uindex_wg
  endif else if idom eq 1 then begin
    figname=ifigdir+'rain_histogram_bob'
    rain=rain_bob
    uindex=uindex_bob
  endif

  ;PLOT SPECS
    csize=0.65
    xsize=2.8 & ysize=1.8
    position = [0.17,0.19,0.87,0.87]

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,39,/silent

  x=indgen(nbin)
  y=[0,1]
  xtitle='[ m/s ]'
  ytitle='[ mm/d ]'
  xrange=[0,nbin]
  yrange=[0,50]
  if idom eq 1 then yrange[1]=90

  xtickname=['0-2','2-4','4-6','6-8','8-10','>10']
  xtickv=findgen(nbin)+0.5

  plot,x,y,/nodata,position=position,charsize=csize,$
    xstyle=9,ystyle=9,xminor=1,$
    xrange=xrange,yrange=yrange,$
    xtitle=xtitle,xticklen=-0.01,yminor=1,ytitle=ytitle,$
    xticks=nbin-1,xtickv=xtickv,xtickname=xtickname

  plots,!x.crange,[0,0]

  for ib=0,nbin-1 do begin

    ;BIN DATA

      if ib lt nbin-1 then $
        it=where(uindex ge bins[ib] and uindex lt bins[ib+1]) $
      else $
        it=where(uindex ge bins[ib])

      irain=rain[it]
      sorted = sort(irain)
      irain=irain[sorted]

      ;NANs
      nans = where(~finite(irain),nnan,complement=good)
      irain=irain[good]
      npts = n_elements(irain)
  
      ind=indgen(npts)
      ind_med = npts/2
      vmed = irain[ind_med]
  
      lowerGroup = ind[0:ind_med]
      higherGroup = ind[ind_med:npts-1]

      q25 = irain[lowerGroup[n_elements(lowerGroup)/2]]
      q75 = irain[higherGroup[n_elements(higherGroup)/2]]
  
      vmax=max(irain[higherGroup])
      vmin=min(irain[lowerGroup])

    ;PLOT IT UP

      thick=1.3
      width=0.18
      ix=ib+0.5
  
      ;MEDIAN
      plots,[ix-width,ix+width],[vmed,vmed],thick=thick,linestyle=0,color=230
  
      ;BOX
      col=60
      plots,[ix-width,ix+width],[q25,q25],thick=thick,linestyle=0,color=col
      plots,[ix-width,ix+width],[q75,q75],thick=thick,linestyle=0,color=col
      plots,[ix-width,ix-width],[q25,q75],thick=thick,linestyle=0,color=col
      plots,[ix+width,ix+width],[q25,q75],thick=thick,linestyle=0,color=col
  
      ;WHISKERS
      col=0
      plots,[ix,ix],[vmin,q25],thick=thick,linestyle=1,color=col
      plots,[ix,ix],[q75,vmax],thick=thick,linestyle=1,color=col
      width*=0.6
      plots,[ix-width,ix+width],[vmax,vmax],thick=thick,linestyle=0,color=col
      plots,[ix-width,ix+width],[vmin,vmin],thick=thick,linestyle=0,color=col

  endfor

  device,/close
  convert_png,figname,res=200;,/remove_eps


endfor ; idom

print,'Done!!'
end
