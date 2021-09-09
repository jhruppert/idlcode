; 
; Calculate indices of DAILY-AVERAGED rainfall based on coastal distance.
;
; Coastal distance downloaded from https://oceancolor.gsfc.nasa.gov/docs/distfromcoast/dist2coast.txt.bz2
;
; James Ruppert
; 12/8/20
; 
pro run_imerg_coast_index

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

iplot_mn_rain=0

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

  maindir=dirs.scdir+'myanmar/'
  npd_imerg=48
  era_dir=dirs.wkdir+'wrf_wps/wrfv4/output_myanmar/JJAS_2013-2017/era5/'
  coast_dist_fil=dirs.home+'idl/code/misc/dist2coast.nc'
  im_fil=dirs.wkdir+'imerg/imerg/data/imerg_3B-HHR.MS.MRG.3IMERG.V06_daily_2000-2021_latlonsubset.nc4'

  figdir=dirs.figdir+'myanmar/imerg/'
  irain_dir=figdir+'rainfall_indices/'

;----TIME ARRAY FOR DAILY AVERAGE TIME SERIES--------------------

  ;SELECTED TIME ARRAY
    time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
      final=julday(mm_plot[1],dd_plot[1],yy_plot[1],23,59,59),step_size=1,units='Days');30,units='Minutes')
    nt=n_elements(time)
    nd=nt;/npd_imerg
    nyr=yy_plot[1]-yy_plot[0]+1

;----READ RAIN--------------------

  rain_sav=read_nc_imerg(time,im_fil,lon=lonim,lat=latim,bounds=bounds) ; already in mm/d
  nx=n_elements(lonim)
  ny=n_elements(latim)

  ;SAVE BASIC TIME-MEAN
  rain_mn_sav=mean(rain_sav,dimension=3,/nan,/double)

;----PLOT TIME-MEAN RAINFALL--------------------

if iplot_mn_rain then begin

  var_str='RAINNC'
  setmax=10 & setmin='0.'
  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figspecs=create_struct(figspecs,'figname',figdir+'imerg_mean_'+dat_str)
  figspecs.cbar_format='(i2)'
  
  wrf_myanmar_map_plot, dirs, rain_mn_sav, lonim, latim, figspecs

endif

;----READ LAND/SEA MASK--------------------

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

;----READ COASTAL DISTANCE--------------------

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

;----PLOT COASTAL DISTANCE--------------------

iplot_coastal=0
if iplot_coastal then begin

  var_str='RAINNC'
  setmax=10 & setmin='0.'
  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figspecs=create_struct(figspecs,'figname',figdir+'coastal_distance')
  ncols=n_elements(figspecs.colors)
  
  max=800.
  min=0
  figspecs.cbar_format='(i3)'
  figspecs.title='Coastal Distance (km)'
  figspecs.cbar_tag=' '
  figspecs.levels=findgen(ncols)/(ncols-1)*(max-min)+min

  wrf_myanmar_map_plot, dirs, coast_dist_im, lonim, latim, figspecs

endif

;----RAINFALL INDICES--------------------

  ;ONLY JJAS FOR MAPS
    caldat,time,mm,dd,yy
    jjas=where((mm ge 6) and (mm le 9))

  ;BOUNDS FOR SELECTING ONSHORE V OFFSHORE RAINFALL
  ;USE COASTAL DISTANCE VARIABLE
;    bounds=[81.5,9.5,97,23.5]
  ;TAKEN FROM run_imerg_wavelet_coastal.pro
    bounds=[85.,7,100,23.5] ; Large BoB
    coast_thresh=100.

  ;SUBSET TO BOUNDING BOX
    ix=where((lonim ge bounds[0]) and (lonim le bounds[2]),nx)
    iy=where((latim ge bounds[1]) and (latim le bounds[3]),ny)
    rain_dm=rain_sav[ix,*,*] & rain_dm=rain_dm[*,iy,*]
    coast_dist2=coast_dist[ix,*] & coast_dist2=coast_dist2[*,iy]
    lsmask2=lsmask[ix,*] & lsmask2=lsmask2[*,iy]

  ;ADD TIME INDICES TO THESE TO SIMPLIFY WHERE FUNCTIONS
      coast_dist2=rebin(coast_dist2,[nx,ny,nd])
      lsmask2=rebin(lsmask2,[nx,ny,nd])

  ;IMPOSE SELECTION RULES

    ;ALL
    irain_all=mean(mean(rain_dm,dimension=1,/nan,/double),dimension=1,/nan,/double)

    ;COAST
    icoast=where(coast_dist2 le coast_thresh,complement=incoast)
    rain_coast=rain_dm
    rain_coast[incoast]=!values.f_nan
    irain_coast=mean(mean(rain_coast,dimension=1,/nan,/double),dimension=1,/nan,/double)

    ;OFFSHORE
    isea=where((coast_dist2 ge coast_thresh) and (lsmask2 le 0.1),complement=inot)
    rain_offshore=rain_dm
    rain_offshore[inot]=!values.f_nan
    irain_offshore=mean(mean(rain_offshore,dimension=1,/nan,/double),dimension=1,/nan,/double)

    ;ONSHORE
;    ionshore=where((coast_dist2 le coast_thresh) and (lsmask2 ge 0.3),complement=ioffshore)
;    rain_onshore=rain_dm
;    rain_onshore[ioffshore]=!values.f_nan
;    irain_onshore=mean(mean(rain_onshore,dimension=1,/nan,/double),dimension=1,/nan,/double)

    ;BOX OVER BANGLADESH
    bounds_bang=[88.,23,93.5,27] ; Bangladesh
    ixb=where((lonim ge bounds_bang[0]) and (lonim le bounds_bang[2]),nx,complement=nix)
    iyb=where((latim ge bounds_bang[1]) and (latim le bounds_bang[3]),ny,complement=niy)
    rain_bang=rain_sav[ixb,*,*] & rain_bang=rain_bang[*,iyb,*]
    irain_bang=mean(mean(rain_bang,dimension=1,/nan,/double),dimension=1,/nan,/double)

    ;ADAMES AND MING SYNOPTIC MONSOON DISTURBANCES (SMD) BOX
;    bounds_smd=[85.,15,90,20.]
;    ixs=where((lonim[ix] ge bounds_smd[0]) and (lonim[ix] le bounds_smd[2]),complement=nix)
;    iys=where((latim[iy] ge bounds_smd[1]) and (latim[iy] le bounds_smd[3]),complement=niy)
;    rain_smd=rain_dm
;    rain_smd[nix,*,*]=!values.f_nan
;    rain_smd[*,niy,*]=!values.f_nan
;    irain_smd=mean(mean(rain_smd,dimension=1,/nan,/double),dimension=1,/nan,/double)

;----PLOT TIME-MEAN RAINFALL WITH RULES--------------------

iplot_rain_rules=1
if iplot_rain_rules then begin

  ifigdir=figdir+'rainfall_indices/'

  var_str='RAINNC'
  setmax=25 & setmin='0.'
  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figspecs=create_struct(figspecs,'figname',' ')
  figspecs.cbar_format='(i2)'
  figspecs.ndivs-=1

  rain_plt=mean(rain_dm[*,*,jjas],dimension=3,/nan,/double)
  figspecs.title='"All"'
  figspecs.figname=ifigdir+'imerg_mean_all_'+dat_str
  wrf_myanmar_map_plot, dirs, rain_plt, lonim[ix], latim[iy], figspecs

  rain_plt=mean(rain_coast[*,*,jjas],dimension=3,/nan,/double)
  figspecs.title='"Coastal" ('+string(coast_thresh,format='(i3)')+' km)'
  figspecs.figname=ifigdir+'imerg_mean_coast_'+dat_str
  wrf_myanmar_map_plot, dirs, rain_plt, lonim[ix], latim[iy], figspecs

  rain_plt=mean(rain_offshore[*,*,jjas],dimension=3,/nan,/double)
  figspecs.title='"Offshore"'
  figspecs.figname=ifigdir+'imerg_mean_offshore_'+dat_str
  wrf_myanmar_map_plot, dirs, rain_plt, lonim[ix], latim[iy], figspecs

;  rain_plt=mean(rain_onshore[*,*,jjas],dimension=3,/nan,/double)
;  figspecs.title='"Onshore"'
;  figspecs.figname=figdir+'imerg_mean_onshore_'+dat_str
;  wrf_myanmar_map_plot, dirs, rain_plt, lonim[ix], latim[iy], figspecs

  rain_plt=mean(rain_bang[*,*,jjas],dimension=3,/nan,/double)
  figspecs.title='"Bangladesh"'
  figspecs.figname=ifigdir+'imerg_mean_bang_'+dat_str
  wrf_myanmar_map_plot, dirs, rain_plt, lonim[ixb], latim[iyb], figspecs

;  rain_plt=mean(rain_smd[*,*,jjas],dimension=3,/nan,/double)
;  figspecs.title='"SMD"'
;  figspecs.figname=ifigdir+'imerg_mean_smd_'+dat_str
;  wrf_myanmar_map_plot, dirs, rain_plt, lonim[ix], latim[iy], figspecs

endif

;----WRITE RESULTS--------------------

  print,'ALL:'
  stats,irain_all
  print,'COAST:'
  stats,irain_coast
  print,'OFFSHORE:'
  stats,irain_offshore
;  print,'ONSHORE:'
;  stats,irain_onshore
  print,'BANGLADESH:'
  stats,irain_bang
;  print,'SMD:'
;  stats,irain_smd

  ascfil=irain_dir+'rain_all_'+dat_str+'.txt'
  openw,1,ascfil
    printf,1,irain_all
  close,1

  ascfil=irain_dir+'rain_coast_'+dat_str+'.txt'
  openw,1,ascfil
    printf,1,irain_coast
  close,1

  ascfil=irain_dir+'rain_offshore_'+dat_str+'.txt'
  openw,1,ascfil
    printf,1,irain_offshore
  close,1

  ascfil=irain_dir+'rain_bangladesh_'+dat_str+'.txt'
  openw,1,ascfil
    printf,1,irain_bang
  close,1

;  ascfil=irain_dir+'rain_onshore_'+dat_str+'.txt'
;  openw,1,ascfil
;    printf,1,irain_onshore
;  close,1
;
;  ascfil=irain_dir+'rain_smd_'+dat_str+'.txt'
;  openw,1,ascfil
;    printf,1,irain_smd
;  close,1

  ;STANDARDIZED RAIN INDEX
;    rain_dom_mn = (temporary(rain_dom_mn) - mean(rain_dom_mn,/nan,/double)) / stddev(rain_dom_mn)



print,'DONE!!'
end
