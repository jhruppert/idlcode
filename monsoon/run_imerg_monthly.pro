; 
; Create large-scale maps of monthly averaged IMERG and ERA5 data.
;
; James Ruppert
; 8/3/21
; 
pro run_imerg_monthly

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

;iseason=1 ; 0 - all-year, 1 - MJJAS , 2 - ONDJFM
          ; Will use only complete samples for each year.

psel_era=850

;SELECT DATE RANGE
  yy_plot=[2000,2020]
  mm_plot=[6,12]

;----OB DIRECTORIES--------------------

;  im_fil=dirs.wkdir+'imerg/imerg/data/3B-MO.MS.MRG.3IMERG.2000-2020_monthly.V06B.HDF5.nc4'
  im_fil=dirs.wkdir+'imerg/imerg/easthem/3B-HHR.MS.MRG.3IMERG.2000-2020_JJAS_allmean.V06B.nc4'
  era_fil_p=dirs.wkdir+'era5/ERA5-20000101-20201231-pl_monthly.nc4'
    ; u,v [ m/s ]
    ; RH [ % ]
  era_fil_s=dirs.wkdir+'era5/ERA5-20000101-20201231-sl_monthly.nc4'
    ; SST [ K ]
    ; PW [ mm ]

  ifigdir=dirs.figdir+'myanmar/monthly/'

;----TIME ARRAY FOR DAILY AVERAGE TIME SERIES--------------------

  ;SELECTED TIME ARRAY
    time=timegen(start=julday(mm_plot[0],1,yy_plot[0],0,0,0),$
      final=julday(mm_plot[1],1,yy_plot[1],23,59,59),step_size=1,units='Months')
    nt=n_elements(time)
    nyr=yy_plot[1]-yy_plot[0]+1

;=====BEGIN READING=========================================================

;----READ RAIN--------------------

  rain=read_nc_imerg(time,im_fil,lon=lonim,lat=latim,bounds=bounds)*24 ; mm/h --> mm/d
  nx=n_elements(lonim)
  ny=n_elements(latim)

;----READ ERA--------------------

  u=read_nc_era5(time,era_fil_p,'var131',plev=psel_era,lon=eralon,lat=eralat,bounds=bounds)
  v=read_nc_era5(time,era_fil_p,'var132',plev=psel_era,lon=eralon,lat=eralat,bounds=bounds)
  nxera=n_elements(eralon)
  nyera=n_elements(eralat)

  ;SST
    cvar=read_nc_era5(time,era_fil_s,'var34',lon=eralon,lat=eralat,bounds=bounds)
  ;PW
;    cvar=read_nc_era5(time,era_fil_s,'var137',lon=eralon,lat=eralat,bounds=bounds)

;=====TIME AVERAGE=========================================================

caldat,time,mm,dd,yy
iy=yy[uniq(yy)]

  ;ALL-SEASON

  it_avg=where(yy ge iy[1])

  rain_all=mean(rain[*,*,it_avg],dimension=3,/nan,/double)
  u_all=mean(u[*,*,it_avg],dimension=3,/nan,/double)
  v_all=mean(v[*,*,it_avg],dimension=3,/nan,/double)

  cvar_all=mean(cvar[*,*,it_avg],dimension=3,/nan,/double)

  ;JJAS

;  it_avg=where((yy ge iy[1]) and (mm ge 6 and mm le 9))
  it_avg=where((mm ge 6 and mm le 9))

  rain_sum=rain;mean(rain[*,*,it_avg],dimension=3,/nan,/double)
  u_sum=mean(u[*,*,it_avg],dimension=3,/nan,/double)
  v_sum=mean(v[*,*,it_avg],dimension=3,/nan,/double)

  cvar_sum=mean(cvar[*,*,it_avg],dimension=3,/nan,/double)

;endif else if iseason eq 2 then begin
;  it_avg=where((yy ge iy[1]) and ((mm ge 10) or (mm le 3)))
;  timtag='_ondjfm'
;endif


;  ;DATE STRING
;    form='(i4)'
;    dat_str=' ('+string(yy_plot[0],format=form)+'-'+string(yy_plot[1],format=form)+')'
;  if do_jjas then $
;    dat_str=' ('+string(yy_plot[0],format=form)+'-'+string(yy_plot[1],format=form)+', JJAS)'

;=====CALCULATE VORTICITY=========================================================

  avor=abs_vorticity(u_sum,v_sum,eralon,eralat)

;=====PLOTTING=========================================================

  var_str='RAINNC'
  setmax=18 & setmin='0'
  setmax=24

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, setndivs=4, set_cint=3
  figspecs=create_struct(figspecs,'figname',' ')
  figspecs.cbar_format='(i2)'
  figspecs.cbar_tag='[ mm/d ]'
  figspecs.ndivs+=2

  figspecs.title='';'Mean Rainfall';+dat_str
  figspecs.figname=ifigdir+'imerg_era5_'+strtrim(psel_era,2)
;  figspecs.figname+=timtag

  if keyword_set(cvar) then begin
    cvar=create_struct('cvar',cvar_sum,'x',eralon,'y',eralat)
    cint=0.5
    figspecs.clevs=(findgen(50)+1)*cint + 290
  endif

;  wrf_monthly_myanmar_plot, dirs, rain_all, rain_sum, u_all, v_all, u_sum, v_sum, $
;    lonim, latim, eralon, eralat, figspecs, cvar=cvar

  ;ZOOMED IN

    bounds=[50,-20,160,30]
    bounds=[66,2,108,25] ; SASM region
    ;SUBSET TO BOUNDING BOX
      ix=where((lonim ge bounds[0]) and (lonim le bounds[2]))
      iy=where((latim ge bounds[1]) and (latim le bounds[3]))
      rain=rain_sum[ix,*] & rain=rain[*,iy]
;      rain=rain_all[ix,*] & rain=rain[*,iy]
      lon=lonim[ix] & lat=latim[iy]

    wind=create_struct('u',u_sum,'v',v_sum,'x',eralon,'y',eralat)
;    wind=create_struct('u',u_all,'v',v_all,'x',eralon,'y',eralat)

    if keyword_set(cvar) then begin
      cvar=create_struct('cvar',cvar_sum,'x',eralon,'y',eralat)
;      cvar=create_struct('cvar',cvar_all,'x',eralon,'y',eralat)
      cint=1.;0.5
      figspecs.clevs=(findgen(50)+1)*cint + 290
    endif

    figspecs.figname=ifigdir+'imerg_era5_'+strtrim(psel_era,2)+'_zoomin'
    wrf_myanmar_map_plot, dirs, rain, lon, lat, figspecs, wind=wind;, cvar=cvar;, /noscalewind

;stop

print,'DONE!!'
end
