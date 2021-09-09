; 
; Program to calculate the diurnal composite from long-term IMERG data.
; 
; James Ruppert
; 8/10/2021
;
pro calc_imerg_dcomp

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

coast_thresh=200 ; km

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

  ;IMERG V06B, JJAS, 2000-2020
  daily_im_fil=dirs.wkdir+'imerg/imerg/easthem/3B-HHR.MS.MRG.3IMERG.2000-2020_JJAS_dayavg.V06B.nc4'
  outdir=dirs.wkdir+'imerg/imerg/easthem/'
  datdir=dirs.scdir+'imerg/'
  spawn,'ls '+datdir+'*JJAS_[0-9]*nc4',imfils

  npd_im=48

  outdir=dirs.wkdir+'imerg/imerg/easthem/'

;----TIME ARRAY FOR DAILY AVERAGE TIME SERIES--------------------

  ;SELECTED TIME ARRAY
    time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
      final=julday(mm_plot[1],dd_plot[1],yy_plot[1],23,59,59),step_size=1,units='Days');30,units='Minutes')
    nt=n_elements(time)
    nd=nt;/npd_imerg
    nyr=yy_plot[1]-yy_plot[0]+1

;=====BEGIN READING=========================================================

;----READ RAIN--------------------

  rain_daily=read_nc_var(daily_im_fil,'precipitationCal') ; [ y, x, nd ]
  lon=read_nc_var(daily_im_fil,'lon')
  lat=read_nc_var(daily_im_fil,'lat')
  nx=n_elements(lon)
  ny=n_elements(lat)

  rain_mean=mean(rain_daily,dimension=3,/nan,/double)
  rain_daily=0
  rain_daily=rebin(rain_mean,[ny,nx,nd])

  rain_dc=fltarr(ny,nx,npd_im)

for it=0,npd_im-1 do begin

  hrstr=string(it,format='(i2.2)')
  print,'Hour: ',hrstr

  imfil=imfils[it]

  if it eq 0 then begin
    lon=read_nc_var(imfil,'lon')
    lat=read_nc_var(imfil,'lat')
  endif

  rain_sav=read_nc_var(imfil,'precipitationCal') ; [ y, x, nd ]

  rain_anom = rain_sav-rain_daily

  rain_dc[*,*,it]=mean(rain_anom,dimension=3,/double,/nan)

endfor

  file=outdir+'idl_detrend_3B-HHR.MS.MRG.3IMERG.2000-2020_JJAS_diurncomp.V06B.nc4'
  write_sing_ncvar,file,rain_dc,'rain_dc',dim2=lon,dim1=lat,dimtag2='lon',dimtag1='lat',dim3=npd_im,dimtag3='npd'

end

