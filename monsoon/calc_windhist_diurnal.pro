; 
; Write out histogram information and binned IMERG based on filtered coast-normal wind.
;
; James Ruppert
; 8/10/21
; 
pro calc_windhist_diurnal

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
;  datdir=dirs.wkdir+'imerg/imerg/data/'
  datdir=dirs.scdir+'imerg/'
  spawn,'ls '+datdir+'*JJAS_[0-9]*nc4',imfils

  outdir=dirs.wkdir+'imerg/imerg/easthem/'
  compdir=outdir+'unorm_composites/'

  maindir=dirs.wkdir
  era_dir=maindir+'era5/'
  era_fil=era_dir+'ERA5-20000101-20201231-pl_dayavg.nc'

  npd_im=48

;----TIME ARRAY FOR DAILY AVERAGE TIME SERIES--------------------

  ;SELECTED TIME ARRAY
    time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
      final=julday(mm_plot[1],dd_plot[1],yy_plot[1],23,59,59),step_size=1,units='Days');30,units='Minutes')
    nt=n_elements(time)
    nd=nt
    nyr=yy_plot[1]-yy_plot[0]+1

  ;SAVE JJAS INDICES
    caldat,time,mm,dd,yy
    jjas=where((mm ge 6) and (mm le 9))
    nd=n_elements(jjas)
    npy=nd/nyr

;=====BEGIN READING=========================================================

;----READ RAIN--------------------

  rain_daily=read_nc_var(daily_im_fil,'precipitationCal') ; [ y, x, nd ]
  lon=read_nc_var(daily_im_fil,'lon')
  lat=read_nc_var(daily_im_fil,'lat')
  nx=n_elements(lon)
  ny=n_elements(lat)

  ;SEASONAL MEAN TO DETREND RAW RAINFALL
  rain_annual=fltarr(ny,nx,nyr)
  ityr=indgen(npy)
  for iyr=0,nyr-1 do $
    rain_annual[*,*,iyr]=mean(rain_daily[*,*,ityr+iyr*npy],dimension=3,/nan,/double)

  allmean=mean(rain_daily,dimension=3,/nan,/double)
  rain_annual -= rebin(allmean,[ny,nx,nyr])
  ;RAIN_ANNUAL NOW REPRESENTS THE SEASONAL-MEAN TREND

;=====SET UP VARIABLES=========================================================

  bins = ['-1','-0.5','-0.25','0.25','0.5','1']
  nbin=n_elements(bins)-1
  nband=2
  ncross=2

  rain_ucomp=fltarr(ny,nx,npd_im,nbin,nband,ncross)

;=====CYCLE OVER ALL CROSS SECTIONS AND BINS=========================================================

for it=0,npd_im-1 do begin
;for it=0,0 do begin

  hrstr=string(it,format='(i2.2)')
  print,'Hour: ',hrstr

  imfil=imfils[it]

  rain_sav=read_nc_var(imfil,'precipitationCal') ; mm/hr [ y, x, nd ]

  ;DETREND
  for iyr=0,nyr-1 do $
    rain_sav[*,*,ityr+iyr*npy] -= rebin(reform(rain_annual[*,*,iyr]),[ny,nx,npy])

for icross=1,2 do begin

;----READ ERA--------------------

  if icross eq 1 then bounds_era=[89.,15.5,95.,22.5] ; Northern Myanmar coastline
  if icross eq 2 then bounds_era=[69.5,8.,77.5,21.] ; Western Ghats
  if icross eq 3 then bounds_era=[78.,14.,89.7,24.5] ; Northwestern BoB coastline

  u=read_nc_era5(time,era_fil,'var131',plev=psel_era,lon=eralon,lat=eralat,bounds=bounds_era)
  v=read_nc_era5(time,era_fil,'var132',plev=psel_era,lon=eralon,lat=eralat,bounds=bounds_era)

;=====CREATE INDEX FOR REGRESSION/BINNING=========================================================

  coastnormal, icross, coast_thresh, u, v, eralon, eralat, uindex=uindex

;----FILTER--------------------

  ;UNORM INDEX
  filter_monsoon, uindex, 'fft', var_bw=u_bw, var_intra=u_intra

;=====GENERATE COMPOSITES=========================================================

;----KEEP JJAS ONLY--------------------

  u_bw=u_bw[jjas]
  u_intra=u_intra[jjas]

;----STANDARDIZE--------------------

  ;UNORM INDEX

    std_bw=stddev(u_bw,/nan,/double)
    std_intra=stddev(u_intra,/nan,/double)

if it eq 0 then begin
    print,'Standard Deviations:'
    print,'BW:',std_bw
    print,'Intra:',std_intra
endif

    u_bw    = (u_bw -    mean(u_bw,/nan,/double))    / std_bw
    u_intra = (u_intra - mean(u_intra,/nan,/double)) / std_intra

    ;N-CASES where |irain| > 1 sigma
    ibw = where(abs(u_bw) ge 1,nbw)
    iintra = where(abs(u_intra) ge 1,nintra)

  ;PRINT STATS

if it eq 0 then begin
    print,'N where |index| >= 1:'
    print,'BW:',nbw
    print,'Intra:',nintra
endif

;----BINNING--------------------

  bandtag=['bw','intra']

  for iband=1,2 do begin ; Biweekly / Intraseasonal
  
    if it eq 0 then print,'iband: ',iband
  
    if iband eq 1 then begin
;      bandtag='bw'
      uband=u_bw
      stdd=std_bw
    endif else if iband eq 2 then begin
;      bandtag='intra'
      uband=u_intra
      stdd=std_intra
    endif
  
    for ibin=0,nbin-1 do begin
  
      bintag=strtrim(ibin+1,2)

      bin_txt = [ bins[ibin] , bins[ibin+1] ]
      bin_p = float(bin_txt) * stdd
      it_sel = where((uband ge bin_p[0]) and (uband le bin_p[1]),np_p)
      if it eq 0 then print,'Count-p:',np_p

;      filetag=bandtag+'_'+bintag
;      title_tag=' ('+bin_txt[0]+' to '+bin_txt[1]+' sigma)'
  
;      for it=0,npd_im-1 do begin
;      
;        hrstr=string(it,format='(i2.2)')
;        print,'Hour: ',hrstr
;      
;        imfil=imfils[it]
;      
;        rain_sav=read_nc_var(imfil,'precipitationCal') ; mm/hr [ y, x, nd ]
;      
;        ;DETREND
;        for iyr=0,nyr-1 do $
;          rain_sav[*,*,ityr+iyr*npy] -= rebin(reform(rain_annual[*,*,iyr]),[ny,nx,npy])

        rain_ucomp[*,*,it,ibin,iband-1,icross-1]=mean(rain_sav[*,*,it_sel],dimension=3,/nan,/double)

;      endfor ; it (hour)

;      rain_file=compdir+'3B-HHR.MS.MRG.3IMERG.2000-2020_JJAS_cross'+strtrim(icross,2)+'_'+hrstr+'H_'+dat_str+'_unormcomp_'+filetag+'.nc4'
;      write_sing_ncvar,rain_file,rain,'rain',dim2=lon,dim1=lat,dimtag2='lon',dimtag1='lat'

    endfor ; ibin
  
  endfor ; iband

endfor ; icross

endfor ; it (hour)

for icross=1,2 do begin
  for iband=1,2 do begin ; Biweekly / Intraseasonal
    filetag=bandtag[iband-1]
    rain_file=compdir+'3B-HHR.MS.MRG.3IMERG.2000-2020_JJAS_cross'+strtrim(icross,2)+'_'+dat_str+'_unormcomp_'+filetag+'.nc4'
    write_sing_ncvar,rain_file,reform(rain_ucomp[*,*,*,*,iband-1,icross-1]),'rain_ucomp',$
      dim2=lon,dim1=lat,dimtag2='lon',dimtag1='lat',dim3=npd_im,dimtag3='npd',dim4=nbin,dimtag4='nbin'
  endfor ; iband
endfor ; icross

print,'Done!!'
end
