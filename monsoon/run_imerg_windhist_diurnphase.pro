; 
; Create maps of diurnal rainfall amplitude and phase (LT of peak)
; as a function of intraseasonal coast-normal wind.
;
; James Ruppert
; 8/8/21
; 
pro run_imerg_windhist_diurnphase

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

;ALL VARIABLES ARE NOW CALCULATED BY calc_windhist_diurnal.pro

;iwrite=1 ; write binned rainfall composites?

;coast_thresh=200 ; km

icross=1 ; Which cross section for unorm calculation?

;BOUNDS FOR PLOTTING
  bounds=[62,5,103,28] ; SASM region

;ERAi SETTINGS
;LEVEL SELECTION
  psel_era=850;500;700;925

;skip=4;12 ; Skip diurnal time steps for testing purposes?

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

  ifigdir=dirs.figdir+'myanmar/imerg/diurnal/'

;  datdir=dirs.wkdir+'imerg/imerg/data/'
  datdir=dirs.scdir+'imerg/'
  outdir=dirs.wkdir+'imerg/imerg/easthem/'
  compdir=outdir+'unorm_composites/'

  maindir=dirs.wkdir
  era_dir=maindir+'era5/'
  era_fil=era_dir+'ERA5-20000101-20201231-pl_dayavg.nc'

  npd_im=48
;  if keyword_set(skip) then npd_im/=skip

;Local SOLAR time conversion
;  local=6;round(mean(dims.lon)/360.*24.) ; deg lon --> hours
;  print,'Adding +'+strtrim(local,2)+' for LT'
;  ltim_test=findgen(npd_im)*24/npd_im+local
;  for it=0,npd_im-1 do ltim_test[it]-=24*(ltim_test[it] ge 24)

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


;=====BEGIN READING=========================================================

;----READ RAIN--------------------

  rain_file=outdir+'3B-HHR.MS.MRG.3IMERG.2000-2020_JJAS_allmean.V06B.nc4'
  lon=read_nc_var(rain_file,'lon')
  lat=read_nc_var(rain_file,'lat')
  nx=n_elements(lon)
  ny=n_elements(lat)

;if iwrite then $
;  rain_sav=read_imerg_dc_jjas(datdir,lon=lon,lat=lat,bounds=bounds,skip=skip) ; [x,y,hr,day] returns mm/day
;  rain_sav=read_imerg_dc_jjas(datdir,lon=lon,lat=lat,skip=skip) ; [x,y,hr,day] returns mm/day

;----READ ERA--------------------

;  if icross eq 1 then bounds_era=[89.,15.5,95.,22.5] ; Northern Myanmar coastline
;  if icross eq 2 then bounds_era=[69.5,8.,77.5,21.] ; Western Ghats
;  if icross eq 3 then bounds_era=[78.,14.,89.7,24.5] ; Northwestern BoB coastline
;
;  u=read_nc_era5(time,era_fil,'var131',plev=psel_era,lon=eralon,lat=eralat,bounds=bounds_era)
;  v=read_nc_era5(time,era_fil,'var132',plev=psel_era,lon=eralon,lat=eralat,bounds=bounds_era)
;  nxera=n_elements(eralon)
;  nyera=n_elements(eralat)
;
;;=====CREATE INDEX FOR REGRESSION/BINNING=========================================================
;
;  coastnormal, icross, coast_thresh, u, v, eralon, eralat, uindex=uindex
;
;;----FILTER--------------------
;
;  ;UNORM INDEX
;  filter_monsoon, uindex, 'fft', var_bw=u_bw, var_intra=u_intra
;;  filter_monsoon, uindex, 'rmean', var_bw=u_bw, var_intra=u_intra
;
;;=====GENERATE COMPOSITES=========================================================
;
;;----KEEP JJAS ONLY--------------------
;
;  nd=n_elements(jjas)
;
;;  uindex=uindex[jjas]
;  u_bw=u_bw[jjas]
;  u_intra=u_intra[jjas]
;
;  u=u[*,*,jjas]
;  v=v[*,*,jjas]
;
;;----STANDARDIZE--------------------
;
;  ;UNORM INDEX
;
;    std_bw=stddev(u_bw,/nan,/double)
;    std_intra=stddev(u_intra,/nan,/double)
;
;    print,'Standard Deviations:'
;    print,'BW:',std_bw
;    print,'Intra:',std_intra
;
;    u_bw    = (u_bw -    mean(u_bw,/nan,/double))    / std_bw
;    u_intra = (u_intra - mean(u_intra,/nan,/double)) / std_intra
;
;    ;N-CASES where |irain| > 1 sigma
;    ibw = where(abs(u_bw) ge 1,nbw)
;    iintra = where(abs(u_intra) ge 1,nintra)
;
;  ;PRINT STATS
;
;    print,'N where |index| >= 1:'
;    print,'BW:',nbw
;    print,'Intra:',nintra

;----BINNING--------------------

  bins = ['-1','-0.5','-0.25','0.25','0.5','1']

  nbin=n_elements(bins)-1

for iband=1,2 do begin ; Biweekly / Intraseasonal

  print,'iband: ',iband

  if iband eq 1 then begin
    bandtag='bw'
;    uband=u_bw
;    stdd=std_bw
  endif else if iband eq 2 then begin
    bandtag='intra'
;    uband=u_intra
;    stdd=std_intra
  endif

  for ibin=0,nbin-1 do begin

    bin_txt = [ bins[ibin] , bins[ibin+1] ]
;    bin_p = float(bin_txt) * stdd
;    it_sel = where((uband ge bin_p[0]) and (uband le bin_p[1]),np_p)
;    print,'Count-p:',np_p
    bintag=strtrim(ibin+1,2)

    filetag=bandtag+'_'+bintag
    title_tag=' ('+bin_txt[0]+' to '+bin_txt[1]+' sigma)'

rain=fltarr(ny,nx,npd_im)
for it=0,npd_im-1 do begin
  hrstr=string(it,format='(i2.2)')
;    rain_file=compdir+'imerg_3B-HHR.MS.MRG.3IMERG.V06_jjas_'+dat_str+'_unormcomp_'+filetag+'.nc4'
    rain_file=compdir+'3B-HHR.MS.MRG.3IMERG.2000-2020_JJAS_cross'+strtrim(icross,2)+'_'+hrstr+'H_'+dat_str+'_unormcomp_'+filetag+'.nc4'
;    if iwrite then begin
;      rain=mean(rain_sav[*,*,*,it_sel],dimension=4,/nan,/double)
;      write_sing_ncvar,rain_file,rain,'rain',dim1=lon,dim2=lat,dimtag1='lon',dimtag2='lat'
;      continue ; SKIP PLOTTING IF WRITING OUT
;    endif else begin
      rain[*,*,it]=read_nc_var(rain_file,'rain')
;      lon=read_nc_var(rain_file,'lon')
;      lat=read_nc_var(rain_file,'lat')
;    endelse
endfor

  rain*=24. ; mm/hr --> mm/d
  rmean=mean(rain,dimension=3,/nan,/double)
  rmean=transpose(temporary(rmean))


;=====SMOOTH RAINFALL=========================================================

  xsmth = 3;9 ; 9-pt running mean in x,y
;  rain = smooth(temporary(rain),[xsmth,xsmth,0],/edge_truncate)
  tsmth = 3 ; 3-pt running mean in hour
;  rain = smooth(temporary(rain),[0,0,tsmth],/edge_wrap)

;=====MAIN VARIABLES=========================================================

  ;REORDER LOCAL TIME TO PUT 0 LT AT FIRST INDEX
;  ishift=1.*local*npd_im/24
;  ltim_test=shift(ltim_test,ishift)
  ltim=round(lon/360.*24)
  for ix=0,nx-1 do rain[*,ix,*]=shift(reform(rain[*,ix,*]),[0,ltim[ix]*npd_im/24])

  phase=fltarr(nx,ny)
  amp=fltarr(nx,ny)
  for iy=0,ny-1 do $
    for ix=0,nx-1 do begin
      vmax=max(reform(rain[iy,ix,*]),ind,/nan)
      phase[ix,iy] = 1.*ind*24/npd_im
      amp[ix,iy] = vmax - min(rain[iy,ix,*],/nan)
    endfor

;=====PLOTS=========================================================

;for iplot=0,3 do begin
for iplot=0,2 do begin

  if iplot eq 0 or iplot eq 1 then begin

    if iplot eq 0 then begin
      var=rmean
      figname=ifigdir+'rain_mean'
      title='Mean'
    endif else if iplot eq 1 then begin
      var=amp
      figname=ifigdir+'rain_diurnal_amp'
      title='Diurnal Range'
    endif

    col_table=75;20;10;75;0;70
    ncols=15
    colors=findgen(ncols)/(ncols-1)*255
    max=24. & min=0
    levels=findgen(ncols)/(ncols-1)*(max-min)+min
    ndivs=6
    cbar_format='(i2)'
    cbar_tag='[ mm/d ]'

  endif else if iplot eq 2 then begin

    var=phase
    figname=ifigdir+'rain_diurnal_phase'
    title='Diurnal Phase'

    col_table=72;70
    levels=indgen(24)
    nlevs=n_elements(levels)
    colors=findgen(nlevs)/(nlevs-1)*254
colors=reverse(colors)
;    colors[12:23]=reverse(colors[12:23])

    ndivs=4
    cbar_format='(i2.2)'
    cbar_tag='[ LT ]'

  endif else if iplot eq 3 then begin

    var=amp / rmean
    figname=ifigdir+'rain_diurnal_ampnorm'
    title='Normalized Diurnal Range'

    col_table=75;20;10;75;0;70
    ncols=15
    colors=findgen(ncols)/(ncols-1)*255
    max=4. & min=0
    levels=findgen(ncols)/(ncols-1)*(max-min)+min
    ndivs=6
    cbar_format='(f3.1)'
    cbar_tag=' '

  endif

  figname+='_'+filetag

  ;PLOT SPECS
    csize=0.65
    position=[0.05,0.00,0.84,0.95]
    xsize=4.2 & ysize=3.8

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  area = [min(lat),min(lon),max(lat),max(lon)]
  if keyword_set(bounds) then $
    area = bounds[[1,0,3,2]]

  map_set,0,0,/mercator,/isotropic,limit=area,xmargin=0,ymargin=0,position=position

  loadct,col_table,/silent,file=dirs.ctfil

  ;FILL SHADING
    for i=0,1 do $
      contour,var,lon,lat,/cell_fill,/overplot,levels=levels,c_colors=colors

  loadct,0,/silent

  ;OVERLAY TITLE
    xyouts,mean(area[[1,3]]),area[2]+2.3,title+title_tag,align=0.5,charsize=csize,/data

  ;LAND
    landcol=0
;    landcol=255
    landthick=2.4
    hires=1;0
    if hires then landthick=0.8
    map_continents,/coasts,/countries,limit=area,color=landcol,mlinethick=landthick*0.35,hires=hires
    ;FILL LAND
    do_fill=0
;    do_fill=1
    if do_fill then begin
      map_continents,/coasts,limit=area,color=254,mlinethick=landthick,hires=hires,fill_continents=1
      map_continents,/coasts,limit=area,color=landcol,mlinethick=landthick,hires=hires
    endif else $
      map_continents,/coasts,limit=area,color=landcol,mlinethick=landthick,hires=hires

      ;LATLON GRID
      map_grid,/box_axes,latdel=10,londel=10,color=0,charsize=csize,glinethick=1.5,glinestyle=1,/no_grid

      ;COLOR BAR
      dy=0.3
      ;dy=0.35
      cpos= [ position[2]+0.048 ,$
      position[1]+dy ,$
      position[2]+0.064 ,$
      position[3]-dy ]
      loadct,col_table,/silent,file=dirs.ctfil
if iplot eq 2 then begin
  colors=[colors,colors[0]]
  setlevels=string(indgen(ndivs+1)*6,format='(i2.2)')
endif else setlevels=0
    colorbar2, colors=colors, range=[min(levels),max(levels)],divisions=ndivs,$
      charsize=csize, position=cpos, /right, /vertical, title=cbar_tag,$
      annotatecolor='black',format=cbar_format,$
      setlevels=setlevels
    loadct,0,/silent

  device,/close
  convert_png,figname,res=200,/remove_eps

endfor ; iplot

  endfor ; ibin

endfor ; iband


print,'Done!!'
end
