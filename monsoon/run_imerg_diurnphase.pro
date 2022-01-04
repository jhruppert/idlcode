;function read_imerg, dirs, lon=lon, lat=lat, bounds=bounds, skip=skip
;
;;datdir=dirs.wkdir+'imerg/imerg/data/'
;;file=datdir+'jjas_2013-2017_diurncomp_3B-HHR.MS.MRG.3IMERG.V06B.nc4'
;datdir=dirs.wkdir+'imerg/imerg/easthem/'
;file_dc=datdir+'idl_detrend_3B-HHR.MS.MRG.3IMERG.2000-2020_JJAS_diurncomp.V06B.nc4'
;
;;iwrite=0
;
;;if iwrite then begin
;;  rain=read_imerg_dc_jjas(datdir,lon=lon,lat=lat,bounds=bounds,skip=skip) ; [x,y,hr,day] returns mm/day
;;  rain=mean(rain,dimension=4,/double,/nan)
;;  write_sing_ncvar,file,rain,'rain',dim1=lon,dim2=lat,dimtag1='lon',dimtag2='lat'
;;endif else begin
;  rain=read_nc_var(file_dc,'rain_dc')
;  lon=read_nc_var(file_dc,'lon')
;  lat=read_nc_var(file_dc,'lat')
;;endelse
;
;return,rain
;
;end
; 
; Create maps of diurnal rainfall amplitude and phase (LT of peak).
;
; James Ruppert
; 7/22/21
; 
pro run_imerg_diurnphase

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

;BOUNDS FOR PLOTTING
  bounds=[66,2,108,25] ; SASM region
;  bounds=[63,5,103,28] ; SASM region
;  bounds=[60,-20,156,34] ; Eastern hemisphere

;PRESSURE LEVEL FOR VORTICITY
;  psel_era=500
  psel_era=850

;skip=4;12 ; Skip diurnal time steps for testing purposes?

;SELECT DATE RANGE
;dat_str='2000-2020'
;dat_str='2013-2017'

;----DIRECTORIES--------------------

  datdir=dirs.wkdir+'imerg/imerg/easthem/'
  ifigdir=dirs.figdir+'myanmar/imerg/diurnal/'

  npd_im=48
;  if keyword_set(skip) then npd_im/=skip

;Local SOLAR time conversion
;  local=6;round(mean(dims.lon)/360.*24.) ; deg lon --> hours
;;  print,'Adding +'+strtrim(local,2)+' for LT'
;  ltim_test=findgen(npd_im)*24/npd_im+local
;  for it=0,npd_im-1 do ltim_test[it]-=24*(ltim_test[it] ge 24)

;=====BEGIN READING=========================================================

  ;DIURNAL COMPOSITE CALCULATED BY FIRST DETRENDING BY SUBTRACTING DAILY MEAN
  file_dc=datdir+'idl_detrend_3B-HHR.MS.MRG.3IMERG.2000-2020_JJAS_diurncomp.V06B.nc4'
  rain_dc=read_nc_var(file_dc,'rain_dc')
  lon=read_nc_var(file_dc,'lon')
  lat=read_nc_var(file_dc,'lat')

;  rain_dc=read_imerg(dirs,lon=lon,lat=lat,bounds=bounds,skip=skip) * 24. ; mm/h --> mm/d

  dims=size(rain_dc,/dimen)
  npd_im=dims[2] ; Overwrite this in case of skip
  nx=dims[1]
  ny=dims[0]

  file_mean=datdir+'3B-HHR.MS.MRG.3IMERG.2000-2020_JJAS_allmean.V06B.nc4'
  rain_mean=reform(read_nc_var(file_mean,'precipitationCal'))

  rain = rain_dc + rebin(rain_mean,dims)
  rain *= 24. ; mm/h --> mm/d

  rmean=transpose(rain_mean*24.)


;=====SMOOTH RAINFALL=========================================================

  xsmth = 3;9 ; 9-pt running mean in x,y
;  rain = smooth(temporary(rain),[xsmth,xsmth,0],/edge_truncate)
  tsmth = 3 ; 3-pt running mean in hour
;  rain = smooth(temporary(rain),[0,0,tsmth],/edge_wrap)

;=====READ ERA=========================================================

  era_fil_p=dirs.wkdir+'era5/ERA5-20000101-20201231-pl_JJAS_allmean.nc4'
    ; u,v [ m/s ]
    ; RH [ % ]

  u=read_nc_era5(time,era_fil_p,'var131',plev=psel_era,lon=eralon,lat=eralat,bounds=bounds)
  v=read_nc_era5(time,era_fil_p,'var132',plev=psel_era,lon=eralon,lat=eralat,bounds=bounds)

  avor=abs_vorticity(u,v,eralon,eralat)

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

;for iplot=0,4 do begin
for iplot=2,2 do begin

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
;max=30
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

  endif else if iplot eq 4 then begin
    var=avor * 1e6 ; 10^-6 /s
var=smooth(var,[3,3],/edge_truncate)
    figname=ifigdir+'mean_avor_'+strtrim(psel_era,2)
    title='Relative Vorticity ('+strtrim(psel_era,2)+' hPa)'

    lon=eralon
    lat=eralat

    col_table=75
    ncols=15
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    max=30. & min=0
    levels=findgen(ncols)/(ncols-1)*(max-min)+min
    ndivs=6
    cbar_format='(i2)'
    cbar_tag='10!U-6!N /s'

  endif

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
;area = [12,85,23,100]

  map_set,0,0,/mercator,/isotropic,limit=area,xmargin=0,ymargin=0,position=position

  loadct,col_table,/silent,file=dirs.ctfil

  ;FILL SHADING
    for i=0,1 do $
      contour,var,lon,lat,/cell_fill,/overplot,levels=levels,c_colors=colors

  loadct,0,/silent

  ;OVERLAY TITLE
    dy=2.3
    if bounds[1] lt 0 then dy=5
    xyouts,mean(area[[1,3]]),area[2]+dy,title,align=0.5,charsize=csize,/data

  if iplot eq 2 then begin
    ;ICROSS = 1
      xcross=[77.8,96.2] ; lon,lat
      ycross=[12.,21.5]   ; lon,lat
      cross=[xcross,ycross]
      plots,cross[0:1],cross[2:3],linestyle=0,thick=3,/data
    ;ICROSS = 3
      xcross=[83.,93.5] ; lon,lat
      ycross=[21.4,10]   ; lon,lat
      cross=[xcross,ycross]
      plots,cross[0:1],cross[2:3],linestyle=1,thick=3,/data
  endif

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
  convert_png,figname,res=300,/remove_eps

endfor ; iplot


print,'Done!!'
end
