; 
; TC tracks.

; TC MOTION IS ALSO CALCULATED HERE

; James Ruppert
; 3/28/19
; 
pro run_tc_track

tcname='maria'
tcname='haiyan'

if tcname eq 'maria' then begin
  tcyear='2017'
  hurdat=read_hurdat(tcname,tcyear)
endif else if tcname eq 'haiyan' then begin
  tcyear='2013'
  hurdat=read_jtwcdat(tcname)
endif

;subdir='moving_nest/'+tcname
;subdir='static_nest/'+tcname
subdir='redux/'+tcname
tc_sim_config, subdir, dirs=dirs, dims=dims, vars=vars;, /verbose
dirs.figdir+=tcname+'/'

underlay=0 ; plot field as underlay?
  undervar_str='SST'
  dom='01'

;VORTEX LEVEL
  psel=700 ; (hPa) level to locate vortex
  izsel=(where(dims.pres eq psel))[0]


;TIME ARRAYS0
  time=dims.time
  nt_full=dims.nt;-1
  npd=dims.npd
  nhrs=1.*nt_full*npd/24.
  hours=indgen(nhrs)
  nd=(nhrs-(nhrs mod 24))/24.

;TIME SUBSET
  it_sel=indgen(nt_full)
;it_sel=indgen(73)
  nt=n_elements(it_sel)

;----READ VARS--------------------


loc=fltarr(dirs.nc,2,nt)
loc[*]=!values.f_nan

;LAND MASK
  vtag='LANDMASK'
  file=dirs.files_raw[0,2,0]
  mask=reform(read_nc_var(file,'LANDMASK'))
;  land=where(mask eq 1,nland)
;  mask=0

for ic=0,dirs.nc-1 do begin
;for ic=0,0 do begin

  print,strupcase(dirs.cases[ic])

    i_nt=nt_full
    if ( strmatch(dirs.cases[ic],'*36h*') or strmatch(dirs.cases[ic],'icrf_*') or strmatch(dirs.cases[ic],'lwcr*') $
      or (dirs.cases[ic] eq 'lwswcrf') or (dirs.cases[ic] eq 'axisym') ) then i_nt-=36
    if strmatch(dirs.cases[ic],'*24h*') then i_nt-=24
    if strmatch(dirs.cases[ic],'*48h*') then i_nt-=48
    if strmatch(dirs.cases[ic],'*60h*') then i_nt-=60
    if strmatch(dirs.cases[ic],'*72h*') then i_nt-=72
    if strmatch(dirs.cases[ic],'*84h*') then i_nt-=84
    if strmatch(dirs.cases[ic],'*96h*') then i_nt-=96
    it_test=indgen(i_nt)+nt_full-i_nt

  ;READ ABSOLUTE VORTICITY
    iv=where(vars.vars eq 'AVOR')
    file=dirs.files_post[ic,iv]
    count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,izsel,0] ; x,y,z,t
    avor=reform(read_nc_var(file,'AVOR',count=count,offset=offset))

  specs=size(avor,/dimensions)
;  i_nt=specs[2]
;  it_test=indgen(i_nt)

  ;SMOOTH
    ixsmooth=round(111./3) ; 1-degree smoothing, run twice
    ismooth=[ixsmooth,ixsmooth,0]
    for i=1,2 do $
      avor=smooth(temporary(avor),ismooth,/edge_truncate,/nan)

  ;MASKING
    imask=fltarr(dims.nx,dims.ny,i_nt)
    for it=0,i_nt-1 do imask[*,*,it]=mask
    land=where(imask eq 1,nland)
    if nland gt 0 then avor[land]=!values.f_nan
    if tcname eq 'maria' then avor=wrf_maria_mask(temporary(avor),time[it_test],hurdat,dims)
    if tcname eq 'haiyan' then avor=wrf_haiyan_mask(temporary(avor),time[it_test],hurdat,dims)

  ;VORTEX TRACKING
    vloc=maria_vortex_locate(avor,dims);,/write)

  ;SMOOTH TRACKS
    for i=0,1 do $
      vloc=smooth(temporary(vloc),[0,3],/edge_truncate,/nan)

  ;SKIP 1ST 1.5 DAYS FOR HAIYAN CONTROL
    if dirs.cases[ic] eq 'ctl' and tcname eq 'haiyan' then $
      vloc[*,0:npd-1+npd/2]=!values.f_nan

;    ipmin = wrf_pres_min(dirs.casedir[ic])
    loc[ic,0,it_test]=vloc[0,*];ipmin.lon
    loc[ic,1,it_test]=vloc[1,*];ipmin.lat

;  ;CALCULATE AND PRINT STORM MOTION
;    loc_x = reform(vloc[0,*])
;    loc_y = reform(vloc[1,*])
;    motion_x = 111d3 * cos(loc_y*!pi/180) * deriv(loc_x) / (deriv(time[it_test]) * 24*3600)
;    motion_y = 111d3 *                      deriv(loc_y) / (deriv(time[it_test]) * 24*3600)
;    motion_file=dirs.casedir[ic]+'storm_motion.txt'
;    openw,1,motion_file
;      printf,1,i_nt
;      printf,1,it_test
;      printf,1,motion_x
;      printf,1,motion_y
;      printf,1,loc_x
;      printf,1,loc_y
;    close,1
;continue

  ;OVERWRITE WITH PRESSURE MIN FOR T=48 AND ONWARD
;  t_loc=where(hours[it_test] ge 48,count)
;  if count gt 0 then begin
;    ipmin = wrf_pres_min(dirs.casedir[ic],time[it_test],hurdat,dims)
;    loc[ic,0,it_test[t_loc]]=ipmin.lon[it_test[t_loc]]
;    loc[ic,1,it_test[t_loc]]=ipmin.lat[it_test[t_loc]]
;  endif

  if underlay then begin
    spawn,'ls '+dirs.casedir[ic]+'wrfout_d'+dom+'_*',rawfils
    undervar=reform(read_nc_var(rawfils[0],undervar_str))
    undervar-=273.15
    lonv=reform(read_nc_var(rawfils[0],'XLONG'))
    latv=reform(read_nc_var(rawfils[0],'XLAT'))
    lmask=reform(read_nc_var(rawfils[0],'LANDMASK'))
    undervar[where(lmask eq 1)]=!values.f_nan
  endif

endfor ; icase


;----CREATE PLOT--------------------


  figname=dirs.figdir+'tracks'

  icplot=indgen(dirs.nc)
  ncplot=n_elements(icplot)

  ;PLOT SPECS
    csize=1.0;0.8
    position=[0.04,0.03,0.92,0.92]
    xsize=4.2 & ysize=3.4

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  if tcname eq 'joaquin' then $
    area=[15,-90,37,-50] $
  else if tcname eq 'maria'  then $
    area=[dims.lat[0],min(dims.lon),max(dims.lat),max(dims.lon)] $
  else $
    area=[dims.lat[0],min(dims.lon),max(dims.lat),max(dims.lon)]
;    area=[0,-80,40,-25]

  map_set,0,0,/mercator,/isotropic,limit=area,xmargin=0,ymargin=0,position=position;,color=254;,$
;    charsize=csize,title=figspecs.title

  ;OVERLAY TITLE
;    xyouts,0.5,0.95,figspecs.title,align=0.5,charsize=csize*1.,/normal

  ;FILL SHADING
  if underlay then begin
    min=20
    max=35
    ncols=150
    colors=findgen(ncols)/(ncols-1)*255
    colors=reverse(colors)
    levels=findgen(ncols)/(ncols-1)*(max-min)+min
    loadct,70,/silent
    for i=0,1 do $
      contour,undervar,reform(lonv[*,0]),reform(latv[0,*]),/cell_fill,/overplot,$
        levels=levels,c_colors=colors
    loadct,0,/silent
    contour,undervar,reform(lonv[*,0]),reform(latv[0,*]),/follow,/overplot,$
      levels=indgen(200)*1,c_colors=0
  endif

  ;LAND
    map_continents,/coasts,limit=area,color=0,mlinethick=1.0,/hires;,/countries
;    map_continents,color=0,/coasts,/hires,mlinethick=0.1

  ;LATLON GRID
    del=10
    map_grid,/box_axes,latdel=del,londel=del,color=0,charsize=csize*0.7,glinethick=1.5,glinestyle=1,/no_grid

  loadct,4,/silent

  ;TC TRACKS
    psav=!P.Font
    !P.Font=-1
    rmax=0.28

    ;HURDAT
      loadct,0,/silent
      hurcol=200
      plots,hurdat.lon,hurdat.lat,thick=4.0,linestyle=0,color=hurcol,/data
      nt_hur=n_elements(hurdat.jultim)
      caldat,hurdat.jultim,mm,dd,yy,hh
      if tcname eq 'maria' then hh_circ=12
      if tcname eq 'haiyan' then hh_circ=0
      i24=where(hh eq hh_circ, n24)
      ;DAILY CIRCLES
      for it=0,n24-1 do $
        polyfill,rmax*cos(findgen(50)/49*2*!pi)/cos(hurdat.lat[i24[it]]*!dtor)+hurdat.lon[i24[it]],$
                 rmax*sin(findgen(50)/49*2*!pi)+hurdat.lat[i24[it]],/data,color=hurcol

    ;MODEL
    loadct,3,/silent
;    loadct,0,/silent
    colors=(findgen(dirs.nc)+1)/dirs.nc*220
    colors[0]=0
    lthick=reverse((findgen(dirs.nc)+3)*4.5/dirs.nc)
lthick[*]=2.
    for ic=0,dirs.nc-1 do begin
;    for ic=1,dirs.nc-1 do begin
;      if ~finite(loc[ic,0,5]) then continue
      ;TRACKS
      plots,loc[ic,0,*],loc[ic,1,*],thick=lthick[ic],linestyle=0,color=colors[ic],/data
;      for it=0,nt-1,24 do $
      cols=[0,160]
;      icol=0
      ;DAILY CIRCLES
      for it=0,nt-1,24 do begin
;        xyouts,loc[ic,0,it],loc[ic,1,it],'!20<!X',charsize=0.5,charthick=2.5,align=0.5,color=colors[ic],/data
;        plots,loc[ic,0,it],loc[ic,1,it],psym=6,color=colors[ic],symsize=0.4,/data
        polyfill,rmax*cos(findgen(50)/49*2*!pi)/cos(loc[ic,1,it]*!dtor)+loc[ic,0,it],rmax*sin(findgen(50)/49*2*!pi)+loc[ic,1,it],$
          /data,color=colors[ic]
;        icol+=1
;        if icol eq 1 then icol=0
      endfor
;print,loc[ic,0,*]
    endfor
    !P.Font=psav

  ;LEGEND
  ileg=1
  if ileg then begin
    csize_fac=0.6;0.7
    margin=0.1
    pspacing=2.0 ; length of lines
    spacing=0.8 ; between lines
    leg_str=strupcase(dirs.cases[icplot])
    leg_style=replicate(0,ncplot)
    leg_thick=lthick;replicate(2,ncplot)
    leg_color=colors[icplot]
    legend2,leg_str,linestyle=leg_style,thick=leg_thick,COLORS=leg_color,$
      charsize=csize*csize_fac,/top_legend,/left_legend,/normal,$
      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.36,0.75]
    ;HURDAT
    loadct,0,/silent
    legend2,['HURDAT'],linestyle=0,thick=4,COLORS=hurcol,$
      charsize=csize*csize_fac,/top_legend,/left_legend,/normal,$
      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.16,0.75]
  endif


  device,/close

  convert_png,figname,res=200;,/remove_eps


print,'DONE!!'
end
