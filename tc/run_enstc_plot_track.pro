; 
; TC tracks for ensemble runs.

; TC MOTION IS ALSO CALCULATED HERE

; James Ruppert
; 12/20/21
; 
pro run_enstc_plot_track

tcname='maria'
case_str='ctl'
tcname='haiyan'
case_str='haiyan'

if tcname eq 'maria' then begin
  tcyear='2017'
  hurdat=read_hurdat(tcname,tcyear)
endif else if tcname eq 'haiyan' then begin
  tcyear='2013'
  hurdat=read_jtwcdat(tcname)
endif

dom='d02'
tc_ens_config, case_str, dom, dirs=dirs, dims=dims, vars=vars;, /verbose
dirs.figdir+=tcname+'/'

underlay=0 ; plot field as underlay?
  undervar_str='SST'
  dom='01'

;VORTEX LEVEL
  psel=700 ; (hPa) level to locate vortex
  izsel=(where(dims.pres eq psel))[0]

;TIME ARRAYS
  time=dims.time
  nt_full=dims.nt;-1
  npd=dims.npd
  nhrs=1.*nt_full*npd/24.
  hours=indgen(nhrs)
  nd=(nhrs-(nhrs mod 24))/24.

;TIME SUBSET
  it_sel=indgen(nt_full)
  nt=n_elements(it_sel)

;----READ VARS--------------------


;locvort=fltarr(dirs.nens,3,nt)
;locvort[*]=!values.f_nan

for ic=0,dirs.nens-1 do begin
;for ic=0,5 do begin
;for ic=0,1 do begin

  print,'Memb: ',ic+1

  ;TC TRACKING
    spawn,'ls '+dirs.ensdir[ic]+'track*',trackfiles,err
    ntrack=n_elements(trackfiles)
    if ~keyword_seT(trackfiles) then continue
;locvort=0

  for ifile=0,ntrack-1 do begin

    temp=read_nc_var(trackfiles[ifile],'locvort')
    temp=reform(temporary(temp),[4,dims.nt,1])

    if ~keyword_set(locvort) then $
      locvort=temp $
    else $
      locvort=[[[locvort]],[[temp]]]

  endfor
vortsize=size(locvort)
if vortsize[0] eq 2 then locvort=reform(temporary(locvort),[4,dims.nt,1])

;  if underlay then begin
;    spawn,'ls '+dirs.casedir[ic]+'wrfout_d'+dom+'_*',rawfils
;    undervar=reform(read_nc_var(rawfils[0],undervar_str))
;    undervar-=273.15
;    lonv=reform(read_nc_var(rawfils[0],'XLONG'))
;    latv=reform(read_nc_var(rawfils[0],'XLAT'))
;    lmask=reform(read_nc_var(rawfils[0],'LANDMASK'))
;    undervar[where(lmask eq 1)]=!values.f_nan
;  endif

endfor ; iens


;----CREATE PLOT--------------------


  figname=dirs.figdir+'tracks'
;figname=dirs.figdir+'tracks_memb'+string(ic+1,format='(i2.2)')

;  icplot=indgen(dirs.nc)
;  ncplot=n_elements(icplot)

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
    area=[dims.lat[0],dims.lon[0],max(dims.lat),max(dims.lon)] $
  else $
    area=[dims.lat[0],dims.lon[0],max(dims.lat),max(dims.lon)]
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
    loadct,13,/silent
;    loadct,0,/silent
    n_track=(size(locvort,/dimen))[2]
    colors=(findgen(n_track)+1)/n_track*254
colors+=30
    colors[0]=0
    lthick=reverse((findgen(n_track)+3)*4.5/n_track)
lthick[*]=1.
;    for itrack=0,dirs.nens-1 do begin
    for itrack=0,n_track-1 do begin
;    for itrack=0,0 do begin
      ;TRACKS
;colors[*]=250
      plots,locvort[0,*,itrack],locvort[1,*,itrack],thick=lthick[itrack],linestyle=0,color=colors[itrack],/data
;      for it=0,nt-1,24 do $
      cols=[0,160]
;      icol=0
      ;DAILY CIRCLES
;      for it=0,nt-1,24 do begin
;        polyfill,rmax*cos(findgen(50)/49*2*!pi)/cos(locvort[itrack,1,it]*!dtor)+locvort[itrack,0,it],rmax*sin(findgen(50)/49*2*!pi)+locvort[itrack,1,it],$
;          /data,color=colors[itrack]
;;        icol+=1
;;        if icol eq 1 then icol=0
;      endfor
    endfor
    !P.Font=psav

  ;LEGEND
  ileg=0
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

  convert_png,figname,res=300,/remove_eps

;endfor

print,'DONE!!'
end
