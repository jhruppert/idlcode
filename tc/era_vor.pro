; 
; Plot basic maps from ERA-I for TC cases.
; 
; James Ruppert
; 1/9/19
; 
pro era_vor

iv=2 ; 1=vo, 2=pv
lev=150;500;150;925 ; pressure level (hPa)
iwind=1 ; plot wind? 0 or 1
isingle=0 ; plot just the first time? 0 or 1

  ;OVERLY WRF MODEL DOMAINS
    mod_dom1=[ -75.7992 , -7.7992 , -14.094 , 32.206 ]
    mod_dom2=[ -62.779  ,  1.8137 , -30.461 , 22.789 ]

  ;DIRS AND FILES
    figdir='/home1/06040/tg853394/idl/figures/tc/era/'
    datdir='/work/06040/tg853394/stampede2/erai/maria_vor/'
    era=datdir+'vor.2017-09-10_2017-09-20.nc'
    eraw=datdir+'wind.2017-09-10_2017-09-20.nc'
    era_daily=datdir+'vor_daily.2017-09-10_2017-09-20.nc'

  ;DATES
    time=read_nc_var(era,'time') ; hours since 1900-01-01
    nt=n_elements(time)
    time=timegen(nt,start=julday(9,10,2017,0,0,0),step_size=6,units='H')

  ;GRID
    lon=read_nc_var(era,'longitude')
    lat=read_nc_var(era,'latitude')
    pres=read_nc_var(era,'level') ; hPa
    nx=n_elements(lon)
    ny=n_elements(lat)
    np=n_elements(pres)

;-- READ DATA -----------------------

  if iv eq 1 then vstr='vo' else vstr='pv'

  ;PRESSURE LEVEL
    ip=(where(pres eq lev))[0]
    if ip eq -1 then message,'Fix level.'

  ;READ
    count=[nx,ny,1,nt] & offset=[0,0,ip,0]
    var=read_nc_var(era,vstr,count=count,offset=offset)
    var=reform(var)
    ;WIND
    u=read_nc_var(eraw,'u',count=count,offset=offset)
    v=read_nc_var(eraw,'v',count=count,offset=offset)
    u=reform(u) & v=reform(v)

  ;SCALE AND CONTOUR SETTINGS
    if iv eq 1 then begin
      title_tag='Rel. Vor.'
      var*=1e5 ; --> 10^-5 /s
      vmax=10.
      ctbl=71
      irev=1
      cbar_tag='[ 10!U-5!N s!U-1!N ]'
    endif
    if iv eq 2 then begin
      title_tag='Pot. Vor.'
      var*=1e6 ; --> PVU ( = 10^-6 K * m**2 / (kg * s ) )
      vmax=2.
      irev=1
      ctbl=71
      cbar_tag='[ PVU ]'
    endif

;-- CREATE PLOT ---------------------

  ;PLOT SPECS
    csize=0.8
    position=[0.05,0.03,0.88,0.95]
    ;LARGE BASIN
      xsize=4.2 & ysize=3.0
;      area=[-1,360-75,35,360-10]
      area=[-14,360-90,40,360-5]
    ;MODEL DOMAIN
;      xsize=4.2 & ysize=3.1
;      area=[2,360-63,23,360-32]

  ;AXES
    x=lon
    y=lat
    if ~keyword_set(xrange) then $
      xrange=[min(x),max(x)]
    if ~keyword_set(xrange) then $
      yrange=[min(y),max(y)]

if isingle then top=0 else top=nt-1
for it=0,top do begin

  caldat,time[it],mm,dd,yy,hh
  timstr=string(dd,format='(i2.2)')+'-'+string(hh,format='(i2.2)')

  levstr=strtrim(lev,2)

  print,'Date: ',timstr

  set_plot,'ps'
  figname=figdir+vstr+'_'+levstr+'mb_'+timstr
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  map_set,0,0,/mercator,/isotropic,limit=area,xmargin=0,ymargin=0,position=position

  loadct,ctbl,/silent

  ;SHADING
    nlevs=25;15
    levs=findgen(nlevs)*2*vmax/(nlevs-1) - vmax
    cols=findgen(nlevs)*255/(nlevs-1)
    cols=reverse(cols)
    for i=0,1 do $
      contour,reform(var[*,*,it]),x,y,/cell_fill,/overplot,$
        levels=levs,c_colors=cols

  loadct,0,/silent

  ;WIND
    if iwind then begin
      nskip=3
      xind=indgen(nx)*nskip & xind=xind[where(xind lt nx)]
      yind=indgen(ny)*nskip & yind=yind[where(yind lt ny)]
      velovect,reform(u[xind,yind,it]),reform(v[xind,yind,it]),x[xind],y[yind],$
        color=0,length=2.0,/overplot
    endif

  ;OVERLAY TITLE
    xyouts,0.5,0.94,title_tag+' ('+levstr+' hPa; '+timstr+')',align=0.5,charsize=csize*1.,/normal

  ;MODEL DOMAIN
    lthick=1
    plots,[mod_dom1[[0,2]]],replicate(mod_dom1[1],2),linestyle=0,thick=lthick,/data
    plots,[mod_dom1[[0,2]]],replicate(mod_dom1[3],2),linestyle=0,thick=lthick,/data
    plots,replicate(mod_dom1[0],2),mod_dom1[[1,3]],linestyle=0,thick=lthick
    plots,replicate(mod_dom1[2],2),mod_dom1[[1,3]],linestyle=0,thick=lthick
    plots,[mod_dom2[[0,2]]],replicate(mod_dom2[1],2),linestyle=0,thick=lthick,/data
    plots,[mod_dom2[[0,2]]],replicate(mod_dom2[3],2),linestyle=0,thick=lthick,/data
    plots,replicate(mod_dom2[0],2),mod_dom2[[1,3]],linestyle=0,thick=lthick
    plots,replicate(mod_dom2[2],2),mod_dom2[[1,3]],linestyle=0,thick=lthick

  ;LAND
    map_continents,/coasts,limit=area,color=0,/hires,mlinethick=0.8
;    map_continents,color=0,/coasts;,mlinethick=0.1

  ;LATLON GRID
    map_grid,/box_axes,limit=area,latdel=10,londel=10,color=0,charsize=csize*0.5,glinethick=0.6,glinestyle=1;,/no_grid
;    map_grid,lons=[indgen((120-50)/10+1)*10+50],lats=[-20,-10,0,10,20],/box_axes,latlab=49,lonlab=-24,label=1,color=0,charsize=csize*0.5,/no_grid

  ;COLOR BAR
    loadct,ctbl,/silent
    cpos= [ position[2]+0.04 ,$
            position[1]+0.3 ,$
            position[2]+0.05 ,$
            position[3]-0.3 ]
    colorbar2, colors=cols, range=[min(levs),max(levs)],$; divisions=ndivs,$
      charsize=csize*0.5, position=cpos, /right, /vertical, title=cbar_tag, $
      divisions=4,annotatecolor='black';,format='(i4)'
    loadct,0,/silent

  device,/close

  convert_png,figname,res=200,/remove_eps

endfor

print,'DONE!!'
end
