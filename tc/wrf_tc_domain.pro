; 
; Plot WRF domain for WRF TC simulation.
; 
; James Ruppert
; 1/10/19
; 
pro wrf_tc_domain

  ;DIRS AND FILES
    wpsdir='/work/06040/tg853394/stampede2/wrfv4/WPS/'
    figdir='/home1/06040/tg853394/idl/figures/tc/'
    dom1fil=wpsdir+'geo_em.d01.nc'
    dom2fil=wpsdir+'geo_em.d02.nc'

  ;GRID
    lon1=read_nc_var(dom1fil,'XLONG_M')
    lat1=read_nc_var(dom1fil,'XLAT_M')
    lon2=read_nc_var(dom2fil,'XLONG_M')
    lat2=read_nc_var(dom2fil,'XLAT_M')
;    nx1=n_elements(lon)
;    ny1=n_elements(lat)

;-- CREATE PLOT ---------------------

  ;PLOT SPECS
    csize=0.8
    position=[0.05,0.03,0.88,0.95]
    xsize=4.2 & ysize=3.1
    area=[-14,360-90,40,360-5]

  set_plot,'ps'
  figname=figdir+'wrf_tc_domain'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  map_set,0,0,/mercator,/isotropic,limit=area,xmargin=0,ymargin=0,position=position

  ;OVERLAY TITLE
    title='Model Domain'
    xyouts,0.5,0.93,title,align=0.5,charsize=csize*1.,/normal

  ;MODEL DOMAIN
    loadct,4,/silent
    plots,lon1,lat1,psym=3,color=200,/data
    plots,lon2,lat2,psym=3,color=70,/data

  loadct,0,/silent

  ;LAND
    map_continents,/coasts,limit=area,color=0,/hires,mlinethick=0.8
;    map_continents,color=0,/coasts;,mlinethick=0.1

  ;LATLON GRID
    map_grid,/box_axes,limit=area,latdel=10,londel=10,color=0,charsize=csize*0.5,glinethick=0.6,glinestyle=1;,/no_grid
;    map_grid,lons=[indgen((120-50)/10+1)*10+50],lats=[-20,-10,0,10,20],/box_axes,latlab=49,lonlab=-24,label=1,color=0,charsize=csize*0.5,/no_grid

  device,/close

  convert_png,figname,res=200,/remove_eps

print,'DONE!!'
end
