; 
; Plot MSE variance maps from WRF TC output.
; 
; James Ruppert
; 3/9/2020
; 
pro wrf_mse_var_maps, figname, lon, lat, var, title

tcname='maria'

;CREATE FIGURE

  ;PLOT SPECS
    csize=0.8
    position=[0.06,0.06,0.89,0.89]
    xsize=3.3 & ysize=3.2

nt=(size(var,/dim))[2]
;for it=12,nt-1,12 do begin
for it=60,60 do begin
;for it=60,96,12 do begin

  print,'IT = ',strtrim(it,2)

  ivar=reform(var[*,*,it])

  ;ITERATIVELY FIND THE PLOTTING DOMAIN BASED ON NANS
  ix=where(finite(reform(ivar[*,20])),nfin)
  ipl=20
  while nfin le 10 do begin
    ix=where(finite(reform(ivar[*,10+ipl])),nfin)
    ipl+=20
  endwhile
  iy=where(finite(reform(ivar[ix[5],*])))

  area = [lat[iy[0]],lon[ix[0]],lat[max(iy)],lon[max(ix)]]

  ititle=title+' ('+strtrim(it,2)+'hr)'
  ifigname=figname+'_'+string(it,format='(i3.3)')+'hr'

  set_plot,'ps'
  epsname=ifigname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  map_set,0,0,/mercator,/isotropic,limit=area,xmargin=0,ymargin=0,position=position,$
    charsize=csize,title=ititle

  col_table=67
  loadct,col_table,/silent

  ;FILL SHADING

    ;MSE VAR
;      scale=1e-14
;      min=0
;      max=30
;      ncols=40
;      irev=1
;      cbar_format='(i2)'

    ;H'LW' and H'SEF' (normalized)
      scale=1.;e2
      max=2.
      min=-1.*max
      ncols=40
      irev=1
      cbar_format='(f4.1)'

    ;H'LW' and H'SEF' (UNnormalized)
      scale=1.e-9
      max=2.
      min=-1.*max
      ncols=40
      irev=1
      cbar_format='(f4.1)'

    colors=findgen(ncols)*255/ncols
    if irev then colors=reverse(colors)
    levels=findgen(ncols)/(ncols-1)*(max-min)+min
    for i=0,1 do $
      contour,ivar*scale,lon,lat,/cell_fill,/overplot,levels=levels,c_colors=colors

  loadct,0,/silent

  ;LAND
    landcol=0
;    landcol=255
    map_continents,/coasts,limit=area,color=landcol,mlinethick=1.8,/hires

  ;LATLON GRID
    map_grid,/box_axes,latdel=2,londel=2,color=0,charsize=csize*0.75,glinethick=1.5,glinestyle=1,/no_grid

  ;COLOR BAR
  icbar=1
  if icbar then begin
    cpos= [ position[2]+0.030 ,$
            position[1]+0.2 ,$
            position[2]+0.054 ,$
            position[3]-0.2 ]
    loadct,col_table,/silent
    colorbar2, colors=colors, range=[min(levels),max(levels)],$;divisions=figspecs.ndivs,$
      charsize=csize*0.6, position=cpos, /right, /vertical, $;title=figspecs.cbar_tag,$
      annotatecolor='black',format=cbar_format;,$
;      setlevels=setlevs
    loadct,0,/silent
  endif

  device,/close

  !p.multi=0

  convert_png,ifigname,/remove_eps,res=150

endfor ; it

end
