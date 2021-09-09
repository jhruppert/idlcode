; 
; Calculate an x,y grid of radius[x,y] (in km) from a center location provided with
;  that location (in deg lon,lat) and a grid of lons[x] and lats[y].
;
; James Ruppert
; 4/29/19
; 
function radius_tc_ll, tcloc, lon, lat

nx=n_elements(lon)
ny=n_elements(lat)

;GENERATE LAT/LON GRID IN CARTESIAN (KM), WITH 0 AT DOMAIN-CENTER
  x=fltarr(nx,ny)
  y=x
  for iy=0,ny-1 do $
    x[*,iy]=111.*(lon-mean(lon)) * cos(lat[iy]*!pi/180)
  for ix=0,nx-1 do $
    y[ix,*]=111.*(lat-mean(lat))

;SHIFT INPUT GRID TO TC CENTER
  
  difx=abs(lon-tcloc[0])
  itcx=(where(difx eq min(difx),count))[0]

  dify=abs(lat-tcloc[1])
  itcy=(where(dify eq min(dify),count))[0]

  x_shift = x-x[itcx,itcy]
  y_shift = y-y[itcx,itcy]
  xy_radius = sqrt(x_shift^2+y_shift^2)

return,xy_radius

end
