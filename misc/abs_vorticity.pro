; Function to calculate absolute vorticity
; 
; Assumes lon and lat are provided in deg, winds in m/s
; 
; Also assumes an x,y slice as var[x,y[,t]], and lon=lon[nx], lat=lat[ny].
;
; James Ruppert
; 12/6/21
function abs_vorticity, u, v, lon, lat

  dims=size(u)
  ndims=dims[0]
  if ndims eq 3 then nt=dims[3]

  nx=n_elements(lon)
  ny=n_elements(lat)

  ;LON TO DISTANCE
  dy=111e3*(lat[1]-lat[0]) ; deg --> m

  if ndims eq 2 then begin

    dudy=fltarr(nx,ny)
    dvdx=fltarr(nx,ny)
  
    for ix=0,nx-1 do $
      dudy[ix,*] = deriv( reform(u[ix,*]) )/dy
    for iy=0,ny-1 do begin
      dx = dy*cos(lat[iy]*!dtor)
      dvdx[*,iy] = deriv( reform(v[*,iy]) )/dx
    endfor
  
    avor = dvdx - dudy
  
  ;  ineg=where(lat le 0)
  ;  avor*=-1.;[*,ineg]*=-1.
  
    ;ADD CORIOLIS
  ;  twoom=2.*7.292e-5 ; /s
  ;  for iy=0,ny-1 do avor[*,iy] += twoom*sin(lat[iy]*!dtor)

  endif else begin

    dudy=fltarr(nx,ny,nt)
    dvdx=fltarr(nx,ny,nt)

    for ix=0,nx-1 do begin
      for it=0,nt-1 do dudy[ix,*,it] = deriv( reform(u[ix,*,it]) )
      dudy[ix,*,*]/=dy
    endfor
    for iy=0,ny-1 do begin
      dx = dy*cos(lat[iy]*!dtor)
      for it=0,nt-1 do dvdx[*,iy,it] = deriv( reform(v[*,iy,it]) )
      dvdx[*,iy,*]/=dx
    endfor

    avor = dvdx - dudy

  ;  ineg=where(lat le 0)
  ;  avor*=-1.;[*,ineg]*=-1.

    ;ADD CORIOLIS
  ;  twoom=2.*7.292e-5 ; /s
  ;  for iy=0,ny-1 do avor[*,iy] += twoom*sin(lat[iy]*!dtor)

  endelse

  return,avor

end
