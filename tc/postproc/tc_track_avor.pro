; 
; Track a TC or precursor vortex via 
;   1 - large-scale smoothing of an input field (e.g., SLP)
;   2 - location prediction based on provided HURDAT data
;
;  Step 2 is important given a large domain where there may be multiple vortices
;    and/or the target vortex is weak.
;
; Set /write and trackfil (=filename) to write out locations to an ascii file.
;
; XX Set kernel to dramatically accelerate process; only needs calculating once.
;
; James Ruppert
; 12/20/21
; 
function tc_track_avor, invar_sav, dims, time, hurdat, ntskip=ntskip, write=write, trackfil=trackfil, kernel=kernel

invar=invar_sav

specs=size(invar,/dimensions)
nt=specs[2]

if ~keyword_set(ntskip) then ntskip=1

nthd=n_elements(hurdat.jultim)

;OUTPUT ARRAY
locvort=fltarr(3,nt)


;2D LON/LAT
  ilon=fltarr(dims.nx,dims.ny)
  ilat=ilon
  for ix=0,dims.nx-1 do ilat[ix,*]=dims.lat
  for iy=0,dims.ny-1 do ilon[*,iy]=dims.lon


;GAUSSIAN SETTINGS

  dx_sigma=111. ; 2 x sigma [km]
  dx = 111.*(dims.lat[1]-dims.lat[0]) ; delta-x in km
  sigma=dx_sigma/dx/2.

;RADIUS DEVIATION SETTINGS

  radmax=5. ; radius limit (in degrees) from Best Track center
  m2deg=1./(111e3)

  for it=0,dims.nt-1,ntskip do begin
;  for it=51,51 do begin

    tmp=reform(invar[*,*,it])

; 1 - IMPOSE GAUSSIAN SMOOTHING

    if ~keyword_set(kernel) then $
      tmp=gauss_smooth(temporary(tmp),sigma,/edge_truncate,/nan,kernel=kernel) $
    else $
      tmp=convol(temporary(tmp),kernel,total(kernel),/edge_truncate)

; 2 - MASK BASED ON TRACK OF REAL STORM

    t_diff=abs(time[it]-hurdat.jultim)
    ithd=(where(t_diff eq min(t_diff)))[0]

    ;PREDICT LOCATION AS X1 = X_0 + CX_0*DT

    dt = (time[it]-hurdat.jultim[ithd])*86400d ; s

    cx = hurdat.motion_x[ithd] ; m/s
    cy = hurdat.motion_y[ithd] ; m/s

    hdlat = hurdat.lat[ithd] + cy*dt*m2deg
    hdlon = hurdat.lon[ithd] + cx*dt*m2deg/(cos(hdlat*!dtor))

    radius=sqrt( (ilon-hdlon)^2 + (ilat-hdlat)^2 )
    inan=where(radius ge radmax)
    tmp[inan]=!values.f_nan

;OVERWRITE INPUT AVOR FIELD
invar_sav[*,*,it]=tmp

    ;VORTEX CENTER LAT/LON
      maxval=max(tmp,loc,/nan)
      ind=array_indices([dims.nx,dims.ny],loc,/dimensions)
      locvort[0,it:it+ntskip-1]=dims.lon[ind[0]]
      locvort[1,it:it+ntskip-1]=dims.lat[ind[1]]
      locvort[2,it:it+ntskip-1]=minval

  endfor ; it


;WRITE TO FILE

  if keyword_set(write) then $
    write_sing_ncvar,trackfil,vmax,'vmax',$
      dim1=indgen(3),dim2=indgen(nt),$
      dimtag1='lonlatvor1e5',dimtag2='time'

return,vmax

end
