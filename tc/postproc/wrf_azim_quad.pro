; 
; Return azimuthally averaged  output for 4 quadrants from Cartesian input.
;
; Returns: quad_var( quadrant, radius, [height,] time )
;
; James Ruppert
; 4/6/19
; 
function wrf_azim_quad, output_dir, invar, var_str, tim_tag, radius=radius


ioverride=0


;RADIUS ARRAY
  nrad=round(800./3) ; 800 km at 3 km spacing
  radius=findgen(nrad)/nrad*800


azim_file=output_dir+var_str+'_'+tim_tag+'_azim.nc'
exists=file_test(azim_file)


;FIRST TRY TO READ
if exists and ~ioverride then begin


  quad_var=read_nc_var(azim_file,var_str)


;OTHERWISE GENERATE
endif else begin


print,'Generating azimuthal var!'


specs=size(invar)
ndims=specs[0]
nx=specs[1]
ny=specs[2]
if ndims eq 3 then i3d=0 $
else if ndims eq 4 then i3d=1 $
else message,'Check dimensions!'


;OUTPUT VAR
  if ~i3d then begin
    nt=specs[3]
    quad_var=fltarr(4,nrad,nt)
  endif else begin
    nz=specs[3]
    nt=specs[4]
    quad_var=fltarr(4,nrad,nz,nt)
  endelse
  quad_var[*]=!values.f_nan


;LOCATE PRESSURE MIN

;-------TAKEN FROM wrf_pres_min.pro------------
;READ AND SMOOTH SLP

  slp=reform(read_nc_var(output_dir+'post/SLP.nc',' ',varid='0'))

  ismooth=[3,3,0]
  for i=0,3 do $
    slp=smooth(temporary(slp),ismooth,/edge_truncate,/nan)

; CHECKED HOW THIS LOOKS; DOES AN ACCURATE JOB, AND EFFECTIVELY REMOVES NOISE AT SMALLEST SCALES

;LOCATE P-MIN

;  pres=fltarr(nt)
;  minloc=fltarr(2,nt)

;  for it=0,nt-1 do begin
;    pres[it]=min(reform(slp[*,*,it]),loc)
;    minloc[*,it]=array_indices([nx,ny],loc,/dimensions)
;  endfor
;----------------------------------------------


;TIME LOOP FOR COORDINATE TRANSFORM

pt_thresh=2 ; minimum n-points to average over

spawn,'ls '+output_dir+'wrfout_d03*',rawfils

for it=0,nt-1 do begin

  lon=reform(read_nc_var(rawfils[it],'XLONG',count=[nx,1,1],offset=[0,0,0]))
  lat=reform(read_nc_var(rawfils[it],'XLAT',count=[1,ny,1],offset=[0,0,0]))

  ;LOCATE PRESSURE MIN
    ipmin=min(reform(slp[*,*,it]),loc)
    minloc=array_indices([nx,ny],loc,/dimensions)

  ;ADJUST COORDINATES TO P-MIN
    center_x=lon[minloc[0]]
    center_y=lat[minloc[1]]
    ilon=lon-center_x
    ilat=lat-center_y

  ;RADIUS FROM CENTER [ KM ]
    dx=fltarr(nx,ny)
    dy=dx
    for iy=0,ny-1 do $
      dx[*,iy]=ilon*cos(lat[iy]*!dtor)
    dx*=111.
    for ix=0,nx-1 do $
      dy[ix,*]=ilat
    dy*=111.
    irad=sqrt(dx^2+dy^2)

  ;DO TRANSFORM

    if ~i3d then ivar=reform(invar[*,*,it]) $
    else ivar=reform(invar[*,*,*,it])

    for ir=1,nrad-1 do begin

      loc1=where(((irad ge radius[ir-1]) and (irad lt radius[ir])), npts)
      if npts lt 2 then continue

      ;UPPER-RIGHT
      loc2=where(((dx[loc1] ge 0) and (dy[loc1] ge 0)), npts)
      if npts ge 2 then begin
        if ~i3d then $
          quad_var[0,ir,it]=mean(ivar[loc1[loc2]],/nan,/double) $
        else $
          for iz=0,nz-1 do $
            quad_var[0,ir,iz,it]=mean((reform(ivar[*,*,iz]))[loc1[loc2]],/nan,/double)
      endif

      ;UPPER-LEFT
      loc2=where(((dx[loc1] lt 0) and (dy[loc1] ge 0)), npts)
      if npts gt 5 then begin
        if ~i3d then $
          quad_var[1,ir,it]=mean(ivar[loc1[loc2]],/nan,/double) $
        else $
          for iz=0,nz-1 do $
            quad_var[1,ir,iz,it]=mean((reform(ivar[*,*,iz]))[loc1[loc2]],/nan,/double)
      endif

      ;LOWER-LEFT
      loc2=where(((dx[loc1] lt 0) and (dy[loc1] lt 0)), npts)
      if npts gt 5 then begin
        if ~i3d then $
          quad_var[2,ir,it]=mean(ivar[loc1[loc2]],/nan,/double) $
        else $
          for iz=0,nz-1 do $
            quad_var[2,ir,iz,it]=mean((reform(ivar[*,*,iz]))[loc1[loc2]],/nan,/double)
      endif

      ;LOWER-RIGHT
      loc2=where(((dx[loc1] ge 0) and (dy[loc1] lt 0)), npts)
      if npts gt 5 then begin
        if ~i3d then $
          quad_var[3,ir,it]=mean(ivar[loc1[loc2]],/nan,/double) $
        else $
          for iz=0,nz-1 do $
            quad_var[3,ir,iz,it]=mean((reform(ivar[*,*,iz]))[loc1[loc2]],/nan,/double)
      endif


    endfor

endfor


  write_sing_ncvar,azim_file,quad_var,var_str


endelse
  

return,quad_var


end
