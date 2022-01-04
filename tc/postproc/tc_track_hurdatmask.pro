;WRITE TO FILE
pro write_tc_track,trackfile,locvort
  write_sing_ncvar,trackfile,locvort,'locvort',$
    dim1=indgen(3),dimtag1='lonlatslpcxcy',dimtag2='time'
end
; 
; Track a TC or precursor vortex
;
; Procedure steps:
;   Large-scale smoothing of SLP in time and space
;   Calculated as anomaly from large-scale and time-average
;   Location prediction based on real TC location from HURDAT
;
;  Step 2 is important given a large domain where there may be multiple vortices
;    and/or the target vortex is weak.
;
; Set /write and trackfile (=filename) to write out locations to an ascii file.
;
; James Ruppert
; 12/20/21
; 
function tc_track_hurdatmask, slp_sav, dims, time, hurdat, write=write, trackfile=trackfile, kernel=kernel

;SETTINGS
  nx_sm=round(0.5*1./(dims.lat[1]-dims.lat[0])) ; Half-degree smoothing in space
  nx_sm=9    ; Following Chavas 2013 (smoothing run x30 times)
  nt_sm=3    ; temporal smoothing (n time steps)
  lmin=1e3   ; Minimum distance threshold [km] between storms
      radmax=7.  ; radius limit (in degrees) from Best Track center
  it_zero=24 ; set all time steps prior to this one to NaN
  c_max=20.  ; a very liberal maximum allowable translation speed for tracking
  pmin_thresh1=-5. ; lifetime-minimum pressure as threshold for identifying a cyclone
  pmin_thresh2=-3. ; after storm ID, sequential time steps must retain center of this strength
  nt_min=24. ; minimum time steps for a cyclone to be retained

  m2deg=1./(111e3)

;SET SOME ARRAYS
  slp=slp_sav
  locvort=fltarr(5,dims.nt) ; lon[deg], lat, SLP[hPa], c_x[m/s], c_y
  locvort[*]=!values.f_nan

;2D LON/LAT
  ilon=fltarr(dims.nx,dims.ny)
  ilat=ilon
  for ix=0,dims.nx-1 do ilat[ix,*]=dims.lat
  for iy=0,dims.ny-1 do ilon[*,iy]=dims.lon


;STEP 1 - MASK OUT FIRST 24 H
  slp[*,*,0:it_zero-1]=!values.f_nan
;MASK OUT REGION WITHIN 8*DX OF DOMAIN EDGE
  iedge=8
  slp[[indgen(iedge),dims.nx-1-indgen(iedge)],*,*]=!values.f_nan
  slp[*,[indgen(iedge),dims.ny-1-indgen(iedge)],*]=!values.f_nan


;STEP 2 - SMOOTHING
  ;SMOOTH IN TIME
  for i=0,3 do $
      slp = smooth(temporary(slp),[0,0,nt_sm],/edge_truncate,/nan)
  ;SMOOTH IN SPACE
  for i=0,29 do $
    slp=smooth(temporary(slp),[nx_sm,nx_sm,0],/edge_truncate,/nan)


;STEP 3 - CALCULATE ANOMALY
; SUBTRACT ALL-DOMAIN-MEAN AT TIME T
  for it=0,dims.nt-1 do $
    slp[*,*,it] -= mean(reform(slp[*,*,it]),/nan,/double)
; SUBTRACT ALL-TIME MEAN SPATIAL PATTERN
  slp -= rebin(mean(slp,dimension=3,/nan,/double),[dims.nx,dims.ny,dims.nt])


;PRINT OUT A COPY OF THE FIELD USED TO ID VORTEX
  write_sing_ncvar,'slp_filter.nc',slp,'slp',$
    dimtag1='lon',dimtag2='lat',dimtag3='time'
exit

  pmin=pmin_thresh1 ; minimum all-time SLP threshold
  it_pmin=-1

  for it=it_zero-1,dims.nt-1 do begin

  ;STEP 4 - MASK OUT POINTS BEYOND RADIUS THRESHOLD OF REAL TC

    t_diff=abs(time[it]-hurdat.jultim)
    ithd=(where(t_diff eq min(t_diff)))[0]

    ;USE HURDAT STORM MOTION AND PREDICT LOCATION FROM:
    ;  X1 = X_0 + CX_0*DT

    dt = (time[it]-hurdat.jultim[ithd])*86400d ; s

    cx = hurdat.motion_x[ithd] ; m/s
    cy = hurdat.motion_y[ithd] ; m/s

    hdlat = hurdat.lat[ithd] + cy*dt*m2deg
    hdlon = hurdat.lon[ithd] + cx*dt*m2deg/(cos(hdlat*!dtor))

    radius=sqrt( (ilon-hdlon)^2 + (ilat-hdlat)^2 )
    inan=where(radius ge radmax)

    islp=reform(slp[*,*,it])
    islp[inan]=!values.f_nan

  ;STEP 4 - IDENTIFY STRONGEST ALL-TIME LOW-PRESSURE CENTER IN MASKED FIELD

    ipmin=min(islp,loc,/nan)

    ifin=where(finite(islp),count)
    if count eq 0 then continue

    if ipmin le pmin then begin
      pmin=ipmin
      ind=array_indices([dims.nx,dims.ny],loc,/dimensions)
      ploc=[dims.lon[ind[0]],dims.lat[ind[1]]]
      it_pmin=it
    endif

  endfor ; it


;KILL IF VORTEX WAS NOT IDENTIFIED
  if it_pmin eq -1 then begin
    ;WRITE OUT AND END
    if keyword_set(write) then $
      write_tc_track,trackfile,locvort
    return,locvort
  endif

;ELSE SET ANCHOR VALUES
  locvort[0:1,it_pmin]=ploc
  locvort[2,it_pmin]=pmin


;STEP 5 - TRACK THE ABOVE-LOCATED CENTER FORWARD AND BACK IN TIME
  dt=(time[1]-time[0])*86400.d ; s
  max_deg = c_max * dt * m2deg ; maximum distance in degrees
  idt=1d/dt

 ;BACKWARDS IN TIME
  for it=it_pmin-1,it_zero,-1 do begin

    islp=reform(slp[*,*,it])

    ;CHECK WITHIN RADIUS CORRESPONDING TO MAX TRANSLATION SPEED
      radius=sqrt( (ilon-locvort[0,it+1])^2 + (ilat-locvort[1,it+1])^2 )
      inan=where(radius gt max_deg)
      islp[inan]=!values.f_nan
      ipmin=min(islp,loc,/nan)
      if ipmin le pmin_thresh2 then begin
        ind=array_indices([dims.nx,dims.ny],loc,/dimensions)
        locvort[0:1,it]=[dims.lon[ind[0]],dims.lat[ind[1]]]
        locvort[2,it]=ipmin 
      endif else begin
        locvort[*,it]=!values.f_nan
        break
      endelse

  endfor

  ;FORWARD IN TIME
  for it=it_pmin+1,dims.nt-1 do begin

    islp=reform(slp[*,*,it])

    ;CHECK WITHIN RADIUS CORRESPONDING TO MAX TRANSLATION SPEED
      radius=sqrt( (ilon-locvort[0,it-1])^2 + (ilat-locvort[1,it-1])^2 )
      inan=where(radius gt max_deg)
      islp[inan]=!values.f_nan
      ipmin=min(islp,loc,/nan)
      if ipmin le pmin_thresh2 then begin
        ind=array_indices([dims.nx,dims.ny],loc,/dimensions)
        locvort[0:1,it]=[dims.lon[ind[0]],dims.lat[ind[1]]]
        locvort[2,it]=ipmin
      endif else begin
        locvort[*,it]=!values.f_nan
        break
      endelse

  endfor

  ifin=where(finite(reform(locvort[0,*])),nt_cyclone)


;EXCLUDE SHORT-LIVED VORTICES
  if nt_cyclone le nt_min then begin
    locvort[*]=!values.f_nan
    if keyword_set(write) then $
      write_tc_track,trackfile,locvort
    return,locvort
  endif


;SMOOTH TRACKS
;  nt_sm=3
;  locvort=smooth(temporary(locvort),[0,nt_sm],/edge_truncate,/nan)


;STORM MOTION (m/s)
  lon=reform(locvort[0,*])
  lat=reform(locvort[1,*])
  motion_x = 111d3 * cos(lat*!pi/180) * deriv(time*24*3600,lon)
  motion_y = 111d3 *                    deriv(time*24*3600,lat)
  locvort[3,*]=motion_x
  locvort[4,*]=motion_y


;WRITE TO FILE
  if keyword_set(write) then $
    write_tc_track,trackfile,locvort

return,locvort

end
