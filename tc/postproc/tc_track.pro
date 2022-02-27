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
;
;  Step 2 is important given a large domain where there may be multiple vortices
;    and/or the target vortex is weak.
;
; Set trackdir to directory where storm track(s) should be written to.
;
; James Ruppert
; 12/20/21
; 
pro tc_track, slp_sav, dims, time, trackdir=trackdir, time_mask=time_mask

;SETTINGS
;  nx_sm=round(0.5*1./(dims.lat[1]-dims.lat[0])) ; Half-degree smoothing in space
  nx_sm=9    ; Following Chavas 2013 (smoothing run x30 times)
  nt_sm=3    ; temporal smoothing (n time steps)
  lmin=1e3   ; Minimum distance threshold [km] between storms
  it_zero=24 ; set all time steps prior to this one to NaN
  c_max=20.  ; a very liberal maximum allowable translation speed for tracking
  pmin_thresh1=-4. ; lifetime-minimum pressure as threshold for identifying a cyclone
;pmin_thresh1=-3.
  pmin_thresh2=-3. ; after storm ID, sequential time steps must retain center of this strength
pmin_thresh2=pmin_thresh1;-2.
  nt_min=24. ; minimum time steps for a cyclone to be retained

;CONSTANTS
  m2deg=1./(111e3)
  ;THE BELOW ARE FOR TRACKING VORTEX SUBJECT TO C_MAX
  dt=(time[1]-time[0])*86400.d ; s
  idt=1d/dt
  rmax_track = c_max * dt * m2deg ; maximum single-time-step displacement [degrees]

;SAVE SLP
  slp=slp_sav

;2D LON/LAT
  ilon=fltarr(dims.nx,dims.ny)
  ilat=ilon
  for ix=0,dims.nx-1 do ilat[ix,*]=dims.lat
  for iy=0,dims.ny-1 do ilon[*,iy]=dims.lon


;STEP 1 - MASK OUT FIRST 24 H
  if keyword_set(time_mask) then slp[*,*,0:it_zero-1]=!values.f_nan
;MASK OUT REGION WITHIN 8*DX OF DOMAIN EDGE
  iedge=81 ; about 2 deg
  slp[[indgen(iedge),dims.nx-1-indgen(iedge)],*,*]=!values.f_nan
  slp[*,[indgen(iedge),dims.ny-1-indgen(iedge)],*]=!values.f_nan


;STEP 2 - SMOOTHING
  ;SMOOTH IN TIME
  for i=0,3 do $
      slp = smooth(temporary(slp),[0,0,nt_sm],/edge_truncate,/nan)
  ;SMOOTH IN SPACE
  for i=0,29 do $
    slp=smooth(temporary(slp),[nx_sm,nx_sm,0],/edge_truncate,/nan)


;;REPEAT: MASK OUT REGION WITHIN 8*DX OF DOMAIN EDGE
;  iedge=81 ; about 2 deg
;  slp[[indgen(iedge),dims.nx-1-indgen(iedge)],*,*]=!values.f_nan
;  slp[*,[indgen(iedge),dims.ny-1-indgen(iedge)],*]=!values.f_nan


;STEP 3 - CALCULATE ANOMALY
; SUBTRACT ALL-DOMAIN-MEAN AT TIME T
  for it=0,dims.nt-1 do $
    slp[*,*,it] -= mean(reform(slp[*,*,it]),/nan,/double)
; SUBTRACT ALL-TIME MEAN SPATIAL PATTERN
  slp -= rebin(mean(slp,dimension=3,/nan,/double),[dims.nx,dims.ny,dims.nt])

;PRINT A COPY
;  write_sing_ncvar,'slp_filter.nc',slp,'slp',$
;    dimtag1='lon',dimtag2='lat',dimtag3='time'


;STEP 4 - IDENTIFY ALL CENTERS THAT MEET 1ST P-THRESHOLD
  ilow=where(slp le pmin_thresh1,complement=inan)
  islp1=slp
  islp1[inan]=!values.f_nan
  ifin=where(finite(islp1),n_finite)


;STEP 5 - TRACK CENTER THROUGH TIME AND MASK OUT SURROUNDINGS, LOOPING OVER N-CENTERS

  ;START LOOP, KILL IF THRESHOLD NOT MET
  if n_finite gt 1 then do_storm_loop=1 else return

  islp2=islp1 ; New array to mask out identified storm tracks

  istorm=1 ; Storm label

  while do_storm_loop do begin

    ipmin=min(islp2,loc,/nan)

    ind=array_indices([dims.nx,dims.ny,dims.nt],loc,/dimensions)
    ploc=[dims.lon[ind[0]],dims.lat[ind[1]]]
    it_pmin=ind[2]

    storm_loc=fltarr(2,dims.nt) ; lon[deg], lat ; temp var to hold location values
    storm_loc[*]=!values.f_nan
    storm_loc[*,it_pmin]=ploc

    ;ISOLATE VORTEX OBJECT BY SETTING NEIGHBORING POINTS TO NAN AT TIME STEP
      radius=sqrt( (ilon-storm_loc[0,it_pmin])^2 + (ilat-storm_loc[1,it_pmin])^2 ) 
      inan=where(radius le lmin)
      slp_it=reform(islp2[*,*,it_pmin])
      slp_it[inan]=!values.f_nan
      islp2[*,*,it_pmin]=slp_it

    ;TRACK BACKWARDS IN TIME
    for it=it_pmin-1,it_zero,-1 do begin

      radius=sqrt( (ilon-storm_loc[0,it+1])^2 + (ilat-storm_loc[1,it+1])^2 ) 
  
      ;CHECK WITHIN RADIUS CORRESPONDING TO MAX TRANSLATION SPEED
        inan=where(radius gt rmax_track)
        slp_it=reform(islp2[*,*,it])
        slp_it[inan]=!values.f_nan

      ipmin=min(slp_it,loc,/nan)
      if ipmin le pmin_thresh2 then begin
        ind=array_indices([dims.nx,dims.ny],loc,/dimensions)
        ;BREAK IF REACHES DOMAIN EDGE
          if (ind[0] le iedge) or (ind[0] ge dims.nx-1-iedge) or $
             (ind[1] le iedge) or (ind[1] ge dims.ny-1-iedge) then break
        storm_loc[*,it]=[dims.lon[ind[0]],dims.lat[ind[1]]]
      ;SET NEIGHBORING POINTS TO NAN AT TIME STEP
        inan=where(radius le lmin)
        slp_it=reform(islp2[*,*,it])
        slp_it[inan]=!values.f_nan
        islp2[*,*,it]=slp_it
      endif else break
  
    endfor
  
    ;TRACK FORWARD IN TIME
    for it=it_pmin+1,dims.nt-1 do begin
  
      radius=sqrt( (ilon-storm_loc[0,it-1])^2 + (ilat-storm_loc[1,it-1])^2 ) 

      ;CHECK WITHIN RADIUS CORRESPONDING TO MAX TRANSLATION SPEED
        inan=where(radius gt rmax_track)
        slp_it=reform(islp2[*,*,it])
        slp_it[inan]=!values.f_nan

      ipmin=min(slp_it,loc,/nan)
      if ipmin le pmin_thresh2 then begin
        ind=array_indices([dims.nx,dims.ny],loc,/dimensions)
        ;BREAK IF REACHES DOMAIN EDGE
          if (ind[0] le iedge) or (ind[0] ge dims.nx-1-iedge) or $
             (ind[1] le iedge) or (ind[1] ge dims.ny-1-iedge) then break
        storm_loc[*,it]=[dims.lon[ind[0]],dims.lat[ind[1]]]
      ;SET NEIGHBORING POINTS TO NAN AT TIME STEP
        inan=where(radius le lmin)
        slp_it=reform(islp2[*,*,it])
        slp_it[inan]=!values.f_nan
        islp2[*,*,it]=slp_it
     endif else break
  
    endfor

    ;CHECK IF MINIMUM DURATION IS REACHED
    ifin=where(finite(reform(storm_loc[0,*])),nt_storm)

    if nt_storm ge nt_min then begin

      locvort=fltarr(4,dims.nt) ; lon[deg], lat[deg], c_x[m/s], c_y[m/s]
      locvort[*]=!values.f_nan

      locvort[[0,1],*]=storm_loc

      ;STORM MOTION (m/s)
        lon=reform(locvort[0,*])
        lat=reform(locvort[1,*])
        motion_x = 111d3 * cos(lat*!pi/180) * deriv(time*24*3600,lon)
        motion_y = 111d3 *                    deriv(time*24*3600,lat)
        locvort[2,*]=motion_x
        locvort[3,*]=motion_y

      trackfile=trackdir+'track_'+strtrim(istorm,2)+'.nc'
      write_tc_track,trackfile,locvort

      print,'ISTORM: ',istorm
      istorm+=1

    endif

    ;CHECK FOR NEW STORM CENTERS
    ilow=where(islp2 le pmin_thresh1,complement=inan)
    islp2[inan]=!values.f_nan
    ifin=where(finite(islp2),n_finite)

    ;CONTINUE LOOP IF THRESHOLD MET
    if n_finite eq 0 then do_storm_loop=0

  endwhile

;SMOOTH TRACKS
;  nt_sm=3
;  locvort=smooth(temporary(locvort),[0,nt_sm],/edge_truncate,/nan)


end
