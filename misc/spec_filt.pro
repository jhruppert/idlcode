;--------------------------------------------------------------------------------------------
;  This routine filters a time-dependent field, as in f(t,x,y), leaving only the desired mode
;  of variability. The Butterwortrh filter is used (IDL routine "butterworth") with a user-
;  specified number of weights.
;  
;  
;  spec_filt, INFIELD, OUTFIELD, NW=NW, LOC=LOC, /DIURNCOMP, NPERDAY=NPERDAY, /NOFILT
;  
;  
;  INFIELD:   Input 1-, 2-, or 3-D array (function of TIME, [X, [and Y,]]: time must be first)
;             to be filtered.
;  OUTFIELD:  Output array of the same format as input, with only the desired mode of
;             variability remaining.
;  
;  
;  OPTIONAL KEYWORDS:
;  
;  NW:        Number of weights to use for the Butterworth filter.
;  LOC:       Location (index) of spectral variability peak that is to be grabbed.
;  DIURNCOMP: Set this keyword when LOC corresponds to the location of the diurnal peak in
;             order to return an array of the form f(n,nx,ny), where n corresponds to the
;             number of timesteps that make up a full 24 hours (see NPERDAY below) and the
;             0,*,* elements correspond with the average pattern at the first diurnal time
;             step, and so on. NPERDAY must be set to use this option.
;  NPERDAY    Also when LOC corresponds to the location of the diurnal peak, set this equal
;             to the number of timesteps per day when T_AVE is set (e.g., set to 4 for data
;             where values correspond to 00, 06, 12, and 18Z throughout the period). NT does
;             not have to complete a full day for the last day (i.e., it is okay if the last
;             few data times correspond to only a portion of a day). The code, however, is
;             designed under the assumption that the data at the beginning of the array
;             correspond to the first time step in the day (e.g., 00Z). Simple code
;             modifications might allow for more flexibility.
;  NOFILT:    Set this keyword to turn off filtering (i.e., if you only want field
;             diurnally composited...set /diurncomp).
;  
;  
;  James Ruppert (ruppert@atmos.colostate.edu)
;  3/10/2011
;--------------------------------------------------------------------------------------------
pro spec_filt,infield,outfield,NW=nw,LOC=loc,DIURNCOMP=diurncomp,NPERDAY=nperday,NOFILT=nofilt,remove_mean=remove_mean
  
  
  on_error,3  ; Return to this program if an error occurs.
  
  
  infield=reform(infield) ; Remove any null dimensions
  
  
  ndims=size(infield,/n_dimensions) & dims=size(infield,/dimensions)
  if ndims eq 0 then message,'Input field must be a vector array with at least 1 dimension!'
  type=size(infield,/type)
  if type ne 4 and type ne 5 then message,'Input field must be either a float or double-precision array!'
  nt=dims[ndims-1]
  if ndims ge 2 then nx=dims[0] & if ndims eq 3 then ny=dims[1]
  
  
  if ~keyword_set(nofilt) then begin
    
    ;CREATE FILTER
      
      if ~keyword_set(nw) then begin
        nw=9 & print,'Used 9 weights for Butterworth filter'
      endif
      
      butter_temp=butterworth(nt,cutoff=2,order=nw,/origin)
      butter_band=shift(butter_temp, -1*((nt-1)/2-loc) )
;      butter_band=fltarr(nt)
;      for i=0,nt-1 do begin
;        j=i-((nt-1)/2-loc)
;        if j ge 0 then butter_band[j]=butter_temp[i]
;      endfor
    
    ;APPLY FILTER
      
      if ndims eq 1 then begin
        outfield=fltarr(nt)
        outfield=fft(fft(infield,-1)*butter_band,1)
      endif else if ndims eq 2 then begin
        outfield=fltarr(nx,nt)
        for ix=0,nx-1 do begin
          iin=reform(infield[ix,*])
          outfield[ix,*]=fft(fft(iin,-1)*butter_band,1)
        endfor
      endif else if ndims eq 3 then begin
        outfield=fltarr(nx,ny,nt)
        for ix=0,nx-1 do $
          for iy=0,ny-1 do begin
            iin=reform(infield[ix,iy,*])
            outfield[ix,iy,*]=fft(fft(iin,-1)*butter_band,1)
          endfor
      endif
    
  endif else outfield=infield
  
  
  ;AVERAGE DAILY CYCLE IN TIME
    
    if keyword_set(diurncomp) then begin

      if ~keyword_set(nperday) then $
        message,'Must set NPERDAY when using the T_AVE option!' $
      else $
        if n_elements(nperday) gt 1 or fix(nperday[0]) ne nperday[0] then $
          message,'NPERDAY must be a single-element whole number value!'
      
      extra = nt mod nperday ; Determine number of data points (at end of array) that don't
                             ; complete a day. Here we assume that these "extra" data points
                             ; correspond to the first daily time, the second daily time (and
                             ; so on...), respectively.
      
      ndays=fix((nt-extra)/nperday) ; number of complete days of data
      
      t_ave=indgen(ndays)*nperday
      
      if ndims eq 1 then $
        diurn_ave=fltarr(nperday) $
      else if ndims eq 2 then $
        diurn_ave=fltarr(nx,nperday) $
      else if ndims eq 3 then $
        diurn_ave=fltarr(nx,ny,nperday)
      
      if keyword_set(remove_mean) then begin
stop
        if ndims eq 1 then $
          mean_field=mean(outfield,/nan) $
        else if ndims ge 2 then $
          mean_field=mean(outfield,/nan,dimension=1)

      endif else mean_field=0
      
      for itime=0,nperday-1 do begin
        
;ONLY USING FULL DAYS
;        if itime lt extra then $
;          t_ind=[t_ave+itime,(ndays)*nperday+itime] $
;        else $
          t_ind=t_ave+itime
        
        if ndims eq 1 then $
          diurn_ave[itime]=mean(outfield[t_ind],/nan)-mean_field $
        else if ndims eq 2 then $
          diurn_ave[*,itime]=mean(outfield[*,t_ind],/nan,dimension=2)-mean_field $
        else if ndims eq 3 then $
          diurn_ave[*,*,itime]=mean(outfield[*,*,t_ind],/nan,dimension=3)-mean_field
        
      endfor
      
      outfield=diurn_ave
      
    endif
  
  
end
