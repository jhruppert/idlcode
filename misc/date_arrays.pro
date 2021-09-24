;---------------------------------------------------------------------------------------------------------------
;  Subroutine to produce a set of time/date arrays based on a start and end date and a frequency (1 day or
;  less).
;  
;  Arrays are fed in in the order they are used below.
;  
;  James Ruppert 5/25/2013
;---------------------------------------------------------------------------------------------------------------
pro date_arrays, date_start, date_end, interval, jul_dates=jul_dates, utc=utc

  ;dates can either be strings of the form MMDDYYYY_HH
  ;  or the Julian start and end dates
  ;interval should be in hours
  ;utc will be a string array of the form YYYYMMDDHH
  
  type=(size(date_start))[1] ; 5 for double, 7 for string
  
  if type eq 5 then begin
    jul_start=date_start & jul_end=date_end
  endif else if type eq 7 then begin
    m1=strmid(date_start,0,2) & d1=strmid(date_start,2,2) & y1=strmid(date_start,4,4) & h1=strmid(date_start,9,2)
    m2=strmid(date_end,0,2) & d2=strmid(date_end,2,2) & y2=strmid(date_end,4,4) & h2=strmid(date_end,9,2)
    jul_start=julday(m1,d1,y1,h1) & jul_end=julday(m2,d2,y2,h2)
  endif else message,'"DATE_START/END": UNKNOWN TYPE'
  
  jul_dates=timegen(start=jul_start,final=jul_end,step_size=interval,units='hours')
  caldat,jul_dates,mm,dd,yy,hh
  utc=strmid(strtrim(yy,2),0,4)+string(mm,format='(i2.2)')+string(dd,format='(i2.2)')+string(hh,format='(i2.2)')

end