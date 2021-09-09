; 
; Identify vortex center for Maria TC WRF output based on ABS VOR.
;
; Set /write to print out hourly locations to an ascii file.
;
; James Ruppert
; 4/22/19
; 
function maria_vortex_locate, avor_sav, dims, write=write

;trackfil='/scratch/06040/tg853394/tc/output/static_nest/maria/ctl/avor700_track.txt'
trackfil='/scratch/06040/tg853394/tc/output/redux/maria/ctl/avor700_track.txt'

specs=size(avor_sav,/dimensions)
nt=specs[2]

;SMOOTHING AND LEVEL SELECTION NOW TAKE PLACE IN PARENT SCRIPT

  avor=avor_sav

;ABSOLUTE VORTICITY SMOOTHED FIRST

;  ixsmooth=round(111./3) ; 1-degree smoothing, run twice
;  ismooth=[ixsmooth,ixsmooth,0]
;  for i=1,2 do $
;    avor=smooth(temporary(avor),ismooth,/edge_truncate,/nan)

;VORTEX CENTER LAT/LON

  locvort=fltarr(2,nt)

  for it=0,nt-1 do begin
    temp=max(reform(avor[*,*,it]),loc,/nan)
    ind=array_indices([dims.nx,dims.ny],loc,/dimensions)
    locvort[0,it]=dims.lon[ind[0]]
    locvort[1,it]=dims.lat[ind[1]]
  endfor

;SMOOTH TRACKS
;  for i=0,1 do $
;    locvort=smooth(temporary(locvort),[0,3],/edge_truncate,/nan)

;WRITE TO FILE

  if keyword_set(write) then begin

spawn,'rm '+trackfil,out,err
openw,1,trackfil

  extra=nt mod 8
  nlin=(nt-extra)/8

  printf,1,'LON (n='+strtrim(nt,2)+'):'

  for il=0,nlin-1 do $
    printf,1,strjoin(strtrim(float(reform(locvort[0,il*8:il*8+7])),2),',')+', & ! JHR'
  if extra eq 1 then $
    printf,1,strtrim(float(reform(locvort[0,nt-1])),2)
  if extra gt 1 then $
    printf,1,strjoin(strtrim(float(reform(locvort[0,nlin*8:nlin*8+extra-1])),2),',')

  printf,1,' '
  printf,1,'LAT (n='+strtrim(nt,2)+'):'
  
  for il=0,nlin-1 do $
    printf,1,strjoin(strtrim(float(reform(locvort[1,il*8:il*8+7])),2),',')+', & ! JHR'
  if extra eq 1 then $
    printf,1,strtrim(float(reform(locvort[1,nt-1])),2)
  if extra gt 1 then $
    printf,1,strjoin(strtrim(float(reform(locvort[1,nlin*8:nlin*8+extra-1])),2),',')

  printf,1,' '
  printf,1,'JULDAY (n='+strtrim(nt,2)+'):'

  t0=julday(1,1,2017,0,0,0)
  juldays=fltarr(nt)
  for it=0,nt-1 do juldays[it] = float( dims.time[it] - t0 )

  for il=0,nlin-1 do $
    printf,1,strjoin(strtrim(float(reform(juldays[il*8:il*8+7])),2),',')+', & ! JHR'
  if extra eq 1 then $
    printf,1,strtrim(float(reform(juldays[nt-1])),2)
  if extra gt 1 then $
    printf,1,strjoin(strtrim(float(reform(juldays[nlin*8:nlin*8+extra-1])),2),',')

close,1

  endif

return,locvort

end
