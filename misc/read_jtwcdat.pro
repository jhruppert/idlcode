;
; Function to read JTWC best tracks a single storm.
;
; Input: storm (string)
; Return: structure with Best Track data (see below).
;
; James Ruppert
; 7/8/2019
function read_jtwcdat, storm

if storm eq 'haiyan' then ifil='bwp312013.dat'
hurdat='/work/06040/tg853394/stampede2/tc_stuff/'+ifil

;READ

;  grep="grep -iE '"+year+".*"+storm+"' "+hurdat
;  spawn,grep,line0,err
;  if keyword_set(err) then message,err
;  
;  split=strsplit(line0,/extract,',')
;  print,'READING HURDAT FOR STORM: ',strtrim(split[0],2),' ',strtrim(split[1],2)
;  nlin=fix(split[2])
;  
;  grep="grep -A "+strtrim(nlin,2)+" -iE '"+year+".*"+storm+"' "+hurdat
;  spawn,grep,lines,err
;  if keyword_set(err) then message,err

  spawn,'cat '+hurdat,lines,err
  nlin=n_elements(lines)

;EXTRACT

  jultim=dblarr(nlin)
  lon=fltarr(nlin)
  lat=lon
  wspd=lon
  pres=lon
  status=strarr(nlin)

  for i=0,nlin-1 do begin
    
    line=strsplit(lines[i],',',/extract)
    
    dat=strtrim(line[2],2) ; YYYYMMDD
    yy=strmid(dat,0,4)
    mm=strmid(dat,4,2)
    dd=strmid(dat,6,2)
;    dat=line[1] ; HHNN
    hh=strmid(dat,8,2)
;    nn=strmid(dat,3,2)
    jultim[i]=julday(mm,dd,yy,hh,00,00)

;    dat=line[2] ; L
    
    dat=line[10] ; TS
    status[i]=strtrim(dat,2)
    
    dat=line[6] ; LAT
    lat[i]=float( strmid(dat,0,5) )/10.
    if strmid(dat,0,1,/reverse) eq 'S' then lat[i]*=-1.
    
    dat=line[7] ; LON
    lon[i]=float( strmid(dat,0,6) )/10.
    if strmid(dat,0,1,/reverse) eq 'W' then lon[i]*=-1.
    
    dat=line[8] ; MAX SUSTAINED WIND
    wspd[i]=float(dat)
    
    dat=line[9] ; MIN PRESSURE (hPa)
    pres[i]=float(dat)
    
  endfor

  ;REMOVE EXTRA TIME STEPS USED FOR RADIUS REPORTING
    iuq=uniq(jultim)
    nlin=n_elements(iuq)
    jultim=jultim[iuq]
    lon=lon[iuq]
    lat=lat[iuq]
    wspd=wspd[iuq]
    pres=pres[iuq]
    status=status[iuq]

;STORM MOTION (m/s)
motion_x = 111d3 * cos(lat*!pi/180) * deriv(lon) / (deriv(jultim) * 24*3600)
motion_y = 111d3 *                    deriv(lat) / (deriv(jultim) * 24*3600)

dat=create_struct('jultim',jultim,'lon',lon,'lat',lat,'wspd',wspd,'pres',pres,'motion_x',motion_x,'motion_y',motion_y,'status',status)

;CHECK
;for i=0,nlin-1 do begin
;  caldat,jultim[i],mm,dd,yy,hh,nn
;  print,i,' ',yy,mm,dd,' ',hh,nn,' ',lon[i],' ',lat[i],' ',wspd[i],' ',pres[i],' ',status[i]
;endfor
;exit

return, dat

end

