;
; Function to read Hurdat data for a single storm.
;
; See bottom for NHC's beta on HURDAT file.
; 
; Input: storm (string), year (string)
; Return: structure with HURDAT data (see below).
;
; James Ruppert
; 3/22/2019
function read_hurdat, storm, year

hurdat='/work/06040/tg853394/stampede2/tc_stuff/hurdat2-1851-2017-050118.txt'

;READ

  grep="grep -iE '"+year+".*"+storm+"' "+hurdat
  spawn,grep,line0,err
  if keyword_set(err) then message,err
  
  split=strsplit(line0,/extract,',')
  print,'READING HURDAT FOR STORM: ',strtrim(split[0],2),' ',strtrim(split[1],2)
  nlin=fix(split[2])
  
  grep="grep -A "+strtrim(nlin,2)+" -iE '"+year+".*"+storm+"' "+hurdat
  spawn,grep,lines,err
  if keyword_set(err) then message,err

;EXTRACT

  jultim=dblarr(nlin)
  lon=fltarr(nlin)
  lat=lon
  wspd=lon
  pres=lon
  status=strarr(nlin)

  for i=0,nlin-1 do begin
    
    line=strsplit(lines[i+1],',',/extract)
    
    dat=line[0] ; YYYYMMDD
    yy=strmid(dat,0,4)
    mm=strmid(dat,4,2)
    dd=strmid(dat,6,2)
    dat=line[1] ; HHNN
    hh=strmid(dat,1,2)
    nn=strmid(dat,3,2)
    jultim[i]=julday(mm,dd,yy,hh,nn,00)
    
    dat=line[2] ; L
    
    dat=line[3] ; TS
    status[i]=strtrim(dat,2)
    
    dat=line[4] ; LAT
    lat[i]=float( strmid(dat,0,5) )
    if strmid(dat,0,1,/reverse) eq 'S' then lat[i]*=-1.
    
    dat=line[5] ; LON
    lon[i]=float( strmid(dat,0,6) )
    if strmid(dat,0,1,/reverse) eq 'W' then lon[i]*=-1.
    
    dat=line[6] ; MAX SUSTAINED WIND
    wspd[i]=float(dat)
    
    dat=line[7] ; MIN PRESSURE (hPa)
    pres[i]=float(dat)
    
  endfor

;STORM MOTION (m/s)
motion_x = 111d3 * cos(lat*!pi/180) * deriv(jultim*24*3600,lon)
motion_y = 111d3 *                    deriv(jultim*24*3600,lat)

dat=create_struct('jultim',jultim,'lon',lon,'lat',lat,'wspd',wspd,'pres',pres,'motion_x',motion_x,'motion_y',motion_y,'status',status)

;CHECK
;for i=0,nlin-1 do begin
;  caldat,jultim[i],mm,dd,yy,hh,nn
;  print,i,' ',yy,mm,dd,' ',hh,nn,' ',lon[i],' ',lat[i],' ',wspd[i],' ',pres[i],' ',status[i]
;endfor

return, dat

end

; Clipped from: https://www.nhc.noaa.gov/data/hurdat/hurdat2-format-atlantic.pdf
;
; There are two types of lines of data in the new format: the header line and the data lines. The format is comma delimited to maximize its ease in use.
; The header line has the following format:

;   AL092011, IRENE, 39,
;   1234567890123456789012345768901234567

; AL (Spaces 1 and 2) – Basin – Atlantic
; 09 (Spaces 3 and 4) – ATCF cyclone number for that year
; 2011 (Spaces 5-8, before first comma) – Year
; IRENE (Spaces 19-28, before second comma) – Name, if available, or else “UNNAMED”
; 39 (Spaces 34-36) – Number of best track entries – rows – to follow

; Notes:

; 1) Cyclone number: In HURDAT2, the order cyclones appear in the file is determined by the date/time of the first tropical or subtropical cyclone record
; in the best track. This sequence may or may not correspond to the ATCF cyclone number. For example, the 2011 unnamed tropical storm AL20 which formed on 1
; September, is sequenced here between AL12 (Katia – formed on 29 Aug) and AL13 (Lee – formed on 2 September). This mismatch between ATCF cyclone
; number and the HURDAT2 sequencing can occur if post-storm analysis alters the relative genesis times between two cyclones. In addition, in 2011 it became
; practice to assign operationally unnamed cyclones ATCF numbers from the end of the list, rather than insert them in sequence and alter the ATCF numbers of
; cyclones previously assigned.
;
; 2) Name: Tropical cyclones were not formally named before 1950 and are thus referred to as “UNNAMED” in the database. Systems that were added into the
; database after the season (such as AL20 in 2011) also are considered “UNNAMED”. Non-developing tropical depressions formally were given names (actually
; numbers, such as “TEN”) that were included into the ATCF b-decks starting in 2003. Non-developing tropical depressions before this year are also
; referred to as “UNNAMED”.

; The remaining rows of data in the new format are the data lines. These have the following format:

; 20110828, 0935, L, TS, 39.4N, 74.4W, 60, 959, 230, 280, 160, 110, 150, 150, 80, 30, 0, 0, 0, 0,
; 123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
; 
; 2011 (Spaces 1-4) – Year
; 08 (Spaces 5-6) – Month
; 28 (Spaces 7-8, before 1st comma) – Day 
; 09 (Spaces 11-12) – Hours in UTC (Universal Time Coordinate)
; 35 (Spaces 13-14, before 2nd comma) – Minutes
; L (Space 17) – Record identifier (see notes below)
;   C – Closest approach to a coast, not followed by a landfall
;   G – Genesis
;   I – An intensity peak in terms of both pressure and wind
;   L – Landfall (center of system crossing a coastline)
;   P – Minimum in central pressure
;   R – Provides additional detail on the intensity of the cyclone when rapid changes are underway
;   S – Change of status of the system
;   T – Provides additional detail on the track (position) of the cyclone
;   W – Maximum sustained wind speed
; TS (Spaces 20-21, before 3rd comma) – Status of system. Options are:
;   TD – Tropical cyclone of tropical depression intensity (< 34 knots)
;   TS – Tropical cyclone of tropical storm intensity (34-63 knots)
;   HU – Tropical cyclone of hurricane intensity (> 64 knots)
;   EX – Extratropical cyclone (of any intensity)
;   SD – Subtropical cyclone of subtropical depression intensity (< 34 knots)
;   SS – Subtropical cyclone of subtropical storm intensity (> 34 knots)
;   LO – A low that is neither a tropical cyclone, a subtropical cyclone, nor an extratropical cyclone (of any intensity)
;   WV – Tropical Wave (of any intensity)
;   DB – Disturbance (of any intensity) 
; 39.4 (Spaces 24-27) – Latitude
; N (Space 28, before 4th comma) – Hemisphere – North or South
; 74.4 (Spaces 31-35) – Longitude
; W (Space 36, before 5th comma) – Hemisphere – West or East
; 60 (Spaces 39-41, before 6th comma) – Maximum sustained wind (in knots)
; 959 (Spaces 44-47, before 7th comma) – Minimum Pressure (in millibars)
; 230 (Spaces 50-53, before 8th comma) – 34 kt wind radii maximum extent in northeastern quadrant (in nautical miles)
; 280 (Spaces 56-59, before 9th comma) – 34 kt wind radii maximum extent in southeastern quadrant (in nautical miles)
; 160 (Spaces 62-65, before 10th comma) – 34 kt wind radii maximum extent in southwestern quadrant (in nautical miles)
; 110 (Spaces 68-71, before 11th comma) – 34 kt wind radii maximum extent in northwestern quadrant (in nautical miles)
; 150 (Spaces 74-77, before 12th comma) – 50 kt wind radii maximum extent in northeastern quadrant (in nautical miles)
; 150 (Spaces 80-83, before 13th comma) – 50 kt wind radii maximum extent in southeastern quadrant (in nautical miles)
; 80 (Spaces 86-89, before 14th comma) – 50 kt wind radii maximum extent in southwestern quadrant (in nautical miles)
; 30 (Spaces 92-95, before 15th comma) – 50 kt wind radii maximum extent in northwestern quadrant (in nautical miles)
; 0 (Spaces 98-101, before 16th comma) – 64 kt wind radii maximum extent in northeastern quadrant (in nautical miles)
; 0 (Spaces 104-107, before 17th comma) – 64 kt wind radii maximum extent in southeastern quadrant (in nautical miles)
; 0 (Spaces 110-113, before 18th comma) – 64 kt wind radii maximum extent in southwestern quadrant (in nautical miles)
; 0 (Spaces 116-119, before 19th comma) – 64 kt wind radii maximum extent in northwestern quadrant (in nautical miles) 
