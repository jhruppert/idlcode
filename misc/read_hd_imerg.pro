;
; Program to read raw IMERG data from individual time files in HDF format.
;
; XXX Writes out a netcdf file with all trmm data.
;
; James Ruppert
; 12 Aug 2020
;
pro read_hd_imerg

  rainvar='precipitationCal'

  ;Reads TRMM 3B42V7 data, units of MM / H

  config_dir,dirs=dirs

  read_dir=dirs.scdir+'myanmar/imerg/data/'
;  outfil=read_dir+'../trmm_3b42v7_01011998H00:00:00-01012014H00:00:00.nc'

;GLOBAL SUBSET

;  sellon=[48.7317,121.268]
;  sellat=[-20.2682,20.2682]

;TIME SPAN

;  d_start=julday(1,1,1998,0,0,0)
;  d_end=julday(1,1,2014,0,0,0)
;  d_start=julday(10,1,2011,0,0,0)
;  d_end=julday(11,16,2011,0,0,0)
;  times=timegen(start=d_start,final=d_end,step_size=3,units='H')
;  nt=n_elements(times)
;  caldat,times,mm,dd,yyall
;  yrs=yyall[uniq(yyall)]
;  nyrs=n_elements(yrs)

;FILE LIST

  spawn,'ls '+read_dir+'3B-HHR*nc4',files,err,count=nt
  if err then message,err

;SPECS FROM FIRST FILE

  file=files[0]
  lat=read_nc_var(file,'lat')
  lon=read_nc_var(file,'lon')
  nx=n_elements(lon)
  ny=n_elements(lat)

;READ AND WRITE

;for iyr=0,nyrs-1 do begin
;for iyr=7,nyrs-1 do begin
;for iyr=0,2 do begin

;  subset=where(yyall eq yrs[iyr],nt)

  rain=fltarr(nx,ny,nt)
  time=dblarr(nt)

print,nt

  for it=0,nt-1 do begin
;  for it=0,5 do begin
if (it mod 100) eq 0 then print,it
;    caldat,times[subset[it]],mm,dd,yy,hh
;    caldat,times[it],mm,dd,yy,hh
  
;    mm=string(mm,format='(i2.2)')
;    dd=string(dd,format='(i2.2)')
;    yy=string(yy,format='(i4.4)')
;    hh=string(hh,format='(i2.2)')
  
;    timtag=yy+mm+dd+'.'+hh
;    print,'YYYYMMDD.HH: ',timtag
  
  ;  file=read_dir+header+yy+mm+dd+'.'+hh+tail
;    loc=where(strmatch(files[subset],+'*'+timtag+'*') eq 1)
;    file=files[subset[loc[0]]]
    file=files[it]
;    if ~strmatch(file,'*'+timtag+'*') then message,'Something wrong with file indexing'

    ;TIME STAMP
    itim=read_nc_var(file,'time')
    time[it]=julday(1,1,1970,0,0,itim)

    ;READ VAR
    tmp=read_nc_var(file,rainvar)
;    nan=where(tmp eq -9999.9,count1)
;    print,'n-NAN = ',count1
;    nan=where(~finite(tmp),count2)
;    print,'n-NAN = ',count2

    ;SET NANS
  ;    nan=where(tmp lt -1000,count)
  ;    if count gt 0 then tmp[nan]=!values.f_nan
  
    rain[*,*,it]=transpose(tmp)
  
  endfor ; it

;outfil_ann=outfil+'_'+yy+'.nc'

;write_sing_ncvar,outfil,rain,'trmm_rain_rate_mmph'

;endfor ; iyr

print,'Done!'
end
