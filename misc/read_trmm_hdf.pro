;
; Program to read raw TRMM 3B42V7 data from individual time files in HDF format.
;
; Writes out a netcdf file with all trmm data.
;
; James Ruppert
; 6 Nov 2018
;
pro read_trmm_hdf

  ;Reads TRMM 3B42V7 data, units of MM / H

  config_dir,dirs=dirs

;  read_dir=dirs.wkdir+'trmm/'
  read_dir=dirs.scdir+'trmm/raw/'
;  outfil=read_dir+'trmm_3b42.nc'
;  outfil=read_dir+'trmm_3b42_global';.nc'
  outfil=read_dir+'../trmm_3b42v7_01011998H00:00:00-01012014H00:00:00.nc'

  header='3B42.'

;GLOBAL SUBSET

;  sellon=[48.7317,121.268]
;  sellat=[-20.2682,20.2682]

;TIME SPAN

  d_start=julday(1,1,1998,0,0,0)
  d_end=julday(1,1,2014,0,0,0)
;  d_start=julday(10,1,2011,0,0,0)
;  d_end=julday(11,16,2011,0,0,0)
  times=timegen(start=d_start,final=d_end,step_size=3,units='H')
  nt=n_elements(times)
  caldat,times,mm,dd,yyall
  yrs=yyall[uniq(yyall)]
  nyrs=n_elements(yrs)
print,yrs

;FILE LIST

  spawn,'ls '+read_dir+header+'*',files,err
  if err then message,err

;SPECS FROM FIRST FILE

  file=files[0];read_dir+header+'20110901.00'+tail
  fid=hdf_sd_start(file,/read)
  vid=hdf_sd_select(fid,0) ; 0 is the index of 'precipitation'
  hdf_sd_getdata,vid,tmp
  hdf_sd_end,fid
  specs=size(tmp,/dimensions)
  nx=specs[1]
  ny=specs[0]


;READ AND WRITE

;for iyr=0,nyrs-1 do begin
;for iyr=7,nyrs-1 do begin
;for iyr=0,2 do begin

;  subset=where(yyall eq yrs[iyr],nt)

  rain=fltarr(nx,ny,nt)

  for it=0,nt-1 do begin
;  for it=0,5 do begin
  
;    caldat,times[subset[it]],mm,dd,yy,hh
    caldat,times[it],mm,dd,yy,hh
  
    mm=string(mm,format='(i2.2)')
    dd=string(dd,format='(i2.2)')
    yy=string(yy,format='(i4.4)')
    hh=string(hh,format='(i2.2)')
  
    timtag=yy+mm+dd+'.'+hh
    print,'YYYYMMDD.HH: ',timtag
  
  ;  file=read_dir+header+yy+mm+dd+'.'+hh+tail
;    loc=where(strmatch(files[subset],+'*'+timtag+'*') eq 1)
;    file=files[subset[loc[0]]]
    file=files[it]
    if ~strmatch(file,'*'+timtag+'*') then message,'Something wrong with file indexing'

    ;READ VAR
      fid=hdf_sd_start(file,/read)
      vid=hdf_sd_select(fid,0) ; 0 is the index of 'precipitation'
      hdf_sd_getdata,vid,tmp
      hdf_sd_end,fid
  
    ;SET NANS
  ;    nan=where(tmp lt -1000,count)
  ;    if count gt 0 then tmp[nan]=!values.f_nan
  
    rain[*,*,it]=transpose(tmp)
  
  endfor ; it

;outfil_ann=outfil+'_'+yy+'.nc'

write_sing_ncvar,outfil,rain,'trmm_rain_rate_mmph'

;endfor ; iyr

print,'Done!'
end
