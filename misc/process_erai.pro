;
; Program to read raw ERA-I data from individual monthly files in NCDF format.
;
; James Ruppert
; 17 Nov 2018
;
pro process_erai

  var_tag='850mb_uwind'
  months=['201110','201111'] ; 'YYYYMM'
  nfil=n_elements(months)

  read_dir='/work/06040/tg853394/stampede2/trmm/'
  outfil=read_dir+'erai_'+var_tag+'_'+strjoin(months[[0,nfil-1]],'-')+'.nc'

;RAW FILES

  files=strarr(nfil)
  for ifil=0,nfil-1 do $
    files[ifil]=read_dir+'ERAI_'+var_tag+'_'+months[ifil]+'.nc'

;LOOP THROUGH FILES

;  time=read_nc_var(files[0],'time')
  var=read_nc_var(files[0],'u') ; [x,y,t]

;READ

rain=fltarr(nx,ny,1,nt)

for it=0,nt-1 do begin

  caldat,times[it],mm,dd,yy,hh

  mm=string(mm,format='(i2.2)')
  dd=string(dd,format='(i2.2)')
  yy=string(yy,format='(i4.4)')
  hh=string(hh,format='(i2.2)')

  print,'YYYYMMDD-HH: ',yy,mm,dd,'-',hh

  file=read_dir+'raw/3B42.2011'+mm+dd+'.'+hh+'.7.HDF'

  ;READ VAR
    fid=hdf_sd_start(file,/read)
    vid=hdf_sd_select(fid,0) ; 0 is the index of 'precipitation'
    hdf_sd_getdata,vid,tmp
    hdf_sd_end,fid

  ;SET NANS
;    nan=where(tmp lt -1000,count)
;    if count gt 0 then tmp[nan]=!values.f_nan

  rain[*,*,0,it]=transpose(tmp)

endfor

write_sing_ncvar,outfil,rain,'trmm_rain_rate_mmph'

print,'Done!'
end
