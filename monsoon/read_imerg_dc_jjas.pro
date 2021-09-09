; 
; Read in diurnally varying IMERG for all of JJAS 2013-2017.
;
; James Ruppert
; 12/31/20
; 
function read_imerg_dc_jjas, datdir, lon=lon, lat=lat, bounds=bounds, skip=skip

  npd_imerg=48

;----READ IMERG RAINFALL--------------------

  if ~keyword_set(skip) then skip=1

;  spawn,'ls '+datdir+'merge_*nc4',imfils
  spawn,'ls '+datdir+'*JJAS_[0-9]*nc4',imfils

  ;TIME
    time=read_nc_var(imfils[0],'time')
    time=julday(1,1,1970,0,0,time)
    nt_read=n_elements(time)

  ;LAT/LON
    im_fil=imfils[0]
    lon=read_nc_var(imfils[0],'lon')
    lat=read_nc_var(imfils[0],'lat')
    nx=n_elements(lon)
    ny=n_elements(lat)

  ;SPACE SUBSET
    ix0=0 & iy0=0
    if keyword_set(bounds) then begin
      ix=where((lon ge bounds[0]) and (lon le bounds[2]),nx)
      iy=where((lat ge bounds[1]) and (lat le bounds[3]),ny)
      ix0=ix[0] & iy0=iy[0]
      lon=lon[ix] & lat=lat[iy]
    endif

  ;DIURNALLY VARYING
    ;READ DATA FROM HOURLY FILES THAT CONTAIN ENTIRE JJAS13-17 DATASET
      rain=fltarr(nx,ny,npd_imerg/skip,nt_read,/nozero)
      ifil=0
      for ih=0,npd_imerg-1,skip do begin
        print,'  Reading: ',string(ih*0.5,format='(f4.1)'),'h'
        ;READ ENTIRE SET OF DAYS
          count=[ny,nx,nt_read] & offset=[iy0,ix0,0] ; y,x,t
;          file=imfils[ih]
          hh=24*ih/48
          hh=string(hh,format='(i2.2)')
          nn=0
          if (24*ih mod 48) ne 0 then nn=30
          nn=string(nn,format='(i2.2)')
          file=datdir+'merge_'+hh+nn+'.nc4'
          tmp=reform(read_nc_var(file,'precipitationCal',count=count,offset=offset)) ; mm/hr
          tmp2=transpose(tmp,[1,0,2])
          rain[*,*,ifil,*]=tmp2
          ifil+=1
      endfor
      rain*=24.    ; mm/hr --> mm/d

  return, rain

end
