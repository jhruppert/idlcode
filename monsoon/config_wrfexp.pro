; 
; Initialize directories and common arrays for a set of WRF simulations.
; 
; James Ruppert
; 8/14/20
; 
pro config_wrfexp, casedir=casedir, cases=cases, dirs=dirs, dims=dims, $
  vars=vars, nfils_set=nfils_set, verbose=verbose, domtag=domtag

;Set nfils_set to number of raw output files to consider in case there's
;an uneven number of output files across directories

;CASE NAMES
  nc=n_elements(cases)
  print,'CASES: ',cases

  casedir+=cases+'/output/'

  ;RAW DOMAIN FILES
  ;if ~keyword_set(domtag) then domtag='d01'
  spawn,'ls '+casedir[0]+'wrfout_'+domtag+'_20*',domfiles
  nt=n_elements(domfiles)
;  files_raw=strarr(nc,3,nt)
  if keyword_set(nfils_set) then begin
    nt=nfils_set
    files_raw=strarr(nc,nfils_set)
  endif else files_raw=strarr(nc,nt)
  for ic=0,nc-1 do begin
;  for ic=0,0 do begin
    dom=domtag
;    for i=3,3 do begin
      spawn,'ls '+casedir[ic]+'wrfout_'+dom+'_20*',fils
      if ~keyword_set(nfils_set) then nfils=n_elements(fils)
      if keyword_set(nfils_set) then $
        files_raw[ic,0:nfils_set-1]=fils[0:nfils_set-1] $
      else $
        files_raw[ic,0:nfils-1]=fils
;    endfor
  endfor

  ;POST FILES
  spawn,'ls '+casedir[0]+'post/'+domtag+'/*.nc',fils
;  spawn,'ls '+casedir[2]+'post/*.nc',fils
  nvar=n_elements(fils)
  files_post=strarr(nc,nvar)
  for ic=0,nc-1 do begin
    spawn,'ls '+casedir[ic]+'post/'+domtag+'/*.nc',fils
    n_fil=n_elements(fils)
    files_post[ic,0:n_fil-1]=fils
  endfor

;VARIABLES

  spawn,'ls '+casedir[0]+'post/'+domtag+'/*.nc',tmp
;  spawn,'ls '+casedir[2]+'/post/*.nc',tmp
  vars=strarr(nvar)
  for iv=0,nvar-1 do begin
    tmp2=strsplit(tmp[iv],'/',/extract,count=nspl)
    tmp2=tmp2[nspl-1]
    tmp3=(strsplit(tmp2,'.',/extract))[0]
    vars[iv]=tmp3
  endfor

  dirs=create_struct(dirs,'nc',nc,'cases',cases,'casedir',casedir,'files_raw',files_raw,'files_post',files_post)
;  dirs=create_struct(dirs,'cases',cases,'casedir',casedir,'files_post',files_post)

if keyword_set(verbose) then begin

  print,''
  print,'DIRS:'
  help,dirs
  print,''

endif

  vars={nvar:nvar,vars:vars}

if keyword_set(verbose) then begin

  print,''
  print,'VARS:'
  help,vars
  print,''
  print,'Post-processed: ',vars.vars
  print,''

endif

;DIMENSIONS

  ;TIME ARRAY
  npd=24
  timstr=strmid(files_raw[0,0],18,19,/reverse_offset)
  yy=strmid(timstr,0,4)
  mm=strmid(timstr,5,2)
  dd=strmid(timstr,8,2)
  hh=strmid(timstr,11,2)
  tim0=julday(mm,dd,yy,hh,0,0)
;  if strmatch(subdir,'*edouard*') then $
;    tim0=julday(9,11,2014,12,0,0) $
;  else if strmatch(subdir,'*maria*') then $
;    tim0=julday(9,14,2017,12,0,0) $
;  else if strmatch(subdir,'*haiyan*') then $
;    tim0=julday(11,1,2013,0,0,0) $
;  else message,'Need to select a time-zero!'
  time=timegen(nt,start=tim0,step_size=24/npd,units='H')

  ifil=where(vars.vars eq 'U',count)
  if count ne 0 then pres=read_nc_var(files_post[0,ifil],'pres') $; deg
;  if count ne 0 then pres=read_nc_var(files_post[1,ifil],'pres') $; deg
  else pres=!values.f_nan
;pres=!values.f_nan

  np=n_elements(pres)

  lon=read_nc_var(files_raw[0,0],'XLONG') ; deg
  lon=reform(lon[*,0])
  nx=n_elements(lon)
  
  lat=read_nc_var(files_raw[0,0],'XLAT') ; deg
  lat=reform(lat[0,*])
  ny=n_elements(lat)
  
  dims={nt:nt,npd:npd,time:time,np:np,pres:pres,nx:nx,ny:ny,lon:lon,lat:lat}

if keyword_set(verbose) then begin

  print,''
  print,'DIMS:'
  help,dims
  print,''

endif

end
