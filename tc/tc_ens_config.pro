; 
; Initialize directories and common arrays for ensemble TC simulations.
; 
; James Ruppert
; 12/20/21
; 
pro tc_ens_config, tcname, case_str, dom, imember=imember, dirs=dirs, dims=dims, vars=vars, verbose=verbose

imembset=0
if keyword_set(imember) then begin
  imembset=1
  imemb=imember-1
endif

;DIRECTORIES

  ;CREATE STRUCTURE WITH MAIN GENERAL DIRECTORIES
  config_dir, dirs=dirs
  figdir=dirs.figdir+'tc/ens/'

;  maindir=dirs.scdir+'wrfenkf/ensemble/'
;  spawn,'ls '+maindir,all_cases
;  nc=n_elements(all_cases)

;for ic=0,nc-1 do begin

  tmpdir=dirs.scdir+'wrfenkf/ensemble/'+tcname+'/'
  spawn,'ls '+tmpdir+' | grep memb',memb,err_result

;print,'MEMBERS: ',memb

  ensdir=tmpdir+memb+'/'+case_str+'/'
  nens=n_elements(memb)

  if imembset then icheck=imemb else icheck=0

  ;RAW DOMAIN FILES
;  dom='d02'
  spawn,'ls '+ensdir[icheck]+'wrfout_'+dom+'_20*',d02,err_result
  nt=n_elements(d02)
  files_raw=strarr(nens,nt)
;  for ic=0,nc-1 do $
  if imembset then begin
      spawn,'ls '+ensdir[imemb]+'wrfout_'+dom+'_20*',fils,err_result
      nfils=n_elements(fils)
      files_raw[imemb,0:nfils-1]=fils
  endif else begin
  for imemb=0,nens-1 do begin
;    for i=3,3 do begin
      spawn,'ls '+ensdir[imemb]+'wrfout_'+dom+'_20*',fils,err_result
      nfils=n_elements(fils)
      files_raw[imemb,0:nfils-1]=fils
;    endfor
  endfor
  endelse

  ;POST FILES
  spawn,'ls '+ensdir[icheck]+'post/'+dom+'/*.nc',fils,err_result
  nvar=n_elements(fils)
  files_post=strarr(nens,nvar)
;  for ic=0,nc-1 do $
  if imembset then begin
    spawn,'ls '+ensdir[imemb]+'post/'+dom+'/*.nc',fils,err_result
    n_fil=n_elements(fils)
    files_post[imemb,0:n_fil-1]=fils
  endif else begin
  for imemb=0,nens-1 do begin
    spawn,'ls '+ensdir[imemb]+'post/'+dom+'/*.nc',fils,err_result
    n_fil=n_elements(fils)
    files_post[imemb,0:n_fil-1]=fils
  endfor
  endelse

;VARIABLES

  spawn,'ls '+ensdir[icheck]+'post/'+dom+'/*.nc',tmp,err_result
  vars=strarr(nvar)
  for iv=0,nvar-1 do begin
    tmp2=strsplit(tmp[iv],'/',/extract)
    nstr=n_elements(tmp2)
    vars[iv]=(strsplit(tmp2[nstr-1],'.',/extract))[0]
  endfor

  dirs.figdir=figdir
  dirs=create_struct(dirs,'nens',nens,'memb',memb,'ensdir',ensdir,'files_raw',files_raw,'files_post',files_post)

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

;MANUAL OVEIRRIDES FOR TIME DIMENSIONALITY
  if tcname eq 'haiyan' then begin
    if case_str eq 'ctl' then begin
      tstr0='2013-11-01_00'
      tstr1='2013-11-08_00'
    endif else if case_str eq 'ncrf' then begin
      tstr0='2013-11-02_12'
      tstr1='2013-11-08_00'
    endif else if case_str eq 'crfon' then begin
      tstr0='2013-11-03_12'
      tstr1='2013-11-04_12'
    endif else message,'Set case name!'
  endif

;  tstr0=strmid(files_raw[0,0],7+11,13,/reverse)
  yy=strmid(tstr0,0,4) & mm=strmid(tstr0,5,2) & dd=strmid(tstr0,8,2)
  hh0=strmid(tstr0,11,2) & nn=strmid(tstr0,14,2)
  tim0=julday(mm,dd,yy,hh0,nn,0)

;  tstr1=strmid(files_raw[0,nt-1],7+11,13,/reverse)
  yy=strmid(tstr1,0,4) & mm=strmid(tstr1,5,2) & dd=strmid(tstr1,8,2)
  hh2=strmid(tstr1,11,2) & nn=strmid(tstr1,14,2)
  timend=julday(mm,dd,yy,hh2,nn,0)

  ;ASSUMING 1-HOURLY OUTPUT
  npd=24
  time=timegen(nt,start=tim0,final=timend,step_size=1,units='H')

  ifil=where(vars.vars eq 'T',count)
  if imembset then icheck=imemb else icheck=0

  if count ne 0 then pres=read_nc_var(files_post[icheck,ifil],'pres') $; deg
  else pres=!values.f_nan

  np=n_elements(pres)

  lon=read_nc_var(files_raw[icheck,0],'XLONG') ; deg
  lon=reform(lon[*,0])
  nx=n_elements(lon)
  
  lat=read_nc_var(files_raw[icheck,0],'XLAT') ; deg
  lat=reform(lat[0,*])
  ny=n_elements(lat)
  
  dims={nt:nt,npd:npd,time:time,np:np,pres:pres,nx:nx,ny:ny,lon:lon,lat:lat}

if keyword_set(verbose) then begin

  print,''
  print,'DIMS:'
  help,dims
  print,''

endif

;endfor ; ICASE

end
