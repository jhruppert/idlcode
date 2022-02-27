; 
; Initialize directories and common arrays for ensemble TC simulations.
; 
; James Ruppert
; 12/20/21
; 
pro tc_ens_config, tcname, case_str, dom, dirs=dirs, dims=dims, vars=vars, verbose=verbose

;DIRECTORIES

  ;CREATE STRUCTURE WITH MAIN GENERAL DIRECTORIES
  config_dir, dirs=dirs
  figdir=dirs.figdir+'tc/ens/'

;  maindir=dirs.scdir+'wrfenkf/ensemble/'
;  spawn,'ls '+maindir,all_cases
;  nc=n_elements(all_cases)

;for ic=0,nc-1 do begin

  tmpdir=dirs.scdir+'wrfenkf/ensemble/'+tcname+'/'
  spawn,'ls '+tmpdir+' | grep memb',memb

print,'MEMBERS: ',memb

  ensdir=tmpdir+memb+'/'+case_str+'/'
  nens=n_elements(memb)

  ;RAW DOMAIN FILES
;  dom='d02'
  spawn,'ls '+ensdir[0]+'wrfout_'+dom+'_20*',d02
  nt=n_elements(d02)
  files_raw=strarr(nens,nt)
;  for ic=0,nc-1 do $
  for imemb=0,nens-1 do begin
;    for i=3,3 do begin
      spawn,'ls '+ensdir[imemb]+'wrfout_'+dom+'_20*',fils
      nfils=n_elements(fils)
      files_raw[imemb,0:nfils-1]=fils
;    endfor
  endfor

  ;POST FILES
  spawn,'ls '+ensdir[0]+'post/'+dom+'/*.nc',fils
  nvar=n_elements(fils)
  files_post=strarr(nens,nvar)
;  for ic=0,nc-1 do $
  for imemb=0,nens-1 do begin
    spawn,'ls '+ensdir[imemb]+'post/'+dom+'/*.nc',fils
    n_fil=n_elements(fils)
    files_post[imemb,0:n_fil-1]=fils
  endfor

;VARIABLES

  spawn,'ls '+ensdir[0]+'post/'+dom+'/*.nc',tmp
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

  tstr=strmid(files_raw[0,0],7+11,13,/reverse)
  yy=strmid(tstr,0,4) & mm=strmid(tstr,5,2) & dd=strmid(tstr,8,2)
  hh0=strmid(tstr,11,2) & nn=strmid(tstr,14,2)
  tim0=julday(mm,dd,yy,hh0,nn,0)

  tstr=strmid(files_raw[0,nt-1],7+11,13,/reverse)
  yy=strmid(tstr,0,4) & mm=strmid(tstr,5,2) & dd=strmid(tstr,8,2)
  hh2=strmid(tstr,11,2) & nn=strmid(tstr,14,2)
  timend=julday(mm,dd,yy,hh2,nn,0)

  ;ASSUMING 1-HOURLY OUTPUT
  npd=24
  time=timegen(nt,start=tim0,final=timend,step_size=1,units='H')

  ifil=where(vars.vars eq 'T',count)
  if count ne 0 then pres=read_nc_var(files_post[0,ifil],'pres') $; deg
  else pres=!values.f_nan

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

;endfor ; ICASE

end
