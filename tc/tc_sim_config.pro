; 
; Initialize directories and common arrays for TC simulations.
; 
; James Ruppert
; 3/12/19
; 
pro tc_sim_config, subdir, dirs=dirs, dims=dims, vars=vars, verbose=verbose

;DIRECTORIES

  ;CREATE STRUCTURE WITH MAIN GENERAL DIRECTORIES
  config_dir, dirs=dirs
  figdir=dirs.figdir+'tc/'

  casedir=dirs.scdir+'tc/output/'+subdir+'/'
  spawn,'ls '+casedir,cases
;  spawn,'ls -I initial '+casedir,cases
;OVERWRITE WITH SELECTION
;cases=['ctl','ctl_nocrf','ctl_qvdry','ctl_nocrf_qvdry','ctl_qvmoist']
;cases=['ctl','ctl_qvdry','ctl_nocrf_qvdry','ctl_qvmoist']
;cases=['ctl','ctl_qvdry','ctl_nocrf_qvdry','ctl_qvmoist','ctl_nocrf_qvmoist']
;cases=['ctl','ctl_nocrf','ctl_nocrf_48h','ctl_nocrf_72h']
;cases=['ctl','ctl_nocrf','nosw','nosw_nocrf']
;cases=['edouard','harvey']
;cases=['joaquin','joaquin_ncrf','harvey','harvey_ncrf']
;cases=['initial','ctl'];,'ncrf_48h']
;cases=cases[[0,4,5,6,7]]
cases=['ctl',reverse(['ncrf_36h','ncrf_48h','ncrf_60h','ncrf_72h','ncrf_84h','ncrf_96h'])]
;cases=['ctl','ncrf_36h','ncrf_48h']
;cases=['ctl','lwcrf2','lwcrf2_pbl','lwcrf2_xpbl','lwcrf2_xcool'];,'ncrf_84h'];,'ncrf_84h','ncrf_36h']
;cases=['ctl','lwcrf2','lwcrf2_pbl','lwcrf2_xpbl','lwcrf2_xcool'];,'lwcrf2_xpbl','lwcrf2_xcool']
cases=['ncrf_84h','ncrf_72h','ncrf_60h','ncrf_48h','ncrf_36h']
cases=['ctl',reverse(['ncrf_36h','ncrf_48h','ncrf_60h','ncrf_72h','ncrf_84h','ncrf_96h','lwcrf'])];,'axisym']
;PAPER SELECTION
  cases=['ctl',reverse(['ncrf_36h','ncrf_60h','ncrf_96h','lwcrf'])];,'axisym']
;FOR RAINFALL PDF
cases=['ctl','lwcrf','ncrf_36h'];,'ncrf_96h']
;HAIYAN
if strmatch(subdir,'*haiy*') then $
  cases=['ctl',reverse(['hncrf_36h','hncrf_60h','hncrf_96h'])]
cases=['ctl','ncrf_36h']
cases=['ctl','axisym']


print,'CASES: ',cases

  casedir+=cases+'/'
  nc=n_elements(cases)

  ;REARRANGE CASES
;  casedir=casedir[[0,1,4,5,2,3]]
;  cases=cases[[0,1,4,5,2,3]]

  ;RAW DOMAIN FILES
  spawn,'ls '+casedir[0]+'wrfout_d01_20*',d01
;  spawn,'ls '+casedir[0]+'wrfout_d03_20*',d01
  nt=n_elements(d01)
;  files_raw=strarr(nc,3,nt)
  files_raw=strarr(nc,nt)
  for ic=0,nc-1 do begin
    dom='d03'
    if (strmatch(subdir,'*maria*') or strmatch(subdir,'*haiyan*')) and (cases[ic] ne 'initial') and (cases[ic] ne '3dom') then dom='d01'
;    for i=3,3 do begin
      spawn,'ls '+casedir[ic]+'wrfout_'+dom+'_20*',fils
      nfils=n_elements(fils)
      files_raw[ic,0:nfils-1]=fils
;    endfor
  endfor

  ;POST FILES
  spawn,'ls '+casedir[0]+'post/*.nc',fils
;  spawn,'ls '+casedir[2]+'post/*.nc',fils
  nvar=n_elements(fils)
  files_post=strarr(nc,nvar)
  for ic=0,nc-1 do begin
    spawn,'ls '+casedir[ic]+'post/*.nc',fils
    n_fil=n_elements(fils)
    files_post[ic,0:n_fil-1]=fils
  endfor

;VARIABLES

  spawn,'ls '+casedir[0]+'post/*.nc',tmp
;  spawn,'ls '+casedir[2]+'/post/*.nc',tmp
  vars=strarr(nvar)
  for iv=0,nvar-1 do vars[iv]=(strsplit((strsplit(tmp[iv],'/',/extract))[9],'.',/extract))[0]

  dirs.figdir=figdir
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

  npd=24
  if strmatch(subdir,'*edouard*') then $
    tim0=julday(9,11,2014,12,0,0) $
  else if strmatch(subdir,'*maria*') then $
    tim0=julday(9,14,2017,12,0,0) $
  else if strmatch(subdir,'*haiyan*') then $
    tim0=julday(11,1,2013,0,0,0) $
  else message,'Need to select a time-zero!'
  time=timegen(nt,start=tim0,step_size=24/npd,units='H')

  ifil=where(vars.vars eq 'T',count)
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
