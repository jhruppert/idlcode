; 
; Call routine to calculate TC tracks for ensemble runs.
;
; James Ruppert
; 12/20/21
; 
pro run_enstc_track

;tcname='maria'
tcname='haiyan'
case_str='ctl'
case_str='haiyan'

if tcname eq 'maria' then begin
  tcyear='2017'
  hurdat=read_hurdat(tcname,tcyear)
endif else if tcname eq 'haiyan' then begin
  tcyear='2013'
  hurdat=read_jtwcdat(tcname)
endif

dom='d02'
tc_ens_config, case_str, dom, dirs=dirs, dims=dims, vars=vars;, /verbose
dirs.figdir+=tcname+'/'

;VORTEX LEVEL
  psel=700 ; (hPa) level to locate vortex
  izsel=(where(dims.pres eq psel))[0]

;TIME ARRAYS
  time=dims.time
  nt_full=dims.nt;-1
  npd=dims.npd
  nhrs=1.*nt_full*npd/24.
  hours=indgen(nhrs)
  nd=(nhrs-(nhrs mod 24))/24.

;TIME SUBSET
;  it_sel=indgen(nt_full)
;  nt=n_elements(it_sel)

;----LOOP OVER SIMULATIONS--------------------

for ic=0,dirs.nens-1 do begin
;for ic=0,4 do begin
;for ic=1,1 do begin

  print,'Memb: ',ic+1

  ;READ ABSOLUTE VORTICITY
;    iv=where(vars.vars eq 'AVOR')
;    file=dirs.files_post[ic,iv]
;    count=[dims.nx,dims.ny,1,nt] & offset=[0,0,izsel,0] ; x,y,z,t
;    avor=reform(read_nc_var(file,'AVOR',count=count,offset=offset))
  ;READ SLP
    iv=where(vars.vars eq 'SLP')
    file=dirs.files_post[ic,iv]
;    count=[dims.nx,dims.ny,1,nt] & offset=[0,0,0,0] ; x,y,z,t
    slp=reform(read_nc_var(file,'SLP'));,count=count,offset=offset))

  ;TC TRACKING
    trackdir=dirs.ensdir[ic]
    tc_track,slp,dims,time,hurdat,trackdir=trackdir

endfor ; icase


print,'DONE!!'
end
