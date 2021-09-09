; 
; Generate azimuthal data from TC output.
;
; James Ruppert
; 4/6/19
; 
pro run_calc_azim

tcname='maria'
;tcname='haiyan'

if tcname eq 'maria' then begin
  tcyear='2017'
  hurdat=read_hurdat(tcname,tcyear)
endif else if tcname eq 'haiyan' then begin
  tcyear='2013'
  hurdat=read_jtwcdat(tcname)
endif

;subdir='moving_nest/'+tcname
;subdir='static_nest/'+tcname
subdir='redux/'+tcname
tc_sim_config, subdir, dirs=dirs, dims=dims, vars=vars;, /verbose
dirs.figdir+=tcname+'/'

;VORTEX TRACKING LEVEL
  psel=700 ; (hPa) level to locate vortex
  izsel=(where(dims.pres eq psel))[0]

;TIME SELECTION
;  hr_sel=[0,200] ; entire simulation
;  hr_sel=[24,120] ; for Maria/CTL "static_nest"
;  hr_sel=[49,72] ; for Maria/other tests
;  hr_sel=[49,120]
if tcname eq 'maria' then $
  hr_sel=[0,144] $; for Maria/CTL "redux" PAPER
else $
  hr_sel=[0,168] ; for Haiyan PAPER

;VAR SELECTION
;  vars_sel=['U10','V10','LH','OLR','OLRC','SLP','PW','GLW','GLWC','rainrate'];,$
;            'QVAPOR','T','H_DIABATIC','RTHRATLWC','RTHRATLW','RTHRATSW','RTHRATSWC','QCLOUD','QICE','W','AVOR']
;            'U','V']
;  vars_sel=['RTHRATLW','QICE']

  ;2D
  vars_sel=['SLP','PW','OLR','OLRC','rainrate','U10','V10']

  ;SET NEEDED FOR PAPER
  vars_sel=['SLP','PW','OLR','OLRC','H_DIABATIC','RTHRATLWC','RTHRATLW'];,'W','U','V','T'];,'QVAPOR','QCLOUD','QICE']
  ;FOR ALL TESTS
;    vars_sel=['U10','V10','RTHRATLW','W']
  ;EXTRA FOR CTL
;    vars_sel=['U','V','QVAPOR','T']

;vars_sel=['SLP','PW','OLR','OLRC','QVAPOR','QCLOUD','QICE','rainrate','U10','V10','H_DIABATIC','RTHRATSW','RTHRATLW']
vars_sel=['H_DIABATIC','RTHRATLWC','RTHRATLW','U','V','W','T']

;----TIME SPECS--------------------

;FULL TIME SERIES
  time=dims.time
  nt=dims.nt;-1
  nt_full=nt
  npd=dims.npd
  nhrs=1.*nt*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nhrs)

;TIME SELECTION
  t_ind=where((time_hrs ge hr_sel[0]) and (time_hrs le hr_sel[1]))
  t_ind=t_ind[where(t_ind le nt-1)]
  nt_sav=n_elements(t_ind)
  ;OVERWRITE THESE
;    nhrs=1.*nt_sav*npd/24.
;    nd=(1.*nhrs-(nhrs mod 24))/24.
;    time_hrs_sav=indgen(nt_sav)+hr_sel[0]

;----READ AND PROCESS VARS--------------------


;LAND MASK
  vtag='LANDMASK'
  file=dirs.files_raw[0,2,0]
  mask=reform(read_nc_var(file,'LANDMASK'))
  land=where(mask eq 1,nland)
  mask=0


;for ic=0,dirs.nc-1 do begin
for ic=1,dirs.nc-1 do begin
;for ic=0,0 do begin

  print,'CASE: ',strupcase(dirs.cases[ic])

  ;VAR SELECTION
    nvsel=n_elements(vars_sel)
    iv_sel=intarr(nvsel)
    for iv=0,nvsel-1 do $
      iv_sel[iv]=(where(strmatch(reform(dirs.files_post[ic,*]),'*/'+vars_sel[iv]+'.*')))[0]

  hr0=0
  if strmatch(dirs.cases[ic],'*36h*') then hr0=36
  if strmatch(dirs.cases[ic],'*24h*') then hr0=24
  if strmatch(dirs.cases[ic],'*48h*') then hr0=48
  if strmatch(dirs.cases[ic],'*60h*') then hr0=60
  if strmatch(dirs.cases[ic],'*72h*') then hr0=72
  if strmatch(dirs.cases[ic],'*84h*') then hr0=84
  if strmatch(dirs.cases[ic],'*96h*') then hr0=96
  if strmatch(dirs.cases[ic],'lwcrf*') then hr0=36
  if strmatch(dirs.cases[ic],'lwcrf2*') then hr0=36
  if strmatch(dirs.cases[ic],'axisym') then hr0=36
  nt_test_full=nt_full - hr0

  ;TIME SELECTION
    t_offset=max([0,hr0-hr_sel[0]])
    ut_offset=max([0,hr_sel[0]-hr0])
    nt_test=nt_sav-t_offset
    t_ind_test=indgen(nt_test)+ut_offset
    time_hrs_test=time_hrs[t_ind_test+hr0]
    hr_tag_test=strtrim(time_hrs_test[0],2)+'-'+strtrim(time_hrs_test[nt_test-1],2)+'hr'

;VORTEX TRACKING

  ;GET TC LOC FOR ENTIRE LENGTH OF TEST; LOC IS INDEXED INSIDE OF wrf_regrid_azim

  ;READ ABSOLUTE VORTICITY
    iv=where(vars.vars eq 'AVOR')
    file=dirs.files_post[ic,iv]
    count=[dims.nx,dims.ny,1,nt_test_full] & offset=[0,0,izsel,0] ; x,y,z,t
    avor=reform(read_nc_var(file,'AVOR',count=count,offset=offset))

  ;SMOOTH
    ixsmooth=round(111./3) ; 1-degree smoothing, run twice
    ismooth=[ixsmooth,ixsmooth,0]
    for i=1,2 do $
      avor=smooth(temporary(avor),ismooth,/edge_truncate,/nan)

  ;MASKING
    if nland gt 0 then avor[land]=!values.f_nan
    if tcname eq 'maria' then avor=wrf_maria_mask(temporary(avor),time[t_ind_test+hr0],hurdat,dims); else stop
    if tcname eq 'haiyan' then avor=wrf_haiyan_mask(temporary(avor),time[t_ind_test+hr0],hurdat,dims); else stop

  ;VORTEX TRACKING
    vloc=maria_vortex_locate(avor,dims);,/write)
    avor=0


;GENERATE AZIMUTHAL DATA FOR VAR SELECTION

for iv=0,nvsel-1 do begin
;for iv=9,nvsel-1 do begin

;SKIP
;if strupcase(dirs.cases[ic]) eq 'NCRF_36H' and iv eq 0 then iv=6

  varfil=dirs.files_post[ic,iv_sel[iv]]
  print,'INFIL: ',varfil

  outfil=dirs.casedir[ic]+'azim_'+vars_sel[iv]+'_'+hr_tag_test+'.nc'
  print,'OUTFIL: ',outfil

  exists=file_test(outfil)
;  if exists then continue

  wrf_regrid_azim,outfil,varfil,vars_sel[iv],t_ind_test,dims,vloc,override=0

endfor ; ifil

endfor ; icase

print,'DONE!!'
end
