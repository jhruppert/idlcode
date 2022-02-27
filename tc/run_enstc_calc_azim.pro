; 
; Generate azimuthal data from ensemble TC output.
;
; Assumes TC tracking has been done by run_enstc_track.pro
;
; Input:
;   imemb - ensemble member (0-19)
;
; James Ruppert
; 12/25/21
; 
pro run_enstc_calc_azim, imemb, iv_calc, tcname, case_str

;tcname='maria'
;case_str='ctl'

if tcname eq 'maria' then begin
  tcyear='2017'
  hurdat=read_hurdat(tcname,tcyear)
endif else if tcname eq 'haiyan' then begin
  tcyear='2013'
  hurdat=read_jtwcdat(tcname)
endif

dom='d02'
tc_ens_config, tcname, case_str, dom, dirs=dirs, dims=dims, vars=vars;, /verbose

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

  ;FOR ML TRAINING
;    vars_sel=['SLP','U','V','W','QVAPOR','T','H_DIABATIC','RTHRATLWC','RTHRATLW','RTHRATSWC','RTHRATSW']
    vars_sel=['SLP','U10','V10','U','V','W','QVAPOR','T','H_DIABATIC','RTHRATLWC','RTHRATLW','RTHRATSWC','RTHRATSW']

;vars_sel=['SLP','PW','OLR','OLRC','QVAPOR','QCLOUD','QICE','rainrate','U10','V10','H_DIABATIC','RTHRATSW','RTHRATLW']
;vars_sel=['RTHRATLWC','RTHRATSWC']

;----TIME SPECS--------------------

;FULL TIME SERIES
  time=dims.time
  nt=dims.nt;-1
  nt_full=nt
  npd=dims.npd
  nhrs=1.*nt*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nhrs)

  nt_test=nt_full
  hr_tag_test=strtrim(time_hrs[0],2)+'-'+strtrim(time_hrs[nt_test-1],2)+'hr'

;----READ AND PROCESS VARS--------------------

;for ic=0,dirs.nens-1 do begin
;for ic=imemb,imemb do begin

ic=imemb-1

  print,'Memb: ',ic+1

  ;TC TRACKING
    spawn,'ls '+dirs.ensdir[ic]+'track*',trackfiles,err
    ntrack=n_elements(trackfiles)
    if ~keyword_seT(trackfiles) then begin
      print,'No track files for this member!'
      return
    endif

for itrack=0,ntrack-1 do begin

  vloc=reform(read_nc_var(trackfiles[itrack],'locvort'))

  t_ind_test=where(finite(reform(vloc[0,*])))

;GENERATE AZIMUTHAL DATA FOR VAR SELECTION

  ;VAR SELECTION
    nvsel=n_elements(vars_sel)
    iv_sel=intarr(nvsel)
    for iv=0,nvsel-1 do $
      iv_sel[iv]=(where(strmatch(reform(dirs.files_post[ic,*]),'*/'+vars_sel[iv]+'.*')))[0]

;  for iv_calc=0,nvsel-1 do begin
  
    varfil=dirs.files_post[ic,iv_sel[iv_calc]]
    print,'INFIL: ',varfil
  
;    outfil=dirs.ensdir[ic]+'azim_'+vars_sel[iv_calc]+'_'+hr_tag_test+'.nc'
    outfil=dirs.ensdir[ic]+'azim_'+vars_sel[iv_calc]+'_'+hr_tag_test+'_track'+strtrim(itrack+1,2)+'.nc'
    print,'OUTFIL: ',outfil
  
    exists=file_test(outfil)
  ;  if exists then continue

    wrf_regrid_azim,outfil,varfil,vars_sel[iv_calc],t_ind_test,dims,vloc
  
;  endfor ; iv_calc

endfor ; itrack

print,'DONE!!'
end
