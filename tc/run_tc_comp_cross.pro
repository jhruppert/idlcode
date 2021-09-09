; 
; Generate PW- or OLR-binned cross section from TC output
;
; Mean water vapor profile used for experiment 'FIXQV' generated here (average over hr_sel=[0,24*2]).
;
; Bin-averaged profiles are also calculated here (see iprof).
;
; James Ruppert
; 3/12/19
; 
pro run_tc_comp_cross

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
subdir='static_nest/'+tcname
subdir='redux/'+tcname
tc_sim_config, subdir, dirs=dirs, dims=dims, vars=vars;, /verbose
dirs.figdir+=tcname+'/'

model_pres_levels_nc=dirs.casedir[0]+'mean_pres_tim0_ctl.nc'

;----PLOT OPTIONS--------------------


binvar_str='pw';'rainrate';'pw';'olr';'pw';
var1_str='RTHRATLW'
  iadd_sw=0 ; Add SW?
;  iadd_sw=1 ; Add SW?
var1_str='RTHRATSW'
;var1_str='RH';'QVAPOR';
;var1_str='W'
;var1_str='H_DIABATIC'

icrf=0
icrf=1 ; Subtract clear sky?
;  olr_thresh=290 ; OLR threshold for clear-sky estimate (W/m2)
;  pw_thresh_crf=[30,40] ; OLR threshold for clear-sky estimate (W/m2)

icont=1 ; 0=stream, 1=cloud, 2=LH, 3=AVOR, 4=W

iprof=1 ; do profiles?
if binvar_str ne 'pw' then iprof=0 ; only do profiles with PW

write_qv=0 ; write qv-profile for dry test?
  if write_qv and ((var1_str ne 'RH') or (binvar_str ne 'pw')) then $
    message,'Need RH, PW to write qv-profile.'
;  pw_thresh_fixqv=[40,55]
  pw_thresh_fixqv=[55,70]
  qv_file=dirs.figdir+'wvapor_profile_edouard_'+string(pw_thresh_fixqv[0],format='(i2.2)')+'-'+$
    string(pw_thresh_fixqv[1],format='(i2.2)')+'mm.txt'

;TIME SELECTION
;  hr_sel=[0,200] ; entire simulation
;  hr_sel=[0,24*2]
  hr_sel0=[49,72]
  hr_sel0=[49,96]
;  hr_sel0=[74,77] ; 3 h centered on local noon
;  hr_sel=[49,60]
;for iadd=-1,3 do begin
for iadd=-0,0 do begin
hr_sel=hr_sel0+iadd*24
;  hr_sel=[49,96]
  hr_tag=string(hr_sel[0],format='(i3.3)')+'-'+string(hr_sel[1],format='(i3.3)')+'hr'

;VORTEX TRACKING LEVEL
  psel=700 ; (hPa) level to locate vortex
  izsel=(where(dims.pres eq psel))[0]

;RADIUS THRESHOLD TO CONSIDER
  radius_thresh=800 ; km


;----TIME SPECS--------------------


;FULL TIME SERIES
  time=dims.time
  nt=dims.nt-1
  npd=dims.npd
  nhrs=1.*nt*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nhrs)

;TIME SELECTION
  t_ind=where((time_hrs ge hr_sel[0]) and (time_hrs le hr_sel[1]))
  t_ind=t_ind[where(t_ind le nt-1)]
  nt_sel=n_elements(t_ind)
  ;OVERWRITE THESE
    nhrs=1.*nt_sel*npd/24.
    nd=(1.*nhrs-(nhrs mod 24))/24.
    time_hrs=indgen(nt_sel)
    time=time[t_ind]
    nt_sav=nt_sel


;----READ VARS--------------------


print,'VAR: ',var1_str

;LAND MASK
  vtag='LANDMASK'
  file=dirs.files_raw[0,2,0]
  mask=reform(read_nc_var(file,'LANDMASK'))
  land=where(mask eq 1,nland)
  mask=0

;NEED LARGE PRES ARRAY FOR RH
  if var1_str eq 'RH' then begin
    prx=fltarr(dims.nx,dims.ny,dims.np,nt_sel)
    for iz=0,dims.np-1 do prx[*,*,iz,*]=dims.pres[iz]*1e2
  endif


;for ic=0,dirs.nc-1 do begin
;for ic=0,3 do begin
for ic=0,0 do begin

  print,'CASE: ',strupcase(dirs.cases[ic])

  hr0=0
  if strmatch(dirs.cases[ic],'*36h*') then hr0=36
  if strmatch(dirs.cases[ic],'*24h*') then hr0=24
  if strmatch(dirs.cases[ic],'*48h*') then hr0=48
  if strmatch(dirs.cases[ic],'*60h*') then hr0=60
  if strmatch(dirs.cases[ic],'*72h*') then hr0=72
  if strmatch(dirs.cases[ic],'*84h*') then hr0=84

  ;TIME SELECTION
    if hr0 gt hr_sel[0] then begin
      print,'Not enough times in test, so skipping...'
      continue
    endif
    t_offset=max([0,hr0-hr_sel[0]])
    ut_offset=max([0,hr_sel[0]-hr0])
    nt_test=nt_sel-t_offset
    t_ind_test=indgen(nt_test)+ut_offset
;    hrs_test=time_hrs[t_ind_test+hr0]

  it0 = t_ind_test[0]

  ;VORTEX TRACKING

    ;READ ABSOLUTE VORTICITY
      iv=where(vars.vars eq 'AVOR')
      file=dirs.files_post[ic,iv]
      count=[dims.nx,dims.ny,1,nt_test] & offset=[0,0,izsel,it0] ; x,y,z,t
      avor=reform(read_nc_var(file,'AVOR',count=count,offset=offset))

    ;SMOOTH
      ixsmooth=round(111./3) ; 1-degree smoothing, run twice
      ismooth=[ixsmooth,ixsmooth,0]
      for i=1,2 do $
        avor=smooth(temporary(avor),ismooth,/edge_truncate,/nan)

    ;MASKING
      if nland gt 0 then avor[land]=!values.f_nan
      if tcname eq 'maria' then avor=wrf_maria_mask(temporary(avor),time,hurdat,dims)
      if tcname eq 'haiyan' then avor=wrf_haiyan_mask(temporary(avor),time,hurdat,dims)

    ;VORTEX TRACKING
      vloc=maria_vortex_locate(avor,dims);,/write)
      avor=0

  ;BIN VAR
    iv=where(strmatch(vars.vars,binvar_str,/fold_case))
;where(vars.vars eq strupcase(binvar_str))
    file=dirs.files_post[ic,iv]
    count=[dims.nx,dims.ny,1,nt_test] & offset=[0,0,0,it0] ; x,y,z,t
    if binvar_str ne 'rainrate' then str=strupcase(binvar_str) else str=binvar_str
    binvar=reform(read_nc_var(file,str,count=count,offset=offset))
    if binvar_str eq 'rainrate' then begin
      binvar*=1./24 ; mm/h
    endif
;    if binvar_str eq 'olr' then begin
;      binvar[*,*,0]=!values.f_nan
;      ix=indgen(3)
;      binvar[[ix,dims.nx-1-ix],*,*]=!values.f_nan
;      binvar[*,[ix,dims.ny-1-ix],*]=!values.f_nan
;    endif

  ;OLR
    iv=where(vars.vars eq 'OLR')
    file=dirs.files_post[ic,iv]
    olr=reform(read_nc_var(file,'OLR',count=count,offset=offset))
    iv=where(vars.vars eq 'OLRC')
    file=dirs.files_post[ic,iv]
    olrc=reform(read_nc_var(file,'OLRC',count=count,offset=offset))

  ;MAIN VARS

    ;VAR1

    if var1_str eq 'RH' then begin

      iv=where(vars.vars eq 'QVAPOR')
      file=dirs.files_post[ic,iv]
      count=[dims.nx,dims.ny,dims.np,nt_test] & offset=[0,0,0,it0] ; x,y,z,t
      qv=read_nc_var(file,'QVAPOR',count=count,offset=offset)

      iv=where(vars.vars eq 'T')
      file=dirs.files_post[ic,iv]
      tmpk=read_nc_var(file,'T',count=count,offset=offset)

      var1=calc_relh(qv,tmpk,prx)

;      if ~write_qv then qv=0
;      tmpk=0

    endif else begin

      iv=where(vars.vars eq var1_str)
      file=dirs.files_post[ic,iv]
      count=[dims.nx,dims.ny,dims.np,nt_test] & offset=[0,0,0,it0] ; x,y,z,t
      var1=read_nc_var(file,var1_str,count=count,offset=offset)

    ;ADD SW
      iv=where(vars.vars eq 'RTHRATSW')
      file=dirs.files_post[ic,iv]
      if iadd_sw and strmatch(var1_str,'*RTHRA*') then var1+=read_nc_var(file,'RTHRATSW',count=count,offset=offset)

    ;CRF
      iv=where(vars.vars eq 'RTHRATLWC')
      file=dirs.files_post[ic,iv]
      if icrf and strmatch(var1_str,'*RTHRA*') then begin
        if strmatch(var1_str,'RTHRATLW') then begin
          var1-=read_nc_var(file,'RTHRATLWC',count=count,offset=offset)
;var1=read_nc_var(file,'RTHRATLWC',count=count,offset=offset)
          if iadd_sw then begin
            iv=where(vars.vars eq 'RTHRATSWC')
            file=dirs.files_post[ic,iv]
            var1-=read_nc_var(file,'RTHRATSWC',count=count,offset=offset) 
          endif
        endif else if strmatch(var1_str,'RTHRATSW') then begin
          iv=where(vars.vars eq 'RTHRATSWC')
          file=dirs.files_post[ic,iv]
          var1-=read_nc_var(file,'RTHRATSWC',count=count,offset=offset)
;var1=read_nc_var(file,'RTHRATSWC',count=count,offset=offset)
        endif else message,'Wrong Rad var for CRF'
      endif

    endelse

    ;TAKE CARE OF BCS FOR RAD
;      if (t_ind[0] eq 0) and strmatch(var1_str,'*RTHRA*') then begin
;        var1[*,*,0,*]=!values.f_nan
;        ix=indgen(3)
;        var1[[ix,dims.nx-1-ix],*,*,*]=!values.f_nan
;        var1[*,[ix,dims.ny-1-ix],*,*]=!values.f_nan
;      endif

    ;PSI
    if icont eq 0 then begin
    ;W
      iv=where(vars.vars eq 'W')
      file=dirs.files_post[ic,iv]
      w=read_nc_var(file,'W',count=count,offset=offset)
      iv=where(vars.vars eq 'T')
      file=dirs.files_post[ic,iv]
      tmpk=read_nc_var(file,'T',count=count,offset=offset)
      iv=where(vars.vars eq 'QVAPOR')
      file=dirs.files_post[ic,iv]
      qvt=read_nc_var(file,'QVAPOR',count=count,offset=offset)
      tvirt = tmpk*(1.+0.61*qvt)
      tmpk=0 & qvt=0
      rho=fltarr(dims.nx,dims.ny,dims.np,nt_sel)
      for iz=0,dims.np-1 do rho[*,*,iz,*] = dims.pres[iz]*1e2 / ( 287. * tvirt[*,*,iz,*] )
      var2=w
      rho=mean(mean(reform(temporary(rho[*,*,*,0])),dimension=1,/nan,/double),dimension=1,/nan,/double)
      w=0
    endif else if icont eq 1 then begin
    ;QC + QI
      iv=where(vars.vars eq 'QCLOUD')
      file=dirs.files_post[ic,iv]
      var2=read_nc_var(file,'QCLOUD',count=count,offset=offset)
      iv=where(vars.vars eq 'QICE')
      file=dirs.files_post[ic,iv]
      var2+=read_nc_var(file,'QICE',count=count,offset=offset)
    endif else if icont eq 2 then begin
    ;Latent Heat
      iv=where(vars.vars eq 'H_DIABATIC')
      file=dirs.files_post[ic,iv]
      var2=read_nc_var(file,'H_DIABATIC',count=count,offset=offset)
    endif else if icont eq 3 then begin
    ;AVOR
      iv=where(vars.vars eq 'AVOR')
      file=dirs.files_post[ic,iv]
      var2=read_nc_var(file,'AVOR',count=count,offset=offset)
    endif else if icont eq 4 then begin
    ;W
      iv=where(vars.vars eq 'W')
      file=dirs.files_post[ic,iv]
      var2=read_nc_var(file,'W',count=count,offset=offset)
    endif


;----SUBTRACT CLEAR SKY--------------------


;if icrf and strmatch(var1_str,'*RTHRA*') then begin
;
;  ix_clear=where(olr ge olr_thresh,count,complement=cloudy)
;  print,'REMOVING CLEAR SKY'
;  print,'...as columns where OLR .le. ',olr_thresh
;  print,'N-clear: ',count
;  print,'N-cloudy: ',n_elements(cloudy)
;
;  clearsky=fltarr(dims.np)
;  for iz=0,dims.np-1 do begin
;    clearsky[iz] = mean((reform(var1[*,*,iz,*]))[ix_clear],/nan,/double)
;    print,'LW clear-sky (K/d): ',strtrim(fix(dims.pres[iz]),2),' hPa, ',clearsky[iz]*3600.*24
;    var1[*,*,iz,*] -= clearsky[iz]
;  endfor
;
;endif


;----TIME BINNING--------------------


;CONCLUSION: OLR has peculiar variation in convective region (i.e., not
;  unimodal). Time averaging likely not a good idea with moving domain.


;  tbin_h=6 ; bin size (hours = time steps)
;  nt_s=1.*(nt_sel-1)/tbin_h
;
;  binvar_s=fltarr(dims.nx,dims.ny,nt_s)
;  var1_s=fltarr(dims.nx,dims.ny,dims.np,nt_s)
;  var2_s=var1_s
;
;  for ih=0,nt_s-1 do begin
;    it_sel=indgen(tbin_h)+1+ih*tbin_h
;    binvar_s[*,*,ih]=mean(binvar[*,*,it_sel],dimension=3,/nan,/double)
;    var1_s[*,*,*,ih]=mean(var1[*,*,*,it_sel],dimension=4,/nan,/double)
;    var2_s[*,*,*,ih]=mean(var2[*,*,*,it_sel],dimension=4,/nan,/double)
;  endfor
;
;  binvar=binvar_s
;  var1=var1_s
;  var2=var2_s


;----SUBSET TO CERTAIN RADIUS--------------------

    radius=binvar
    radius[*]=!values.f_nan
    for it=0,nt_test-1 do begin

      ;TC LOCATION
        ivloc = [ vloc[0,it] , vloc[1,it] ]
        tcrad = radius_tc_ll(ivloc,dims.lon,dims.lat)
        irad = where((tcrad le radius_thresh),count_rad,complement=nan)

      tcrad[nan]=!values.f_nan
      radius[*,*,it]=tcrad

    endfor

;      ibv=reform(binvar[*,*,it])
;      ibv[nan]=!values.f_nan
;      binvar[*,*,it]=ibv
;
;      iolr=reform(olr[*,*,it])
;      iolr[nan]=!values.f_nan
;      olr[*,*,it]=iolr
;
;      iolrc=reform(olrc[*,*,it])
;      iolrc[nan]=!values.f_nan
;      olrc[*,*,it]=iolrc
;
;      for iz=0,dims.np-1 do begin
;
;        iv1=reform(var1[*,*,iz,it])
;        iv1[nan]=!values.f_nan
;        var1[*,*,iz,it]=iv1
;
;        iv2=reform(var2[*,*,iz,it])
;        iv2[nan]=!values.f_nan
;        var2[*,*,iz,it]=iv2
;
;      endfor
;
;    endfor


;----BINNING--------------------


  if binvar_str eq 'pw' then begin
    min=25
    max=80
    xtitle='PW [ mm ]'
  endif else if binvar_str eq 'olr' then begin
    min=70
    max=330
    xtitle='OLR [ W m!U-2!N ]'
  endif else if binvar_str eq 'rainrate' then begin
    min=0
    max=75
    xtitle='Rainfall [ mm h!U-1!N ]'
  endif else message,'add bin var here'

  nbin=40
  delta=1.*(max-min)/nbin
  bins=findgen(nbin)*delta+min ; BIN VARIABLE

  bin1=fltarr(nbin,dims.np)
  bin1[*]=!values.f_nan
  bin2=bin1
;  bin_tmpk=bin1
;  bin_qv=bin1
  bin_olr=fltarr(nbin)
  bin_olr[*]=!values.f_nan
  bin_olrc=bin_olr
  bin_count=bin_olr

  nthresh=10;25
  for ibin=0,nbin-1 do begin
    x0 = min + delta*ibin
;    ix=where((binvar ge x0) and (binvar lt x0+delta),count)
    ix=where( ( (binvar ge x0) and (binvar lt x0+delta) and finite(radius) ) ,count)
    if count le nthresh then begin
      bins[ibin]=!values.f_nan
      continue
    endif

    bin_count[ibin]=count

    bin_olr[ibin]=mean(olr[ix],/nan,/double)
    bin_olrc[ibin]=mean(olrc[ix],/nan,/double)
    for iz=0,dims.np-1 do begin
      bin1[ibin,iz]=mean((reform(var1[*,*,iz,*]))[ix],/nan,/double)
      bin2[ibin,iz]=mean((reform(var2[*,*,iz,*]))[ix],/nan,/double)
;      bin_tmpk[ibin,iz]=mean((reform(tmpk[*,*,iz,*]))[ix],/nan,/double)
;      bin_qv[ibin,iz]=mean((reform(qv[*,*,iz,*]))[ix],/nan,/double)
    endfor
  endfor

  if binvar_str eq 'olr' then begin
    bins=reverse(temporary(bins))
    bin_count=reverse(temporary(bin_count))
    bin_olr=reverse(temporary(bin_olr))
    bin_olrc=reverse(temporary(bin_olrc))
    bin1=reverse(temporary(bin1),1)
    bin2=reverse(temporary(bin2),1)
  endif


;----PROFILES--------------------


if iprof then begin

  pw_sel1=[30,50]  ; Clear
  pw_sel2=[60,150] ; Cloudy

  ;Clear
  isel=where((bins ge pw_sel1[0]) and (bins le pw_sel1[1]))
  prof1=mean(bin1[isel,*],dimension=1,/nan,/double)

  ;Cloudy
  isel=where((bins gt pw_sel2[0]) and (bins le pw_sel2[1]))
  prof2=mean(bin1[isel,*],dimension=1,/nan,/double)

prof2-=prof1
prof1[*]=0
endif


;----SUBTRACT CLEAR SKY--------------------


;if icrf and strmatch(var1_str,'*RTHRA*') then begin
;
;  ix_clear=where((bins ge pw_thresh_crf[0]) and (bins le pw_thresh_crf[1]),count,complement=cloudy)
;  print,'REMOVING CLEAR SKY'
;  print,'...as columns where PW is b/t ',pw_thresh_crf
;  print,'N-clear: ',count
;  print,'N-cloudy: ',n_elements(cloudy)
;
;  clearsky=fltarr(dims.np)
;  for iz=0,dims.np-1 do begin
;    clearsky[iz] = mean((reform(bin1[*,iz]))[ix_clear],/nan,/double)
;    print,'LW clear-sky (K/d): ',strtrim(fix(dims.pres[iz]),2),' hPa, ',clearsky[iz]*3600.*24
;    bin1[*,iz] -= clearsky[iz]
;  endfor
;
;endif


;----STREAMFUNCTION--------------------


  ;CALCULATE SIMPLE STREAMFUNCTION IN PW-BINNED SPACE AS IN BRETHERTON ET AL. (2005)
  if icont eq 0 then begin
    psi=fltarr(nbin,dims.np)
    x0=1 & xend=nbin-1 & dir=1
    if binvar_str eq 'olr' then begin
      x0=nbin-2 & xend=0 & dir=-1
    endif
    for ix=x0,xend,dir do $
      for iz=0,dims.np-1 do $
        psi[ix,iz] = psi[ix-1,iz] + rho[iz]*bin2[ix,iz] ; kg / m2 / s
;        psi[ix,iz] = psi[ix+1,iz] + rho[iz]*bin2[ix,iz] ; kg / m2 / s
;        psi[ix,iz] = psi[ix-1,iz] + (bin_w[ix,iz]*bin_rho[ix,iz]) ; kg / m2 / s
        ;RESULTS ARE VIRTUALLY IDENTICAL FOR TAKING W*RHO BEFORE OR AFTER BINNING
;    psi *= 1e2 ; 10^-2 kg / m2 / s
    bin2=psi & psi=0
  endif


;----WRITE QV-PROFILE--------------------


  if write_qv then begin

    print,'WRITING QV-PROFILE TO ',qv_file

    if dirs.cases[ic] ne 'ctl' then message,"Sure you don't want to use CTL?"

    ;AVERAGE FROM BINNED ARRAY
;    loc_sel=where((bins ge pw_thresh_fixqv[0]) and (bins le pw_thresh_fixqv[1]),npts,NCOMPLEMENT=ncomplement)
;    print,'Count of selection: ',npts,', %: ',1d*npts/ncomplement*1e2
;    qv_prof=fltarr(dims.np)
;    for iz=0,dims.np-1 do begin
;      qv_prof[iz]=mean(reform(bin_qv[loc_sel,iz]),/nan,/double)
;;      print,'QV: ',dims.pres[iz],' hPa, ',qv_prof[iz]*1e3,' g/kg'
;    endfor

    ;AVERAGE FROM RAW 2D CARTESIAN VARIABLE
;      loc_sel=where((binvar ge pw_thresh_fixqv[0]) and (binvar le pw_thresh_fixqv[1]),npts,NCOMPLEMENT=ncomplement)
;      print,'Count of selection: ',npts,', %: ',1d*npts/ncomplement*1e2
      ;TIME AVERAGE FIRST
        qvavg=mean(qv,dimension=4,/nan,/double)
        pwavg=mean(binvar,dimension=3,/nan,/double)
      loc_sel=where((pwavg ge pw_thresh_fixqv[0]) and (pwavg le pw_thresh_fixqv[1]),npts,NCOMPLEMENT=ncomplement)
      print,'Count of selection: ',npts,', %: ',1.*npts/(npts+ncomplement)*1e2
      qv_prof=fltarr(dims.np)
      for iz=0,dims.np-1 do begin
        iqv=reform(qvavg[*,*,iz])
        qv_prof[iz]=mean(iqv[loc_sel],/nan,/double)
;        print,'QV1, QV2: ',dims.pres[iz],' hPa, ',qv_prof[iz]*1e3,' ',qv_prof2[iz]*1e3,' g/kg'
      endfor
;qv_file+='2'

    ;INTERPOLATE ONTO MODEL LEVELS
      pres_model=reform(read_nc_var(model_pres_levels_nc,'',varid='3'))
      np_mod=n_elements(pres_model)
      qv_prof_mod = interpol([qv_prof,min(qv_prof)*0.05],[dims.pres*1e2,1000],pres_model)
      for iz=0,np_mod-1 do print,'QV: ',pres_model[iz]*1e-2,' hPa, ',qv_prof_mod[iz]*1e3,' g/kg'

    openw,1,qv_file
      forstr='(f10.8)'
      printf,1,'/ '+strjoin(string(qv_prof_mod[0:7],format=forstr),',',/single)+', &'
      for i=0,3 do printf,1,strjoin(string(qv_prof_mod[8+i*8:15+i*8],format=forstr),',',/single)+', &'
      printf,1,strjoin(string(qv_prof_mod[16+3*8:np_mod-1],format=forstr),',',/single)+'/'
    close,1

    return

  ;OVERWRITE RH WITH DRY WATER VAPOR PROFILE
    overwrite_rh=0
    if overwrite_rh then begin
      for iz=0,dims.np-1 do $
        bin_qv[*,iz]=qv_prof[iz]
      for ibin=0,ibin-1 do $
        bin1[ibin,*]=calc_relh(bin_qv[ibin,*],bin_tmpk[ibin,*],dims.pres*1e2)
    endif

  endif


;----CREATE PLOTS--------------------


;if ic eq 0 then begin
;  ctl=
;endif else begin
;
;endelse

;  if strmatch(var1_str,'*RTHRA*') then begin
;    setmax=4
;    setmin=-4
;  endif

  tc_figspecs, var1_str, figspecs, setmax=setmax, setmin=setmin

  figdir=dirs.figdir+'/comp_cross/'
  spawn,'mkdir '+figdir,tmp,tmpe
  figname=figdir+dirs.cases[ic]+'_'+var1_str

  if strmatch(var1_str,'*RTHRA*') then begin
    if icrf then begin
      figname+='_crf'
;figname+='_csk'
;      radtag='CRF'
      radtag='LW-CRF'
;radtag='LW-CSK'
    endif else $
      if iadd_sw then radtag='LW+SW' else radtag=strmid(var1_str,1,2,/reverse)
    figspecs.cbar_tag=radtag+' '+figspecs.cbar_tag
    irad=1
  endif else irad=0
;irad=0

  if icont eq 0 then begin
  ;PSI
    cvar_str='psi'
    cvar=bin2*1e2
    cint=0.5
    clevs=indgen(200)*cint+cint
  endif else if icont eq 1 then begin
  ;CLOUD
    cvar_str='qc'
    cvar=bin2*1e6 ; kg/kg --> mg/kg
;    clevs=[1,5,50,200,500]
    clevs=[5,10,25,50,200,500]
  endif else if icont eq 2 then begin
  ;Latent heat
    cvar_str='lh'
    cvar=bin2*3600d*24 ; K/d
;for i=0,0 do cvar=gauss_smooth(temporary(cvar),[3,0],/edge_truncate)
    cint=4
;    clevs=indgen(50)*cint+cint
    clevs=[1,2,5,10,50,100,500,1000,2000,3000]
  endif else if icont eq 3 then begin
  ;AVOR
    cvar_str='avor'
    cvar=bin2 ; 10^-5 /s
;for i=0,0 do cvar=gauss_smooth(temporary(cvar),[3,0],/edge_truncate)
    clevs=[1,2,5,10,50,100,500,1000,2000,3000]
  endif else if icont eq 4 then begin
  ;W
    cvar_str='w'
    cvar=bin2*1e2 ; cm/s
;for i=0,0 do cvar=gauss_smooth(temporary(cvar),[3,0],/edge_truncate)
    clevs=[1,2,5,10,50,100,500,1000,2000,3000]
  endif

;PROFILES
if iprof then begin
  figname_prof=figname+'_'+binvar_str+'_prof_'+hr_tag
  figspecs=create_struct(figspecs,'figname',figname_prof)
  prof1*=3600d*24
  prof2*=3600d*24
;stop
  wrf_tc_comp_prof, prof1, prof2, dims.pres, figspecs, xtitle=var1_str
  print,'Done!!'
  exit
endif


;CROSS SECTIONS

  figname_cross=figname+'_'+binvar_str+'_'+cvar_str+'_'+hr_tag
  figspecs=create_struct(figspecs,'figname',figname_cross)

  wrf_tc_comp_cross, bin1, bins, dims.pres, figspecs, cvar=cvar, xtitle=xtitle, $
;    clevs=clevs, olr=bin_olr, csol=bin_olrc, irad=irad
    clevs=clevs, olr=bin_count, irad=irad ; replace OLR with COUNT

endfor ; icase

print,'DONE!!'
endfor ; iadd (12 h increments)
end
