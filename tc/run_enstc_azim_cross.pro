; 
; Plot azimuthal cross sections from ensemble TC output.
;
; James Ruppert
; 12/26/21
; 
pro run_enstc_azim_cross

;tcname='maria'
tcname='haiyan'
case_str='ctl';'ncrf';'ctl'

;if tcname eq 'maria' then begin
;  tcyear='2017'
;  hurdat=read_hurdat(tcname,tcyear)
;endif else if tcname eq 'haiyan' then begin
;  tcyear='2013'
;  hurdat=read_jtwcdat(tcname)
;endif

dom='d02'
tc_ens_config, tcname, case_str, dom, dirs=dirs, dims=dims, vars=vars;, /verbose
dirs.figdir+=tcname+'/'

;TIME SPECS
  nt_full=dims.nt
  npd=dims.npd
  nhrs=1.*nt_full*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nhrs)

;RADIUS SUBSAMPLE
  rad_sel=[0,1300];1800];800] ; km , max=1300 km for redux
  ;USE ENTIRE MODEL DOMAIN FOR REMOVING MEANS

;for iex=0,dirs.nens-1 do begin
for iex=0,0 do begin

;AZIM FILES
;  if tcname eq 'maria'  then hr_sel=[24,120]
;  if tcname eq 'haiyan' then hr_sel=[24,144]
;  if strmatch(subdir,'redux*') then begin
;    if tcname eq 'maria'  then hr_sel=[0,144]
;    if tcname eq 'haiyan' then hr_sel=[0,168]
;  endif
;  t_ind=where((time_hrs ge hr_sel[0]) and (time_hrs le hr_sel[1]))
;  t_ind=t_ind[where(t_ind le nt_full-1)]
;  nt_sav=n_elements(t_ind)


;----PLOT OPTIONS--------------------


var1_str='RTHRATLW'
;var1_str='RTHRATTOT'
;var1_str='RTHRATSW'
;var1_str='RH';'QVAPOR';
;var1_str='W'
;var1_str='T'
;var1_str='theta_e'
;var1_str='LR'
;var1_str='AVOR'
;var1_str='H_DIABATIC'
;var1_str='madv'
;var1_str='u_rad'
;var1_str='v_tan'
;var1_str='mfr' ; Vertical mass flux
;var1_str='mfz' ; Radial mass flux

iremove_ctl=0 ; remove CTL?

remove_azmn=0
remove_azmn=1

do_smooth=1 ; smooth var1?

ihght=1 ; convert vertical dimension to height?
  ztop=15. ; km

icrf=0
;icrf=1

;for irun=0,1 do begin
;if irun eq 0 then begin
;  remove_azmn=1
;  var1_str='RTHRATLW'
;endif else begin
;  remove_azmn=0
;  var1_str='H_DIABATIC'
;endelse

cont_var=1 ; Contour var: 0-psi, 1-cld, 2-u_rad, 3-w, 4-m, 5-LH, 6-AVOR, 7-v_tan, 8-inertial stability
var_1d=-1 ; 1D var in sub-panel: -1: no plot, 0-OLR, 1-PW, 2-v_tan

;AZIMUTHAL SUBSET?
nazim=8 ; octants
nazim=1 ; full 360


;----SPECIAL OPTIONS--------------------


overlay_seq=0 ; Overlay fields from HEW's SEQ model?
  seq_subdir='';'total';'radanom_cool';ftwarm';xpbl';_pbl';
;  seq_subdir='mp'
;  seq_cvar='UUU'
;  seq_cvar='VG_OUT'
  seq_cvar='WWW'
;  seq_cvar='DVT'

;if overlay_seq=1
plot_forcing_v1=1 ; Overwrite var1 with SEQ forcing (read in from SEQ directory)?

;TIME SCALE TO SATURATION (USING MADV INSTEAD OF LCL NOW)
  ;if overlay_seq=1
  ilcl=0 ; calculate time to saturation from SEQ motion?

write_force_seq=0 ; write var out for HEW's SEQ code?
  snd_file='/work/06040/tg853394/stampede2/tc_stuff/FCIRC_19_local/DATA/WRF_SND.txt'
;  seq_file='/work/06040/tg853394/stampede2/tc_stuff/FCIRC_19_local/DATA/Test_0/seq_force_lwcrf.txt'
  seq_file='/work/06040/tg853394/stampede2/tc_stuff/FCIRC_19_local/DATA/Test_0/seq_force_lwradanom.txt'
;  seq_file='/work/06040/tg853394/stampede2/tc_stuff/FCIRC_19_local/DATA/Test_0/seq_force_lwradanom_'+hr_tag_check+'.txt'
;  seq_file='/work/06040/tg853394/stampede2/tc_stuff/FCIRC_19_local/DATA/Test_0/seq_force_mpheat.txt'

write_crf=0 ; write profile of CRF?
  add_icrf=1 ; Replace CRF with ICRF (not to CTL)?
;  rad_thresh=[0,250] ; Radius range to average over (km) [ Used 0-250 for ICRF_* ]
;;  rad_thresh=[0,300] ; Radius range to average over (km)
;  crf_file=dirs.figdir+'crf_profile_maria_ctl_'+string(rad_thresh[0],format='(i3.3)')+'-'+$
;    string(rad_thresh[1],format='(i3.3)')+'km_radius_'+string(hr_plot[0],format='(i3.3)')+'-'+string(hr_plot[1],format='(i3.3)')+'hr.txt'
;  model_pres_levels_nc='/scratch/06040/tg853394/tc/output/static_nest/maria/ctl/p_total_16-1200-17-1200_tavg.nc'
;
write_qv=0 ; write qv-profile for dry test?
;;  if write_qv and ((var1_str ne 'RH') or (binvar_str ne 'pw')) then $
;;    message,'Need RH, PW to write qv-profile.'
;;;  pw_thresh_fixqv=[40,55]
;;  pw_thresh_fixqv=[55,70]
;;  qv_file=dirs.figdir+'wvapor_profile_edouard_'+string(pw_thresh_fixqv[0],format='(i2.2)')+'-'+$
;;    string(pw_thresh_fixqv[1],format='(i2.2)')+'mm.txt'


;----SETUP FOR PAPER FIGURE 5--------------------

;OPTIONS ARE ALL OVERWRITTEN HERE

single=1 ; single panel case

;FOR ipanel=2,3 DO BEGIN
FOR ipanel=single,single DO BEGIN

  do_smooth=1 ; smooth var1?
  ihght=1 & ztop=15. ; km
  icrf=0
  nazim=1 ; full 360
  var_1d=-1

  seq_subdir=''
  write_force_seq=0 ; switch on to write SEQ using IR'
;  write_force_seq=1

  if write_force_seq then begin
    ipanel=2
    ihght=0
    print,'WRITING SEQ FORCING!'
  endif else begin
    print,'DOING PANEL ',strtrim(ipanel,2)
  endelse

;case ipanel of
;  1: begin ; V_TAN, CLOUD
;       var1_str='v_tan'
;       remove_azmn=0
;       cont_var=1 ; cloud
;       overlay_seq=0
;     end
;  2: begin ; IR', W_CRF
;       var1_str='RTHRATLW'
;       remove_azmn=1
;       cont_var=7 ; vtan, will overlay RMW
;       overlay_seq=1
;         plot_forcing_v1=1
;         seq_cvar='WWW'
;         ilcl=0
;     end
;  3: begin ; RH, TAU_SATURATION
;       var1_str='RH'
;       remove_azmn=0
;       cont_var=1 ; cloud, will be overwritten by tau_sat
;       overlay_seq=1
;         plot_forcing_v1=0
;         seq_cvar='WWW'
;         ilcl=1
;     end
;endcase


;----READ VARS--------------------


;AZIMUTHAL SPECS
  nrad=267
  if tcname eq 'haiyan' then nrad=433
  naz=360

print,'VAR: ',var1_str

if (cont_var eq 0 or strmatch(var1_str,'mf*') or ihght or overlay_seq) then dopsi=1 else dopsi=0

;for ic=0,dirs.nc-1 do begin
;for ic=0,dirs.nc-1,2 do begin
for ic=0,0 do begin

  print,'CASE: ',strupcase(dirs.cases[ic])

  hr0=0
  if strmatch(dirs.cases[ic],'*36h*') then hr0=36
  if strmatch(dirs.cases[ic],'*48h*') then hr0=48
  if strmatch(dirs.cases[ic],'*60h*') then hr0=60
  if strmatch(dirs.cases[ic],'*72h*') then hr0=72
  if strmatch(dirs.cases[ic],'*84h*') then hr0=84
  if strmatch(dirs.cases[ic],'*96h*') then hr0=96
  if strmatch(dirs.cases[ic],'lwcrf*') then hr0=36
;  if strmatch(dirs.cases[ic],'lwcrf') then hr0=49;36
;  if strmatch(dirs.cases[ic],'axisym') then hr0=49;36

;if strmatch(subdir,'redux*') and $
;  total(strmatch(['ncrf_36h','lwcrf','axisym'],dirs.cases[ic])) then begin
;
;;    hr0=0
;;    nt_sav=24
;    hr_fil='49-72hr'
;    hr_tag_plot='049-072hr'
;    nt_plot=24
;    t_ind_plot=indgen(nt_plot)+0
;;  endif else begin
;;  ;CTL
;;    ;REDUX
;;    hr_fil='0-144hr'
;;    hr_tag_plot='049-096hr'
;;    nt_plot=48;24
;;    t_ind_plot=indgen(nt_plot)+(fix(strmid(hr_tag_plot,0,3)) - fix(strmid(hr_fil,0,2)))
;;  endelse
;
;endif else begin

  ;TIME SELECTION
    t_offset=max([0,hr0-hr_sel[0]])
    ut_offset=max([0,hr_sel[0]-hr0])
    nt_test=nt_sav-t_offset
    t_ind_test=indgen(nt_test)+ut_offset
    hrs_test=time_hrs[t_ind_test+hr0]
    hr_fil=strtrim(hrs_test[0],2)+'-'+strtrim(hrs_test[nt_test-1],2)+'hr'
    t_ind_plot=where(hrs_test ge hr_plot[0] and hrs_test le hr_plot[1],nt_plot)
    hrs_plot=hrs_test[t_ind_plot]
    hr_tag_plot=string(hrs_plot[0],format='(i3.3)')+'-'+string(hrs_plot[nt_plot-1],format='(i3.3)')+'hr'

;endelse

  if hr_tag_plot ne hr_tag_check then begin
    print,'Not enough times in test, so skipping...'
    continue
  endif

  it0 = t_ind_plot[0]

  ;RADIUS SUBSAMPLE
    file=dirs.casedir[ic]+'azim_SLP_'+hr_fil+'.nc'
    radius=read_nc_var(file,'radius')
    azimuth=read_nc_var(file,'azmiuth')
    xrad=where(radius ge rad_sel[0] and radius le rad_sel[1],nrad)
    radius=radius[xrad]

  ;READ STORM MOTION (CALCULATED IN RUN_WRF_TC_TRACKS)
    motion_file=dirs.casedir[ic]+'storm_motion.txt'
    openr,1,motion_file
      readf,1,tmp_nt
      tmp_it_test=intarr(tmp_nt)
      readf,1,tmp_it_test
      motion_x=fltarr(tmp_nt) & motion_y=motion_x
      tc_lon=motion_x & tc_lat=motion_x
      readf,1,motion_x
      readf,1,motion_y
      readf,1,tc_lon
      readf,1,tc_lat
    close,1
    motion_x=motion_x[t_ind_plot]
    motion_y=motion_y[t_ind_plot]
    tc_lon=tc_lon[t_ind_plot]
    tc_lat=tc_lat[t_ind_plot]

  ;2D VARS
    count=[nrad,naz,1,nt_plot] & offset=[xrad[0],0,0,it0] ; x,y,z,t

  ;SLP
    file=dirs.casedir[ic]+'azim_SLP_'+hr_fil+'.nc'
    slp=reform(read_nc_var(file,'SLP',count=count,offset=offset))

  var1d=slp & var1d[*]=!values.f_nan
  var1d2=var1d

  if var_1d eq 0 then begin
  ;OLR
    file=dirs.casedir[ic]+'azim_OLR_'+hr_fil+'.nc'
    var1d=reform(read_nc_var(file,'OLR',count=count,offset=offset))
  ;OLRC
    file=dirs.casedir[ic]+'azim_OLRC_'+hr_fil+'.nc'
    var1d2=reform(read_nc_var(file,'OLRC',count=count,offset=offset))
  endif else if var_1d eq 1 then begin
  ;PW
    file=dirs.casedir[ic]+'azim_PW_'+hr_fil+'.nc'
    var1d=reform(read_nc_var(file,'PW',count=count,offset=offset))
  endif else if var_1d eq 2 then begin
  ;v_tan - Tangential Wind
      file=dirs.casedir[ic]+'azim_U10_'+hr_fil+'.nc'
      u=reform(read_nc_var(file,'U10',count=count,offset=offset))
      file=dirs.casedir[ic]+'azim_V10_'+hr_fil+'.nc'
      v=reform(read_nc_var(file,'V10',count=count,offset=offset))
      ;SUBTRACT STORM MOTION
      for it=0,nt_plot-1 do begin
        u[*,*,it] -= motion_x[it]
        v[*,*,it] -= motion_y[it]
      endfor
      wnd_azim=azim_wind_conv(u,v,azimuth) & u=0 & v=0
      var1d=wnd_azim.v_tan
  endif

  ;3D VARS

    count=[nrad,naz,dims.np,nt_plot] & offset=[xrad[0],0,0,it0] ; x,y,z,t

  ;VAR1

    dummy=fltarr(count)
    dummy2=fltarr(count)
    dummy3=fltarr(count)

    var2=dummy

    if var1_str eq 'RH' then begin

      ;NEED LARGE PRES ARRAY FOR RH
        prx=fltarr(nrad,naz,dims.np,nt_plot)
        for iz=0,dims.np-1 do prx[*,*,iz,*]=dims.pres[iz]*1e2

      file=dirs.casedir[ic]+'azim_QVAPOR_'+hr_fil+'.nc'
      qv=read_nc_var(file,'QVAPOR',count=count,offset=offset)
      file=dirs.casedir[ic]+'azim_T_'+hr_fil+'.nc'
      tmpk=read_nc_var(file,'T',count=count,offset=offset)

      var1=calc_relh(qv,tmpk,prx)

      if ~write_qv then qv=0
      tmpk=0

    endif else if var1_str eq 'RTHRATTOT' then begin

      file=dirs.casedir[ic]+'azim_RTHRATLW_'+hr_fil+'.nc'
      var1=read_nc_var(file,'RTHRATLW',count=count,offset=offset)
      file=dirs.casedir[ic]+'azim_RTHRATSW_'+hr_fil+'.nc'
      var1+=read_nc_var(file,'RTHRATSW',count=count,offset=offset)

    endif else if var1_str eq 'LR' then begin

      file=dirs.casedir[ic]+'azim_T_'+hr_fil+'.nc'
      tmpk=read_nc_var(file,'T',count=count,offset=offset)

      pr=tmpk
      for iz=0,dims.np-1 do pr[*,*,iz,*]=dims.pres[iz]
      var1=theta(tmpk-273.15,pr) & tmpk=0 & pr=0

    endif else if var1_str eq 'theta_e' then begin

      file=dirs.casedir[ic]+'azim_T_'+hr_fil+'.nc'
      tmpk=read_nc_var(file,'T',count=count,offset=offset)
      file=dirs.casedir[ic]+'azim_QVAPOR_'+hr_fil+'.nc'
      qv=read_nc_var(file,'QVAPOR',count=count,offset=offset)
      file=dirs.casedir[ic]+'azim_QCLOUD_'+hr_fil+'.nc'
      cld=read_nc_var(file,'QCLOUD',count=count,offset=offset)
      file=dirs.casedir[ic]+'azim_QICE_'+hr_fil+'.nc'
      cld+=read_nc_var(file,'QICE',count=count,offset=offset)

      pr=tmpk
      for iz=0,dims.np-1 do pr[*,*,iz,*]=dims.pres[iz]

      var1=theta_e(tmpk,pr*1e2,qv,qv+cld) & tmpk=0 & pr=0 & qv=0 & cld=0

    endif else if var1_str eq 'u_rad' or var1_str eq 'v_tan' or strmatch(var1_str,'mf*') or var1_str eq 'madv' then begin

      file=dirs.casedir[ic]+'azim_U_'+hr_fil+'.nc'
      u=read_nc_var(file,'U',count=count,offset=offset)
      file=dirs.casedir[ic]+'azim_V_'+hr_fil+'.nc'
      v=read_nc_var(file,'V',count=count,offset=offset)

      ;SUBTRACT STORM MOTION
      for it=0,nt_plot-1 do begin
        u[*,*,*,it] -= motion_x[it]
        v[*,*,*,it] -= motion_y[it]
      endfor

      wnd_azim=azim_wind_conv(u,v,azimuth) & u=0 & v=0
      if var1_str eq 'u_rad' then var1=wnd_azim.u_rad
      if var1_str eq 'v_tan' then var1=wnd_azim.v_tan

      if var1_str eq 'madv' then begin
        file=dirs.casedir[ic]+'azim_QVAPOR_'+hr_fil+'.nc'
        qv=read_nc_var(file,'QVAPOR',count=count,offset=offset)
        var1=qv
      endif

    endif else begin

      file=dirs.casedir[ic]+'azim_'+var1_str+'_'+hr_fil+'.nc'
      var1=read_nc_var(file,var1_str,count=count,offset=offset)

    endelse

    ;CONTOURS VARIABLE
    if cont_var eq 2 or cont_var eq 4 or cont_var eq 7 or cont_var eq 8 then begin

      file=dirs.casedir[ic]+'azim_U_'+hr_fil+'.nc'
      u=read_nc_var(file,'U',count=count,offset=offset)
      file=dirs.casedir[ic]+'azim_V_'+hr_fil+'.nc'
      v=read_nc_var(file,'V',count=count,offset=offset)

      ;SUBTRACT STORM MOTION
      for it=0,nt_plot-1 do begin
        u[*,*,*,it] -= motion_x[it]
        v[*,*,*,it] -= motion_y[it]
      endfor

      wnd_azim=azim_wind_conv(u,v,azimuth) & u=0 & v=0
      if cont_var eq 2 then var2=wnd_azim.u_rad
      if cont_var eq 7 or cont_var eq 8 then var2=wnd_azim.v_tan
      if cont_var eq 4 then begin
        ;fcor=2.*7.292e-5*sin(13.*!pi/180)
        var2=wnd_azim.u_rad
        for it=0,nt_plot-1 do begin
          fcor=2.*7.292e-5*sin(tc_lat[it]*!pi/180)
          for irad=0,nrad-1 do var2[irad,*,*,it] = (radius[irad]*1e3)*wnd_azim.v_tan[irad,*,*,it] + fcor*((radius[irad]*1e3)^2)/2
        endfor
      endif

    endif else if cont_var eq 3 then begin

      file=dirs.casedir[ic]+'azim_W_'+hr_fil+'.nc'
      var2=read_nc_var(file,'W',count=count,offset=offset)

    endif else if cont_var eq 5 then begin

      file=dirs.casedir[ic]+'azim_H_DIABATIC_'+hr_fil+'.nc'
      var2=read_nc_var(file,'H_DIABATIC',count=count,offset=offset)

    endif else if cont_var eq 6 then begin

      file=dirs.casedir[ic]+'azim_AVOR_'+hr_fil+'.nc'
      var2=read_nc_var(file,'AVOR',count=count,offset=offset)

    endif

    ;LCL HEIGHT FOR SEQ
    ;See Romps 2017, JAS for formulae
    if overlay_seq and ilcl then begin
      ;NEED LARGE PRES ARRAY FOR RH
        prx=fltarr(nrad,naz,dims.np,nt_plot)
        for iz=0,dims.np-1 do prx[*,*,iz,*]=dims.pres[iz]*1e2
      file=dirs.casedir[ic]+'azim_QVAPOR_'+hr_fil+'.nc'
      qv=read_nc_var(file,'QVAPOR',count=count,offset=offset)
      file=dirs.casedir[ic]+'azim_T_'+hr_fil+'.nc'
      tmpk=read_nc_var(file,'T',count=count,offset=offset)
      rh=calc_relh(qv,tmpk,prx);,/noice)
;;var1=rh
;;var1_str='RH'
;      prx=0 ;& qv=0
;      ;LCL rough estimate from Espy's formula
;      tdew=rh2tdew(tmpk,rh)
;;      dummy2=125.*(tmpk-tdew)
;      ;LCL according to Bolton (1980)
;      cpm=(1.-qv)*1004.+qv*(1418.+461) ; CPM according to Romps
;      dummy2 = cpm/9.81 * ( tmpk - 55. - (((1./(tmpk-55)) - (alog10(rh*1e-2)/2840.))^(-1)) )
;USE VERTICAL ADVECTION INSTEAD
      qvs=mixr_sat(tmpk,prx*1e-2)*1e-3
      satdef=qvs-qv
      dummy2=satdef
      dummy3=qv
      rh=0 & tmpk=0 & tdew=0 & qv=0
    endif

    ;NEED THETA0 FOR SEQ
    if overlay_seq or write_force_seq then begin
      file=dirs.casedir[ic]+'azim_T_'+hr_fil+'.nc'
      tmpk=read_nc_var(file,'T',count=count,offset=offset)
      pr=tmpk
      for iz=0,dims.np-1 do pr[*,*,iz,*]=dims.pres[iz]
      dummy=theta(tmpk-273.15,pr) & tmpk=0 & pr=0
    endif

    ;REMOVE CLEAR SKY
    if strmatch(var1_str,'RTHRA*') and icrf then begin

      if var1_str eq 'RTHRATTOT' then begin
        file=dirs.casedir[ic]+'azim_RTHRATLWC_'+hr_fil+'.nc'
        var1-=read_nc_var(file,'RTHRATLWC',count=count,offset=offset)
;var1=read_nc_var(file,'RTHRATLWC',count=count,offset=offset)
        file=dirs.casedir[ic]+'azim_RTHRATSWC_'+hr_fil+'.nc'
        var1-=read_nc_var(file,'RTHRATSWC',count=count,offset=offset)
;var1+=read_nc_var(file,'RTHRATSWC',count=count,offset=offset)
      endif else begin
        file=dirs.casedir[ic]+'azim_'+var1_str+'C_'+hr_fil+'.nc'
        var1-=read_nc_var(file,var1_str+'C',count=count,offset=offset)
;var1=read_nc_var(file,var1_str+'C',count=count,offset=offset)
      endelse

    endif

    ;DENSITY FOR PSI
    if dopsi then begin
    ;W
      file=dirs.casedir[ic]+'azim_W_'+hr_fil+'.nc'
      w=read_nc_var(file,'W',count=count,offset=offset)
    ;U_RAD
      file=dirs.casedir[ic]+'azim_U_'+hr_fil+'.nc'
      u=read_nc_var(file,'U',count=count,offset=offset)
      file=dirs.casedir[ic]+'azim_V_'+hr_fil+'.nc'
      v=read_nc_var(file,'V',count=count,offset=offset)
      wnd_azim=azim_wind_conv(u,v,azimuth) & u=0 & v=0
      u_rad=wnd_azim.u_rad
    ;DENSITY
      file=dirs.casedir[ic]+'azim_T_'+hr_fil+'.nc'
      tmpk=read_nc_var(file,'T',count=count,offset=offset)
      file=dirs.casedir[ic]+'azim_QVAPOR_'+hr_fil+'.nc'
      qvt=read_nc_var(file,'QVAPOR',count=count,offset=offset)
      tvirt = tmpk*(1.+0.61*qvt)
      tmpk=0 & qvt=0
      rho=w
      for iz=0,dims.np-1 do rho[*,*,iz,*] = dims.pres[iz]*1e2 / ( 287. * tvirt[*,*,iz,*] )
      if strmatch(var1_str,'mf*') then var1=var2
    endif 

    if cont_var eq 1 then begin
    ;QC + QI
      file=dirs.casedir[ic]+'azim_QCLOUD_'+hr_fil+'.nc'
      var2=read_nc_var(file,'QCLOUD',count=count,offset=offset)
      file=dirs.casedir[ic]+'azim_QICE_'+hr_fil+'.nc'
      var2+=read_nc_var(file,'QICE',count=count,offset=offset)
    endif

;SAV VARS FOR AZIMUTHAL LOOP
  slpsav=slp
  var1dsav=var1d
  var1d2sav=var1d2
  var1sav=var1
  var2sav=var2
  dummysav=dummy
  dummy2sav=dummy2
  dummy3sav=dummy3
  if dopsi then begin
    wsav=w
    u_radsav=u_rad
    rhosav=rho
  endif

;SUBSET AS A FUNCTION OF AZIMUTH
for iazim=0,nazim-1 do begin

  if nazim gt 1 then print,'IAZIM = ',iazim
  azimsel = indgen(naz/nazim) + iazim*naz/nazim

  ;AZIMUTHALLY AVERAGE
    slp=mean(slpsav[*,azimsel,*],dimension=2,/nan,/double)
;    pw=mean(temporary(pw),dimension=2,/nan,/double)
    var1d=mean(var1dsav[*,azimsel,*],dimension=2,/nan,/double)
    var1d2=mean(var1d2sav[*,azimsel,*],dimension=2,/nan,/double)
    var1=mean(var1sav[*,azimsel,*,*],dimension=2,/nan,/double)
    var2=mean(var2sav[*,azimsel,*,*],dimension=2,/nan,/double)
    dummy=mean(dummysav[*,azimsel,*,*],dimension=2,/nan,/double)
    dummy2=mean(dummy2sav[*,azimsel,*,*],dimension=2,/nan,/double)
    dummy3=mean(dummy3sav[*,azimsel,*,*],dimension=2,/nan,/double)
;    cld=mean(temporary(cld),dimension=2,/nan,/double)
    if dopsi then begin
      w=mean(wsav[*,azimsel,*,*],dimension=2,/nan,/double)
      u_rad=mean(u_radsav[*,azimsel,*,*],dimension=2,/nan,/double)
      rho=mean(rhosav[*,azimsel,*,*],dimension=2,/nan,/double)
    endif

  ;TIME AVERAGE
    slp=mean(temporary(slp),dimension=2,/nan,/double)
;    pw=mean(temporary(pw),dimension=2,/nan,/double)
    var1d=mean(temporary(var1d),dimension=2,/nan,/double)
    var1d2=mean(temporary(var1d2),dimension=2,/nan,/double)
    var1=mean(temporary(var1),dimension=3,/nan,/double)
    var2=mean(temporary(var2),dimension=3,/nan,/double)
    dummy=mean(temporary(dummy),dimension=3,/nan,/double)
    dummy2=mean(temporary(dummy2),dimension=3,/nan,/double)
    dummy3=mean(temporary(dummy3),dimension=3,/nan,/double)
;    cld=mean(temporary(cld),dimension=3,/nan,/double)
    if dopsi then begin
      w=mean(temporary(w),dimension=3,/nan,/double)
      u_rad=mean(temporary(u_rad),dimension=3,/nan,/double)
      rho=mean(temporary(rho),dimension=3,/nan,/double)
    endif

  ;REMOVE AZIMUTHAL MEAN
  if remove_azmn then begin
    print,'REMOVING AZIM MEAN!!'
    radall=fltarr(nrad,dims.np)
    for ip=0,dims.np-1 do radall[*,ip]=radius
    totrad=total(radius,/double,/nan)
    var_mn=total(var1*radall,1,/double,/nan)/totrad
    for ip=0,dims.np-1 do var1[*,ip]-=var_mn[ip]
  endif

  ;REMOVE BELOW-GROUND POINTS
  for ip=0,dims.np-1 do begin
    nan=where(slp lt dims.pres[ip],count)
    if count gt 0 then begin
      var1[nan,ip]=!values.f_nan
      var2[nan,ip]=!values.f_nan
;      cld[nan,ip]=!values.f_nan
      dummy[nan,ip]=!values.f_nan
      dummy2[nan,ip]=!values.f_nan
      dummy3[nan,ip]=!values.f_nan
      if dopsi then begin
        w[nan,ip]=!values.f_nan
        u_rad[nan,ip]=!values.f_nan
        rho[nan,ip]=!values.f_nan
      endif
    endif
  endfor

  ;LAPSE RATE
  if var1_str eq 'LR' then begin
    th=var1
    for ir=0,nrad-1 do $
      var1[ir,*]=deriv(dims.pres,reform(dummy[ir,*]))
    var1*=-1.e2 ; 10^-2 K / hPa
  endif

  if var1_str eq 'madv' then begin
    for ip=0,dims.np-1 do $
      var1[*,ip] = -1. * var2[*,ip] * deriv(radius,reform(var1[*,ip]))
  endif

  ;REMOVE CTL?
  if iremove_ctl then begin
;    print,'ASSUMING IC=0 IS CTL!'
;    if dirs.cases[ic] eq 'ctl' then $
    if dirs.cases[ic] eq 'lwcrf2' then $
;    if ic eq 0 then $
      varctl=var1 else $
      var1=varctl-var1
;      var1-=varctl
  endif


;----STREAMFUNCTION--------------------
;----ALSO PRODUCING 2D HEIGHT HERE--------------------


  if dopsi then begin

    ;HEIGHT FROM HYDROSTATIC
      z = fltarr(nrad,dims.np) ; m
      z[*]=!values.f_nan
      dz=z
      for ir=0,nrad-1 do begin
        ;TAKE CARE OF BELOW-GROUND POINTS
          kdo=where(dims.pres lt slp[ir],nkdo)
          ik=kdo[0]
          z[ir,ik] = (slp[ir]-dims.pres[ik])*1e2 / 9.81 / rho[ir,ik]
        for ik=kdo[1],kdo[nkdo-1] do $
          z[ir,ik] = z[ir,ik-1] + (dims.pres[ik-1]-dims.pres[ik])*1e2 / 9.81 / mean(rho[ir,ik-1:ik],/double,/nan)
        dz[ir,kdo]=deriv(reform(z[ir,kdo]))
      endfor

    ;INTEGRAL OVER R
    psir=fltarr(nrad,dims.np)
;    psir[*]=!values.f_nan
    mfr=psir
    idr=(radius[1]-radius[0])*1e3 ; m
    for ir=1,nrad-1 do $
      for iz=0,dims.np-1 do begin
        mfr[ir,iz]  = (radius[ir]*1e3) * rho[ir,iz] * w[ir,iz] ; JUST MASS FLUX [ kg / m s ]
        psir[ir,iz] = psir[ir-1,iz] + ( mfr[ir,iz]*idr ) ; kg / s
      endfor

    ;INTEGRAL OVER Z
    psiz=fltarr(nrad,dims.np)
    mfz=psiz
    ;PSI_Z should be zero at lowest model level
    for ir=0,nrad-1 do $
      for iz=1,dims.np-2 do begin
        mfz[ir,iz]  = (radius[ir]*1e3) * rho[ir,iz] * u_rad[ir,iz] ; JUST MASS FLUX [ kg / m s ]
        psiz[ir,iz] = psiz[ir,iz-1] - ( mfz[ir,iz] * dz[ir,iz] ) ; kg / s
      endfor

    psi = 0.5*(psir + psiz)
    psi *= 1e-8 ; 10^8 kg / s

    if cont_var eq 0 then var2=psi

    psi=0 & psir=0 & psiz=0 & u_rad=0

    if strmatch(var1_str,'mf*') then begin
      if var1_str eq 'mfr' then begin
        var1=mfr
        for i=0,1 do var1=smooth(var1,[3,0],/edge_truncate,/nan)
      endif else if var1_str eq 'mfz' then $
        var1=mfz
      if iremove_ctl then begin
        if dirs.cases[ic] eq 'ctl' then $
          varctl=var1 else $
          var1-=varctl
      endif
    endif

  ;NOW DOING THE VERTICAL REGRIDDING AT THE BOTTOM OF SCRIPT TO INCORPORATE SEQ VARIABLES

;    ;CONVERT VERTICAL AXIS TO HEIGHT
;    if ihght then begin
;      dz=0.25
;      nz=ztop/dz;+1
;      hght=findgen(nz)*dz+dz
;      var1_z=fltarr(nrad,nz) & var1_z[*]=!values.f_nan
;      var2_z=var1_z
;      for ir=0,nrad-1 do begin
;        kdo=where(finite(reform(z[ir,*])))
;        p_out=interpol(dims.pres[kdo],reform(z[ir,kdo])*1e-3,hght,/nan,/quadratic)
;        var1_z[ir,*]=interpol(reform(var1[ir,kdo]),dims.pres[kdo],p_out,/nan,/quadratic)
;        var2_z[ir,*]=interpol(reform(var2[ir,kdo]),dims.pres[kdo],p_out,/nan,/quadratic)
;      endfor
;      var1=var1_z
;      var2=var2_z
;    endif

  endif


;----WRITE QV-PROFILE--------------------


;  if write_qv then begin
;
;    print,'WRITING QV-PROFILE TO ',qv_file
;
;    if dirs.cases[ic] ne 'ctl' then message,"Sure you don't want to use CTL?"
;
;    ;AVERAGE FROM BINNED ARRAY
;;    loc_sel=where((bins ge pw_thresh_fixqv[0]) and (bins le pw_thresh_fixqv[1]),npts,NCOMPLEMENT=ncomplement)
;;    print,'Count of selection: ',npts,', %: ',1d*npts/ncomplement*1e2
;;    qv_prof=fltarr(dims.np)
;;    for iz=0,dims.np-1 do begin
;;      qv_prof[iz]=mean(reform(bin_qv[loc_sel,iz]),/nan,/double)
;;;      print,'QV: ',dims.pres[iz],' hPa, ',qv_prof[iz]*1e3,' g/kg'
;;    endfor
;
;    ;AVERAGE FROM RAW 2D CARTESIAN VARIABLE
;;      loc_sel=where((binvar ge pw_thresh_fixqv[0]) and (binvar le pw_thresh_fixqv[1]),npts,NCOMPLEMENT=ncomplement)
;;      print,'Count of selection: ',npts,', %: ',1d*npts/ncomplement*1e2
;      ;TIME AVERAGE FIRST
;        qvavg=mean(qv,dimension=4,/nan,/double)
;        pwavg=mean(binvar,dimension=3,/nan,/double)
;      loc_sel=where((pwavg ge pw_thresh_fixqv[0]) and (pwavg le pw_thresh_fixqv[1]),npts,NCOMPLEMENT=ncomplement)
;      print,'Count of selection: ',npts,', %: ',1.*npts/(npts+ncomplement)*1e2
;      qv_prof=fltarr(dims.np)
;      for iz=0,dims.np-1 do begin
;        iqv=reform(qvavg[*,*,iz])
;        qv_prof[iz]=mean(iqv[loc_sel],/nan,/double)
;;        print,'QV1, QV2: ',dims.pres[iz],' hPa, ',qv_prof[iz]*1e3,' ',qv_prof2[iz]*1e3,' g/kg'
;      endfor
;;qv_file+='2'
;
;    ;INTERPOLATE ONTO MODEL LEVELS
;      pres_model=reform(read_nc_var(model_pres_levels_nc,'',varid='3'))
;      np_mod=n_elements(pres_model)
;      qv_prof_mod = interpol([qv_prof,min(qv_prof)*0.05],[dims.pres*1e2,1000],pres_model)
;      for iz=0,np_mod-1 do print,'QV: ',pres_model[iz]*1e-2,' hPa, ',qv_prof_mod[iz]*1e3,' g/kg'
;
;    openw,1,qv_file
;      forstr='(f10.8)'
;      printf,1,'/ '+strjoin(string(qv_prof_mod[0:7],format=forstr),',',/single)+', &'
;      for i=0,3 do printf,1,strjoin(string(qv_prof_mod[8+i*8:15+i*8],format=forstr),',',/single)+', &'
;      printf,1,strjoin(string(qv_prof_mod[16+3*8:np_mod-1],format=forstr),',',/single)+'/'
;    close,1
;
;    return
;
;  ;OVERWRITE RH WITH DRY WATER VAPOR PROFILE
;    overwrite_rh=0
;    if overwrite_rh then begin
;      for iz=0,dims.np-1 do $
;        bin_qv[*,iz]=qv_prof[iz]
;      for ibin=0,ibin-1 do $
;        bin1[ibin,*]=calc_relh(bin_qv[ibin,*],bin_tmpk[ibin,*],dims.pres*1e2)
;    endif
;
;  endif


;----ROUTINES FOR HEW'S SEQ CODE--------------------


  if write_force_seq or overlay_seq then begin

    ;--> SEQ CODE OF PENDERGRASS AND WILLOUGHBY (2009, MWR)

    ;HEIGHT FROM HYDROSTATIC
      z = fltarr(nrad,dims.np) ; m
      z[*]=!values.f_nan
      dz=z
      for ir=0,nrad-1 do begin
        ;TAKE CARE OF BELOW-GROUND POINTS
          kdo=where(dims.pres lt slp[ir],nkdo)
          ik=kdo[0]
          z[ir,ik] = (slp[ir]-dims.pres[ik])*1e2 / 9.81 / rho[ir,ik]
        for ik=kdo[1],kdo[nkdo-1] do $
          z[ir,ik] = z[ir,ik-1] + (dims.pres[ik-1]-dims.pres[ik])*1e2 / 9.81 / mean(rho[ir,ik-1:ik],/double,/nan)
        dz[ir,kdo]=deriv(reform(z[ir,kdo]))
      endfor
    ;HEIGHT FROM HYDROSTATIC
;      z = fltarr(dims.np) ; m
;      z[0]=25.
;      for iz=1,dims.np-1 do $
;        z[iz] = z[iz-1] + (dims.pres[iz-1]-dims.pres[iz])*1e2 / 9.81 / mean(rho[*,iz-1:iz],/double,/nan)
;      dz=deriv(z)

    ;SEQ GRID
      nr_seq=501;601;751
      nz_seq=21
      rad_seq=findgen(nr_seq)*2
      z_seq=findgen(nz_seq)*1e3

    ;SHOULD NOT EXTRAPOLATE ONTO SEQ GRID
      if max(rad_seq) gt max(radius) then message,'SHOULD NOT EXTRAPOLATE ONTO SEQ GRID!'

    ;WRITE OUT MEAN PROFILE
      ex=exner(dims.pres)
      mth=fltarr(dims.np)
      mz=fltarr(dims.np)
      irad=where(radius le max(rad_seq))
      totrad=total(radius[irad],/double)
      for iz=0,dims.np-1 do begin
        mth[iz]=total(reform(dummy[irad,iz])*radius[irad],/double)/totrad
        mz[iz]=total(reform(z[irad,iz])*radius[irad],/double)/totrad
      endfor
      mtmpk=theta(mth,dims.pres,/reverse)+273.15
      ddz=deriv(mz,mtmpk)
;      for iz=0,dims.np-1 do print,iz,dims.pres[iz],mtmpk[iz]-273.15,ddz[iz]*1e3
      kstrat=36
      tstrat=mtmpk[kstrat]-273.15
      if write_force_seq then begin
        openw,1,snd_file
          ;SPECS
          printf,1,dims.np,' ',kstrat,' ',tstrat,' T'
          ;SURFACE
          printf,1,1013.0,' ',exner(1013.),' ',mth[0],' ',0.
          ;PROFILE
          for iz=0,dims.np-1 do $
            printf,1,dims.pres[iz],' ',ex[iz],' ',mth[iz],' ',mz[iz]
        close,1
      endif

    ;CONVERT TO BUOYANCY SOURCE (from q to Q)
      cp=1004. & g=9.81
      ;fac=g;/cp ; SMALL q SHOULD BE THETA-DOT * CP, SO CP'S CANCEL OUT
      q = var1 * g / dummy ; K / s --> m / s^3

    ;INTERPOLATE ONTO SEQ GRID
      var_seq=fltarr(nz_seq,nrad)
      var1_seq=fltarr(nz_seq,nr_seq)
      for irad=0,nrad-1 do $
        var_seq[*,irad]=interpol(reform(q[irad,*]),reform(z[irad,*]),z_seq)
      for iz=0,nz_seq-1 do $
        var1_seq[iz,*]=interpol(reform(var_seq[iz,*]),radius,rad_seq)

    ;REPLACE NANS WITH ZERO
    locnan=where(~finite(var1_seq),count)
    if count gt 0 then var1_seq[locnan]=0.

    ;FOR MP_HEAT, RAMP DOWN TO ZERO AT OUTER RADII
    ;FOR RAD VARS, REPLACE EXTRAPOLATED COLUMNS WITH OUTER-RADIUS MEAN
    if var1_str eq 'H_DIABATIC' then begin
      r1=700.
      r2=800.
      for irad=0,nr_seq-1 do begin
        fac=(r2-rad_seq[irad])/(r2-r1)
        fac=min([1,fac])
        fac=max([0,fac])
        var1_seq[*,irad] *= fac
      endfor
    endif; else if strmatch(var1_str,'RTHRA*') then begin
;      r1=800.;600.
;      r2=1000.
;      rloc=where((rad_seq ge r1) and (rad_seq le r2))
;      vmean=fltarr(nz_seq)
;      rtot=total(rad_seq[rloc],/double)
;      for iz=0,nz_seq-1 do vmean[iz]=total(reform(var1_seq[iz,rloc])*rad_seq[rloc],/double)/rtot
;      loc_replace=where(rad_seq ge r1)
;      ;IMPOSE THIS MEAN PROFILE FROM r >= R1
;        for iz=0,nz_seq-1 do var1_seq[iz,loc_replace]=vmean[iz]
;    endif

    ;NEED TO REMOVE AZIMUTHAL MEAN AFTER INTERPOLATION
    if remove_azmn then begin
      radall=fltarr(nr_seq,nz_seq)
      for ip=0,nz_seq-1 do radall[*,ip]=rad_seq
      totrad=total(rad_seq,/double,/nan)
      var_mn=total(var1_seq*radall,1,/double,/nan)/totrad
      for ip=0,nz_seq-1 do var1_seq[*,ip]-=var_mn[ip]
    endif

    ;ZERO-OUT CERTAIN LAYERS
    zero_lay=0
    if zero_lay then begin
      ;PRESSURE ON SEQ GRID
        pr_seq=interpol(dims.pres,z,z_seq,/quadratic)
;      loc_zero=where(pr_seq ge 800.)
      loc_zero=where(pr_seq ge 300.)
;      loc_zero=where((pr_seq ge 800.) or (pr_seq lt 300.))
      var1_seq[loc_zero,*]=0.
    endif

    ;WRITE OUT
      if write_force_seq then begin
        openw,1,seq_file
          printf,1,var1_seq
        close,1
        print,'Done writing out!'
        exit
      endif

  if overlay_seq then begin

    ;PRESSURE ON SEQ GRID
      pr_seq=interpol(dims.pres,mz,z_seq,/quadratic)

    ;VAR SETTINGS
    doread=1
    if seq_cvar eq 'VG_OUT' or seq_cvar eq 'VVV' then begin;var1_str eq 'v_tan' then begin
      cint_seq=1
      fact=1.
    endif else if seq_cvar eq 'UUU' then begin;var1_str eq 'u_rad' then begin
      cint_seq=0.1;1;.1
      fact=1.
    endif else if seq_cvar eq 'WWW' then begin;var1_str eq 'W' then begin
      cint_seq=.1;2
      fact=1e2 ; m/s --> cm/s
    endif else if seq_cvar eq 'DVT' then begin ; Delta(v_tan)
      cint_seq=0.01;1
      fact=1. ; already in m/s / hour from HEW's code
    endif

    if doread then begin

    ;READ
      var_seq=fltarr(nr_seq,nz_seq)
      seq_read_file='/work/06040/tg853394/stampede2/tc_stuff/FCIRC_19_local/DATA/Test_0/'+seq_subdir+'/'+seq_cvar+'.txt'
      line=''
      openr,1,seq_read_file
        readf,1,line
        if seq_cvar eq 'VVV' then readf,1,line
        for il=0,nz_seq-1 do begin
          readf,1,line
          var_seq[*,il]=strsplit(line,/extract)
        endfor
      close,1

      var_seq*=fact

  ;PLOT HEAT FORCING
    if plot_forcing_v1 then begin
    ;READ
      force_seq=fltarr(nr_seq,nz_seq)
      seq_read_file='/work/06040/tg853394/stampede2/tc_stuff/FCIRC_19_local/DATA/Test_0/'+seq_subdir+'/FQQ.txt'
      line=''
      openr,1,seq_read_file
        readf,1,line
        for il=0,nz_seq-1 do begin
          readf,1,line
          force_seq[*,il]=strsplit(line,/extract)
        endfor
      close,1
    ;INTERPOLATE
      ;TH1_SEQ IS A DUMMY FOR INTERPOLATING SEQ HEAT FORCING FROM SEQ GRID ONTO WRF AZIM GRID
      th_seq=fltarr(nrad,nz_seq)
      ;th_seq=fltarr(nr_seq,dims.np)
      th1_seq=fltarr(nrad,dims.np)
      for iz=0,nz_seq-1 do $
        th_seq[*,iz]=interpol(reform(force_seq[*,iz]),rad_seq,radius)
      for irad=0,nrad-1 do $
        th1_seq[irad,*]=interpol(reform(th_seq[irad,*]),z_seq,reform(z[irad,*]))
      var1=th1_seq*dummy/9.81 ; K/s --> K/d
    ;SET TO NANS IF EXTRAPOLATED BEYOND SEQ_RADIUS
      iradnan=where(radius gt max(rad_seq),nnan)
      if nnan gt 0 then var1[iradnan,*]=!values.f_nan
    endif

  ;CALCULATE TIME SCALE TO SATURATION
    if ilcl then begin
      ;VAR_SEQ IS VERTICAL MOTION W_CRF
      ;DUMMY2 IS SATURATION DEFICIT (QS - QV)
      ;DUMMY3 IS QV
      seq1=fltarr(nrad,nz_seq)
      seq2=fltarr(nrad,dims.np)
      for iz=0,nz_seq-1 do seq1[*,iz]=interpol(reform(var_seq[*,iz]),rad_seq,radius)
      for ir=0,nrad-1 do   seq2[ir,*]=interpol(reform(seq1[ir,*]),z_seq,reform(z[ir,*])) ; VERTICAL MOTION (M/S)
      ;SET TO NANS IF EXTRAPOLATED BEYOND SEQ_RADIUS
        iradnan=where(radius gt max(rad_seq),nnan)
        if nnan gt 0 then seq2[iradnan,*]=!values.f_nan
      ;SMOOTH W_CRF
        seq2=smooth(temporary(seq2),[3,3],/edge_truncate,/nan)
      ;dq / dz
      for ir=0,nrad-1 do var2[ir,*]=deriv(reform(z[ir,*]),reform(dummy3[ir,*])) ; kg/kg / m
    ;CALCULATE TIME TO SATURATION VIA LIFTING AS
      ; XX H_LCL / W
      ; tau = -1. * Delta_qv / madv, where MADV = W*dqdz
      seq2*=1./fact ; w_crf from cm/s --> m/s
      seq2[where(seq2 le 1e-4)]=!values.f_nan
;      seq2*=24.*3600 ; -> m/day
;      seqsat = dummy2 / seq2 ; --> units of days
      madv = seq2*var2 ; kg/kg / s
      seqsat = -1. * dummy2 / madv ; seconds
      seq2=0 & dummy2=0 & dummy3=0
    endif

    endif

  endif
  endif


;----WRITE LW CRF-PROFILE--------------------


  ;CALCULATES THE RADIUS-NORMALIZED RADIAL AVERAGE OF LW-CRF TO BE IMPOSED AS FORCING IN WRF
  ;INTERPOLATES ONTO WRF VERTICAL GRID

  if write_crf and (var1_str eq 'RTHRATLW') and (dirs.cases[ic] eq 'ctl') then begin

    if ~icrf then message,'Turn on CRF to write profile out!'
    if remove_azmn then message,'Should not be removing azimuthal mean for CRF forcing profile!'

    print,'WRITING CRF-PROFILE TO ',crf_file

    ;AVERAGE FROM AXISYMMETRIC OUTPUT
      rad_sel=where((radius ge rad_thresh[0]) and (radius le rad_thresh[1]),nrad_sel)
      rada=fltarr(nrad_sel,dims.np)
      for iz=0,dims.np-1 do rada[*,iz]=radius[rad_sel]
      radtot=total(radius[rad_sel],/nan,/double)
      crf_force = total( (var1[rad_sel,*]*rada) , 1, /double,/nan) / radtot
;      for iz=0,dims.np-1 do print,'LW-CRF: ',dims.pres[iz],' hPa, ',crf_force[iz]*3600*24,' K/d'

    ;INTERPOLATE ONTO MODEL LEVELS
      pres_model=reform(read_nc_var(model_pres_levels_nc,'P'))
      lon_sel=where(dims.lon ge -55 and dims.lon le -45)
      lat_sel=where(dims.lat ge 7 and dims.lat le 17)
      pres_model=mean(mean(temporary(pres_model[lon_sel,lat_sel,*]),dimension=1,/nan,/double),dimension=1,/nan,/double)
      np_mod=n_elements(pres_model)
      p_top=1000
      crf_force_model = interpol([crf_force,0],[dims.pres*1e2,p_top],pres_model)

    ;IMPOSE POSITIVE ONLY
;      ineg=where(crf_force_model lt 0 and pres_model le 400*1e2)
;      crf_force_model[ineg[0]:np_mod-1]=0.d

    ;VERTICALLY SMOOTH

    ;PRINT INTERPOLATED
      for iz=0,np_mod-1 do print,'LW-CRF: ',pres_model[iz]*1e-2,' hPa, ',crf_force_model[iz]*3600*24,' K/d'
;      for iz=0,np_mod-1 do print,'LW-CRF: ',pres_model[iz]*1e-2,' hPa, ',crf_force_model[iz]

    openw,1,crf_file
      forstr='(f10.8)'
;      printf,1,'/ '+strjoin(string(crf_force_model[0:7],format=forstr),',',/single)+', &'
;      for i=0,3 do printf,1,strjoin(string(crf_force_model[8+i*8:15+i*8],format=forstr),',',/single)+', &'
;      printf,1,strjoin(string(crf_force_model[16+3*8:np_mod-1],format=forstr),',',/single)+'/'
      extra = np_mod mod 8
      n8=(np_mod-extra)/8
      printf,1,'/ '+strjoin(strtrim(crf_force_model[0:7],2),',',/single)+', &'
      for i=1,n8-1 do printf,1,strjoin(strtrim(crf_force_model[i*8:7+i*8],2),',',/single)+', &'
      printf,1,strjoin(strtrim(crf_force_model[n8*8:np_mod-1],2),',',/single)+'/'
    close,1

    continue
;    return

  endif

  ;ADD FORCING TO RADHEAT
  if write_crf and add_icrf and (dirs.cases[ic] ne 'ctl') then begin
    fact=1.
    rad_lin=[200,300] ; USED 200-300 FOR ICRF_*
;    rad_lin=[250,350]
    for irad=0,nrad-1 do begin
      fact = (rad_lin[1] - radius[irad]) / (rad_lin[1]-rad_lin[0])
      fact=max([0.,fact])
      fact=min([1.,fact])
      for iz=0,dims.np-1 do $
        var1[irad,iz] += crf_force[iz] * fact
    endfor
  endif


;----CREATE PLOTS--------------------


;  if strmatch(var1_str,'*RTHRA*') then begin
;    setmax=3;4
;    setmin=-1.*setmax;4
;  endif

;  stats,var1

  tc_figspecs, var1_str, figspecs, setmax=setmax, setmin=setmin

  figdir=dirs.figdir+'/azim_cross/'
  spawn,'mkdir '+figdir,tmp,tmpe
  figname=figdir+dirs.cases[ic]+'_'+var1_str

  if strmatch(var1_str,'*RTHRA*') then begin
    if var1_str eq 'RTHRATTOT' then radtag='LW+SW' else radtag=strmid(var1_str,1,2,/reverse)
    if icrf then begin
      radtag='LW-CRF'
;radtag='SW-CRF'
radtag='CRF'
;radtag='CLEAR'
      figname+='_crf'
    endif
    figspecs.cbar_tag=radtag+' '+figspecs.cbar_tag
    irad=1
  endif else irad=0
;irad=0

;  if icrf then olr-=olrc

  if remove_azmn then figname+='_azprm'

  figname+='_'+hr_tag_plot
  if nazim gt 1 then figname+='_az'+strtrim(iazim,2)
  figspecs=create_struct(figspecs,'figname',figname)

;OVERLAY VAR

;  if ihght then inz=nz else 
inz=dims.np

  if cont_var eq 0 then begin
  ;PSI
    cvar=var2
    cint=3
    clevs=indgen(200)*cint+cint
  endif else if cont_var eq 1 then begin
  ;CLOUD
;    for i=0,1 do var2=smooth(var2,[3,0],/edge_truncate,/nan)
    cvar=var2*1e6 ; kg/kg --> mg/kg
;     cvar=cld*1e6 ; kg/kg --> mg/kg
;    clevs=[1,5,10,25,50,200,300]
     clevs=[5,10,25,50]
;     clevs=[5,10,25,50,200,500]
     clevs=[10,20,50,100]
  endif else if cont_var eq 2 then begin
  ;URAD
    cvar=var2 ; m/s
    cint=1;0.5;1;3
    clevs=indgen(50)*cint+cint
  endif else if cont_var eq 3 then begin
  ;W
    cvar=var2*1e2 ; cm/s
    cint=2
    clevs=indgen(50)*cint+cint
  endif else if cont_var eq 4 then begin
  ;M
    cvar=var2*1e-5 ; 10^5 m2/s
    cint=5
    clevs=indgen(50)*cint+cint
  endif else if cont_var eq 5 then begin
  ;LH
    cvar=var2*3600d*24 ; K/d
for i=0,0 do cvar=gauss_smooth(temporary(cvar),[3,0],/edge_truncate)
    cint=4
;    clevs=indgen(50)*cint+cint
    clevs=[1,2,5,10,50,100,500,1000,2000,3000]
  endif else if cont_var eq 6 then begin
  ;AVOR
    cvar=var2 ; 10^-5 /s
for i=0,0 do cvar=gauss_smooth(temporary(cvar),[3,0],/edge_truncate)
    clevs=[1,2,5,10,50,100,500,1000,2000,3000]
  endif else if cont_var eq 7 then begin
  ;VTAN
;    for i=0,1 do $
;      var2=smooth(temporary(var2),[3,3],/edge_truncate,/nan)
    cvar=var2 ; m/s
    cint=3
    clevs=indgen(50)*cint+cint
  endif else if cont_var eq 8 then begin
  ;INERTIAL STABILITY
    fcor=2.*7.292e-5*sin(mean(tc_lat)*!pi/180)
    for iz=0,inz-1 do $
      var2[*,iz] = (fcor + 2*var2[*,iz]/(radius*1e3)) * (fcor + deriv(radius*1e3,radius*1e3*var2[*,iz])/(radius*1e3)) ; /s
    cvar = var2*1e7 ; 10^-7/s
    clevs=10^(findgen(20)-5)
  endif

;TIME TO SATURATION VIA SEQ UPWARD MOTION
  if overlay_seq and ilcl then begin
    var_seq=0 ; so that it doesn't contour this instead
    cvar=seqsat/3600 ; s --> hours
    cint=30
    ;clevs=indgen(50)*cint+cint
    clevs=[1,2,4,indgen(4)*6+6,indgen(4)*12+36,96]
  endif

;1D VAR

  ;-1=no plot, 0=OLR, 1=PW, 2=v_tan
  if var_1d eq 0 then begin
    ytitle1d='OLR!C[ W m!U-2!N ]'
;    yrange1d=[60,310]
    yrange1d=[100,310]
  endif else if var_1d eq 1 then begin
    ytitle1d='PW!C[ mm ]'
    yrange1d=[40,70];[25,80]
  endif

  if var_1d ne -1 then $
    figspecs=create_struct(figspecs,'ytitle1d',ytitle1d,'yrange1d',yrange1d)

;FOR IRUN TO ADD TWO DIFFERENT HEAT SOURCES
;  if irun eq 0 then var1sav=var1
;  if irun eq 1 then var1+=var1sav

  ;SMOOTH VAR1
  if do_smooth then begin
;    for i=0,1 do $
      var1=smooth(temporary(var1),[3,3],/edge_truncate,/nan)
      cvar=smooth(temporary(cvar),[3,3],/edge_truncate,/nan)
    if keyword_set(var_seq) then $
      var_seq=smooth(temporary(var_seq),[3,3],/edge_truncate,/nan)
  endif

  ;CONVERT VERTICAL AXIS TO HEIGHT
  if ihght then begin
    dz=0.25
    nz=ztop/dz;+1
    hght=findgen(nz)*dz+dz
    y=hght
    var1_z=fltarr(nrad,nz) & var1_z[*]=!values.f_nan
    cvar_z=var1_z
    for ir=0,nrad-1 do begin
      kdo=where(finite(reform(z[ir,*])))
      p_out=interpol(dims.pres[kdo],reform(z[ir,kdo])*1e-3,hght,/quadratic)
      var1_z[ir,*]=interpol(reform(var1[ir,kdo]),dims.pres[kdo],p_out,/quadratic)
      cvar_z[ir,*]=interpol(reform(cvar[ir,kdo]),dims.pres[kdo],p_out,/quadratic)
    endfor
    var1=var1_z
    cvar=cvar_z
    if keyword_set(var_seq) then begin
      var_seq_z=fltarr(nr_seq,nz) & var_seq_z[*]=!values.f_nan
      for ir=0,nr_seq-1 do begin
        p_out=interpol(pr_seq,z_seq*1e-3,hght,/quadratic)
        var_seq_z[ir,*]=interpol(reform(var_seq[ir,*]),pr_seq,p_out,/quadratic)
      endfor
      var_seq=var_seq_z
      pr_seq=hght
    endif
  endif else y=dims.pres

  ;CALCULATE RMW
    if var1_str eq 'v_tan' or cont_var eq 7 then begin
      if ihght then inz=nz else inz=dims.np
      rmw=fltarr(inz)
      if var1_str eq 'v_tan' then vtn=var1 else vtn=cvar
      for iz=0,inz-1 do begin
        max=max(reform(vtn[*,iz]),loc,/nan)
        rmw[iz]=radius[xrad[loc]]
      endfor
      rmw=smooth(temporary(rmw),3,/nan,/edge_truncate)
    endif

  wrf_tc_azim_cross, var1, radius, y, figspecs, irad=irad, cvar=cvar, clevs=clevs, $
    var1d=var1d, var2_1d=var1d2, $;, cv2=cvar2, cl2=clevs2
    pr_seq=pr_seq, rad_seq=rad_seq, var_seq=var_seq, cint_seq=cint_seq, rmw=rmw

endfor ; iazim

endfor ; icase
;endfor ; irun

endfor ; ipanel ; for Paper Fig 4

endfor ; iex (extra hours)

print,'DONE!!'
end
