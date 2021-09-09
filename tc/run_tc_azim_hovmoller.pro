; 
; Plot Hovmoller from azimuthally averaged TC output (generated by run_calc_azim.pro).
;
; James Ruppert
; 1/29/20
; 
pro run_tc_azim_hovmoller

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

;TIME SPECS
  nt_full=dims.nt-1
  if strmatch(subdir,'redux*') and tcname eq 'maria' then nt_full=145
  if strmatch(subdir,'redux*') and tcname eq 'haiyan' then nt_full=169
  npd=dims.npd
  nhrs=1.*nt_full*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nhrs)

;AZIM FILES
  hr_sel=[24,120]
  if tcname eq 'haiyan' then hr_sel=[24,144]
  if strmatch(subdir,'redux*') then begin
    if tcname eq 'maria'  then hr_sel=[0,144]
    if tcname eq 'haiyan' then hr_sel=[0,168]
  endif
  t_ind=where((time_hrs ge hr_sel[0]) and (time_hrs le hr_sel[1]))
  t_ind=t_ind[where(t_ind le nt_full-1)]
  nt_sav=n_elements(t_ind)

;----PLOT OPTIONS--------------------

ismooth=1 ; smooth time series?

radius_range=[0,800];600];800];200] ; km, max=1300 km for redux

var1_str='RTHRATLW'
;var1_str='RTHRATTOT'
;var1_str='RTHRATSW'
;var1_str='W'
;var1_str='T'
;var1_str='AVOR'
;var1_str='H_DIABATIC'
;var1_str='rainrate'
var1_str='CWV'
;var1_str='RH';'QVAPOR';
;var1_str='SLP'
;var1_str='shear'
;var1_str='SEF'
;var1_str='wind'
;  wvar='v_tan'
;  wvar='u_rad'
;  wvar='wspd'

remove_azmn=0
remove_azmn=1

cont_var=0 ; Contour var: -1: none, 0: v_tan

icrf=0

idiff=0 ; Difference from CTL?

;----READ VARS--------------------


;AZIMUTHAL SPECS
  nrad=267
  naz=360
  if tcname eq 'haiyan' or strmatch(subdir,'redux*') then nrad=433

print,'VAR: ',var1_str

;for ic=0,dirs.nc-1 do begin
;for ic=0,4,4 do begin
for ic=0,0 do begin
;for ic=1,dirs.nc-1 do begin

  print,'CASE: ',strupcase(dirs.cases[ic])

  hr0=0
  if strmatch(dirs.cases[ic],'*36h*') then hr0=36
  if strmatch(dirs.cases[ic],'*24h*') then hr0=24
  if strmatch(dirs.cases[ic],'*48h*') then hr0=48
  if strmatch(dirs.cases[ic],'*60h*') then hr0=60
  if strmatch(dirs.cases[ic],'*72h*') then hr0=72
  if strmatch(dirs.cases[ic],'*84h*') then hr0=84
  if strmatch(dirs.cases[ic],'*96h*') then hr0=96
  if strmatch(dirs.cases[ic],'lwcrf2*') then hr0=36
  if strmatch(dirs.cases[ic],'lwcrf*') then hr0=36
  nt_test_full=nt_full - hr0

if strmatch(subdir,'redux*') then begin

  if total(strmatch(['lwcrf','axisym'],dirs.cases[ic])) then begin
;    hr_fil='49-72hr'
;    nt_test=24
    hr0='36'
    nt_test=hr_sel[1]-fix(hr0)+1
    hr_fil=hr0+'-'+strtrim(hr_sel[1],2)+'hr'
  endif else if strmatch(dirs.cases[ic],'*ncrf_*') then begin
    hr0=strmid(dirs.cases[ic],2,2,/reverse_offset)
    nt_test=hr_sel[1]-fix(hr0)+1
    hr_fil=hr0+'-'+strtrim(hr_sel[1],2)+'hr'
  endif else begin
  ;CTL
    ;REDUX
  if tcname eq 'maria' then begin
    hr_fil='0-144hr'
;    hr_tag_plot='049-096hr'
    nt_test=145
;    t_ind_plot=indgen(nt_plot)+(fix(strmid(hr_tag_plot,0,3)) - fix(strmid(hr_fil,0,2)))
  endif else begin
    hr_fil='0-168hr'
    nt_test=169
  endelse
  endelse

  t_ind_final=indgen(nt_test)+hr0

endif else begin

  ;TIME SELECTION
    t_offset=max([0,hr0-hr_sel[0]])
    ut_offset=max([0,hr_sel[0]-hr0])
    nt_test=nt_sav-t_offset
    t_ind_test=indgen(nt_test)+ut_offset
    hrs_test=time_hrs[t_ind_test+hr0]
    hr_fil=strtrim(hrs_test[0],2)+'-'+strtrim(hrs_test[nt_test-1],2)+'hr'
    t_ind_final=indgen(nt_test)+hr0;-24*(hr0 ge 24)

endelse

  it0 = 0;t_ind_test[0]

  ;READ STORM MOTION (CALCULATED IN RUN_WRF_TC_TRACKS)
    motion_file=dirs.casedir[ic]+'storm_motion.txt'
    openr,1,motion_file
      readf,1,tmp_nt
      tmp_it_test=intarr(tmp_nt)
      readf,1,tmp_it_test
      motion_x=fltarr(tmp_nt) & motion_y=motion_x
      readf,1,motion_x
      readf,1,motion_y
    close,1

  ;2D VARS

    count=[nrad,naz,1,nt_test] & offset=[0,0,0,it0] ; x,y,z,t

  ;SLP
    file=dirs.casedir[ic]+'azim_SLP_'+hr_fil+'.nc'
    slp=reform(read_nc_var(file,'SLP',count=count,offset=offset))
  ;OLR
;    file=dirs.casedir[ic]+'azim_OLR_'+hr_fil+'.nc'
;    olr=reform(read_nc_var(file,'OLR',count=count,offset=offset))
  ;OLRC
;    file=dirs.casedir[ic]+'azim_OLRC_'+hr_fil+'.nc'
;    olrc=reform(read_nc_var(file,'OLRC',count=count,offset=offset))

  radius=read_nc_var(file,'radius')
  azimuth=read_nc_var(file,'azmiuth')

  ;RMW
    file=dirs.casedir[ic]+'azim_U10_'+hr_fil+'.nc'
    u=reform(read_nc_var(file,'U10',count=[nrad,naz,1,nt_test],offset=[0,0,0,it0]))
    file=dirs.casedir[ic]+'azim_V10_'+hr_fil+'.nc'
    v=reform(read_nc_var(file,'V10',count=[nrad,naz,1,nt_test],offset=[0,0,0,it0]))
;       ;WINDS ALOFT
;      psel=850;950;850;925 ; hPa
;      iz=where(dims.pres eq psel)
;      file=dirs.casedir[ic]+'azim_U_'+hr_fil+'.nc'
;      u=reform(read_nc_var(file,'U',count=[nrad,naz,1,nt_test],offset=[0,0,iz,it0]))
;      file=dirs.casedir[ic]+'azim_V_'+hr_fil+'.nc'
;      v=reform(read_nc_var(file,'V',count=[nrad,naz,1,nt_test],offset=[0,0,iz,it0]))
    ;SUBTRACT STORM MOTION
    for it=0,nt_test-1 do begin
      u[*,*,it] -= motion_x[it+it0];+hr_plot0[0]]
      v[*,*,it] -= motion_y[it+it0];+hr_plot0[0]]
    endfor
    wnd_azim=azim_wind_conv(u,v,azimuth) & u=0 & v=0
    vtan=wnd_azim.v_tan
    vtan=mean(temporary(vtan),dimension=2,/nan,/double)
    rmw=fltarr(nt_test)
    for it=0,nt_test-1 do begin
      max=max(reform(vtan[*,it]),loc,/nan)
      rmw[it]=radius[loc]
    endfor
    if ismooth then begin
      ntsmooth=3 ; = n * 1 hr
      rmw=gauss_smooth(temporary(rmw),ntsmooth,/edge_truncate)
    endif
    if cont_var eq 0 then cvar=vtan

  ;3D VARS

    count=[nrad,naz,dims.np,nt_test] & offset=[0,0,0,it0] ; x,y,z,t

  ;VAR1

    if var1_str eq 'RH' or var1_str eq 'CWV' then begin

      ;NEED LARGE PRES ARRAY FOR RH
        prx=fltarr(nrad,naz,dims.np,nt_test)
        for iz=0,dims.np-1 do prx[*,*,iz,*]=dims.pres[iz]*1e2

      file=dirs.casedir[ic]+'azim_QVAPOR_'+hr_fil+'.nc'
      qv=read_nc_var(file,'QVAPOR',count=count,offset=offset)
qv[where(qv lt 0)]=0.
;qv[where(abs(qv) gt 50e-3)]=!values.f_nan
      file=dirs.casedir[ic]+'azim_T_'+hr_fil+'.nc'
      tmpk=read_nc_var(file,'T',count=count,offset=offset)

      psel=500 ; hPa
      iz=where(dims.pres eq psel)

      if var1_str eq 'RH' then begin

        rh=calc_relh(qv,tmpk,prx,/noice) ; problem with ice, probably due to huge array
        var1=reform(rh[*,*,iz,*])

      endif else begin

;var1=calc_relh(qv,tmpk,prx)
;print,reform(var1[100,50,*,50])

        ;SAVE QV AND QV*

        var1=qv ; kg/kg

        iice=where(tmpk lt 273.15,complement=iliq)
        var11=tmpk & var11[*]=0.
        Mw=18.0160 ; molec mass of water
        Md=28.9660 ; molec mass of dry air
        X = ESAT(1d*tmpk[iliq])  ; Returns e in hPa; ESAT accepts temperatures in Celsius or K
        var11[iliq] = 1d*Mw/Md* X / (prx[iliq]*1e-2 - X)
        ;X = EICE(1d*tmpk[iice])  ; Returns e in hPa; ESAT accepts temperatures in Celsius or K
            ; From EICE function of Dominik Brunner
                A=-2663.5D
                B=12.537
                logp=A/(1d*tmpk[iice])+B
                X=10^logp/100. ; conversion to hPa
        var11[iice] = Mw/Md* X / (prx[iice]*1e-2 - X)
        ;THE ABOVE YIELDS SATURATION MIXR IN KG/KG
        X=0
;print,reform((1e2*var1/var11)[100,50,*,50])
;CHECKS OUT THAT Q/QS YIELDS THE PROPER REL HUM

      endelse

      qv=0
      tmpk=0

    endif else if var1_str eq 'RTHRATTOT' then begin

      file=dirs.casedir[ic]+'azim_RTHRATLW_'+hr_fil+'.nc'
      var1=read_nc_var(file,'RTHRATLW',count=count,offset=offset)
      file=dirs.casedir[ic]+'azim_RTHRATSW_'+hr_fil+'.nc'
      var1+=read_nc_var(file,'RTHRATSW',count=count,offset=offset)

    endif else if var1_str eq 'rainrate' then begin

      file=dirs.casedir[ic]+'azim_'+var1_str+'_'+hr_fil+'.nc'
      var1=reform(read_nc_var(file,var1_str,count=[nrad,naz,1,nt_test],offset=offset))
      var1*=1./24 ; mm/d --> mm/h

    endif else if var1_str eq 'shear' then begin

      iz850=(where(dims.pres eq '850'))[0]
      iz200=(where(dims.pres eq '200'))[0]
      file=dirs.casedir[ic]+'azim_U_'+hr_fil+'.nc'
      u850=reform(read_nc_var(file,'U',count=[nrad,naz,1,nt_test],offset=[0,0,iz850,it0]))
      u200=reform(read_nc_var(file,'U',count=[nrad,naz,1,nt_test],offset=[0,0,iz200,it0]))
      file=dirs.casedir[ic]+'azim_V_'+hr_fil+'.nc'
      v850=reform(read_nc_var(file,'V',count=[nrad,naz,1,nt_test],offset=[0,0,iz850,it0]))
      v200=reform(read_nc_var(file,'V',count=[nrad,naz,1,nt_test],offset=[0,0,iz200,it0]))
stats,u850
stats,v850
v850[where(abs(v850) ge 5e2)]=!values.f_nan
exit
      ;CALCULATE SHEAR AS MAGNITUDE OF THE VECTOR DIFFERENCE, AS IN AIYYER AND THORNCROFT (2006)
;      udiff=u200-u850
;      vdiff=v200-v850
      var1=u200;sqrt(udiff^2 + vdiff^2)

    endif else if var1_str eq 'wind' or var1_str eq 'rmw' then begin

      file=dirs.casedir[ic]+'azim_U10_'+hr_fil+'.nc'
      u=reform(read_nc_var(file,'U10',count=[nrad,naz,1,nt_test],offset=[0,0,0,it0]))
      file=dirs.casedir[ic]+'azim_V10_'+hr_fil+'.nc'
      v=reform(read_nc_var(file,'V10',count=[nrad,naz,1,nt_test],offset=[0,0,0,it0]))

;      psel=850;950;850;925 ; hPa
;      iz=where(dims.pres eq psel)
;
;      file=dirs.casedir[ic]+'azim_U_'+hr_fil+'.nc'
;      u=reform(read_nc_var(file,'U',count=[nrad,naz,1,nt_test],offset=[0,0,iz,it0]))
;      file=dirs.casedir[ic]+'azim_V_'+hr_fil+'.nc'
;      v=reform(read_nc_var(file,'V',count=[nrad,naz,1,nt_test],offset=[0,0,iz,it0]))

      ;SUBTRACT STORM MOTION
      for it=0,nt_test-1 do begin
        u[*,*,it] -= motion_x[it+it0];+hr_plot0[0]]
        v[*,*,it] -= motion_y[it+it0];+hr_plot0[0]]
      endfor

      wnd_azim=azim_wind_conv(u,v,azimuth)

      if wvar eq 'v_tan' then $
        var1=wnd_azim.v_tan $
      else if wvar eq 'u_rad' then $
        var1=wnd_azim.u_rad $
      else if wvar eq 'wspd' then $
        var1=sqrt(u^2 + v^2)

      u=0 & v=0

    endif else if var1_str eq 'SLP' then begin

      file=dirs.casedir[ic]+'azim_SLP_'+hr_fil+'.nc'
      var1=reform(read_nc_var(file,'SLP',count=[nrad,naz,1,nt_test],offset=[0,0,0,it0]))

    endif else if var1_str eq 'T' then begin

      psel=700;950;850;925 ; hPa
      iz=where(dims.pres eq psel)

      file=dirs.casedir[ic]+'azim_'+var1_str+'_'+hr_fil+'.nc'
      var1=reform(read_nc_var(file,var1_str,count=[nrad,naz,1,nt_test],offset=[0,0,iz,it0]))

    endif else if var1_str eq 'SEF' then begin

      file=dirs.casedir[ic]+'azim_HFX_'+hr_fil+'.nc'
      var1=reform(read_nc_var(file,'HFX',count=[nrad,naz,1,nt_test],offset=[0,0,0,it0])) ; W/m2
      file=dirs.casedir[ic]+'azim_LH_'+hr_fil+'.nc'
      var1+=reform(read_nc_var(file,'LH',count=[nrad,naz,1,nt_test],offset=[0,0,0,it0])) ; W/m2

    endif else begin

      file=dirs.casedir[ic]+'azim_'+var1_str+'_'+hr_fil+'.nc'
      var1=read_nc_var(file,var1_str,count=count,offset=offset)

      ;DON'T MULTIPLY W BY RHO BECAUSE VERTICAL P-INTEGRAL IS ALREADY MASS-WEIGHTED
;      if var1_str eq 'W' then begin
;        file=dirs.casedir[ic]+'azim_T_'+hr_fil+'.nc'
;        tmpk=read_nc_var(file,'T',count=count,offset=offset)
;        file=dirs.casedir[ic]+'azim_QVAPOR_'+hr_fil+'.nc'
;        qvt=read_nc_var(file,'QVAPOR',count=count,offset=offset)
;        tvirt = tmpk*(1.+0.61*qvt)
;        tmpk=0 & qvt=0
;        rho=var1
;        for iz=0,dims.np-1 do rho[*,*,iz,*] = dims.pres[iz]*1e2 / ( 287. * tvirt[*,*,iz,*] )
;        var1*=rho & rho=0 ; W --> rhoW
;      endif

    endelse

    ;REMOVE CLEAR SKY
    if strmatch(var1_str,'RTHRA*') and icrf then begin

      print,'REMOVING CLEAR SKY!!'

      if var1_str eq 'RTHRATTOT' then begin
        file=dirs.casedir[ic]+'azim_RTHRATLWC_'+hr_fil+'.nc'
        var1-=read_nc_var(file,'RTHRATLWC',count=count,offset=offset)
        file=dirs.casedir[ic]+'azim_RTHRATSWC_'+hr_fil+'.nc'
        var1-=read_nc_var(file,'RTHRATSWC',count=count,offset=offset)
      endif else begin
        file=dirs.casedir[ic]+'azim_'+var1_str+'C_'+hr_fil+'.nc'
        var1-=read_nc_var(file,var1_str+'C',count=count,offset=offset)
      endelse

    endif

  ;AZIMUTHALLY AVERAGE
    slp=mean(temporary(slp),dimension=2,/nan,/double)
;    olr=mean(temporary(olr),dimension=2,/nan,/double)
;    olrc=mean(temporary(olrc),dimension=2,/nan,/double)
    var1=mean(temporary(var1),dimension=2,/nan,/double)
    if var1_str eq 'CWV' then var11=mean(temporary(var11),dimension=2,/nan,/double)
;    var2=mean(temporary(var2),dimension=2,/nan,/double)

ndims=size(var1,/N_DIMENSIONS)

  ;REMOVE AZIMUTHAL MEAN
  if remove_azmn then begin
    print,'REMOVING AZIM MEAN!!'
    totrad=total(radius,/double,/nan)
    if ndims eq 3 then begin
      radall=fltarr(nrad,dims.np,nt_test)
      for ir=0,nrad-1 do radall[ir,*,*]=radius[ir]
      var_mn=total(var1*radall,1,/double,/nan)/totrad
      for ir=0,nrad-1 do var1[ir,*,*]-=var_mn
    endif else begin
      ;totrad=total(radius,/double,/nan)
      for it=0,nt_test-1 do begin
        imean=total(reform(var1[*,it])*radius,/nan,/double)/totrad
        var1[*,it]-=imean
      endfor
    endelse
  endif

  ;REMOVE BELOW-GROUND POINTS
  if ndims eq 3 then begin
    for ip=0,dims.np-1 do begin
      nan=where(slp lt dims.pres[ip],count)
      iv=reform(var1[*,ip,*])
      if count gt 0 then iv[nan]=!values.f_nan
      var1[*,ip,*]=iv
    endfor
  endif


;----PROCESS VAR--------------------


;if var1_str eq 'H_DIABATIC' then begin
if ndims eq 3 then begin

  ;DP FOR VERTICAL INTEGRATION
    dp=fltarr(nrad,dims.np,nt_test)
    idp=deriv(dims.pres)
    idp*=-1e2
    for iz=0,dims.np-1 do dp[*,iz,*]=idp[iz]

  if var1_str eq 'H_DIABATIC' or strmatch(var1_str,'RTHRAT*') then var1*=1004. ; K units --> W/m2 (once integrated)
  if var1_str eq 'H_DIABATIC' then var1*=3600d/(2.5e6) ; W/m2 --> mm/s --> mm/hr

  ;VERTICALLY INTEGRATE
    ip = where(dims.pres ge 50)
;ip = where(dims.pres ge 100)
    if var1_str eq 'CWV' then ip = where((dims.pres le 850) and (dims.pres ge 100))
    sum1 = (1./9.81) * total( var1[*,ip,*] * dp[*,ip,*]  ,2,/nan,/double)
    var1=sum1 & sum1=0

  if var1_str eq 'CWV' then begin
    sum1 = (1./9.81) * total( var11[*,ip,*] * dp[*,ip,*] ,2,/nan,/double)
    var11=sum1 & sum1=0
    var1=temporary(var1)/var11
  endif

endif

;FOR RADIATION, SET 1ST TIME STEP TO NAN, SINCE THIS IS CTL
  if strmatch(var1_str,'RTHRAT*') then var1[*,0]=!values.f_nan

if ismooth then begin
  nrsmooth=1 ; = n * 3km
  ntsmooth=3;1 ; = n * 1 hr
  var1=gauss_smooth(var1,[nrsmooth,ntsmooth],/nan,/edge_truncate)
  if keyword_set(cvar) then cvar=gauss_smooth(cvar,[nrsmooth,ntsmooth],/nan,/edge_truncate)
endif

;HAIYAN NCRF-36H GOES OFF MAP - NEED TO CUT IT OFF
if dirs.cases[ic] eq 'hncrf_36h' then begin
  nnan=11
  var1[*,nt_test-1-nnan:nt_test-1]=!values.f_nan
  rmw[nt_test-1-nnan:nt_test-1]=!values.f_nan
endif

;----CREATE PLOT--------------------

var1_str_tmp=var1_str

  if var1_str eq 'H_DIABATIC' then begin
    figtag=var1_str
  endif else if var1_str eq 'W' then begin
    setmax=15
    setmin=-1.*setmax;'0'
    figtag='wmean'
  endif else if var1_str eq 'rainrate' then begin
    setmax=20;40
    setmin=0.1
    figtag='rain'
  endif else if var1_str eq 'T' then begin
    setmax=10;5
    setmin='0';-1.*setmax
    figtag='T'
  endif else if var1_str eq 'SEF' then begin
    setmax=500
    setmin=-100
    figtag='sef'
  endif else if strmatch(var1_str,'RTHRAT*') then begin
    setmax=100
    setmin=-1.*setmax
    figtag=strlowcase(var1_str)
  endif else if var1_str eq 'RH' then begin
    setmax=100
    setmin=30
    figtag='rh'
  endif else if var1_str eq 'CWV' then begin
    setmax=100
    setmin=60;50
    if idiff and dirs.cases[ic] ne 'ctl' then begin
      setmax=10
      setmin=-1.*setmax
    endif
    figtag='crh'
  endif else if var1_str eq 'shear' then begin
    figtag='shear'
  endif else if var1_str eq 'wind' then begin
    var1_str_tmp=wvar
    figtag=wvar;'wind'
    if wvar eq 'v_tan' then begin
      setmax=50
      setmin='0'
      ;ZERO OUT NEGATIVE VALUES
      var1[where(var1 lt 0)]=0.
    endif else if wvar eq 'u_rad' then begin
      setmax='0'
      setmin=-5;10;20
    endif else if wvar eq 'wspd' then begin
      setmax=20
      setmin='0'
      var1_str_tmp='v_tan' ; Just to use the same plot/figure settings
    endif
    if idiff and dirs.cases[ic] ne 'ctl' then begin
      setmax=2.5
      setmin=-1.*setmax
    endif
  endif else if var1_str eq 'rmw' then begin
    figtag='rmw'
  endif else if var1_str eq 'SLP' then begin
    setmax=1015
    setmin=995
    figtag='slp'
  endif

;CONTOUR VARIABLE
  if cont_var eq 0 then begin
  ;VTAN
    cint=5;10
    clevs=indgen(30)*cint+cint
  endif

;DIFFERENCE FROM CTL
  if idiff then begin
    if dirs.cases[0] eq 'ctl' then begin
      if ic eq 0 then var_ctl=var1 else $
        var1=var_ctl[*,t_ind_final]-var1
    endif
  endif

  tc_figspecs, var1_str_tmp, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figdir=dirs.figdir
  figname=figdir+'azim_hov/tser_hov_'+figtag+'_'+strlowcase(dirs.cases[ic])

  if strmatch(var1_str,'*RTHRA*') then begin
    if icrf then begin
      figname+='_crf'
    endif
  endif

;  if icrf then olr-=olrc

  if remove_azmn then figname+='_azprm'

  ;PLOT SPECS
    csize=0.8
    position=[0.21,0.15,0.82,0.97]
    xsize=3.0 & ysize=3.2
;    xtitle='Time [ hr ]'
    ytitle='Day'; in Sept [ UTC ]'
    ;if tcname eq 'haiyan' then ytitle='Date in Nov [ UTC ]'
    xtitle='Radius [ km ]'

  ;AXES
    y=findgen(nt_sav)/24
    x=radius
    xrange=radius_range
    if ~keyword_set(xrange) then $
      xrange=[min(x),max(x)]
    if ~keyword_set(yrange) then $
      yrange=[min(y),max(y)]

  if tcname eq 'maria' then begin
    y+=14.5
    yrange=[16,20.5]
    yticks=4
    ytickv=indgen(yticks+1)+16
;yrange=[14.5,20.5]
;yticks=5
;ytickv=indgen(yticks+1)+15
  endif else begin
    y+=1
    yrange=[2,8]
    yticks=6
    ytickv=indgen(yticks+1)+2
;yrange=[1,8]
;yticks=7
;ytickv=indgen(yticks+1)+1
  endelse

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  plot,x,y,/nodata,position=position,$
    xstyle=9,ystyle=9,$
    yticklen=-0.022,xticklen=-0.018,$
    yticks=yticks,ytickv=ytickv,$;ytickname=ytickname,$
    xrange=xrange,xminor=2,yrange=yrange,yminor=4,$
    xtitle=xtitle,ytitle=ytitle,$
    charsize=csize,$
    title=title

  loadct,figspecs.col_table,/silent,file=dirs.ctfil

  ;FILL SHADING
    for i=0,1 do $
      contour,var1*figspecs.scale,x,y[t_ind_final],/cell_fill,/overplot,$
        levels=figspecs.levels,c_colors=figspecs.colors

  loadct,0,/silent

  ;CONTOURS
    if cont_var ge 0 then begin
      contour,cvar,x,y[t_ind_final],/follow,/overplot,$
        levels=clevs,c_colors=0,c_thick=1.5,c_charsize=0.9*csize;,c_labels=replicate(1,200)
      contour,cvar,x,y[t_ind_final],/follow,/overplot,c_linestyle=1,$
        levels=-1.*reverse(clevs),c_colors=0,c_thick=1.5,c_charsize=0.9*csize;,c_labels=replicate(1,200)
    endif

  ;MARKER FOR ONSET OF RI
  ionset=0
  if ionset and dirs.cases[ic] eq 'ctl' then begin
    if tcname eq 'maria' then t_onset=18.12
    if tcname eq 'haiyan' then t_onset=5.5
    plots,!x.crange,replicate(t_onset,2),linestyle=2,thick=1
  endif

  ;BOX AROUND PLOT
    for i=0,1 do begin
      plots,xrange,replicate(yrange[i],2),linestyle=0,color=0,thick=1,/data
      plots,replicate(xrange[i],2),yrange,linestyle=0,color=0,thick=1,/data
    endfor

  ;PLOT RMW
  if keyword_set(rmw) then begin
;    loadct,7,/silent
    oplot,rmw,y[t_ind_final],color=0,linestyle=0,thick=4
;    loadct,0,/silent
  endif

  ;COLOR BAR
  if figspecs.icbar then begin
    ybuff=0.21
    cpos= [ position[2]+0.020 ,$
            position[1]+ybuff ,$
            position[2]+0.039 ,$
            position[3]-ybuff ]
    loadct,figspecs.col_table,/silent,file=dirs.ctfil
    if var1_str eq 'W' then begin
      ndivs=(n_elements(figspecs.levels)-1)/2+1
      figspecs.ndivs=ndivs
      setlevs=strtrim(figspecs.levels[indgen(ndivs)*2],2)
      figspecs.ndivs-=1
    endif else setlevs=''
    colorbar2, colors=figspecs.colors, range=[min(figspecs.levels),max(figspecs.levels)],divisions=figspecs.ndivs,$
      charsize=csize*0.9, position=cpos, /right, /vertical, title=figspecs.cbar_tag,$
      annotatecolor='black',format=figspecs.cbar_format,$
      setlevels=setlevs
    loadct,0,/silent
  endif

  device,/close

  convert_png,figname,res=400,/remove_eps


endfor ; icase

print,'DONE!!'
end
