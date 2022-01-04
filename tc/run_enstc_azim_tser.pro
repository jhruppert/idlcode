; 
; Plot time series from azimuthally averaged ensemble TC output.
;
; James Ruppert
; 12/26/21
; 
pro run_enstc_azim_tser

tcname='maria'
case_str='ctl'

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

;AZIM FILE INFO
  hr_plot=[0,144]
  hr_tag_full=strtrim(hr_plot[0],2)+'-'+string(hr_plot[1],format='(i3.3)')+'hr'
  hr_fil=hr_tag_full
  print,'Plotting: ',hr_tag_full

;EXTRA TIME SPECS
  nt_full=dims.nt
  npd=dims.npd
  nhrs=1.*nt_full*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nhrs)

;AZIM FILES
  hr_sel=[0,144]
  t_ind=where((time_hrs ge hr_sel[0]) and (time_hrs le hr_sel[1]))
  t_ind=t_ind[where(t_ind le nt_full-1)]
  t_ind_test=t_ind
  nt_sav=n_elements(t_ind)
  nt_test=nt_sav

;----PLOT OPTIONS--------------------

ismooth=0;1 ; smooth time series?

radius_range=[300,500];800] ; max is 800 km
radius_range=[0,600] ; max is 1400 km xx800 km

var1_str='RTHRATLW'
;var1_str='RTHRATTOT'
;var1_str='RTHRATSW'
;var1_str='RH';'QVAPOR';
;var1_str='W'
;var1_str='T'
;var1_str='AVOR'
;var1_str='H_DIABATIC'
;var1_str='rainrate'
;var1_str='psi'
;var1_str='shear'
var1_str='wind'
;var1_str='rmw'
;var1_str='sef'

remove_azmn=0
;remove_azmn=1

icrf=0

;----READ VARS--------------------


print,'VAR: ',var1_str

var1_tser=fltarr(dirs.nens,nt_sav)
var1_tser[*]=!values.f_nan

;for ic=0,dirs.nens-1 do begin
for ic=0,0 do begin

  print,'Memb: ',ic+1

  ;2D VARS

  ;SLP
    file=dirs.ensdir[ic]+'azim_SLP_'+hr_fil+'.nc'
    slp=reform(read_nc_var(file,'SLP',count=count,offset=offset))

  t_ind_final=reform(read_nc_var(file,'time'))
  nt_test=n_elements(t_ind_final)
  it0=0

  radius=read_nc_var(file,'radius')
  azimuth=read_nc_var(file,'azmiuth')
  nrad=n_elements(radius)
  naz=n_elements(azimuth)

  ;3D VARS

    count=[nrad,naz,dims.np,nt_test] & offset=[0,0,0,it0] ; x,y,z,t

  ;VAR1

    if var1_str eq 'RH' then begin

      ;NEED LARGE PRES ARRAY FOR RH
        prx=fltarr(nrad,naz,dims.np,nt_test)
        for iz=0,dims.np-1 do prx[*,*,iz,*]=dims.pres[iz]*1e2

      file=dirs.ensdir[ic]+'azim_QVAPOR_'+hr_fil+'.nc'
      qv=read_nc_var(file,'QVAPOR',count=count,offset=offset)
      file=dirs.ensdir[ic]+'azim_T_'+hr_fil+'.nc'
      tmpk=read_nc_var(file,'T',count=count,offset=offset)

      var1=calc_relh(qv,tmpk,prx)

      if ~write_qv then qv=0
      tmpk=0

    endif else if var1_str eq 'RTHRATTOT' then begin

      file=dirs.ensdir[ic]+'azim_RTHRATLW_'+hr_fil+'.nc'
      var1=read_nc_var(file,'RTHRATLW',count=count,offset=offset)
      file=dirs.ensdir[ic]+'azim_RTHRATSW_'+hr_fil+'.nc'
      var1+=read_nc_var(file,'RTHRATSW',count=count,offset=offset)

    endif else if var1_str eq 'rainrate' then begin

      file=dirs.ensdir[ic]+'azim_'+var1_str+'_'+hr_fil+'.nc'
      var1=reform(read_nc_var(file,var1_str,count=[nrad,naz,1,nt_test],offset=offset))

    endif else if var1_str eq 'sef' then begin

      file=dirs.ensdir[ic]+'azim_LH_'+hr_fil+'.nc'
      var1=reform(read_nc_var(file,'LH',count=[nrad,naz,1,nt_test],offset=offset))
      file=dirs.ensdir[ic]+'azim_HFX_'+hr_fil+'.nc'
      var1+=reform(read_nc_var(file,'HFX',count=[nrad,naz,1,nt_test],offset=offset))

    endif else if var1_str eq 'shear' then begin

      iz850=(where(dims.pres eq '850'))[0]
      iz200=(where(dims.pres eq '200'))[0]
      file=dirs.ensdir[ic]+'azim_U_'+hr_fil+'.nc'
      u850=reform(read_nc_var(file,'U',count=[nrad,naz,1,nt_test],offset=[0,0,iz850,it0]))
      u200=reform(read_nc_var(file,'U',count=[nrad,naz,1,nt_test],offset=[0,0,iz200,it0]))
      file=dirs.ensdir[ic]+'azim_V_'+hr_fil+'.nc'
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

      psel=1000;950;700 ; (hPa)
      izsel=(where(dims.pres eq psel))[0]

      file=dirs.ensdir[ic]+'azim_U10_'+hr_fil+'.nc'
      u=reform(read_nc_var(file,'U10',count=[nrad,naz,1,nt_test],offset=[0,0,0,it0]))
      file=dirs.ensdir[ic]+'azim_V10_'+hr_fil+'.nc'
      v=reform(read_nc_var(file,'V10',count=[nrad,naz,1,nt_test],offset=[0,0,0,it0]))
;      file=dirs.ensdir[ic]+'azim_U_'+hr_fil+'.nc'
;      u=reform(read_nc_var(file,'U',count=[nrad,naz,1,nt_test],offset=[0,0,izsel,it0]))
;      file=dirs.ensdir[ic]+'azim_V_'+hr_fil+'.nc'
;      v=reform(read_nc_var(file,'V',count=[nrad,naz,1,nt_test],offset=[0,0,izsel,it0]))

      ;SUBTRACT STORM MOTION
      ;(doesn't visibly change anything)
      for it=0,nt_test-1 do begin
        u[*,*,it] -= motion_x[it+it0];+hr_plot0[0]]
        v[*,*,it] -= motion_y[it+it0];+hr_plot0[0]]
      endfor

      wnd_azim=azim_wind_conv(u,v,azimuth) & u=0 & v=0

;      var1=sqrt(u^2 + v^2)
;      var1=wnd_azim.u_rad
      var1=wnd_azim.v_tan

    endif else if var1_str eq 'psi' then begin

      ;U_RAD
        file=dirs.ensdir[ic]+'azim_U_'+hr_fil+'.nc'
        u=reform(read_nc_var(file,'U',count=count,offset=offset))
        file=dirs.ensdir[ic]+'azim_V_'+hr_fil+'.nc'
        v=reform(read_nc_var(file,'V',count=count,offset=offset))
  
        wnd_azim=azim_wind_conv(u,v,azimuth) & u=0 & v=0
        var1=wnd_azim.u_rad
  
      ;W
        file=dirs.ensdir[ic]+'azim_W_'+hr_fil+'.nc'
        var2=read_nc_var(file,'W',count=count,offset=offset)
  
      ;DENSITY
        file=dirs.ensdir[ic]+'azim_T_'+hr_fil+'.nc'
        tmpk=read_nc_var(file,'T',count=count,offset=offset)
        file=dirs.ensdir[ic]+'azim_QVAPOR_'+hr_fil+'.nc'
        qvt=read_nc_var(file,'QVAPOR',count=count,offset=offset)
        tvirt = tmpk*(1.+0.61*qvt)
        tmpk=0 & qvt=0
        var3=var2
        for iz=0,dims.np-1 do var3[*,*,iz,*] = dims.pres[iz]*1e2 / ( 287. * tvirt[*,*,iz,*] )
        tvirt=0

    endif else begin

      file=dirs.ensdir[ic]+'azim_'+var1_str+'_'+hr_fil+'.nc'
      var1=read_nc_var(file,var1_str,count=count,offset=offset)

      ;DON'T MULTIPLY W BY RHO BECAUSE VERTICAL P-INTEGRAL IS ALREADY MASS-WEIGHTED
;      if var1_str eq 'W' then begin
;        file=dirs.ensdir[ic]+'azim_T_'+hr_fil+'.nc'
;        tmpk=read_nc_var(file,'T',count=count,offset=offset)
;        file=dirs.ensdir[ic]+'azim_QVAPOR_'+hr_fil+'.nc'
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

      if var1_str eq 'RTHRATTOT' then begin
        file=dirs.ensdir[ic]+'azim_RTHRATLWC_'+hr_fil+'.nc'
        var1-=read_nc_var(file,'RTHRATLWC',count=count,offset=offset)
        file=dirs.ensdir[ic]+'azim_RTHRATSWC_'+hr_fil+'.nc'
        var1-=read_nc_var(file,'RTHRATSWC',count=count,offset=offset)
      endif else begin
        file=dirs.ensdir[ic]+'azim_'+var1_str+'C_'+hr_fil+'.nc'
        var1-=read_nc_var(file,var1_str+'C',count=count,offset=offset)
      endelse

    endif

  ;AZIMUTHALLY AVERAGE
    slp=mean(temporary(slp),dimension=2,/nan,/double)
;    olr=mean(temporary(olr),dimension=2,/nan,/double)
;    olrc=mean(temporary(olrc),dimension=2,/nan,/double)
    var1=mean(temporary(var1),dimension=2,/nan,/double)
    if var1_str eq 'psi' then var2=mean(temporary(var2),dimension=2,/nan,/double)
    if var1_str eq 'psi' then var3=mean(temporary(var3),dimension=2,/nan,/double)

  ;REMOVE AZIMUTHAL MEAN
  if remove_azmn then begin
    print,'REMOVING AZIM MEAN!!'
    radall=fltarr(nrad,dims.np)
    for ip=0,dims.np-1 do radall[*,ip]=radius
    totrad=total(radall,1,/double,/nan)
    var_mn=total(var1*radall,1,/double,/nan)/totrad
    for ip=0,dims.np-1 do var1[*,ip]-=var_mn[ip]
  endif

  ;REMOVE BELOW-GROUND POINTS
ndims=size(var1,/N_DIMENSIONS)
if ndims eq 3 then begin
  for ip=0,dims.np-1 do begin
    nan=where(slp lt dims.pres[ip],count)
    iv=reform(var1[*,ip,*])
    if count gt 0 then iv[nan]=!values.f_nan
    var1[*,ip,*]=iv
    if var1_str eq 'psi' then begin
      iv=reform(var2[*,ip,*])
      if count gt 0 then iv[nan]=!values.f_nan
      var2[*,ip,*]=iv
      iv=reform(var3[*,ip,*])
      if count gt 0 then iv[nan]=!values.f_nan
      var3[*,ip,*]=iv
    endif
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

  if var1_str eq 'H_DIABATIC' then var1*=1004. ; K units --> W/m2 (once integrated)

  if var1_str eq 'psi' then begin

    psi=var1
    psi[*]=0.

    for it=0,nt_test-1 do begin

      ;HEIGHT FROM HYDROSTATIC
      z = fltarr(dims.np) ; m
      for iz=1,dims.np-1 do $
        z[iz] = z[iz-1] + (dims.pres[iz-1]-dims.pres[iz])*1e2 / 9.81 / mean(var3[*,iz-1:iz,it],/double,/nan)
      dz=deriv(z)

      ;INTEGRAL OVER R
      psir=fltarr(nrad,dims.np)
      idr=(radius[1]-radius[0])*1e3 ; m
      for ir=1,nrad-1 do $
        for iz=0,dims.np-1 do begin
          mfr  = (radius[ir]*1e3) * var3[ir,iz,it] * var2[ir,iz,it] ; VERTICAL MASS FLUX [ kg / m s ]
          psir[ir,iz] = psir[ir-1,iz] + ( mfr*idr ) ; kg / s
        endfor
  
      ;INTEGRAL OVER Z
      psiz=fltarr(nrad,dims.np)
      ;PSI_Z should be zero at lowest model level
      for ir=0,nrad-1 do $
        for iz=1,dims.np-2 do begin
          mfz  = (radius[ir]*1e3) * var3[ir,iz,it] * var1[ir,iz,it] ; RADIAL MASS FLUX [ kg / m s ]
          psiz[ir,iz] = psiz[ir,iz-1] - ( mfz * dz[iz] ) ; kg / s
        endfor
  
      psi[*,*,it] = 0.5*(psir + psiz) ; units of kg / s

    endfor

    var1=psi
    psi=0 & var2=0 & var3=0

  endif

  ;VERTICALLY INTEGRATE
;    ip = where(dims.pres ge 50)
;    sum1 = (1./9.81) * total( var1[*,ip,*] * dp[*,ip,*] ,2,/nan,/double)

  ;VERTICALLY AVERAGE
    ip = where(dims.pres ge 50)
    sum1 = total( var1[*,ip,*] * dp[*,ip,*] ,2,/nan,/double) / total(dp[*,ip,*],2,/nan,/double)

  var1=sum1 & sum1=0

endif


  ;RADIUS SUBSET
    irad=where((radius ge radius_range[0]) and (radius le radius_range[1]) , nradtmp)
    rad_tmp=fltarr(nradtmp,nt_test)
    for it=0,nt_test-1 do rad_tmp[*,it]=radius[irad]
    tot_rad = total(radius[irad],/double)

if var1_str eq 'wind' or var1_str eq 'rmw' then begin

  var1_tser[ic,t_ind_final] = max(var1[irad,*],dimension=1,/nan)
;  var1_tser[ic,t_ind_final] = max(-1.*var1[irad,*],dimension=1,/nan)

  if var1_str eq 'rmw' then begin
    for it=0,nt_test-1 do begin
      max=max(reform(var1[irad,it]),loc,/nan)
      var1_tser[ic,t_ind_final[it]]=radius[loc]
    endfor
  endif

endif else if var1_str eq 'shear' then begin

  ;RADIALLY AVERAGE
    au200=total( u200[irad,*] * rad_tmp ,1,/double,/nan) / tot_rad
    au850=total( u850[irad,*] * rad_tmp ,1,/double,/nan) / tot_rad
    av200=total( v200[irad,*] * rad_tmp ,1,/double,/nan) / tot_rad
    av850=total( v850[irad,*] * rad_tmp ,1,/double,/nan) / tot_rad

  ;CALCULATE SHEAR
    udif=au200-au850
    vdif=av200-av850

  ;SMOOTH
;    nsmooth=3
;    if ismooth then $
;      tot=smooth(temporary(tot),nsmooth,/edge_truncate)

  var1_tser[ic,t_ind_final] = sqrt(udif^2+vdif^2);tot

endif else if var1_str eq 'psi' then begin

  ;RADIAL AVERAGE
    tot=total( var1[irad,*] * rad_tmp ,1,/double,/nan) / tot_rad

  ;PSI IS NOW VERTICALLY INTEGRATED, RADIALLY AVERAGED
  ; kg/s * s2/m * kg/m/s2 --> kg^2 / m2 / s

  ;NOW:
  ;PSI IS NOW VERTICALLY AVERAGED, RADIALLY AVERAGED
  ; kg/s

;  var1_tser[ic,t_ind_final] = tot*1e-12 ; 10^12 kg / m2 / s
  var1_tser[ic,t_ind_final] = tot*1e-8 ; 10^8 kg / s

  if ismooth then begin
  ;FILL IN EACH CASE WITH CTL, THEN OVERWRITE LATER TIME STEPS
    if dirs.cases[ic] eq 'ctl' then $
      for icc=1,dirs.nc-1 do begin
        var1_tser[icc,t_ind_final]=var1_tser[0,t_ind_final]
      endfor
  ;SMOOTH TIME SERIES
    nsmooth=[0,6]
    ;for i=0,1 do $
      var1_tser[ic,*]=smooth(var1_tser[ic,*],nsmooth,/edge_truncate,/nan)
  ;NOW REPLACE NANS FOR TIMES PRIOR TO TEST START
    if dirs.cases[ic] ne 'ctl' then $
      var1_tser[ic,0:t_ind_final[0]-1]=!values.f_nan
  endif

endif else begin

  ;RADIALLY AVERAGE
    tot=total( var1[irad,*] * rad_tmp ,1,/double,/nan) / tot_rad

  ;SMOOTH
  if ismooth then begin
    ;FILL IN EACH CASE WITH CTL, THEN OVERWRITE LATER TIME STEPS
      if dirs.cases[ic] eq 'ctl' then $
        for icc=1,dirs.nc-1 do begin
          var1_tser[icc,t_ind_final]=var1_tser[0,t_ind_final]
        endfor
    ;SMOOTH TIME SERIES
      nsmooth=[0,12];6]
      ;for i=0,1 do $
        var1_tser[ic,*]=smooth(var1_tser[ic,*],nsmooth,/edge_truncate,/nan)
    ;NOW REPLACE NANS FOR TIMES PRIOR TO TEST START
      if dirs.cases[ic] ne 'ctl' then $
        var1_tser[ic,0:t_ind_final[0]-1]=!values.f_nan
  endif

endelse

;endif else if var1_str eq 'W' then begin

;  var1_tser[ic,t_ind_final] = tot

;endif


endfor ; icase


;----CREATE PLOT--------------------


  if var1_str eq 'H_DIABATIC' then begin
    ytitle1='MP Heating [ W m!U-2!N ]'
    figtag=var1_str
;    yrange=[600,2400] ; For 0-300 radius
;    yrange=[300,1900] ; 0-400
    yrange=[300,1400] ; 0-500
;    yrange=[200,1020] ; 0-600
;    yrange=[0,800] ; 0-800
  endif else if var1_str eq 'W' then begin
    ytitle1='VMF [ kg m!U-1!N s!U-1!N ]'
;    ytitle1='W" [ m!U2!N s!U-2!N ]'
    figtag='wmean'
    yrange=[0,500]
  endif else if var1_str eq 'rainrate' then begin
    ytitle1='Rain rate [ mm/h ]'
    figtag='rain'
    yrange=[0,60]
  endif else if var1_str eq 'shear' then begin
    ytitle1='Shear [ m/s ]'
    figtag='shear'
    yrange=[0,30]
  endif else if var1_str eq 'wind' then begin
    ytitle1='Max(v) [ m/s ]'
    figtag='wind'
    yrange=[0,80]
;    ytitle1='Max(-u) [ m/s ]'
;    yrange=[0,10]
  endif else if var1_str eq 'sef' then begin
    ytitle1='SEF [ W m!U-2!N ]'
    figtag='sef'
    yrange=[0,1000]
  endif else if var1_str eq 'rmw' then begin
    ytitle1='RMW [ km ]'
    figtag='rmw'
    yrange=[0,500]
  endif else if var1_str eq 'psi' then begin
;    ytitle1='Psi [ 10!U12!N kg m!U-2!N s!U-1!N ]'
    ytitle1='Psi [ 10!U8!N kg/s ]'
    figtag='psi'
    yrange=[0,15]
  endif


  figdir=dirs.figdir
  rad_str='rad'+string(radius_range[0],format='(i3.3)')+'-'+string(radius_range[1],format='(i3.3)')+'km'
  figname=figdir+'azim_tser/tser_azim_'+figtag+'_'+rad_str

  if strmatch(var1_str,'*RTHRA*') then begin
    if icrf then begin
      figname+='_crf'
    endif
  endif

;  if icrf then olr-=olrc

  if remove_azmn then figname+='_azprm'

  figname+='_'+hr_tag_full

  ;PLOT SPECS
    csize=0.8
    position=[0.14,0.20,0.89,0.89]
    xsize=4.2 & ysize=2
;    xtitle='Time [ hr ]'
    xtitle='Date in Sept [ UTC ]'
    ytitle2='OLR [ W m!U-2!N ]'
    yrange2=[60,310]

  ;AXES
    x=indgen(nt_sav)
    y=indgen(10)
    if ~keyword_set(xrange) then $
      xrange=[min(x),max(x)]
    if ~keyword_set(yrange) then $
      yrange=[min(y),max(y)]

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  xticks=5;3
  xtickv=indgen(xticks+1)*24+12
;  xtickname=strtrim(indgen(xticks+1)+16,2)
  xtickname=strtrim(indgen(xticks+1)+15,2)

  plot,x,y,/nodata,position=position,$
    xstyle=9,ystyle=9,$
;    yticklen=yticklen,xticklen=-0.033,$
    xticks=xticks,xtickv=xtickv,xtickname=xtickname,$
    xrange=xrange,xminor=2,yrange=yrange,yminor=2,$
    xtitle=xtitle,ytitle=ytitle1,$
    charsize=csize,$
    title=title

  cols=(findgen(dirs.nc)+1)/dirs.nc*220
  lstyle=replicate(0,dirs.nc);indgen(dirs.nc)
;  lstyle[[4,5]]=[1,2]
  lthick=reverse((findgen(dirs.nc)+1)/dirs.nc*3.5)

  loadct,3,/silent

  ;MAIN VAR
  for ic=0,dirs.nc-1 do $
    oplot,x,var1_tser[ic,*],linestyle=lstyle[ic],thick=lthick[ic],color=cols[ic]

;  loadct,0,/silent

  ;OLR
;  axis,yaxis=1,ystyle=9,ytitle=ytitle2,charsize=csize,yrange=yrange2,yminor=2,yticklen=yticklen,/save
;  oplot,x,olr,linestyle=2,thick=2,color=0

  ;LEGEND
  ileg=1
  if ileg then begin
;  if ileg then begin
    csize_fac=0.7
    margin=0.1
    pspacing=2.0 ; length of lines
    spacing=0.8 ; between lines
    ncplot=dirs.nc
    icplot=indgen(ncplot)
    leg_str=strupcase(dirs.cases[icplot])
    leg_style=lstyle[icplot]
    leg_thick=lthick
    leg_color=cols;replicate(0,ncplot)
    legend2,leg_str,linestyle=leg_style,thick=leg_thick,COLORS=leg_color,$
      charsize=csize*csize_fac,/top_legend,/left_legend,/normal,$
      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.16,0.75]
;      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.56,0.35]
  endif

  device,/close

  convert_png,figname,res=200;,/remove_eps


print,'DONE!!'
end
