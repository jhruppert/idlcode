;
; Retrieve variables and plotting specs for ya-averaged cross sections
; from ICON 3D output.
; 
; Contains: yavg_dcomp
;
; James Ruppert
; 16.12.17
;
pro yavg_call_plots,itimeheight,exp_name,datfils,fil_ind,t_ind,plot_var,hght,hghtw,$ ; INPUT
  figname=figname,pres_plot=pres_plot

npd=24

;READ DATA
  yavg_read,datfils,fil_ind,t_ind,struct_vars
  x_yavg=struct_vars.x_yavg
  nx=n_elements(x_yavg)
  dx=x_yavg[1]-x_yavg[0] ; [ m ]
  specs=size(struct_vars.u_yavg,/dimensions)
  nt=specs[2]
  nz=specs[1]

;  pcp = struct_vars.pcp_yavg ; mm/h
;  pres = struct_vars.pres_yavg ; mm/h

;DIURNAL COMPOSITE
idcomp=1
if idcomp then begin
  yavg_dcomp,npd,struct_vars,struct_vars_dcomp
  struct_vars=struct_vars_dcomp
endif

;NEW TIME SPECS (IN CASE OF DCOMP)
  pcp = struct_vars.pcp_yavg ; mm/h
  specs=size(struct_vars.u_yavg,/dimensions)
  nt=specs[2]

;CALCULATE STREAMFUNCTION (PSI)

  tmpk=struct_vars.tmpk_yavg
  ex=struct_vars.ex_yavg
  qv=struct_vars.qv_yavg

  reps=461.5/287.04
  t_v=tmpk*(1.0+reps*qv)/(1.0+qv)
  rho = exner(ex,/reverse) / (287.*t_v) ; kg/m3

  w=struct_vars.w_yavg ; m/s
;  w=reverse(w,1)
  u=struct_vars.u_yavg ; m/s
;  u=reverse(u,1)
;  rho=reverse(rho,1)

  psi=fltarr(nx,nz,nt)
  psir=psi
  psiz=psi
  for it=0,nt-1 do begin

;  for irad=1,nrad-1 do begin
;  for irad=nrad-2,1,-1 do $
;    for iz=0,nz-1 do $
      ;psi[irad,iz,it] = psi[ix-1,iz,it] + ( rho[ix,iz,it] * w[ix,iz,it] ) ; kg / m2 / s
;      psi[irad,iz] = psi[irad-1,iz] + ( rho[irad,iz] * $
;                                       ( w[irad,iz] + w[irad,iz+1] )*0.5 ) ; kg / m2 / s
;      psi[irad,iz] = psi[irad-1,iz] + x_yavg[irad]*rho[irad,iz] * $
;        ( 0.5*(w[irad,iz]+w[irad,iz+1])*dr - $
;        u[irad,iz]*(hghtw[iz]-hghtw[iz+1]) )

    wi = reform(w[*,*,it])
    ui = reform(u[*,*,it])

    ;SMOOTH VELOCITIES IN R
    for i=0,2 do begin
      wi=gauss_smooth(wi,[5,0],/edge_truncate,/nan)
      ui=gauss_smooth(ui,[5,0],/edge_truncate,/nan)
    endfor

    ;INTEGRAL OVER R
    for iz=0,nz-1 do $
      for irad=nx-2,0,-1 do $
;      for irad=1,nrad-1 do $
        psir[irad,iz,it] = psir[irad+1,iz,it] - x_yavg[irad] * rho[irad,iz,it] * $
          ( wi[irad,iz] + wi[irad,iz+1] )*0.5*dx

    ;INTEGRAL OVER Z
    for irad=0,nx-1 do $
      for iz=nz-2,0,-1 do $
        psiz[irad,iz,it] = psiz[irad,iz+1,it] - x_yavg[irad] * rho[irad,iz,it] * $
          ui[irad,iz] * (hghtw[iz]-hghtw[iz+1])

  endfor

  psi=(psir+psiz)*1e-7 ; 10^8 kg/s
  psi[*]=0
;  psi=psir*1e-8
;  psi=psiz*1e-8

;  psi=reverse(psi,1)
;  psi *= 10 ; 10^-1 kg / m2 / s
;  psi *= 1e-5 ; 10^-5 m2 / s


if plot_var eq 'pres' then begin

  title=strupcase(exp_name)

  if pres_plot eq 'flux' then begin
    v1 = struct_vars.sh_yavg ; W/m2
    v2 = struct_vars.lh_yavg ; W/m2
    vtag1='SH' & vtag2='LH'
    axrange=[-50,400]
    axtitle='Surface Flux [ W m!U-2!N ]'
  endif else if pres_plot eq 'pres' then begin
    v1 = struct_vars.u_yavg ; radial m/s
    v2 = struct_vars.v_yavg ; tangential m/s
    vtag1='u (radial)' & vtag2='v (tangen)'
    axrange=[-20,60]
    axtitle='Wind Speed [ m s!U-1!N ]'
  endif else if pres_plot eq 'theta' then begin
    tmpk = struct_vars.tmpk_yavg
    ex = struct_vars.ex_yavg
    th = tmpk / ex
;    reps=461.5/287.04
;    th_virt = th * (1.0+reps*qv)/(1.0+qv)
    v1 = th - mean(th,/nan,/double)
    vtag1="theta'"
    axrange=[-3,9]
    axtitle='Theta [ K ]'
  endif

  pres = struct_vars.pres_yavg * 1e-2 ; hPa

  azimuth_cross_plot2d,figname,title,x_yavg,pres,pcp,$
    v1=v1,v2=v2,$
    axrange=axrange,axtitle=axtitle,vtag1=vtag1,vtag2=vtag2

endif else if plot_var eq 'lw' or plot_var eq 'sw' or plot_var eq 'rad' then begin

  if plot_var eq 'lw' then begin
    title='LW'
    maxv=8 ; max shading value
    minv=-8 ; max shading value
    col_table=1
    var1=struct_vars.lw_yavg
  endif else if plot_var eq 'sw' then begin
    title='SW'
    maxv=10 ; max shading value
    minv='0' ; max shading value
    col_table=3
    var1=struct_vars.sw_yavg
  endif else if plot_var eq 'rad' then begin
    title="!8Q!DR!N'!X";+cgSymbol('prime')
    maxv=12 ; max shading value
    minv=-12 ; max shading value
    col_table=70
    var1=struct_vars.lw_yavg
    var1+=struct_vars.sw_yavg
  endif

  ;FOR TIME-HEIGHT
  xind_var1=where((x_yavg ge 50e3) and (x_yavg le 300e3))
  xind_var2=xind_var1
  xind_pcp=xind_var1
  xind_psi=where((x_yavg ge 200e3) and (x_yavg le 500e3))

  ;FOR TIME-HEIGHT
  ;SUBTRACT DOMAIN AVERAGE
  if itimeheight then begin
    maxv=6 ; max shading value
    minv=-6 ; max shading value
    clevels2=[5,25,50,100,150,200]
    v1m=mean(var1,dimension=1,/double)
    for irad=0,nrad-1 do $
      var1[irad,*,*] -= v1m
  endif

  var1*=3600.*24 ; K/s --> K/d
  var2=var1 & var2[*]=0

  ;CLOUDS
  icloud=1
  if icloud then begin
    var2=struct_vars.qc_yavg * 1e6
    var2+=struct_vars.qi_yavg * 1e6 ; [ mg/kg ]
  endif

  ;FOR FILL
  nlevs=255;25
  colors=findgen(nlevs)*254/(nlevs-1)
  colors=reverse(colors)
  levels=findgen(nlevs)/(nlevs-1)*(maxv-minv)+minv

  ;CONTOURS
  icont=[0,icloud]
  cstyle=[0,0]
  clevels2=[5,25,50,100,125,150,175]

  icbar=1
  cbar_tag='[ K d!U-1!N ]'
  cbar_format='(i3)'

  title+=' ('+strupcase(exp_name)+')'

endif else if plot_var eq 'u' or plot_var eq 'w' then begin

  if plot_var eq 'u' then begin
    title='u'
    maxv=10 ; max shading value
    minv=-10 ; max shading value
    col_table=70
    var1=struct_vars.u_yavg ; [ m/s ]
    cbar_tag='[ m s!U-1!N ]'
  endif else if plot_var eq 'w' then begin
    title='w'
    maxv=15 ; max shading value
    minv=-15 ; max shading value
    col_table=70
    var1=struct_vars.w_yavg*1e2 ; [ cm/s ]
    var12=struct_vars.u_yavg
    for it=0,nt-1 do $
      for ix=0,nx-1 do $
        var12[ix,*,it]=interpol(reform(var1[ix,*,it]),findgen(nz+1),findgen(nz)+0.5)
    var1=var12
    cbar_tag='[ cm s!U-1!N ]'
  endif

  ;FOR TIME-HEIGHT
  xind_var1=where((x_yavg ge 50e3) and (x_yavg le 300e3))
  xind_var2=xind_var1
  xind_pcp=xind_var1
  xind_psi=where((x_yavg ge 200e3) and (x_yavg le 500e3))

  var2=var1 & var2[*]=0

  ;FOR FILL
  nlevs=255;25
  colors=findgen(nlevs)*254/(nlevs-1)
  colors=reverse(colors)
  levels=findgen(nlevs)/(nlevs-1)*(maxv-minv)+minv

  ;CONTOURS
  icont=[0,0]
  cstyle=[0,0]
  clevels2=[5,25,50]

  icbar=1
  cbar_format='(i3)'

  title+=' ('+strupcase(exp_name)+')'

endif else if plot_var eq 'tmp' then begin

  title="!8T'!X"
  maxv=10 ; max shading value
  minv=-10 ; max shading value
  col_table=70
  var1=struct_vars.tmpk_yavg ; [ K ]
  cbar_tag='[ K ]'

  ;FOR TIME-HEIGHT
  xind_var1=where((x_yavg ge 50e3) and (x_yavg le 300e3))
  xind_var2=xind_var1
  xind_pcp=xind_var1
  xind_psi=where((x_yavg ge 200e3) and (x_yavg le 500e3))

  var2=var1 & var2[*]=0

  ;FOR FILL
  nlevs=255;25
  colors=findgen(nlevs)*254/(nlevs-1)
  colors=reverse(colors)
  levels=findgen(nlevs)/(nlevs-1)*(maxv-minv)+minv

  ;CONTOURS
  icont=[1,0]
  cstyle=[0,0]
  cint=1
  clevels1=findgen(200)*cint+cint
  clevels2=[5,25,50]

  icbar=1
  cbar_format='(i3)'

  title+=' ('+strupcase(exp_name)+')'

endif else if strmatch(plot_var,'q1_*') then begin

  title='Theta-dot'
  maxv=10 ; max shading value
  minv=-10 ; max shading value
  if itimeheight then begin
    maxv=5 ; max shading value
    minv=-5 ; max shading value
  endif
  col_table=70
  cbar_tag='[ K h!U-1!N ]'

  ;HEAT BUDGET
    tmpk=struct_vars.tmpk_yavg ; K
    ex=struct_vars.ex_yavg ; -
    th=tmpk/ex
    u=struct_vars.u_yavg ; m/s
  ;  v=struct_vars.v_yavg ; m/s
    w=struct_vars.w_yavg ; m/s
    azimuth_heat_budget,npd,x_yavg,hght,th,u,w,struct_hbudg

  tag2=strmid(plot_var,3)
  title+=' ('+tag2+')'

  q1 = struct_hbudg.ddt + struct_hbudg.uddx + struct_hbudg.wddz

  if tag2 eq 'total' then begin
    var1=q1
  endif else if tag2 eq 'ddt' then begin
    var1=struct_hbudg.ddt
  endif else if tag2 eq 'uddx' then begin
    var1=struct_hbudg.uddx
  endif else if tag2 eq 'wddz' then begin
    var1=struct_hbudg.wddz
  endif else if tag2 eq 'qc' then begin
    qr = struct_vars.lw_yavg + struct_vars.sw_yavg
    qr *= 3600. ; K/s --> K/h
    var1 = q1 - qr
  endif

  ;FOR TIME-HEIGHT
  xind_var1=where((x_yavg ge 50e3) and (x_yavg le 300e3))
  xind_var2=xind_var1
  xind_pcp=xind_var1
  xind_psi=where((x_yavg ge 200e3) and (x_yavg le 500e3))

  var2=var1 & var2[*]=0

  ;FOR FILL
  nlevs=255;25
  colors=findgen(nlevs)*254/(nlevs-1)
  colors=reverse(colors)
  levels=findgen(nlevs)/(nlevs-1)*(maxv-minv)+minv

  ;CONTOURS
  icont=[1,0]
  cstyle=[0,0]
  cint=3
  clevels1=findgen(200)*cint+cint
  clevels2=[5,25,50]

  icbar=1
  cbar_format='(i3)'

  title+=' ('+strupcase(exp_name)+')'

endif else message,'Pick a valid variable name!'


;CREATE PLOTS

  if plot_var ne 'pres' then begin

  sh = struct_vars.sh_yavg
  lh = struct_vars.lh_yavg
  t_g = struct_vars.t_g_yavg; - 273.15

  ;CROSS SECTIONS
  if ~itimeheight then $
  
    yavg_cross_plot3d,figname,title,x_yavg,hght,pcp,$
    levels,colors,col_table,$
    icbar,cbar_tag,cbar_format,$
    icont,cstyle,$
    clevels1,clevels2,$
    sh=sh,lh=lh,t_g=t_g,$
    v1s=var1,v2s=var2,psi=psi $
  
    ;TIME-HEIGHT PLOT
  else begin
  
    ;AVERAGE VARIABLES
    var1m = mean(var1[xind_var1,*,*],dimension=1,/double,/nan)
    var2m = mean(var2[xind_var2,*,*],dimension=1,/double,/nan)
    pcpm = mean(pcp[xind_pcp,*],dimension=1,/double,/nan)
  
    ;AVERAGE PSI
    psim = mean(psi[xind_psi,*,*],dimension=1,/double,/nan)
stop
    azimuth_cross_timeheight,figname,title,hght,pcpm,$
      levels,colors,col_table,$
      icbar,cbar_tag,cbar_format,$
      icont,cstyle,$
      clevels1,clevels2,$
      v1s=var1m,v2s=var2m,psi=psim
  
  endelse

  endif

end

; 
; CREATE DIURNAL COMPOSITE OF AZIMUTHALLY AVERAGED VARIABLES
; 
pro yavg_dcomp,npd,struct_vars,struct_vars_dcomp

  specs=size(struct_vars.u_yavg,/dimensions)
  nrad=specs[0]
  nz=specs[1]
  nt=specs[2]

;DIURNAL COMPOSITE

  nd=nt/npd
  if nd ne 1.*nt/npd then stop
  pcp_dc=fltarr(nrad,npd)
  pw_dc=pcp_dc & lh_dc=pcp_dc & sh_dc=pcp_dc
  t_g_dc=pcp_dc
  w_so_dc=fltarr(nrad,8,npd)
  t_so_dc=fltarr(nrad,9,npd)
  w_dc=fltarr(nrad,nz+1,npd)
  ex_dc=fltarr(nrad,nz,npd)
  tmpk_dc=ex_dc & u_dc=ex_dc & v_dc=ex_dc & lw_dc=ex_dc & sw_dc=ex_dc
  qv_dc=ex_dc & qc_dc=ex_dc & qi_dc=ex_dc
  for it=0,npd-1 do begin

    t_ind=indgen(nd)*npd+it

    pcp_dc[*,it]=mean(struct_vars.pcp_yavg[*,t_ind],dimension=2,/double,/nan)
    pw_dc[*,it]=mean(struct_vars.pw_yavg[*,t_ind],dimension=2,/double,/nan)
    lh_dc[*,it]=mean(struct_vars.lh_yavg[*,t_ind],dimension=2,/double,/nan)
    sh_dc[*,it]=mean(struct_vars.sh_yavg[*,t_ind],dimension=2,/double,/nan)
    t_g_dc[*,it]=mean(struct_vars.t_g_yavg[*,t_ind],dimension=2,/double,/nan)

    w_so_dc[*,*,it]=mean(struct_vars.w_so_yavg[*,*,t_ind],dimension=3,/double,/nan)
    t_so_dc[*,*,it]=mean(struct_vars.t_so_yavg[*,*,t_ind],dimension=3,/double,/nan)

    w_dc[*,*,it]=mean(struct_vars.w_yavg[*,*,t_ind],dimension=3,/double,/nan)
    ex_dc[*,*,it]=mean(struct_vars.ex_yavg[*,*,t_ind],dimension=3,/double,/nan)
    tmpk_dc[*,*,it]=mean(struct_vars.tmpk_yavg[*,*,t_ind],dimension=3,/double,/nan)
    u_dc[*,*,it]=mean(struct_vars.u_yavg[*,*,t_ind],dimension=3,/double,/nan)
    v_dc[*,*,it]=mean(struct_vars.v_yavg[*,*,t_ind],dimension=3,/double,/nan)
    lw_dc[*,*,it]=mean(struct_vars.lw_yavg[*,*,t_ind],dimension=3,/double,/nan)
    sw_dc[*,*,it]=mean(struct_vars.sw_yavg[*,*,t_ind],dimension=3,/double,/nan)
    qv_dc[*,*,it]=mean(struct_vars.qv_yavg[*,*,t_ind],dimension=3,/double,/nan)
    qc_dc[*,*,it]=mean(struct_vars.qc_yavg[*,*,t_ind],dimension=3,/double,/nan)
    qi_dc[*,*,it]=mean(struct_vars.qi_yavg[*,*,t_ind],dimension=3,/double,/nan)

  endfor


struct_vars_dcomp = { pcp_yavg:pcp_dc, $
  lh_yavg:lh_dc, sh_yavg:sh_dc, pw_yavg:pw_dc, t_g_yavg:t_g_dc, $
  w_so_yavg:w_so_dc, t_so_yavg:t_so_dc, $
  u_yavg:u_dc, v_yavg:v_dc, w_yavg:w_dc, $
  lw_yavg:lw_dc, sw_yavg:sw_dc, $
  ex_yavg:ex_dc, tmpk_yavg:tmpk_dc, qv_yavg:qv_dc, $
  qi_yavg:qi_dc, qc_yavg:qc_dc }


end
