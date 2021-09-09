; 
; Plot integrated mass flux for all simulations.
;
; James Ruppert
; 5/30/19
; 
pro run_tc_massflx

tcname='maria';'edouard';
tcyear='2017'
hurdat=read_hurdat(tcname,tcyear)

;subdir='moving_nest/'+tcname
subdir='static_nest/'+tcname
tc_sim_config, subdir, dirs=dirs, dims=dims, vars=vars;, /verbose
dirs.figdir+=tcname+'/'

;AZIM FILE INFO
  hr_plot=[49,72]
;  hr_plot=[61,84]
;  hr_plot=[61,84]+24;12
;  hr_plot=[97,120]
  hr_tag_check=string(hr_plot[0],format='(i3.3)')+'-'+string(hr_plot[1],format='(i3.3)')+'hr'
  print,'Plotting: ',hr_tag_check

;EXTRA TIME SPECS
  nt_full=dims.nt-1
  npd=dims.npd
  nhrs=1.*nt_full*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nhrs)

;AZIM FILES
  hr_sel=[24,120]
  t_ind=where((time_hrs ge hr_sel[0]) and (time_hrs le hr_sel[1]))
  t_ind=t_ind[where(t_ind le nt_full-1)]
  nt_sav=n_elements(t_ind)


;----READ VARS--------------------


;AZIMUTHAL SPECS
  nrad=267
  naz=360

inflow=fltarr(dirs.nc,nrad)
inflow[*]=!values.f_nan
outflow=inflow
vertmf=inflow
vertmfz=fltarr(dirs.nc,dims.np)
vertmfz[*]=!values.f_nan

for ic=0,dirs.nc-1 do begin
;for ic=0,0 do begin

  print,'CASE: ',strupcase(dirs.cases[ic])

  hr0=0
  if strmatch(dirs.cases[ic],'*36h*') then hr0=36
  if strmatch(dirs.cases[ic],'*24h*') then hr0=24
  if strmatch(dirs.cases[ic],'*48h*') then hr0=48
  if strmatch(dirs.cases[ic],'*60h*') then hr0=60
  if strmatch(dirs.cases[ic],'*72h*') then hr0=72
  if strmatch(dirs.cases[ic],'*84h*') then hr0=84
  nt_test_full=nt_full - hr0

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

  if hr_tag_plot ne hr_tag_check then begin
    print,'Not enough times in test, so skipping...'
    continue
  endif

  it0 = t_ind_plot[0]

  ;2D VARS

    count=[nrad,naz,1,nt_plot] & offset=[0,0,0,it0] ; x,y,z,t

  ;SLP
    file=dirs.casedir[ic]+'azim_SLP_'+hr_fil+'.nc'
    slp=reform(read_nc_var(file,'SLP',count=count,offset=offset))

  radius=read_nc_var(file,'radius')
  azimuth=read_nc_var(file,'azmiuth')

  ;3D VARS

    count=[nrad,naz,dims.np,nt_plot] & offset=[0,0,0,it0] ; x,y,z,t

  ;VAR1

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
      var2=w & var2[*]=!values.f_nan

  ;AZIMUTHALLY AVERAGE
    slp=mean(temporary(slp),dimension=2,/nan,/double)
    w=mean(temporary(w),dimension=2,/nan,/double)
    u_rad=mean(temporary(u_rad),dimension=2,/nan,/double)
    rho=mean(temporary(rho),dimension=2,/nan,/double)

  ;TIME AVERAGE
    slp=mean(temporary(slp),dimension=2,/nan,/double)
    w=mean(temporary(w),dimension=3,/nan,/double)
    u_rad=mean(temporary(u_rad),dimension=3,/nan,/double)
    rho=mean(temporary(rho),dimension=3,/nan,/double)

  ;REMOVE BELOW-GROUND POINTS
  for ip=0,dims.np-1 do begin
    nan=where(slp lt dims.pres[ip],count)
    if count gt 0 then begin
      w[nan,ip]=!values.f_nan
      u_rad[nan,ip]=!values.f_nan
      rho[nan,ip]=!values.f_nan
    endif
  endfor


;----STREAMFUNCTION--------------------


    ;HEIGHT FROM HYDROSTATIC
    z = fltarr(dims.np) ; m
    for iz=1,dims.np-1 do $
      z[iz] = z[iz-1] + (dims.pres[iz-1]-dims.pres[iz])*1e2 / 9.81 / mean(rho[*,iz-1:iz],/double,/nan)
    dz=deriv(z)

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
        psiz[ir,iz] = psiz[ir,iz-1] - ( mfz[ir,iz] * dz[iz] ) ; kg / s
      endfor

;    psi = psiz;0.5*(psir + psiz)
;    psi *= 1e-8 ; 10^8 kg / s

    ;INTEGRATE OVER LAST DIMENSION

    ;INFLOW
    locp=where(dims.pres ge 400)
    for ir=0,nrad-1 do inflow[ic,ir]=total(mfz[ir,locp]*dz[locp],/nan,/double) ; kg / s

    ;OUTFLOW
    locp=where(dims.pres lt 400 and dims.pres ge 125)
    for ir=0,nrad-1 do outflow[ic,ir]=total(mfz[ir,locp]*dz[locp],/nan,/double) ; kg / s

    ;UPDRAFT RADIAL STRUCTURE
    locp=where(dims.pres ge 125)
    for ir=0,nrad-1 do vertmf[ic,ir]=total(mfr[ir,locp]*dz[locp],/nan,/double) ; kg / s

    ;UPDRAFT VERTICAL STRUCTURE
    locr=where(radius le 500)
    for iz=0,dims.np-1 do vertmfz[ic,iz]=total(mfr[locr,iz]*idr) ; kg / s

endfor ; icase


;----CREATE PLOTS--------------------


  tc_figspecs, 'mfr', figspecs, setmax=setmax, setmin=setmin

  figdir=dirs.figdir+'/mf_cross/'
  spawn,'mkdir '+figdir,tmp,tmpe

for ifig=0,3 do begin

  if ifig eq 0 then begin
    figname=figdir+'inflow_mf'
    title='Inflow Mass Flux'
    yrange=[-25,2]
    var=inflow*1e-8
    scale='10!U8!N'
  endif else if ifig eq 1 then begin
    figname=figdir+'outflow_mf'
    title='Outflow Mass Flux'
    yrange=[-5,25]
    var=outflow*1e-8
    scale='10!U8!N'
  endif else if ifig eq 2 then begin
    figname=figdir+'vert_mf'
    title='Vertical Mass Flux'
    yrange=[-4,10]
    var=vertmf*1e-7
for i=0,1 do var=smooth(var,[0,3],/edge_truncate)
    scale='10!U7!N'
  endif else if ifig eq 3 then begin
    figname=figdir+'vert_mfz'
    title='Vertical Mass Flux'
    yrange=[-2,20]
    var=vertmfz*1e-8
    scale='10!U8!N'
  endif

  figname+='_'+hr_tag_check

  if ifig le 2 then begin
    x=radius
    xlog=0
    xrange=[max(x),min(x)]
    xtitle='Radius [ km ]'
  endif else begin
    x=dims.pres
    xlog=1
    xtitle='Pressure [ hPa ]'
    xtickv=[1000,700,500,400,300,200,150,100];,50]
    xrange=[max(x),min(xtickv)]
    xticks=n_elements(xtickv)-1
  endelse

  ;PLOT SPECS
    csize=0.8
    position=[0.14,0.18,0.89,0.89]
    xsize=4.2 & ysize=2
    ytitle='[ '+scale+' kg s!U-1!N ]'

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  plot,x,var[0,*],/nodata,position=position,xlog=xlog,$
    xstyle=9,ystyle=9,$
    xrange=xrange,xminor=2,yrange=yrange,yminor=2,$
    xticks=xticks,xtickv=xtickv,$;xtickname=xtickname,xticklen=0.034,$
    xtitle=xtitle,ytitle=ytitle,$
    charsize=csize,title=title

  plots,!x.crange,[0,0],linestyle=0,thick=0.5,/data

  loadct,3,/silent

  cols=(findgen(dirs.nc)+1)/dirs.nc*220
  cols[0]=[0]
  lstyle=replicate(0,dirs.nc)
  lthick=reverse((findgen(dirs.nc)+1)*4.5/dirs.nc)

  for ic=0,dirs.nc-1 do $
    oplot,x,var[ic,*],linestyle=lstyle[ic],thick=lthick[ic],color=cols[ic]

  device,/close
  convert_png,figname,res=200;,/remove_eps

endfor ; ifig

print,'DONE!!'
end
