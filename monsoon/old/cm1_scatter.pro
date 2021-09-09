; 
; Create a cross section from CM1 netCDF output.
; 
; James Ruppert
; 08.12.17
; 
pro cm1_scatter


;SOME SETTINGS

  nx_island=100

  iz_sel = 29 ; 4.96 km

  nd_avg=2 ; number of days to average over (at simulation end)


;DIRECTORIES
  home='/home/mpim/m300462/'
  ;  datdir=home+'cm1r19/run/'
  datdir=home+'cm1r19/output/'
  figdir=home+'misthome/IDLWorkspace85/Default/figures/monsoon/cm1/'


;exp_all=['24h_30w','24h_60w','24h_90w','24h_120w']
;hpd_all=replicate(12,4)
exp_all=['12h_120w','24h_60w','24h_120w','48h_60w']
hpd_all=[6,12,12,24]
;exp_all=['12h_60w','24h_60w','48h_60w','72h_60w']
;hpd_all=[6,12,24,36]
exp_all=['24h_60w','24h_abs_60w','const_60w']
hpd_all=[12,12,12]
nexp=n_elements(exp_all)

;READ DIMS
  datfil=datdir+exp_all[0]+'/'+'cm1out.nc'
  z=read_nc_var(datfil,'z')*1e3 ; [ m ]
  nz=n_elements(z)
  x=read_nc_var(datfil,'xh') ; [ km ]
  nx=n_elements(x)
  time=read_nc_var(datfil,'xh') ; [ km ]
  nx=n_elements(x)
  time=reform(read_nc_var(datfil,'time'))/3600./24. ; [ d ]
  nt=n_elements(time)
  npd=48
  nd=(nt-1)/npd;10

;QUANTITIES AVERAGED OVER THE ISLAND AND AT SELECTED HEIGHT

  nt_avg=nd_avg*npd
  t_ind_avg = indgen(nt_avg)+nt-1-nt_avg
  x_ind_avg=indgen(nx_island)+nx/2-nx_island/2

  sh_int=fltarr(nt_avg,nexp)
  w_int=sh_int
  eddyh_int=sh_int

for iexp=0,nexp-1 do begin

  exp=exp_all[iexp];'24h_60w';'48h_60w';'12h_60w';'24h_120w';
  print,strupcase(exp)

  ;TIME EXPERIMENT
    time_exp=findgen(npd*20)/(1.*npd)*(12./hpd_all[iexp])

  ;OUTPUT FILE
    datfil=datdir+exp+'/'+'cm1out.nc'

  ;READ VARS
    
    count=[nx,1,nt_avg] & offset=[0,0,t_ind_avg[0]]
    thflux=reform(read_nc_var(datfil,'thflux',count=count,offset=offset)) ; [ K m / s ]

    count=[nx,1,nz,nt_avg] & offset=[0,0,0,t_ind_avg[0]]
    th=reform(read_nc_var(datfil,'th',count=count,offset=offset)) ; [ K ]
    p=reform(read_nc_var(datfil,'prs',count=count,offset=offset)) ; [ Pa ]
    ip = mean(mean(p,dimension=1,/double),dimension=2,/double)
    dp=deriv(ip)
    tmpk=theta(th,p*1e-2,/reverse)+273.15
    w=reform(read_nc_var(datfil,'winterp',count=count,offset=offset)) ; [ m/s ]
    w_int[*,iexp]=mean(reform(w[x_ind_avg,iz_sel,*]),dimension=1,/double)

  ;NEW VARS

    ;SENSIBLE HEAT FLUX
      rd=287. & cp=1005.7
      rho=p/(rd*tmpk)
      sh=cp*reform(rho[*,0,*])*thflux
      sh_int[*,iexp]=mean(sh[x_ind_avg,*],dimension=1,/double)

    ;DRY STATIC ENERGY
      zall=fltarr(nx,nz,nt_avg)
      for iz=0,nz-1 do zall[*,iz,*]=z[iz] ; m
      dse=1d*tmpk*cp + 1d*9.81*zall ; J/kg

    ;BASE-STATE DEVIATIONS
;      spert=dse
;      for it=0,nt_avg-1 do $
;        spert[*,*,it] -= dse[*,*,0]

    ;LOCAL DEVIATIONS
      sprim=dse
      wprim=w
      sxmean=mean(dse,dimension=1,/double)
      wxmean=mean(w,dimension=1,/double)
      for ix=0,nx-1 do begin
        sprim[ix,*,*] -= sxmean
        wprim[ix,*,*] -= wxmean
      endfor
      wpsp=sprim*wprim ; J/kg * m/s

    ;VERTICAL EDDY HEAT CONVERGENCE
      eddyh = fltarr(nz,nt_avg)
      wpspb = mean(reform(wpsp[x_ind_avg,*,*]),dimension=1) & wpsp=0
      z_sel=where(z gt 1e3 and z lt 10e3)
      for it=0,nt_avg-1 do begin
        eddyh[*,it] = -1.*deriv(z, reform(wpspb[*,it])) ; J/kg/s
        ;eddyh_int[*,iexp]=reform(eddyh[iz_sel,*])/cp ; --> K/s
        eddyh_int[it,iexp]=-1.*total(eddyh[z_sel,it]*dp[z_sel])/9.81 ; W/m2
      endfor

endfor ; iexp


;AVERAGE IN TIME
  eddyh = mean(eddyh_int,dimension=1,/double)
  w = mean(w_int,dimension=1,/double)
  sh = max(abs(sh_int),dimension=1)

  ;DAILY-INTEGRATED SH
  sh_tint=fltarr(nexp)
  for iexp=0,nexp-1 do begin
    ish=reform(sh_int[*,iexp])
    loc=where(ish gt 0)
    loc = loc[indgen(hpd_all[iexp]*npd/24)]
;    loc = nt - indgen(hpd_all[iexp]*2*npd/24)
    sh_tint[iexp] = total(ish[loc],/double) ;$
                    ;/ (1d*nd_avg*3600*24)
  endfor


;CREATE FIGURE


;vars=['th','u','w']
;var_str='th';'u';'w';
;nvar=n_elements(vars)

;for ivar=0,nvar-1 do begin

var_str='w_sh';'w_shmax';'w_wpthp';
;vars[ivar]

  figname = figdir+var_str+'_scatter'

  title=var_str;+' ('+strupcase(exp)+')'

  if var_str eq 'w_sh' then begin
    title+=' (SH = abs. max)'
    xvar=sh
    xtitle='SH [ W m!U-2!N ]'
    xrange=[20,150]
    yvar=w*1e2 ; cm/s
    ytitle='w [ cm s!U-1!N ]'
    yrange=[5,35]
  endif else if var_str eq 'w_shmax' then begin
    title+=' (SH = daytime accum)'
    xvar=sh_tint
    xtitle='SH [ W m!U-2!N ]'
    xrange=[0,2500]
    yvar=w*1e2 ; cm/s
    ytitle='w [ cm s!U-1!N ]'
    yrange=[5,25]
  endif else if var_str eq 'w_wpthp' then begin
    ;xvar=eddyh*3600d*24
    ;xvar*=1e2
    xvar=eddyh
    ;xtitle='!8Q!Deddy!N!X [ 10!U-2!N K d!U-1!N ]'
    xtitle='!8Q!Deddy!N!X [ W m!U-2!N ]'
    xrange=[0,15]
    yvar=w*1e2 ; cm/s
    ytitle='w [ cm s!U-1!N ]'
    yrange=[5,25]
  endif

  ;PLOT SETTINGS
    icbar=0 ; color bar?
    ileg=1 ; legend?
    icont=0 ; line contours?
    linethick=2
    csize=0.7
    position=[0.15,0.15,0.97,0.91]
;    xticklen=-0.03 & yticklen=-0.012
    xsize=2.5 & ysize=2.5

;TIME SERIES

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  plot,xvar,yvar,/nodata,position=position,$
    xstyle=9,ystyle=9,color=0,$
    xtick_get=xtickv,ytick_get=ytickv,$
    xticklen=xticklen,yticklen=yticklen,yminor=4,$
    xrange=xrange,yrange=yrange,$;yticks=4,ytickname=strtrim([0,12,0,12,24],2),
;    xrange=[x03,x03+4],yrange=[0,60],$
    xtitle=xtitle,ytitle=ytitle,$
    charsize=csize,title=title

;  loadct,col_table,/silent

  psym=[1,4,5,6]
  psym=psym[0:nexp-1]
  size=0.8
  symsize=replicate(size,4);[2,2,2,2]
  for iexp=0,nexp-1 do $
    plots,xvar[iexp],yvar[iexp],psym=psym[iexp],symsize=symsize[iexp],$
      thick=linethick,color=0,/data

;  loadct,0,/silent

  ;LEGEND
  if ileg then begin
;  if ileg then begin
    csize_fac=0.7
    margin=0.1
    pspacing=2.0 ; length of lines
    spacing=0.8 ; between lines
    leg_str=strupcase(exp_all)
    leg_style=replicate(0,nexp)
    symsize=symsize
    leg_thick=replicate(linethick,nexp)
    leg_color=replicate(0,nexp)
    psym=psym
    legend2,leg_str,psym=psym,symsize=symsize,thick=leg_thick,COLORS=leg_color,$
      charsize=csize*csize_fac,/bottom_legend,/right_legend,/normal,$
      margin=margin,pspacing=pspacing,spacing=spacing,box=0
  endif

  device,/close
  convert_png,figname,res=200,/remove_eps

;endfor ; ivar

print,'Done with all!'

end
