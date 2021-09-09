; 
; Create a cross section from CM1 netCDF output.
; 
; James Ruppert
; 05.12.17
; 
pro cm1_cross

;SOME SETTINGS

exp_all=['24h_60w','48h_60w','12h_60w','24h_120w']
nexp=n_elements(exp_all)

;for iexp=0,nexp-1 do begin
for iexp=0,0 do begin

  exp=exp_all[iexp];'24h_60w';'48h_60w';'12h_60w';'24h_120w';
  print,strupcase(exp)

  iflux=0 ; just plot flux time series?

  home='/home/mpim/m300462/'
;  datdir=home+'cm1r19/run/'
  datdir=home+'cm1r19/output/'+exp+'/'
  figdir=home+'misthome/IDLWorkspace85/Default/figures/monsoon/cm1/'+exp+'/'

  ;READ TIME STATS
;    nd=10 & npd=48
;    nt=48*nd

  ;OUTPUT FILE
;    datfil=datdir+'cm1out_davg.nc'
    datfil=datdir+'cm1out.nc'

  ;READ DIMS
    z=read_nc_var(datfil,'z') ; [ km ]
    nz=n_elements(z)
    x=read_nc_var(datfil,'xh') ; [ km ]
    nx=n_elements(x)

  ;READ VARS
    th=reform(read_nc_var(datfil,'th')) ; [ K ]
    p=reform(read_nc_var(datfil,'prs')) ; [ Pa ]
    if iflux then $
      thflux=reform(read_nc_var(datfil,'thflux')) $ ; [ K m / s ]
    else begin
      u=reform(read_nc_var(datfil,'uinterp')) ; [ m/s ]
      w=reform(read_nc_var(datfil,'winterp')) ; [ m/s ]
    endelse

  ;TIME
    time=reform(read_nc_var(datfil,'time'))/3600./24. ; [ d ]
    nt=n_elements(time)
    npd=48
    nd=(nt-1)/npd;10
    ;nt=48*nd
    ;if nt gt n_elements(time) then stop
;    specs=size(th)
;    nt2=specs[3]
;    time2=indgen(nt2)

  ;FLUX TIME SERIES
    rd=287. & cp=1005.7
    tmpk=theta(th,p*1e-2,/reverse)+273.15
    if iflux then begin
      rho=p/(rd*tmpk)
      sh=cp*reform(rho[*,0,*])*thflux
      figname = figdir+'flux_tser_'+exp
      cm1_island_flux_tser,figname,sh,x,time,$
        exp=exp
      return
    endif

  ;DRY STATIC ENERGY
    zall=fltarr(nx,nz,nt)
    for iz=0,nz-1 do zall[*,iz,*]=z[iz]*1e3 ; m
    dse=1d*tmpk*cp + 1d*9.81*zall ; J/kg

  ;BASE-STATE DEVIATION
    thp=th
    for it=0,nt-1 do $
      thp[*,*,it] -= th[*,*,0]

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
    nx_island=100
    x_ind_avg=indgen(nx_island)+nx/2-nx_island/2
    eddyh = fltarr(nz,nt)
    wpspb = mean(reform(wpsp[*,*,*]),dimension=1) & wpsp=0
    for it=0,nt-1 do $
      eddyh[*,it] = -1.*deriv(z*1e3, reform(wpspb[*,it])) ; J/kg/s
    eddyh*=1d/cp ; --> K/s

  ;DAILY AVERAGES
    nd/=2
    npd*=2
    wd=fltarr(nx,nz,nd)
    ud=wd
    thd=wd
    for id=0,nd-1 do begin
      t_ind=indgen(npd)+id*npd
      wd[*,*,id]=mean(w[*,*,t_ind],dimension=3,/double)
      ud[*,*,id]=mean(u[*,*,t_ind],dimension=3,/double)
      thd[*,*,id]=mean(thp[*,*,t_ind],dimension=3,/double)
    endfor


;CREATE FIGURE


vars=['th','u','w','eddy']
;var_str='th';'u';'w';
nvar=n_elements(vars)

;for ivar=0,nvar-1 do begin
for ivar=3,3 do begin

var_str=vars[ivar]

;for id=0,nd-1 do begin
for id=nd-1,nd-1 do begin

  figname = figdir+var_str;+'_'+strtrim(id,2)
;  figname = figdir+var_str+'_prof_'+strtrim(id,2)
;  figname = figdir+var_str;+'_'+strtrim(id,2)
  title=var_str+' ('+strupcase(exp)+')'

  if var_str eq 'w' then begin
;    iz=29 ; ~5 km
;    var=reform(wd[*,iz,*])
    var=reform(wd[*,*,id])
;    var=smooth(var,[15,0],/edge_mirror )
    maxv=3.0
    cbar_tag='[ m s!U-1!N ]'
  endif else if var_str eq 'u' then begin
    var=reform(ud[*,*,id])
    maxv=6.0
    cbar_tag='[ m s!U-1!N ]'
  endif else if var_str eq 'th' then begin
    var=reform(thd[*,*,id])
;    x=mean(th[*,*,0],dimension=1,/double)
;    x2=mean(reform(thd[*,*,nd-1]),dimension=1,/double)
;    var=reform(th[*,*,id*npd+npd/4.]) - reform(th[*,*,0])
    maxv=1.5
    cbar_tag='[ K ]'
  endif else if var_str eq 'eddy' then begin
    var=eddyh*1d2*3600*24
;    var=gauss_smooth(var,[0,3])
    maxv=100
    cbar_tag='[ 10!U-2!N K d!U-1!N ]'
  endif


  ;SHADING AND CONTOURS

    col_table=70;13

    levels=findgen(255)/254*2*maxv-maxv
    colors=indgen(255)
    colors=reverse(colors)

    cint=0.5 & nlevs=30
    clevels=(findgen(nlevs)+1)*cint

  ;PLOT SETTINGS
    icbar=1 ; color bar?
    ileg=0 ; legend?
    icont=0 ; line contours?
    linethick=3
    csize=0.7
    position=[0.12,0.17,0.96,0.91]
;    position=[0.15,0.17,0.93,0.91]
    xticklen=-0.03 & yticklen=-0.012
    xsize=4.5 & ysize=2.5
;    xsize=2.5
    xtitle='x [ km ]'
;    xtitle='Th [ K ]'
    xtitle='Time [ d ]'
    ytitle='z [ km ]'

  ;AXES
    if ~keyword_set(xrange) then $
      xrange=[min(x),max(x)]
;    xrange=[299.0,301.5]
    x=time
    xrange=[0,4];20]
    ;y=time2
    y=z
    if ~keyword_set(yrange) then $
      yrange=[min(y),max(y)]
      yrange=[0,12]

;TIME SERIES

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  plot,x,y,/nodata,position=position,$
    xstyle=9,ystyle=9,color=0,$
    xtick_get=xtickv,ytick_get=ytickv,$
    xticklen=xticklen,yticklen=yticklen,yminor=4,$
    xrange=xrange,yrange=yrange,$;yticks=4,ytickname=strtrim([0,12,0,12,24],2),
;    xrange=[x03,x03+4],yrange=[0,60],$
    xtitle=xtitle,ytitle=ytitle,$
    charsize=csize,title=title

  loadct,col_table,/silent

  ;FILL SHADING
;  for i=0,1 do $
    contour,transpose(var),x,y,/cell_fill,/overplot,$
      levels=levels,c_colors=colors

;  oplot,x,y,linestyle=2,thick=2.0,color=0
;  oplot,x2,y,linestyle=2,thick=2.0,color=220

  loadct,0,/silent

  ;CONTOURS
  if icont then $
    contour,pdf,x,y,/follow,/overplot,$
      levels=clevels,c_charsize=csize*0.7

  ;OVERPLOT AXIS LINES
  for i=0,1 do $
    plots,!x.crange[[i,i]],!y.crange
  plots,!x.crange,!y.crange[[0,0]]

  ;PRECIP RATE AND AREA FRAC
  plot_agg_size=0
  if plot_agg_size then begin

;    axis,yaxis=1,ystyle=1,charsize=csize,yticklen=yticklen,$
;      yrange=[0,10],ytitle='Rainfall [ mm d!U-1!N ]',/save,/data
;    oplot,times_ave,precip,linestyle=2,thick=2.0
;
;    ;DIURNAL AMPLITUDE
;    if diurn_pcp_amp then $
;      oplot,time_days,d_amp,color=0,linestyle=3,thick=thick_2

    xloc=103
;    axis,xloc,yaxis=1,ystyle=1,charsize=csize,yticklen=yticklen,$
;      yrange=[0,100],ytitle='Ascent Area Fraction [ % ]',/save,/data
    axis,yaxis=1,ystyle=1,charsize=csize,yticklen=yticklen,$
      yrange=[0,100],ytitle='Ascent Area Fraction [ % ]',/save,/data
    oplot,times_ave,frac_up,linestyle=2,thick=2.0

  endif

  ;COLOR BAR
  if icbar then begin
    loadct,col_table,/silent
;    cpos= [ position[2]+0.03 ,$
;    cpos= [ 0.91 ,$
;      position[1] ,$
;;      position[2]+0.04 ,$
;      0.922 ,$
;      position[3] ]
;    if plot_agg_size then $
;      cpos[[0,2]]+=0.02
    cpos= [ position[0]+0.03 , 0.54 , position[0]+0.03+0.01 , position[3]-0.04 ]
    colorbar2, colors=colors, range=[min(levels),max(levels)], divisions=5, $
      charsize=csize*0.7, position=cpos, /right, /vertical, title=cbar_tag, $
      annotatecolor='black',format='(f5.1)'
    loadct,0,/silent
  endif

  ;LEGEND
  if ileg then begin
;  if ileg then begin
    csize_fac=0.7
    margin=0.1
    pspacing=2.0 ; length of lines
    spacing=0.8 ; between lines
    if plot_agg_size then begin
      ;leg_str=['Mean','Standard Dev.','Area Frac.']
      leg_str=['Mean','Precip','Area Frac.']
      leg_str=['Mean','Area Frac.']
      leg_style=[0,2]
      leg_thick=[2,2]
      leg_color=[0,0]
    endif else begin
      leg_str=['Mean','Standard Dev.']
      leg_style=[0,2]
      leg_thick=[2,2]
      leg_color=[0,0]
    endelse
    legend2,leg_str,linestyle=leg_style,thick=leg_thick,COLORS=leg_color,$
      charsize=csize*csize_fac,/top_legend,/right_legend,/normal,$
      margin=margin,pspacing=pspacing,spacing=spacing,box=0
  endif

  device,/close
  convert_png,figname,res=200,/remove_eps

endfor ; id

endfor ; ivar

endfor ; iexp

print,'Done with all!'

end
