; 
; Time series of IMERG rainfall over the northern BoB region, along with other monsoon indices.
;
; James Ruppert
; 12/13/20
; 
pro run_imerg_tser

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

var_plot='rain';
;var_plot='pw';'rain';

;ERAi SETTINGS
;LEVEL SELECTION
  psel_era=925;500;700;925

;USE FULL JJAS 2013-2017
  yy_plot=[2013,2017]
  mm_plot=[6,9]
  dd_plot=[1,30] ; inclusive

  ;DATE STRING
    form2='(i2.2)'
    form4='(i4)'
    dat_str=string(mm_plot[0],format=form2)+string(dd_plot[0],format=form2)+strmid(strtrim(yy_plot[0],2),2,2)+'-'+$
            string(mm_plot[1],format=form2)+string(dd_plot[1],format=form2)+strmid(strtrim(yy_plot[1],2),2,2)

;USE WRF MET_EM DATA FOR LAND/SEA MASK
expname='myanmar'
cases=['ctl']
idir='2dom'
if idir eq '2dom' then $
  domtag='d01'
casedir=dirs.scdir+expname+'/WRF/experiment/'+idir+'/'
dirs.figdir+=expname+'/'
config_wrfexp, casedir=casedir,cases=cases,dirs=dirs,$
  dims=dims, vars=vars, nfils_set=nfils_set, domtag=domtag;, /verbose


;----OB DIRECTORIES--------------------

  maindir=dirs.scdir+'myanmar/'
  im_fil=maindir+'imerg/data/jjas_2013-2017_daymean_3B-HHR.MS.MRG.3IMERG.V06B.nc4'
;  npd_imerg=48

  era_dir=maindir+'erai/jjas_2013-2017/';ERAi-JJAS13-17-pl.nc4'
  era_fil=era_dir+'ERAi-JJAS13-17-pl_dayavg.nc4'
  era_sfil=era_dir+'ERAi-JJAS13-17-sl_dayavg.nc4'
;  npd_era=4

  irain_dir=dirs.figdir+'imerg/'

;Local SOLAR time conversion
;local=6;round(mean(dims.lon)/360.*24.) ; deg lon --> hours
;print,'Adding +'+strtrim(local,2)+' for LT'


;----ONE TIME ARRAY FOR DAILY AVERAGE TIME SERIES--------------------

  ;ACCOMMODATE MULTIPLE YEARS
    if yy_plot[0] eq yy_plot[1] then $
      time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
        final=julday(mm_plot[1],dd_plot[1],yy_plot[1],23,59,59),step_size=1,units='Day') $
    else begin
      time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
        final=julday(mm_plot[1],dd_plot[1],yy_plot[0],23,59,59),step_size=1,units='Day')
      for iyy=yy_plot[0]+1,yy_plot[1] do begin
        itime=timegen(start=julday(mm_plot[0],dd_plot[0],iyy,0,0,0),$
          final=julday(mm_plot[1],dd_plot[1],iyy,23,59,59),step_size=1,units='Day')
        time=[time,itime]
      endfor
    endelse
    nt=n_elements(time)
    nd=nt


;----READ OBS--------------------

  ;IMERG RAINFALL

    rain=reform(read_nc_var(im_fil,'precipitationCal'))*24. ; change to mm/d
    rain=transpose(temporary(rain),[1,0,2])

  ;LAT/LON
    lonim=read_nc_var(im_fil,'lon')
    latim=read_nc_var(im_fil,'lat')
    nx=n_elements(lonim)
    ny=n_elements(latim)

;;READ TOPOGRAPHY AND LAND MASK FROM ERA5 MET_EM FILES
;  spawn,'ls '+era_dir+'met_em.'+domtag+'*',era_fil
;  fil=era_fil[0]
;  topo=reform(read_nc_var(fil,'HGT_M'))*1e3 ; m --> km
;  lsmask=reform(read_nc_var(fil,'LANDMASK'))


;----READ ERA--------------------

  ;DIMENSIONS
  eralon=read_nc_var(era_fil,'lon')
  eralat=read_nc_var(era_fil,'lat')
  nxera=n_elements(eralon)
  nyera=n_elements(eralat)

  ;PLEVS
  p_era=reform(read_nc_var(era_fil,'plev'))*1d-2 ; Pa --> hPa
  nzera=n_elements(p_era)

  ;SELECTED LEVEL
  izlev_era=(where(p_era eq psel_era,count))[0]
  if count eq 0 then stop

  ;WINDS

    count=[nxera,nyera,1,nt] & offset=[0,0,izlev_era,0] ; x, y, p, t
    u=reform(read_nc_var(era_fil,'var131',count=count,offset=offset))
    v=reform(read_nc_var(era_fil,'var132',count=count,offset=offset))

  icalc_pw=0
  if icalc_pw then begin

    qv=reform(read_nc_var(era_fil,'var133')) ; kg/kg
    pr=fltarr(nxera,nyera,nzera,nt) & for iz=0,nzera-1 do pr[*,*,iz,*]=p_era[iz]*1d2 ; --> Pa
    ;REVERSE VERTICAL
      qv=reverse(temporary(qv),3)
      pr=reverse(temporary(pr),3)

    ;FILL IN SURFACE VALUES
      prsfc=reform(read_nc_var(era_sfil,'var134')) ; Pa
      pr[*,*,0,*]=prsfc

    pw=fltarr(nxera,nyera,nt)
    ipsel=indgen(nzera-1)+1
    for it=0,nt-1 do begin
    for ix=0,nxera-1 do begin
    for iy=0,nyera-1 do begin
      ip=where(pr[ix,iy,ipsel,it] lt pr[ix,iy,0,it])
      ip=[0,ip+1]
      dp=deriv(reform(pr[ix,iy,ip,it]))*(-1.)
      pw[ix,iy,it]=total(reform(qv[ix,iy,ip,it])*dp,/double)
    endfor
    endfor
    endfor
    pw/=9.81

  endif


;----MEAN MAPS--------------------

iplot_mn_rain=0
if iplot_mn_rain then begin

if var_plot eq 'rain' then begin
  var_str='RAINNC'
  setmax=30 & setmin='0.'
  cbform='(i2)'
  var_mn=mean(rain,dimension=3,/nan,/double)
  lon=lonim
  lat=latim
endif else if var_plot eq 'pw' then begin
  var_str=var_plot
  setmax=70 & setmin=30
  cbform='(i2)'
  var_mn=mean(pw,dimension=3,/nan,/double)
  lon=eralon
  lat=eralat
endif

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=3
  figdir=dirs.figdir+'imerg/'
  figspecs=create_struct(figspecs,'figname',figdir+'imerg_mean_'+var_plot+'_'+dat_str)
  figspecs.cbar_format=cbform

  u_mn=mean(u,dimension=3,/nan,/double)
  v_mn=mean(v,dimension=3,/nan,/double)
  wind=create_struct('u',u_mn,'v',v_mn,'x',eralon,'y',eralat)
  wspd=sqrt(u_mn^2+v_mn^2)
  cvar=create_struct('cvar',wspd,'x',eralon,'y',eralat)

  wrf_myanmar_map_plot, dirs, var_mn, lon, lat, figspecs, wind=wind, cvar=cvar

endif


;----READ IN RAINFALL TIME SERIES--------------------

  irain_all=fltarr(nd)
  irain_coast=irain_all
  irain_offshore=irain_all
  irain_onshore=irain_all
  irain_smd=irain_all
  openr,1,irain_dir+'rain_all_'+dat_str+'.txt'
  readf,1,irain_all & close,1
  openr,1,irain_dir+'rain_coast_'+dat_str+'.txt'
  readf,1,irain_coast & close,1
  openr,1,irain_dir+'rain_offshore_'+dat_str+'.txt'
  readf,1,irain_offshore & close,1
  openr,1,irain_dir+'rain_onshore_'+dat_str+'.txt'
  readf,1,irain_onshore & close,1
  openr,1,irain_dir+'rain_smd_'+dat_str+'.txt'
  readf,1,irain_smd & close,1

  ;OFFSHORE - ONSHORE; COAST
  irain_cdiff = irain_offshore - irain_onshore

  ;STANDARDIZE

  std_all=stddev(irain_all,/nan,/double)
  std_coast=stddev(irain_coast,/nan,/double)
  std_offshore=stddev(irain_offshore,/nan,/double)
  std_onshore=stddev(irain_onshore,/nan,/double)
  std_smd=stddev(irain_smd,/nan,/double)
  std_cdiff=stddev(irain_cdiff,/nan,/double)

;  print,'Standard Deviations:'
;  print,'All:',std_all
;  print,'Coast:',std_coast
;  print,'Offshore:',std_offshore
;  print,'Onshore:',std_onshore

  irain_all = (irain_all - mean(irain_all,/nan,/double)) / std_all
  irain_coast = (irain_coast - mean(irain_coast,/nan,/double)) / std_coast
  irain_offshore = (irain_offshore - mean(irain_offshore,/nan,/double)) / std_offshore
  irain_onshore = (irain_onshore - mean(irain_onshore,/nan,/double)) / std_onshore
  irain_smd = (irain_smd - mean(irain_smd,/nan,/double)) / std_smd
  irain_cdiff = (irain_cdiff - mean(irain_cdiff,/nan,/double)) / std_cdiff


;----READ IN IMI TIME SERIES--------------------

  imi_fil=dirs.home+'idl/code/monsoon/ind_monsoon_index_full.ascii'

  grep='grep "2013    '+strtrim(6,2)+'" '+imi_fil
  spawn,grep,lines
  for im=7,9 do begin
    grep='grep "2013    '+strtrim(im,2)+'" '+imi_fil
    spawn,grep,tmp
    lines=[lines,tmp]
  endfor
  for n=4,7 do begin
  for im=6,9 do begin
    grep='grep "201'+strtrim(n,2)+'    '+strtrim(im,2)+'" '+imi_fil
    spawn,grep,tmp
    lines=[lines,tmp]
  endfor
  endfor
  nlin=n_elements(lines)

  imi_all=fltarr(nd)
  imi_all2=imi_all
  for id=0,nd-1 do begin
    tmp1=strsplit(lines[id],/extract)
    tmp=float(tmp1[0:4])
    imi_all[id]=tmp[3]  ; Index (stand)
    imi_all2[id]=tmp[4] ; Index (non-stand, mm/d)
  endfor

  ;STANDARDIZE

  std_all=stddev(imi_all,/nan,/double)
  std_all2=stddev(imi_all2,/nan,/double)

;  print,'Standard Deviation:'
;  print,'All:',std_all
;  print,'All:',std_all2

  imi_all  = (imi_all -  mean(imi_all,/nan,/double))  / std_all
  imi_all2 = (imi_all2 - mean(imi_all2,/nan,/double)) / std_all2

;----;SMOOTHING--------------------------

ism=3;5

irain_all=smooth(temporary(irain_all),ism,/edge_truncate)
irain_coast=smooth(temporary(irain_coast),ism,/edge_truncate)
irain_offshore=smooth(temporary(irain_offshore),ism,/edge_truncate)
irain_onshore=smooth(temporary(irain_onshore),ism,/edge_truncate)
irain_smd=smooth(temporary(irain_smd),ism,/edge_truncate)
irain_cdiff=smooth(temporary(irain_cdiff),ism,/edge_truncate)

imi_all=smooth(temporary(imi_all),ism,/edge_truncate)
imi_all2=smooth(temporary(imi_all2),ism,/edge_truncate)


;----PLOT TIME SERIES--------------------

for iy=2013,2017 do begin

  figtag=strtrim(iy,2)
  caldat,time,mm,dd,yy,nn
  itsel=where(yy eq iy)

;------

  ;PLOT SPECS
    csize=0.75
    position=[0.10,0.15,0.96,0.95]
    xsize=6.0 & ysize=2.4
    xtitle='Date ('+strtrim(iy,2)+')'
    ytitle='Sigma'

  var_str='RAINNC'
  myan_figspecs, var_str, figspecs
  figdir=dirs.figdir+'imerg/'
  figspecs=create_struct(figspecs,'figname',' ')

  figspecs.title=' '
  figspecs.figname=figdir+'imerg_tser_'+figtag;strtrim(psel_era,2)

  ;AXES
    x=time[itsel]
    xrange=[min(x),max(x)]
    y=indgen(5)
    yrange=[-3,3]
    if ~keyword_set(xrange) then $
      xrange=[min(x),max(x)]
    if ~keyword_set(yrange) then $
      yrange=[min(y),max(y)]

  labels=label_date(date_format=['%n/%d'])
;  labeldays=timegen(start=x[0],final=max(x),step_size=10,units='days')
;  xticks=n_elements(labeldays)-1

  set_plot,'ps'
  epsname=figspecs.figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  plot,x,y,/nodata,position=position,$
    xstyle=9,ystyle=9,$
    ;xticks=xticks,xtickv=labeldays,
    xtickunits='time',xtickformat='label_date',$
    yticklen=-0.022,xticklen=-0.018,$
    xrange=xrange,yrange=yrange,$;xminor=5,$
    xtitle=xtitle,ytitle=ytitle,$
    charsize=csize;,$
    ;title=figspecs.title

  plots,!x.crange,[0,0]

  loadct,39,/silent

  leg_thick=[1,1,1];,1]
  leg_style=[0,0,1]
  leg_color=[230,80,0];,230]
  leg_str=['Coast','Offshore','Offshore-Onshore']

  ip=0
  oplot,x,irain_coast[itsel],linestyle=leg_style[ip],thick=leg_thick[ip],color=leg_color[ip]
  ip=1
  oplot,x,irain_offshore[itsel],linestyle=leg_style[ip],thick=leg_thick[ip],color=leg_color[ip]
;  ip=2
;  oplot,x,imi_all2[itsel],linestyle=leg_style[ip],thick=leg_thick[ip],color=leg_color[ip]
  ip=2
  oplot,x,irain_cdiff[itsel],linestyle=leg_style[ip],thick=leg_thick[ip],color=leg_color[ip]

  ;LEGEND
  ileg=1
  if ileg then begin
;  if ileg then begin
    csize_fac=0.7
    margin=0.1
    pspacing=2.3 ; length of lines
    spacing=0.8 ; between lines
    legend2,leg_str,linestyle=leg_style,thick=leg_thick,COLORS=leg_color,$
      charsize=csize*csize_fac,/top_legend,/left_legend,/normal,$
      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.16,0.85]
;      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.56,0.35]
;  loadct,0,/silent
  endif

  device,/close

  convert_png,figspecs.figname,res=200,/remove_eps

endfor ; iy


print,'DONE!!'
end
