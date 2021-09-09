; 
; Regression analysis of IMERG and ERAi data using Indian Monsoon index.
;
; James Ruppert
; 12/13/20
; 
pro run_imerg_imi_regress

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

  icalc_pw=1
;  if var_plot eq 'pw' then icalc_pw=1 else icalc_pw=0
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
;    time_ini=julday(tmp[1],tmp[2],tmp[0],0,0,0)
    imi_all[id]=tmp[3]  ; Index (stand)
    imi_all2[id]=tmp[4] ; Index (non-stand, mm/d)
  endfor

  ;STANDARDIZE

  std_all=stddev(imi_all,/nan,/double)
  std_all2=stddev(imi_all2,/nan,/double)

  print,'Standard Deviation:'
  print,'All:',std_all
  print,'All:',std_all2

  imi_all  = (imi_all -  mean(imi_all,/nan,/double))  / std_all
  imi_all2 = (imi_all2 - mean(imi_all2,/nan,/double)) / std_all2


;----CALCULATE REGRESSIONS--------------------

ido_reg=0
if ido_reg then begin

  ;REGRESS RAIN
  rain_reg=fltarr(nx,ny)
  rain_reg2=rain_reg
  for ix=0,nx-1 do begin
  for iy=0,ny-1 do begin
    rain_reg[ix,iy]  = (regress(imi_all, reform(rain[ix,iy,*])))[0]
    rain_reg2[ix,iy] = (regress(imi_all2,reform(rain[ix,iy,*])))[0]
  endfor
  endfor

  ;REGRESS ERAi
  u_reg=fltarr(nxera,nyera) & v_reg=u_reg
  u_reg2=u_reg
  v_reg=u_reg & v_reg2=u_reg
  pw_reg=u_reg
  pw_reg2=u_reg
  for ix=0,nxera-1 do begin
  for iy=0,nyera-1 do begin
    u_reg[ix,iy]=(regress(imi_all,reform(u[ix,iy,*])))[0]
    u_reg2[ix,iy]=(regress(imi_all2,reform(u[ix,iy,*])))[0]
    v_reg[ix,iy]=(regress(imi_all,reform(v[ix,iy,*])))[0]
    v_reg2[ix,iy]=(regress(imi_all2,reform(v[ix,iy,*])))[0]
    pw_reg[ix,iy]=(regress(imi_all,reform(pw[ix,iy,*])))[0]
    pw_reg2[ix,iy]=(regress(imi_all2,reform(pw[ix,iy,*])))[0]
  endfor
  endfor

  rain_reg*=std_all
  rain_reg2*=std_all2

  u_reg*=std_all
  u_reg2*=std_all2
  v_reg*=std_all
  v_reg2*=std_all2

  pw_reg*=std_all
  pw_reg2*=std_all2


;----REGRESSION PLOTS--------------------

;iplot_reg=0
;if iplot_reg then begin

  if var_plot eq 'rain' then begin
    var_str='RAINNC'
    setmax=20 & setmin='0.'
    cbform='(i3)'
    cbtag='Rain [ mm d!U-1!N ]'
    lon=lonim
    lat=latim
  endif else if var_plot eq 'pw' then begin
    var_str=var_plot
    setmax=15 & setmin=-15
    cbform='(i3)'
    cbtag='PW [ mm ]'
    lon=eralon
    lat=eralat
  endif

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, setndivs=4, set_cint=3
  figdir=dirs.figdir+'imerg/'
  figspecs=create_struct(figspecs,'figname',' ')
  figspecs.cbar_format=cbform
  figspecs.cbar_tag=cbtag

for itest=0,1 do begin

  if itest eq 0 then begin
    if var_plot eq 'rain' then $
      var_plt=rain_reg $
    else if var_plot eq 'pw' then $
      var_plt=pw_reg
    u_plt=u_reg
    v_plt=v_reg
    itag='Ini1'
  endif else if itest eq 1 then begin
    if var_plot eq 'rain' then $
      var_plt=rain_reg2 $
    else if var_plot eq 'pw' then $
      var_plt=pw_reg2
    u_plt=u_reg2
    v_plt=v_reg2
    itag='Ini2'
  endif

  wind=create_struct('u',u_plt,'v',v_plt,'x',eralon,'y',eralat)
  wspd=sqrt(u_plt^2+v_plt^2);*10.
  cvar=create_struct('cvar',wspd,'x',eralon,'y',eralat)

  figspecs.title='Regression: "'+itag+'"'
  figspecs.figname=figdir+'imerg_reg_'+var_plot+'_'+strlowcase(itag)+'_'+strtrim(psel_era,2)

  wrf_myanmar_map_plot, dirs, var_plt, lon, lat, figspecs, wind=wind, cvar=cvar

endfor
;endif

endif ; All regression analysis (ido_reg)


;----THRESHOLD AVERAGES--------------------

ido_thresh=1
if ido_thresh then begin

  thresh=1.5

  it_all=where(imi_all ge thresh,count)
  print,'N-all: ',count
  it_allm=where(imi_all le -1.*thresh,count)
  print,'N-all-neg: ',count
  it_all2=where(imi_all2 ge thresh,count)
  print,'N-all2: ',count
  it_allm2=where(imi_all2 le -1.*thresh,count)
  print,'N-all2-neg: ',count

  if var_plot eq 'rain' then begin
    var_str='RAINNC'
    var=rain
    setmax=60 & setmin='0.'
    cbform='(i3)'
    lon=lonim
    lat=latim
  endif else if var_plot eq 'pw' then begin
    var_str=var_plot
    var=pw
    setmax=70 & setmin=30
    cbform='(i2)'
    lon=eralon
    lat=eralat
  endif

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, setndivs=6, set_cint=3
  figdir=dirs.figdir+'imerg/'
  figspecs=create_struct(figspecs,'figname',' ')
  figspecs.cbar_format=cbform

for itest=0,3 do begin

  if itest eq 0 then begin
    it_sel=it_all
    itag='Ini1p'
  endif else if itest eq 1 then begin
    it_sel=it_allm
    itag='Ini1m'
  endif else if itest eq 2 then begin
    it_sel=it_all2
    itag='Ini2p'
  endif else if itest eq 3 then begin
    it_sel=it_allm2
    itag='Ini2m'
  endif

  var_plt=mean(var[*,*,it_sel],dimension=3,/nan,/double)
  u_plt=mean(u[*,*,it_sel],dimension=3,/nan,/double)
  v_plt=mean(v[*,*,it_sel],dimension=3,/nan,/double)

  wind=create_struct('u',u_plt,'v',v_plt,'x',eralon,'y',eralat)
  wspd=sqrt(u_plt^2+v_plt^2)
  cvar=create_struct('cvar',wspd,'x',eralon,'y',eralat)

  figspecs.title='Average: "'+itag+'"'
  figspecs.figname=figdir+'imerg_thresh_'+var_plot+'_'+strlowcase(itag)+'_'+strtrim(psel_era,2)

  wrf_myanmar_map_plot, dirs, var_plt, lon, lat, figspecs, wind=wind, cvar=cvar

endfor
  

endif ; ido_thresh

print,'DONE!!'
end
