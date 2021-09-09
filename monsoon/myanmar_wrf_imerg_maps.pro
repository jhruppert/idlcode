; 
; Maps of Myanmar WRF output side-by-side with IMERG

; James Ruppert
; 8/14/20
; 
pro myanmar_wrf_imerg_maps

config_dir,dirs=dirs

;EXPERIMENT SETTINGS
  expname='myanmar'
  cases=['ctl','morr_mynn','morr_ysu','morr_acm2']
  cases=['kfcu','kfcu_nudge','raincomp_kfcu']

;idir='9km'
idir='2dom'
if idir eq '2dom' then begin
  domtag='d01'
;  domtag='d02'
endif

;Set this if uneven n-output files
;  nfils_set=6*24+1;10*24+1
casedir=dirs.scdir+expname+'/WRF/experiment/'+idir+'/'
dirs.figdir+='experiments/'+expname+'/'
config_wrfexp, casedir=casedir,cases=cases,dirs=dirs,$
  dims=dims, vars=vars, nfils_set=nfils_set, domtag=domtag;, /verbose

;OB DIRECTORIES
maindir=dirs.scdir+'myanmar/'
imergfil=maindir+'imerg/data/jjas_2013-2017_diurncomp_3B-HHR.MS.MRG.3IMERG.V06B.nc4';'imerg/20150620-20150629.nc4'
era_fil=maindir+'era5/jjas_2013-2017/ERA5-JJAS_2013-2017_mean-pl.nc';'era5/ERA5-20150623-20150630-uv.nc'

;----PLOT OPTIONS--------------------

nd_avg=5 ; If doing daily average, use this many days (tail-end of simulation)

iwind=1 ; Plot 10m wind vectors?

allvars=['RAINNC']
nvsel=n_elements(allvars)

;LEVEL SELECTION
  psel=925;700;800;700 700 is the winner of 600, 700, 800, 850, 925
  izlev=(where(dims.pres eq psel))[0]

;----TIME SPECS--------------------

;FULL TIME SERIES
  time=dims.time
  nt_full=dims.nt
  npd=dims.npd
  nhrs=1.*nt_full*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nhrs)

;----READ OBS--------------------

  ;IMERG RAINFALL
    imerg=read_nc_var(imergfil,'precipitationCal') * 24 ; mm/hr --> mm/d
    imerg=transpose(temporary(imerg))

  ;TIME AVERAGE
    imerg=mean(temporary(imerg),dimension=1,/nan,/double)

;    timeim=read_nc_var(imergfil,'time')
;    timeim=julday(1,1,1970,0,0,temporary(timeim))
;    ;FIND OFFSET (IMERG TIME SERIES STARTS EARLIER)
;      tdiff=abs(timeim-time[0])
;      it_offset=(where(tdiff eq min(tdiff)))[0]

;    specs=size(rain,/dimensions)
;    ntim=specs[0] & nxim=specs[1] & nyim=specs[2]

    ;IMERG LAT/LON
    lonim=read_nc_var(imergfil,'lon')
    latim=read_nc_var(imergfil,'lat')

;    npdim=48

  ;ERA5 WINDS

    ;FROM FILES DOWNLOADED AS NETCDF
;    eralon=read_nc_var(era_fil,'longitude')
;    eralat=read_nc_var(era_fil,'latitude')
;    nxera=n_elements(eralon) & nyera=n_elements(eralat)
;    erahh=read_nc_var(era_fil,'time')
;    eratime=julday(1,1,1900,0+erahh,0,0)
;    ntera=n_elements(eratime)
;    npdera=24 & ndera=ntera/npdera
;    ;FIND OFFSET (IMERG TIME SERIES STARTS EARLIER)
;      tdiff=abs(eratime-time[0])
;      it_offset_era=(where(tdiff eq min(tdiff)))[0]

    ;FROM COMPLETE FILES DOWNLOADED AS GRIB
    eralon=read_nc_var(era_fil,'lon')
    eralat=read_nc_var(era_fil,'lat')
    nxera=n_elements(eralon) & nyera=n_elements(eralat)

  if iwind then begin
    ;PLEVS
;    ilevera=[925, 950, 975, 1000]
    p_era=reform(read_nc_var(era_fil,'plev'))
    izlev_era=(where(p_era*1d-2 eq psel))[0]
    count=[nxera,nyera,1,1] & offset=[0,0,izlev_era,0] ; x, y, p, t
;    uera=reform(read_nc_var(era_fil,'u',count=count,offset=offset))
;    vera=reform(read_nc_var(era_fil,'v',count=count,offset=offset))
    varstr='var131'
    uera=reform(read_nc_var(era_fil,varstr,count=count,offset=offset))
    varstr='var132'
    vera=reform(read_nc_var(era_fil,varstr,count=count,offset=offset))
  endif

;----WRF OUTPUT--------------------

for ic=0,dirs.nc-1 do begin
;for ic=0,0 do begin

  print,'CASE: ',strupcase(dirs.cases[ic])

    i_nt=nt_full
    if ( strmatch(dirs.cases[ic],'*36h*') or strmatch(dirs.cases[ic],'icrf_*') or strmatch(dirs.cases[ic],'lwcr*') $
      or (dirs.cases[ic] eq 'lwswcrf') or (dirs.cases[ic] eq 'axisym') ) then i_nt-=36
    if strmatch(dirs.cases[ic],'*24h*') then i_nt-=24
    if strmatch(dirs.cases[ic],'*48h*') then i_nt-=48
    it_test=indgen(i_nt)+nt_full-i_nt

nd=i_nt/npd

var_str=allvars[0]
print,'VAR: ',var_str

;  if var_str eq 'rainrate' or var_str eq 'RAINNC' then begin
    ;ACCUMULATED RAIN OR RAIN RATE
    iv_str=var_str
    iv=where(vars.vars eq iv_str)
    file=dirs.files_post[ic,iv]

    count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
    rain=reform(read_nc_var(file,iv_str,count=count,offset=offset)) ; Accumulated rain, mm

  ;ADD CUMULUS PARAM RAIN IF USED
;    if strmatch(dirs.cases[ic],'kfcu*') then begin
;      iv_str='RAINC'
;      iv=where(vars.vars eq iv_str)
;      file=dirs.files_post[ic,iv]
;      rain+=reform(read_nc_var(file,iv_str,count=count,offset=offset)) ; Accumulated rain, mm
;      iv_str='RAINSH'
;      iv=where(vars.vars eq iv_str)
;      file=dirs.files_post[ic,iv]
;      rain+=reform(read_nc_var(file,iv_str,count=count,offset=offset)) ; Accumulated rain, mm
;    endif

  ;WIND
    if iwind then begin
    ;U
      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,izlev,0] ; x,y,z,t
;      iv=where(vars.vars eq 'U10')
      iv=where(vars.vars eq 'U')
      file=dirs.files_post[ic,iv]
;      u=reform(read_nc_var(file,'',varid='0',count=count,offset=offset))
      u=reform(read_nc_var(file,'U',count=count,offset=offset))
    ;V
;      iv=where(vars.vars eq 'V10')
      iv=where(vars.vars eq 'V')
      file=dirs.files_post[ic,iv]
;      v=reform(read_nc_var(file,'',varid='0',count=count,offset=offset))
      v=reform(read_nc_var(file,'V',count=count,offset=offset))
;      wind=create_struct('u',u,'v',v)
    endif

;----TIME LOOP--------------------

iu=0 & iv=0
iuera=0 & ivera=0

;for iday=0,nd-1 do begin
;for iday=1,1 do begin
for iday=nd-1,nd-1 do begin

  ;AVERAGE OVER N DAYS
    t_ind1=npd + npd*iday
    t_ind0=t_ind1-nd_avg*npd;0+npd*iday-1*(iday gt 0)
    if t_ind0 lt 0 then t_ind0=0

    irain = reform( (rain[*,*,t_ind1] - rain[*,*,t_ind0]) ) / (1.*nd_avg) ; mm --> mm/d

  ;WRF WIND
    if iwind then begin
      ;t_ind=indgen(npd)+npd*iday
      t_ind=indgen(10*npd)+288-10*npd
    ;U
      iu=mean(u[*,*,t_ind],dimension=3,/nan,/double)
    ;V
      iv=mean(v[*,*,t_ind],dimension=3,/nan,/double)
    endif

  ;IMERG DAILY ACCUMULATED RAINFALL
;    t_ind=indgen(npdim)+iday*npdim + it_offset
    rainim=imerg;mean(imerg[t_ind,*,*],dimension=1,/nan,/double)

  ;ERA5 WIND
    if iwind then begin
;      t_ind=fltarr(npdera)+npdera*iday + it_offset_era
;    ;U
;      iuera=mean(uera[*,*,t_ind],dimension=3,/nan,/double)
;    ;V
;      ivera=mean(vera[*,*,t_ind],dimension=3,/nan,/double)
      iuera=uera
      ivera=vera
    endif

;----CREATE PLOTS--------------------

;  if var_str eq 'RAINNC' then begin
    setmin=0
    setmax=36;50;200
;  endif

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figdir=dirs.figdir+dirs.cases[ic]+'/'+var_str+'/'
  spawn,'mkdir '+figdir,tmp,tmpe
  figspecs=create_struct(figspecs,'figname',' ')
;  figspecs.title+=' ('+strupcase(dirs.cases[ic]);+')'
  titlesav=figspecs.title+' ('+strupcase(dirs.cases[ic]);+')'

  ;TIME INFO
    caldat,time[iday*npd],mm,dd,yy,hh,mn
    form='(i2.2)'
    t_stamp_fig=string(mm,format=form)+string(dd,format=form)
    t_stamp=string(mm,format=form)+'/'+string(dd,format=form)
t_stamp_fig='jjas2013-2017'
t_stamp=t_stamp_fig

    print,t_stamp

    figspecs.title=titlesav+'; '+t_stamp+')'

    figname=figdir+t_stamp_fig+'_'+var_str+'_'+dirs.cases[ic];+'_dayavg'
    if var_str eq 'avor' then figname+='_'+string(psel,format='(i3.3)')
    figspecs.figname=figname

    s_wrf=create_struct('lon',dims.lon,'lat',dims.lat,'var',irain,'u',iu,'v',iv)
    s_obs=create_struct('lon',lonim,'lat',latim,'var',rainim,'lonwind',eralon,'latwind',eralat,'u',iuera,'v',ivera)

    plot_myanmar_obcomp_maps, dirs, figname, s_wrf, s_obs, figspecs

endfor ; icase

;endfor ; ivar

endfor ; iday

print,'DONE!!'
end
