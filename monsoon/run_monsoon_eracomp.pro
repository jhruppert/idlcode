; 
; Maps of Myanmar WRF output side-by-side with climatological mean ERA5
;
; For mean rainfall use myanmar_wrf_imerg_maps.pro
;
; James Ruppert
; 10/17/20
; 
pro run_monsoon_eracomp

config_dir,dirs=dirs

;EXPERIMENT SETTINGS
  expname='myanmar'
  cases=['ctl','morr_mynn','morr_ysu','morr_acm2']
  cases=['nudge_xpbl','morr_acm2']

;idir='9km'
idir='2dom'
if idir eq '2dom' then begin
  domtag='d01'
;  domtag='d02'
endif

;Set this if uneven n-output files
;  nfils_set=240-24
casedir=dirs.scdir+expname+'/WRF/experiment/'+idir+'/'
dirs.figdir+=expname+'/'
config_wrfexp, casedir=casedir,cases=cases,dirs=dirs,$
  dims=dims, vars=vars, nfils_set=nfils_set, domtag=domtag;, /verbose

;OB DIRECTORIES
maindir=dirs.scdir+'myanmar/'
imergfil=maindir+'imerg/data/jjas_2013-2017_diurncomp_3B-HHR.MS.MRG.3IMERG.V06B.nc4';'imerg/20150620-20150629.nc4'
;era_fil=maindir+'era5/jjas_2013-2017/ERA5-JJAS_2013-2017_mean-pl.nc';'era5/ERA5-20150623-20150630-uv.nc'
era_dir=dirs.wkdir+'wrf_wps/wrfv4/output_myanmar/JJAS_2013-2017/era5/'
;era_dir=dirs.wkdir+'wrf_wps/wrfv4/output_myanmar/JJAS_2013-2017/erai/'
era_fil=era_dir+'met_em.d01.2013-07-20_avg.nc'

;----PLOT OPTIONS--------------------

nd_avg=5 ; Average over this many days (tail-end of simulation)

iwind=1 ; Plot 10m wind vectors?

allvars=['wspd','RAINNC','rainrate','HFX','avor','slp','lh','th_e','olr','pw']
;allvars=['avor']
allvars=['wspd']
;allvars=['HFX']
;allvars=['lh']
allvars=['pw']
nvsel=n_elements(allvars)
nvsel=n_elements(allvars)

;LEVEL SELECTION
  psel=925;700;800;700 700 is the winner of 600, 700, 800, 850, 925
  izlev=(where(dims.pres eq psel))[0]

;IVAR LOOP
;for ivar_sel=0,nvsel-1 do begin
for ivar_sel=0,0 do begin

var_str=allvars[ivar_sel]
print,'VAR: ',var_str

;----TIME SPECS--------------------

;FULL TIME SERIES
  time=dims.time
  nt_full=dims.nt
  npd=dims.npd
  nhrs=1.*nt_full*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nhrs)

;----READ OBS--------------------

;  ;IMERG RAINFALL
;    imerg=read_nc_var(imergfil,'precipitationCal') * 24 ; mm/hr --> mm/d
;    imerg=transpose(temporary(imerg))
;
;  ;TIME AVERAGE
;    imerg=mean(temporary(imerg),dimension=1,/nan,/double)
;
;;    timeim=read_nc_var(imergfil,'time')
;;    timeim=julday(1,1,1970,0,0,temporary(timeim))
;;    ;FIND OFFSET (IMERG TIME SERIES STARTS EARLIER)
;;      tdiff=abs(timeim-time[0])
;;      it_offset=(where(tdiff eq min(tdiff)))[0]
;
;;    specs=size(rain,/dimensions)
;;    ntim=specs[0] & nxim=specs[1] & nyim=specs[2]
;
;    ;IMERG LAT/LON
;    lonim=read_nc_var(imergfil,'lon')
;    latim=read_nc_var(imergfil,'lat')

;    npdim=48

  ;ERA5 (NOW USING MET_EM FILES)

  eralon=dims.lon
  eralat=dims.lat

  ;PLEVS
  ilevera=925;[925, 950, 975, 1000]
  nz_era=38
  p_era=reform(read_nc_var(era_fil,'PRES',count=[1,1,nz_era,1],offset=[0,0,0,0])) ; Pa
  izlev_era=(where(p_era*1d-2 eq psel))[0]

  if var_str eq 'pw' then begin

    count=[dims.nx,dims.ny,nz_era,1] & offset=[0,0,0,0] ; x,y,z,t
    rh=reform(read_nc_var(era_fil,'RH',count=count,offset=offset))
    tmpk=reform(read_nc_var(era_fil,'TT',count=count,offset=offset))
    pr=reform(read_nc_var(era_fil,'PRES',count=count,offset=offset))
    iice=where(tmpk lt 273.,complement=noice)
    qv=tmpk & qv[*]=0.
    qv[iice]=rh2mixr(rh[iice],pr[iice]*1e-2,tmpk[iice],/ice)
    qv[noice]=rh2mixr(rh[noice],pr[noice]*1e-2,tmpk[noice],/ice)
    var_era=fltarr(dims.nx,dims.ny)
    ipsel=indgen(nz_era-1)+1
    for ix=0,dims.nx-1 do begin
    for iy=0,dims.ny-1 do begin
      ip=where(pr[ix,iy,ipsel] lt pr[ix,iy,0])
      ip=[0,ip+1]
      dp=deriv(reform(pr[ix,iy,ip]))*(-1.)
      var_era[ix,iy]=total(qv[ix,iy,ip]*dp,/double)
    endfor
    endfor
    var_era/=9.81

  endif else if var_str eq 'wspd' then begin

    count=[dims.nx+1,dims.ny,1,1] & offset=[0,0,izlev_era,0] ; x, y, p, t
    iu=reform(read_nc_var(era_fil,'UU',count=count,offset=offset))
    count=[dims.nx,dims.ny+1,1,1] & offset=[0,0,izlev_era,0] ; x, y, p, t
    iv=reform(read_nc_var(era_fil,'VV',count=count,offset=offset))

    uera=fltarr(dims.nx,dims.ny)
    vera=uera

    ;DESTAGGER
      ixin=indgen(dims.nx+1) & ixout=findgen(dims.nx)+0.5
      for iy=0,dims.ny-1 do uera[*,iy]=interpol(reform(iu[*,iy]),ixin,ixout)
      iyin=indgen(dims.ny+1) & iyout=findgen(dims.ny)+0.5
      for ix=0,dims.nx-1 do vera[ix,*]=interpol(reform(iv[ix,*]),iyin,iyout)

    var_era=sqrt(uera^2+vera^2)

  endif

  if iwind then begin

      count=[dims.nx+1,dims.ny,1,1] & offset=[0,0,izlev_era,0] ; x, y, p, t
      iu=reform(read_nc_var(era_fil,'UU',count=count,offset=offset))
      count=[dims.nx,dims.ny+1,1,1] & offset=[0,0,izlev_era,0] ; x, y, p, t
      iv=reform(read_nc_var(era_fil,'VV',count=count,offset=offset))

      uera=fltarr(dims.nx,dims.ny)
      vera=uera

      ;DESTAGGER
        ixin=indgen(dims.nx+1) & ixout=findgen(dims.nx)+0.5
        for iy=0,dims.ny-1 do uera[*,iy]=interpol(reform(iu[*,iy]),ixin,ixout)
        iyin=indgen(dims.ny+1) & iyout=findgen(dims.ny)+0.5
        for ix=0,dims.nx-1 do vera[ix,*]=interpol(reform(iv[ix,*]),iyin,iyout)

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

  if var_str eq 'rainrate' or var_str eq 'RAINNC' then begin

    ;ACCUMULATED RAIN OR RAIN RATE
    iv_str=var_str
    iv=where(vars.vars eq iv_str)
    file=dirs.files_post[ic,iv]

    count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
    rain=reform(read_nc_var(file,iv_str,count=count,offset=offset)) ; Accumulated rain, mm

  endif else if var_str eq 'wspd' then begin

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
    var = sqrt(u^2+v^2)
    u=0 & v=0

  endif else begin

    iv_str=strupcase(var_str)
    iv=where(vars.vars eq iv_str)
    file=dirs.files_post[ic,iv]

    if var_str eq 'avor' then klev=izlev else klev=0

    count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,klev,0] ; x,y,z,t
    var=reform(read_nc_var(file,iv_str,count=count,offset=offset))

  endelse

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

  ;AVERAGE OVER ND_AVG DAYS
    t_ind=indgen(npd*nd_avg)+nd_avg*npd

    ivar = mean(var[*,*,t_ind],dimension=3,/nan,/double)

  ;WRF WIND
    if iwind then begin
    ;U
      iu=mean(u[*,*,t_ind],dimension=3,/nan,/double)
    ;V
      iv=mean(v[*,*,t_ind],dimension=3,/nan,/double)
    endif

  ;ERA5 WIND
    if iwind then begin
      iuera=uera
      ivera=vera
    endif

;----CREATE PLOTS--------------------

  if var_str eq 'RAINNC' then begin
    setmin=0
    setmax=36;50;200
  endif else if var_str eq 'pw' then begin
    setmin=30
    setmax=70
  endif else if var_str eq 'wspd' then begin
    setmin=0
    setmax=15
  endif

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

    s_wrf=create_struct('lon',dims.lon,'lat',dims.lat,'var',ivar,'u',iu,'v',iv)
    s_obs=create_struct('lon',eralon,'lat',eralat,'var',var_era,'lonwind',eralon,'latwind',eralat,'u',iuera,'v',ivera)

    plot_myanmar_obcomp_maps, dirs, figname, s_wrf, s_obs, figspecs

endfor ; icase

;endfor ; ivar

endfor ; iday

endfor ; ivar

print,'DONE!!'
end
