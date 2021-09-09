; 
; Maps of Myanmar WRF output, with time-averaging.

; James Ruppert
; 8/14/20
; 
pro run_avg_myanmar_maps

config_dir,dirs=dirs

;EXPERIMENT SETTINGS
  expname='myanmar'
  cases=['ctl','morr_mynn']

;idir='9km'
idir='2dom'
if idir eq '2dom' then begin
  domtag='d01'
;  domtag='d02'
endif

;Set this if uneven n-output files
  nfils_set=10*24;240-24
casedir=dirs.scdir+expname+'/WRF/experiment/'+idir+'/'
dirs.figdir+=expname+'/'
config_wrfexp, casedir=casedir,cases=cases,dirs=dirs,$
  dims=dims, vars=vars, nfils_set=nfils_set, domtag=domtag;, /verbose

;----PLOT OPTIONS--------------------

idaily=1 ; daily-mean?
iwind=1 ; Plot 10m wind vectors?
ismooth=0 ; smooth output vars?

allvars=['wspd','RAINNC','rainrate','HFX','avor','slp','lh','th_e','olr','pw']
allvars=['RAINNC']
;allvars=['avor']
;allvars=['wspd']
;allvars=['HFX']
;allvars=['lh']
;allvars=['sef']
nvsel=n_elements(allvars)

;LEVEL SELECTION
  psel=925;700 700 is the winner of 600, 700, 800, 850, 925
  izlev=(where(dims.pres eq psel))[0]
  psel_wind=psel;500;700;925
  izlev_wind=(where(dims.pres eq psel_wind))[0]

;----TIME SPECS--------------------

;FULL TIME SERIES
  time=dims.time
  nt_full=dims.nt
  npd=dims.npd
  nhrs=1.*nt_full*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nhrs)

;----READ VARS--------------------

for ic=0,dirs.nc-1 do begin
;for ic=1,1 do begin

  print,'CASE: ',strupcase(dirs.cases[ic])

    i_nt=nt_full
    if ( strmatch(dirs.cases[ic],'*36h*') or strmatch(dirs.cases[ic],'icrf_*') or strmatch(dirs.cases[ic],'lwcr*') $
      or (dirs.cases[ic] eq 'lwswcrf') or (dirs.cases[ic] eq 'axisym') ) then i_nt-=36
    if strmatch(dirs.cases[ic],'*24h*') then i_nt-=24
    if strmatch(dirs.cases[ic],'*48h*') then i_nt-=48
    it_test=indgen(i_nt)+nt_full-i_nt

;IVAR LOOP
;for ivar_sel=0,nvsel-1 do begin
for ivar_sel=0,0 do begin

var_str=allvars[ivar_sel]
print,'VAR: ',var_str


  iwmax=0

  if var_str eq 'th_e' then begin
    ;T
      iv=where(vars.vars eq 'T2')
      file=dirs.files_post[ic,iv]
      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
      t=reform(read_nc_var(file,'',varid='0',count=count,offset=offset))
    ;q
      iv=where(vars.vars eq 'Q2')
      file=dirs.files_post[ic,iv]
      q=reform(read_nc_var(file,'',varid='0',count=count,offset=offset))
    ;PSFC
      iv=where(vars.vars eq 'PSFC')
      file=dirs.files_post[ic,iv]
      p=reform(read_nc_var(file,'',varid='0',count=count,offset=offset))
    ;THETA-E (ASSUME R_TOT = R_V; EQN SOURCE: BRYAN AND FRITSCH 2002, P.2921)
      var = theta_e(t,p,q,q)
  endif else if var_str eq 'wspd' then begin
    ;U
      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,izlev_wind,0] ; x,y,z,t
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
    iwmax=1
  endif else if var_str eq 'lw' then begin
    ;VERTICALLY INTEGRATED LONGWAVE RADIATIVE HEATING
      iv=where(vars.vars eq 'RTHRATLW')
      file=dirs.files_post[ic,iv]
      count=[dims.nx,dims.ny,dims.np,i_nt] & offset=[0,0,0,0] ; x,y,z,t
      lw=reform(read_nc_var(file,'RTHRATLW',count=count,offset=offset))*1004. ; K/s --> J/kg/s = W/kg
      var=total(lw,3,/double)*2500./9.81 ; W/kg --> W/m2
      i_rm_mean=1
      if i_rm_mean then begin
        for it=0,i_nt-1 do begin
          varm=mean(var[*,*,it],/nan,/double)
          var[*,*,it]-=varm
        endfor
      endif
  endif else if var_str eq 'rainrate' or var_str eq 'RAINNC' then begin
    ;ACCUMULATED RAIN OR RAIN RATE
    iv_str=var_str
    iv=where(vars.vars eq iv_str)
    file=dirs.files_post[ic,iv]

    count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
    var=reform(read_nc_var(file,iv_str,count=count,offset=offset)) ; mm/d

    ;ACCUMULATED RAIN OVER LAST n DAYS OF TIME SERIES
    ndavg=3
;if idir eq '2dom' then ndavg=4
;    var=reform(var[*,*,i_nt-1] - var[*,*,i_nt-1-ndavg*npd])
;var/=(1.*ndavg)
;    t_ind=indgen(npd*ndavg)+i_nt-npd*ndavg - 1
;    varm=total(var[*,*,t_ind],3,/double)

    ;DAILY AVG
    ndavg=fix(nd)
    vard=fltarr(dims.nx,dims.ny,ndavg)
    for iday=0,ndavg-1 do begin
      t_ind=indgen(npd)+iday*npd
;      t_ind=indgen(npd)+(i_nt-1-npd) - (ndavg-1-iday)*npd
;      vard[*,*,iday]=mean(var[*,*,t_ind],dimension=3,/nan,/double) / (1.*npd)
      ind0=t_ind[0]-1
      if ind0 lt 0 then ind0=0
      ind1=t_ind[npd-1]
      vard[*,*,iday] = (var[*,*,ind1] - var[*,*,ind0]) ; mm/d
    endfor

  endif else if var_str eq 'sef' then begin

    iv=where(vars.vars eq 'HFX')
    file=dirs.files_post[ic,iv]
    count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
    var=reform(read_nc_var(file,'HFX',count=count,offset=offset))
    iv=where(vars.vars eq 'LH')
    file=dirs.files_post[ic,iv]
    count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
    var+=reform(read_nc_var(file,'LH',count=count,offset=offset))

  endif else begin

    iv_str=strupcase(var_str)
    iv=where(vars.vars eq iv_str)
    file=dirs.files_post[ic,iv]

    if var_str eq 'avor' then klev=izlev else klev=0

    count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,klev,0] ; x,y,z,t
    var=reform(read_nc_var(file,iv_str,count=count,offset=offset))

  endelse

  ;RAIN RATE
;  if var_str eq 'rainrate' then var*=1./24 ; mm/d --> mm/h

  ;WIND
    if iwind then begin
    ;U
      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,izlev_wind,0] ; x,y,z,t
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

  if var_str eq 'rainrate' or var_str eq 'RAINNC' then irain=1 else irain=0

  ;DAILY AVERAGE
  if idaily then begin
    ndavg=fix(nd)
    if ~irain then vard=fltarr(dims.nx,dims.ny,ndavg)
    ud=fltarr(dims.nx,dims.ny,ndavg)
    vd=ud
    for id=0,ndavg-1 do begin
      t_ind=indgen(npd)+id*npd
      if ~irain then vard[*,*,id]=mean(var[*,*,t_ind],dimension=3,/nan,/double)
      ud[*,*,id]=mean(u[*,*,t_ind],dimension=3,/nan,/double)
      vd[*,*,id]=mean(v[*,*,t_ind],dimension=3,/nan,/double)
    endfor
  endif

;----CREATE PLOTS--------------------

  if var_str eq 'th_e' then begin
    setmin=340
    setmax=370
  endif else if var_str eq 'wspd' then begin
    setmin=0
    setmax=20
  endif else if var_str eq 'slp' then begin
    setmin=1000
    setmax=1015
  endif else if var_str eq 'pw' then begin
    setmin=30
    setmax=70
  endif else if var_str eq 'olr' then begin
    setmin=90
    setmax=290
  endif else if var_str eq 'lh' then begin
    setmin='0'
    setmax=300
  endif else if var_str eq 'HFX' then begin
    setmin='0'
    setmax=100
  endif else if var_str eq 'sef' then begin
    setmin='0'
    setmax=300
  endif else if var_str eq 'avor' then begin
    setmin=0
    setmax=20;40
  endif else if var_str eq 'rainrate' then begin
    setmin=0
    setmax=60
  endif else if var_str eq 'RAINNC' then begin
    setmin=0
    setmax=80
  endif else if var_str eq 'SST' then begin
    setmin=26
    setmax=32
  endif

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figdir=dirs.figdir+dirs.cases[ic]+'/'+var_str+'/'
  spawn,'mkdir '+figdir,tmp,tmpe
  figspecs=create_struct(figspecs,'figname',' ')
;  figspecs.title+=' ('+strupcase(dirs.cases[ic]);+')'
  titlesav=figspecs.title+' ('+strupcase(dirs.cases[ic]);+')'

  ;LOOP THROUGH TIMES
;  for it=0,ndavg-1 do begin
;  for it=4,4 do begin
  for it=ndavg-1,ndavg-1 do begin

    ;TIME INFO

    t_stamp='d'+strtrim(it,2);+'-d avg'
    t_stamp_fig=string(it,format='(i2.2)')+'d'
    ivar=reform(vard[*,*,it])
    if iwind then begin
      iu = reform(ud[*,*,it])
      iv = reform(vd[*,*,it])
    endif

    ;AVERAGE OVER N-DAYS
    iavg_all=1
    if iavg_all then begin
      all_ndavg=4;10
      t_ind=indgen(all_ndavg)+ndavg-all_ndavg
      ivar=mean(vard[*,*,t_ind],dimension=3,/double)
      iu=mean(ud[*,*,t_ind],dimension=3,/double)
      iv=mean(vd[*,*,t_ind],dimension=3,/double)
      t_stamp=strtrim(all_ndavg,2)+'-d avg'
      t_stamp_fig='avg_'+string(all_ndavg,format='(i2.2)')+'d'
    endif

    print,t_stamp_fig

    figspecs.title=titlesav+'; '+t_stamp+')'

    ;SMOOTH VARIABLES
    if ismooth then begin
      ndeg=0.5 ; output grid (degrees)
      nskip=round(ndeg/(dims.lon[1]-dims.lon[0]))
      ivar=gauss_smooth(temporary(ivar),replicate(nskip,2),/edge_truncate)
      iu=gauss_smooth(temporary(iu),replicate(nskip,2),/edge_truncate)
      iv=gauss_smooth(temporary(iv),replicate(nskip,2),/edge_truncate)
    endif

    if iwind then begin
      wind=create_struct('u',iu,'v',iv)
      cvar=sqrt(iu^2+iv^2)
    endif

    figname=figdir+t_stamp_fig+'_'+var_str+'_'+dirs.cases[ic]+'_'+domtag
    if var_str eq 'avor' then figname+='_'+string(psel,format='(i3.3)')
    figspecs.figname=figname

    stats,ivar

    wrf_myanmar_map_plot, dirs, ivar, dims.lon, dims.lat, figspecs, wind=wind;, cvar=cvar

  endfor

endfor ; icase

;endfor ; itfil (plotting via raw files)

endfor ; ivar

print,'DONE!!'
end
