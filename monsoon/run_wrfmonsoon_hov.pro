; 
; Hovmoller diagrams from Myanmar/Meiyu WRF output
;
; James Ruppert
; 10/13/20
; 
pro run_wrfmonsoon_hov

config_dir,dirs=dirs

;EXPERIMENT SETTINGS
;  expname='meiyu'
;  cases=['ctl']
  expname='myanmar'
  cases=['ctl','morr_mynn','morr_ysu','morr_acm2']
  cases=['kfcu']

;idir='9km'
idir='2dom'
if idir eq '2dom' then begin
  domtag='d01'
;  domtag='d02'
endif

;Set this if uneven n-output files
;  nfils_set=6.*24+1;10*24;240-24
casedir=dirs.scdir+expname+'/WRF/experiment/'+idir+'/'
dirs.figdir+='experiments/'+expname+'/'
config_wrfexp, casedir=casedir,cases=cases,dirs=dirs,$
  dims=dims, vars=vars, nfils_set=nfils_set, domtag=domtag;, /verbose

;Local SOLAR time conversion
local=round(mean(dims.lon)/360.*24.) ; deg lon --> hours
print,'Adding +'+strtrim(local,2)+' for LT'

;----PLOT OPTIONS--------------------

if expname eq 'myanmar' then begin
  ;CROSS SECTION
    xcross=[77.8,97.4] ; lon,lat
    ycross=[12.,20.]   ; lon,lat
    width=3.5 ; width of cross section (degrees)
endif else if expname eq 'meiyu' then begin
endif

ismooth=0 ; smooth output vars?

dcomp=0 ; diurnally composite?
nd_comp=5 ; days to composite over (from end of simulation)

iwind=0 ; overlay wind?

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

;TERRAIN HEIGHT
  file=dirs.files_raw[0,0]
  topo=reform(read_nc_var(file,'HGT'))*1e-3 ; m --> km
  lsmask=reform(read_nc_var(file,'LANDMASK'))

;INTERPOLATE ONTO HOV GRID
  cross_topo=cross_diag(topo,dims.lon,dims.lat,width,x_bounds=xcross,y_bounds=ycross,lonout=xlon,latout=xlat)
  topo=mean(cross_topo,dimension=2,/nan,/double) & cross_topo=0
  cross_lsmask=cross_diag(lsmask,dims.lon,dims.lat,width,x_bounds=xcross,y_bounds=ycross)
  lsmask=mean(cross_lsmask,dimension=2,/nan,/double) & cross_lsmask=0

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
    var=reform(read_nc_var(file,iv_str,count=count,offset=offset))

    ;ADD CUMULUS PARAM RAIN IF USED
;      if strmatch(dirs.cases[ic],'kfcu*') then begin
;        iv_str='RAINC'
;        iv=where(vars.vars eq iv_str)
;        file=dirs.files_post[ic,iv]
;        var+=reform(read_nc_var(file,iv_str,count=count,offset=offset)) ; Accumulated rain, mm
;        iv_str='RAINSH'
;        iv=where(vars.vars eq iv_str)
;        file=dirs.files_post[ic,iv]
;        var+=reform(read_nc_var(file,iv_str,count=count,offset=offset)) ; Accumulated rain, mm
;      endif

    rr=var
    rr[*]=!values.f_nan
    for it=1,i_nt-2 do rr[*,*,it]=0.5*(var[*,*,it+1]-var[*,*,it-1]) ; mm/h
    var=rr & rr=0

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


;----DIURNALLY COMPOSITE--------------------

if dcomp then begin

  vard=fltarr(dims.nx,dims.ny,npd)

  nd=i_nt/npd
  iday_comp=indgen(nd_comp)+(nd-nd_comp)
  it_hr=iday_comp*npd

  for it=0,npd-1 do vard[*,*,it]=mean(var[*,*,it_hr+it],dimension=3,/nan,/double)
  var=vard & vard=0

  if iwind then begin
    ud=fltarr(dims.nx,dims.ny,npd)
    vd=ud
    for it=0,npd-1 do begin
      ud[*,*,it]=mean(u[*,*,it_hr+it],dimension=3,/nan,/double)
      vd[*,*,it]=mean(v[*,*,it_hr+it],dimension=3,/nan,/double)
    endfor
    u=ud & v=vd & ud=0 & vd=0
  endif

endif

;----AVERAGE IN Y TO MAKE HOVMOLLER--------------------

;ZONAL HOV
  ;x_ind=where(dims.lon ge lonsub[0] and dims.lon le lonsub[1])
  ;y_ind=where(dims.lat ge latavg[0] and dims.lat le latavg[1])
  ;lonhov=dims.lon[x_ind]
;varavg=mean(var[x_ind,y_ind,*],dimension=2,/double,/nan)
;var=varavg & varavg=0

;INTERPOLATE VARS ONTO HOV GRID
  cross=cross_diag(var,dims.lon,dims.lat,width,x_bounds=xcross,y_bounds=ycross)
  var=mean(cross,dimension=2,/nan,/double) & cross=0

if iwind then begin
;  um=mean(u[x_ind,y_ind,*],dimension=2,/double,/nan)
;  u=um & um=0
;  vm=mean(v[x_ind,y_ind,*],dimension=2,/double,/nan)
;  v=vm & vm=0
  ucross=cross_diag(u,dims.lon,dims.lat,width,x_bounds=xcross,y_bounds=ycross)
  u=mean(ucross,dimension=2,/nan,/double) & ucross=0
  vcross=cross_diag(v,dims.lon,dims.lat,width,x_bounds=xcross,y_bounds=ycross)
  v=mean(vcross,dimension=2,/nan,/double) & vcross=0
endif

;DISTANCE RELATIVE TO COASTLINE (KM)
  ;DISTANCE FROM LAT/LON
    nxhov=n_elements(xlon)
    xhov=fltarr(nxhov)
    for ix=0,nxhov-1 do $
      xhov[ix]=111.*sqrt( (xlon[ix]*cos(xlat[ix]*!dtor))^2 + xlat[ix]^2 ) ; Converts to km
  ;ADJUST BY COASTLINE LOCATION
    ixcheck=indgen(nxhov/2)+nxhov/2
    ixcoast=(where(lsmask[ixcheck] gt 0))[0]
    xhov=reverse(xhov)
    xhov-=xhov[ixcheck[ixcoast]]

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
    setmax=3.;1.5
  endif else if var_str eq 'SST' then begin
    setmin=26
    setmax=32
  endif

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figdir=dirs.figdir+'hov/'
  spawn,'mkdir '+figdir,tmp,tmpe
  figspecs=create_struct(figspecs,'figname',' ')
;  figspecs.title+=' ('+strupcase(dirs.cases[ic]);+')'
  titlesav=figspecs.title+' ('+strupcase(dirs.cases[ic])+')'

  if var_str eq 'RAINNC' then figspecs.cbar_format='(f3.1)'

    figspecs.title=titlesav;+'; '+t_stamp+')'

    ;SMOOTH VARIABLES
    if ismooth then begin
      ndeg=0.5 ; output grid (degrees)
      nskip=round(ndeg/(dims.lon[1]-dims.lon[0]))
      ivar=gauss_smooth(temporary(ivar),replicate(nskip,2),/edge_truncate)
      iu=gauss_smooth(temporary(iu),replicate(nskip,2),/edge_truncate)
      iv=gauss_smooth(temporary(iv),replicate(nskip,2),/edge_truncate)
    endif

    if dcomp then begin

      figname=figdir+var_str+'_dcomp_'+dirs.cases[ic]+'_'+domtag

      ltim=indgen(npd)+local
      for it=0,npd-1 do ltim[it]-=24*(ltim[it] ge 24)
      
      ;REORDER TO PUT 0 LT AT FIRST INDEX
      ltim=shift(temporary(ltim),local)
      var=shift(temporary(var),[0,local])

      if iwind then begin
        u=shift(temporary(u),[0,local])
        v=shift(temporary(v),[0,local])
      endif

    endif else $
      figname=figdir+var_str+'_'+dirs.cases[ic]+'_'+domtag

    if var_str eq 'avor' then figname+='_'+string(psel,format='(i3.3)')
    figspecs.figname=figname

    stats,var

    if iwind then wind=create_struct('u',iu,'v',iv)

    if dcomp then $
      wrf_monsoon_dcomp_hov, dirs, figspecs, var, ltim, xhov, wind=wind, cvar=cvar $
    else $
      wrf_monsoon_hov, dirs, figspecs, var, time, npd, xhov, wind=wind, cvar=cvar

endfor ; icase

;endfor ; itfil (plotting via raw files)

endfor ; ivar

print,'DONE!!'
end
