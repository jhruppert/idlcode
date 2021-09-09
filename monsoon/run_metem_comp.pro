; 
; Maps of MET_EM WPS output files side-by-side.
;
; James Ruppert
; 10/29/20
; 
pro run_metem_comp

config_dir,dirs=dirs

;EXPERIMENT SETTINGS
  expname='myanmar'
  cases=['ctl']

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

;MET_EM DIRECTORIES
dir1=dirs.wkdir+'wrf_wps/wrfv4/output_myanmar/JJAS_2013-2017/erai/'
met_fil1=dir1+'met_em.d01.2013-07-20_avg.nc'
dir2=dirs.wkdir+'wrf_wps/wrfv4/output_myanmar/raindays_800km_12/erai/'
met_fil2=dir2+'met_em.d01.2013-07-20_avg.nc'

t_stamp_fig='jjas13-17_comp'

;----PLOT OPTIONS--------------------

iwind=1 ; Plot wind vectors?

allvars=['wspd','slp','pw']
;allvars=['slp']
;allvars=['wspd']
;allvars=['pw']
nvsel=n_elements(allvars)
nvsel=n_elements(allvars)

;LEVEL SELECTION
  psel=925;700;800;700 700 is the winner of 600, 700, 800, 850, 925
;  izlev=(where(dims.pres eq psel))[0]

;IVAR LOOP
for ivar_sel=0,nvsel-1 do begin
;for ivar_sel=0,0 do begin

var_str=allvars[ivar_sel]
print,'VAR: ',var_str

;----READ MET_EM FILES--------------------

  lon=dims.lon
  lat=dims.lat

  ;PLEVS
  nz_era=38
  p_era=reform(read_nc_var(met_fil1,'PRES',count=[1,1,nz_era,1],offset=[0,0,0,0])) ; Pa
  izlev_era=(where(p_era*1d-2 eq psel))[0]

for imet=0,1 do begin

  if imet eq 0 then era_fil=met_fil1 else era_fil=met_fil2

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

  endif else if var_str eq 'slp' then begin

    count=[dims.nx,dims.ny,1] & offset=[0,0,0] ; x, y, p, t
    var_era=reform(read_nc_var(era_fil,'PMSL',count=count,offset=offset))

  endif else begin

    count=[dims.nx+1,dims.ny,1,1] & offset=[0,0,izlev_era,0] ; x, y, p, t
    var_era=reform(read_nc_var(era_fil,var_str,count=count,offset=offset))

  endelse

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

  if imet eq 0 then begin
    var1=var_era
    u1=uera
    v1=vera
  endif else begin
    var2=var_era
    u2=uera
    v2=vera
  endelse

endfor ; imet


;----CREATE PLOTS--------------------

  if var_str eq 'RAINNC' then begin
    setmin=0
    setmax=36;50;200
  endif else if var_str eq 'pw' then begin
    setmin=40
    setmax=75
  endif else if var_str eq 'wspd' then begin
    setmin=0
    setmax=12;15
  endif else if var_str eq 'slp' then begin
    setmin=995
    setmax=1015
  endif

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figdir=dirs.figdir+'/eracomp/'
  spawn,'mkdir '+figdir,tmp,tmpe
  figspecs=create_struct(figspecs,'figname',' ')
  titlesav=figspecs.title

  t_stamp=t_stamp_fig

    print,t_stamp

    figspecs.title=titlesav

    figname=figdir+t_stamp_fig+'_'+var_str
    if var_str eq 'avor' then figname+='_'+string(psel,format='(i3.3)')
    figspecs.figname=figname

    s_wrf=create_struct('lon',lon,'lat',lat,'var',var1,'u',u1,'v',v1)
    s_obs=create_struct('lon',lon,'lat',lat,'var',var2,'lonwind',lon,'latwind',lat,'u',u2,'v',v2)

    plot_myanmar_obcomp_maps, dirs, figname, s_wrf, s_obs, figspecs

endfor ; ivar

print,'DONE!!'
end
