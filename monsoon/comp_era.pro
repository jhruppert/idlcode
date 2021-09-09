; 
; Compare boundary conditions from ERA5 and ERA-i
; 
; James Ruppert
; 9/25/20
; 
pro comp_era

config_dir,dirs=dirs

bcdir=dirs.wkdir+'wrf_wps/wrfv4/output_myanmar/JJAS_2013-2017/'

for iera=0,1 do begin

if iera eq 0 then begin
  eratag='ERA5'
  eraf=bcdir+'era5_avg.nc'
  pres=reform(read_nc_var(bcdir+'era5_pres.nc','PRES'))*1e-2
endif else begin
  eratag='ERAi'
  eraf=bcdir+'erai_avg.nc'
  pres=reform(read_nc_var(bcdir+'erai_pres.nc','PRES'))*1e-2
endelse

;----PLOT OPTIONS--------------------

iwind=1 ; Plot 10m wind vectors?

ismooth=0 ; smooth output vars?

var_str='wspd'
;var_str='U'
;var_str='V'
;var_str='RH'
var_str='AVOR'
;var_str='pw'

;LEVEL SELECTION
  psel=500;700
  izlev=(where(pres eq psel))[0]
  psel_wind=psel
  izlev_wind=izlev

;----READ WINDS--------------------

  ;GET DIMS
    lon=reform((reform(read_nc_var(eraf,'XLONG_M')))[*,20])
    lat=reform((reform(read_nc_var(eraf,'XLAT_M')))[20,*])
    nx=n_elements(lon) & ny=n_elements(lat)
    np=n_elements(pres)

  ;GET WINDS
    ;U
      count=[nx+1,ny,1,1] & offset=[0,0,izlev,0] ; x,y,z,t
      u=reform(read_nc_var(eraf,'UU',count=count,offset=offset))
    ;V
      count=[nx,ny+1,1,1] & offset=[0,0,izlev,0] ; x,y,z,t
      v=reform(read_nc_var(eraf,'VV',count=count,offset=offset))

  ;UNSTAGGER WINDS
    unst=fltarr(nx,ny)
    for iy=0,ny-1 do unst[*,iy]=interpol(reform(u[*,iy]),indgen(nx+1),findgen(nx)+0.5)
    u=unst & unst[*]=0.
    for ix=0,nx-1 do unst[ix,*]=interpol(reform(v[ix,*]),indgen(ny+1),findgen(ny)+0.5)
    v=unst & unst=0.

;----READ VARS--------------------

  if var_str eq 'wspd' then begin
    var = sqrt(u^2+v^2)
    ;u=0 & v=0
  endif else if var_str eq 'U' then begin
    var = u
  endif else if var_str eq 'V' then begin
    var = v
  endif else if var_str eq 'AVOR' then begin
    fo = fltarr(nx,ny)
    dvdx=fo & dudy=fo
    for iy=0,ny-1 do begin
      fo[*,iy]=2.*7.292e-5*sin(lat[iy]*!dtor)
      dvdx[*,iy] = deriv((lon*111.e3*cos(lat[iy]*!dtor)),reform(v[*,iy]))
    endfor
    for ix=0,nx-1 do $
      dudy[ix,*] = deriv((lat*111.e3),reform(u[ix,*]))
    var = (fo + dvdx - dudy) * 1e5 ; 10^-5 /s
  endif else if var_str eq 'pw' then begin

    count=[nx,ny,np,1] & offset=[0,0,0,0] ; x,y,z,t
    rh=reform(read_nc_var(eraf,'RH',count=count,offset=offset))
    tmpk=reform(read_nc_var(eraf,'TT',count=count,offset=offset))
    pr=reform(read_nc_var(eraf,'PRES',count=count,offset=offset))
    iice=where(tmpk lt 273.,complement=noice)
    qv=tmpk & qv[*]=0.
;    p=qv & for iz=0,np-1 do p[*,*,iz]=pres[iz]
    qv[iice]=rh2mixr(rh[iice],pr[iice]*1e-2,tmpk[iice],/ice)
    qv[noice]=rh2mixr(rh[noice],pr[noice]*1e-2,tmpk[noice],/ice)
    var=fltarr(nx,ny)
    ipsel=indgen(np-1)+1
    for ix=0,nx-1 do begin
    for iy=0,ny-1 do begin
      ip=where(pr[ix,iy,ipsel] lt pr[ix,iy,0])
      ip=[0,ip+1]
      dp=deriv(reform(pr[ix,iy,ip]))*(-1.)
      var[ix,iy]=total(qv[ix,iy,ip]*dp,/double)
    endfor
    endfor
    var/=9.81

  endif else begin
    iv_str=strupcase(var_str)
    count=[nx,ny,1,1] & offset=[0,0,izlev,0] ; x,y,z,t
    var=reform(read_nc_var(eraf,iv_str,count=count,offset=offset))
  endelse

  ;RAIN RATE
;  if var_str eq 'rainrate' then var*=1./24 ; mm/d --> mm/h

;----CREATE PLOTS--------------------

  if var_str eq 'wspd' then begin
    setmin=0
    setmax=12
  endif else if var_str eq 'U' or var_str eq 'V' then begin
    setmax=10
    setmin=-1.*setmax
  endif else if var_str eq 'slp' then begin
    setmin=1000
    setmax=1015
  endif else if var_str eq 'SST' then begin
    setmin=26
    setmax=32
  endif else if var_str eq 'RH' then begin
    setmin=50
    setmax=100
  endif else if var_str eq 'AVOR' then begin
    setmin='0'
    setmax=20
  endif else if var_str eq 'pw' then begin
    setmin=20
    setmax=80
  endif

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figdir=dirs.figdir+'myanmar/eracomp/'
  spawn,'mkdir '+figdir,tmp,tmpe
  figspecs=create_struct(figspecs,'figname',' ')
  if var_str eq 'pw' then $
    figspecs.title+=' ('+eratag+')' else $
    figspecs.title+=' ('+string(psel,format='(i3.3)')+' hPa; '+eratag+')'

    ivar=var

    ;SMOOTH VARIABLES
    if ismooth then begin
      ndeg=0.5 ; output grid (degrees)
      nskip=round(ndeg/(dims.lon[1]-dims.lon[0]))
      ivar=gauss_smooth(temporary(ivar),replicate(nskip,2),/edge_truncate)
      iu=gauss_smooth(temporary(iu),replicate(nskip,2),/edge_truncate)
      iv=gauss_smooth(temporary(iv),replicate(nskip,2),/edge_truncate)
    endif

    if iwind then begin
      iu = u
      iv = v
      wind=create_struct('u',iu,'v',iv)
      cvar=sqrt(iu^2+iv^2)
    endif

    figname=figdir+var_str+'_'+eratag
    ;if var_str eq 'avor' then $
    figname+='_'+string(psel,format='(i3.3)')
    figspecs.figname=figname

    stats,ivar

    wrf_myanmar_map_plot, dirs, ivar, lon, lat, figspecs, wind=wind;, cvar=cvar

endfor ; iera

print,'DONE!!'
end
