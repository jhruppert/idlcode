; 
; Maps from GFS
;
; James Ruppert
; 3/21/19
; 
pro run_tc_gfs_maps

config_dir, dirs=dirs
figdir=dirs.figdir+'tc/gfs/maria/'

datdir=dirs.wkdir+'tc_stuff/gfs/maria/analysis_ncdf/'
spawn,'ls '+datdir+'met_em*',datfils

nt=n_elements(datfils)

tcname='maria'
tcyear='2017'
hurdat=read_hurdat(tcname,tcyear)


;----PLOT OPTIONS--------------------


iwind=1

allvars=['PMSL','VOR']
nvsel=n_elements(allvars)

;for ivar_sel=0,nvsel-1 do begin
for ivar_sel=1,1 do begin

var_str=allvars[ivar_sel]

print,'VAR: ',var_str

if var_str eq 'PMSL' then $
  z_sel=0 $
else $
  z_sel=13;9;19;9 ; 0=sfc, 9=700, 13=500, 19=200 hPa
;(ncdump on work/tc_stuff/gfs/maria/analysis_ncdf files)

plev=reform(read_nc_var(datfils[0],'PRES',count=[1,1,1,1],offset=[0,0,z_sel,0]))

pstr=strtrim(string(plev*1e-2,format='(i4)'),2)
print,'PRES: ',pstr

;for ifil=0,nt-1 do begin
for ifil=0,nt-12 do begin
;for ifil=12,12 do begin

file=datfils[ifil]


;----TIME SPECS--------------------


tim_str=strmid(file,21,13,/reverse);reform(read_nc_var(file,'Times'))
print,'TIME: ',tim_str

yy=strmid(tim_str,0,4)
mm=strmid(tim_str,5,2)
dd=strmid(tim_str,8,2)
hh=strmid(tim_str,11,2)
jultim=julday(mm,dd,yy,hh,0,0)

it_hur=where(hurdat.jultim eq jultim,count)
if count then $
  tcloc=[hurdat.lon[it_hur],hurdat.lat[it_hur]]


;----DOMAIN SPECS--------------------


lon=reform(read_nc_var(file,'XLONG_M'))
lat=reform(read_nc_var(file,'XLAT_M'))
lonu=reform(read_nc_var(file,'XLONG_U'))
latu=reform(read_nc_var(file,'XLAT_U'))
lonv=reform(read_nc_var(file,'XLONG_V'))
latv=reform(read_nc_var(file,'XLAT_V'))

dims=size(lon,/dimensions)
nx=dims[0]
ny=dims[1]

triangulate,lonu,latu,triu
triangulate,lonv,latv,triv


;----READ VARS--------------------


  if var_str eq 'VOR' then begin

    ;RELATIVE VORTICITY
    ;U
      count=[nx+1,ny,1,1] & offset=[0,0,z_sel,0] ; x,y,z,t
      u=reform(read_nc_var(file,'UU',count=count,offset=offset))
      u=griddata(lonu,latu,temporary(u),/linear,triangles=triu,xout=lon,yout=lat)
      u=reform(temporary(u),nx,ny)
    ;V
      count=[nx,ny+1,1,1] & offset=[0,0,z_sel,0] ; x,y,z,t
      v=reform(read_nc_var(file,'VV',count=count,offset=offset))
      v=griddata(lonv,latv,temporary(v),/linear,triangles=triv,xout=lon,yout=lat)
      v=reform(temporary(v),nx,ny)
    var=fltarr(nx,ny)
    ;Zeta = dvdx - dudy
      x=reform(lon[*,0])
      y=reform(lat[0,*])
      mpdeg=111e3 ; m
      ;-DUDY
      for ix=0,nx-1 do $
        var[ix,*]-=deriv(y*mpdeg,reform(u[ix,*]))
      ;DVDX
      for iy=0,ny-1 do begin
        dx=mpdeg*cos(y[iy]*!dtor)
        var[*,iy]+=deriv(x*dx,reform(v[*,iy]))
      endfor

  endif else $
    var=reform(read_nc_var(file,var_str))

  ;WIND
    if iwind then begin
      ;U
      count=[nx+1,ny,1,1] & offset=[0,0,z_sel,0] ; x,y,z,t
      u=reform(read_nc_var(file,'UU',count=count,offset=offset))
      u=griddata(lonu,latu,temporary(u),/linear,triangles=triu,xout=lon,yout=lat)
      u=reform(temporary(u),nx,ny)
      ;V
      count=[nx,ny+1,1,1] & offset=[0,0,z_sel,0] ; x,y,z,t
      v=reform(read_nc_var(file,'VV',count=count,offset=offset))
      v=griddata(lonv,latv,temporary(v),/linear,triangles=triv,xout=lon,yout=lat)
      v=reform(temporary(v),nx,ny)
      wind=create_struct('u',u,'v',v)
      cvar=sqrt(u^2+v^2)
    endif
 

;----CREATE PLOTS--------------------


  if var_str eq 'PMSL' then begin
    setmin=980
    setmax=1020
    landcol=255
  endif else if var_str eq 'VOR' then begin
    setmax=10.
    setmin=-1.*setmax
    landcol=0
  endif

  gfs_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=2.5
  spawn,'mkdir '+figdir,tmp,tmpe

  figname=figdir+var_str+'_'+pstr+'hpa_'+tim_str

  figspecs=create_struct(figspecs,'figname',figname)
  figspecs=create_struct(figspecs,'landcol',landcol)
  figspecs.title=tim_str

  gfs_map_plot, var, lon, lat, figspecs, cvar=cvar, wind=wind, tcloc=tcloc


endfor ; ivar

endfor ; ifil

print,'DONE!!'
end
