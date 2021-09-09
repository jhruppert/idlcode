; 
; Generate azimuthal dataset (resolved in azimuth) from Cartesian input.
;
; Writes out to file: azim_var( radius, azimuth, height, time )
;
; James Ruppert
; 4/25/19
; 
pro wrf_regrid_azim, azim_file, infil, var_tag, t_ind_read, dims, tcloc, override=override


;UP-LEVEL IF EXISTS
exists=file_test(azim_file)
;if exists and ~override then begin
;  print,azim_file+' already there. Moving onto next var...'
;  return
;endif


;OTHERWISE GENERATE
print,'Generating azimuthal var!'


;LOCATE PRESSURE MIN

;Location info now provided (via 700-hPa AVOR)


;SET UP GRIDS

  ;VERTICAL LEVELS
    fid=ncdf_open(infil)
    did=ncdf_dimid(fid,'level')
    ncdf_diminq,fid,did,name,nz
    ncdf_close,fid

  ;RADIUS ARRAY
;    rmax=800. ; km
    rmax=1300. ; km <-- Files used for paper
    drad=3.   ; km
    nrad=round(rmax/drad) ; 800 km at 3 km spacing
    radius=findgen(nrad)/nrad*rmax
  ;AZMIUTH ARRAY
    nazim=360
    azim=indgen(nazim) ; 0-359 deg
  
  ;INPUT LAT/LON GRID IN CARTESIAN (KM), WITH 0 AT DOMAIN-CENTER
    x=fltarr(dims.nx,dims.ny)
    y=x
    for iy=0,dims.ny-1 do $
      x[*,iy]=111.*(dims.lon-mean(dims.lon)) * cos(dims.lat[iy]*!pi/180)
    for ix=0,dims.nx-1 do $
      y[ix,*]=111.*(dims.lat-mean(dims.lat))


;LOOP THROUGH TIME

nt=n_elements(t_ind_read)
azim_var=fltarr(nrad,nazim,nz,nt)

for it=0,nt-1 do begin

print,'IT: ',it,' of ',nt-1

;SHIFT INPUT GRID TO TC CENTER
  
  tclon0=tcloc[0,t_ind_read[it]]
  itcx=(where(dims.lon eq tclon0))[0]

  tclat0=tcloc[1,t_ind_read[it]]
  itcy=(where(dims.lat eq tclat0))[0]

  x_shift = x-x[itcx,itcy]
  y_shift = y-y[itcx,itcy]
  xy_radius = sqrt(x_shift^2+y_shift^2)

;itcloc = [ tcloc[0,t_ind_read[it]] , tcloc[1,t_ind_read[it]] ]
;xy_radius2=radius_tc_ll(itcloc,dims.lon,dims.lat)

  triangulate,x_shift,y_shift,tri

;GENERATE LATS/LONS OF CYLINDRICAL GRID

  newx=fltarr(nrad,nazim)
  newy=newx

  for ir=0,nrad-1 do begin
    for iaz=0,nazim-1 do begin
      newx[ir,iaz] = radius[ir] * cos(azim[iaz]*!pi/180)
      newy[ir,iaz] = radius[ir] * sin(azim[iaz]*!pi/180)
    endfor
  endfor

;READ AND REMAP VAR

  count=[dims.nx,dims.ny,nz,1] & offset=[0,0,0,t_ind_read[it]] ; x,y,z,t
  invar=read_nc_var(infil,var_tag,count=count,offset=offset)
  invar=reform(temporary(invar),dims.nx,dims.ny,nz)

;LOOP THROUGH HEIGHT AND GRID
for iz=0,nz-1 do $
  azim_var[*,*,iz,it] = griddata(x_shift,y_shift,reform(invar[*,*,iz]),triangles=tri,xout=newx,yout=newy,/linear,missing=!values.f_nan)

endfor ; it

  write_sing_ncvar,azim_file,azim_var,var_tag,$
    dim1=radius,dim2=azim,$
    dimtag1='radius',dimtag2='azmiuth',dimtag3='level',dimtag4='time'


end
