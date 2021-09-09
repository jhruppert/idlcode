; 
; Calculate the coast-normal wind speed and mask, average, and filter it.
;
; Called by run_imerg_windhist.pro
;
; James Ruppert
; 8/6/21
; 
function read_lsmask, lonout, latout

config_dir,dirs=dirs

nxout=n_elements(lonout)
nyout=n_elements(latout)

;----DIRECTORIES--------------------

  metfil_dir=dirs.wkdir+'wrf_wps/wrfv4/output_myanmar/JJAS_2013-2017/era5/'

;----READ LAND/SEA MASK--------------------

  spawn,'ls '+metfil_dir+'met_em.d01*',era_fil
  fil=era_fil[0]
;  topo=reform(read_nc_var(fil,'HGT_M'))*1e3 ; m --> km
  lsmask=reform(read_nc_var(fil,'LANDMASK'))
  dimsw=size(lsmask,/dimen)
  nxw=dimsw[0] & nyw=dimsw[1]

  lonw=read_nc_var(fil,'XLONG_M') ; deg
  lonw=reform(lonw[*,0])
  latw=read_nc_var(fil,'XLAT_M') ; deg
  latw=reform(latw[0,*])

  ;INTERPOLATE ONTO IMERG GRID
    lsmaskx=fltarr(nxout,nyw)
    for iy=0,nyw-1 do lsmaskx[*,iy]=interpol(reform(lsmask[*,iy]),lonw,lonout)
    lsm2=fltarr(nxout,nyout)
    for ix=0,nxout-1 do lsm2[ix,*]=interpol(reform(lsmaskx[ix,*]),latw,latout)
    lsmask=lsm2

return,lsmask

end
