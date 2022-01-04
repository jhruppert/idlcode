; 
; Calculate the coast-normal wind speed and mask, average, and filter it.
;
; Called by run_imerg_windhist.pro
;
; James Ruppert
; 8/6/21
; 
function read_coastal_dist, lonout, latout

config_dir,dirs=dirs

nxout=n_elements(lonout)
nyout=n_elements(latout)

;----DIRECTORIES--------------------

  coast_dist_fil=dirs.wkdir+'dist2coast.nc'

;----READ COASTAL DISTANCE--------------------

  dist_fil=read_nc_var(coast_dist_fil,'coast_dist')
  londist=reform(dist_fil[0,*])
  latdist=reform(dist_fil[1,*])
  coast_dist=reform(dist_fil[2,*])

  ;KEEP ONLY ERA5 DOMAIN
  ixkeep=where((londist ge min(lonout)) and (londist le max(lonout)) and (latdist ge min(latout)) and (latdist le max(latout)))
  londist=londist[ixkeep]
  latdist=latdist[ixkeep]
  coast_dist=coast_dist[ixkeep]

  ;INTERPOLATE ONTO ERA5 GRID
    triangulate,londist,latdist,tri
    coast_dist_im=griddata(londist,latdist,coast_dist,/linear,triangles=tri,xout=lonout,yout=latout,/grid)
    coast_dist=coast_dist_im

return,coast_dist

end
