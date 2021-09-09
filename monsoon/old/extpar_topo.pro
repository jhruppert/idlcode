; 
; Write topography into EXTPAR file for ICON.
; 
; Assumes a Gaussian 2D ridge.
; 
; James Ruppert
; 15.12.17
; 
pro extpar_topo

;SOME SETTINGS

  hmax = 1.5e3 ; max topography height
  xscale = 0.14 ; scale for standard deviation

;DIRECTORIES / FILES

  homedir='/home/mpim/m300462/'
  extpardir=homedir+'mistwork/extpar_lem/modify/'
  extparfil_in=extpardir+'extpar_modified_50x58_3000m.nc'
  extparfil_out=extpardir+'extpar_modified_topo_50x58_3000m.nc'
  figdir=homedir+'misthome/IDLWorkspace85/Default/figures/monsoon/'

;OPEN READ FILE
  fid=ncdf_open(extparfil_in,/nowrite)

  ;LAT and LON
    vid=ncdf_varid(fid,'lon')
    ncdf_varget,fid,vid,lon
    vid=ncdf_varid(fid,'lat')
    ncdf_varget,fid,vid,lat
    nx=n_elements(lon)

  ;LAND/SEA
    vid=ncdf_varid(fid,'FR_LAND')
    ncdf_varget,fid,vid,iland
    xland = where(iland eq 1,nland,complement=xocean)

;CLOSE FILE
  ncdf_close,fid

;OPEN WRITE FILE
  fid=ncdf_open(extparfil_out,/write)

  ;TOPO
    vid_topo=ncdf_varid(fid,'topography_c')
;    ncdf_varget,fid,vid_topo,topog

  ;EDIT TOPO
    topog=fltarr(nx)
    topog[*]=0.
    xswath = max(lon[xland]) - min(lon[xland])
    xmid =0.5*(max(lon[xland]) + min(lon[xland]))
    c = xswath*xscale ; standard dev of Gaussian
    for ii=0,nland-1 do begin
      ix=xland[ii]
    ;ALL POINTS
;    for ix=0,nx-1 do begin
      topog[ix] = hmax * exp( -(lon[ix] - xmid)^2 / (2.*(c^2)) )
    endfor

  ;ADJUST COASTLINE AT SEA LEVEL
;    zdaj = min(topo)

  ;CHECK
    print,'land'
    stats,topog[xland]
    print,'ocean'
    stats,topog[xocean]

  ;WRITE TOPO TO NEW FILE
    ncdf_varput,fid,vid_topo,topog

;CLOSE FILE
  ncdf_close,fid


;PLOT MAP

  x0=min(lon) & x1=max(lon)
  x=findgen(100)/99 * (x1-x0) + x0
  xrange=[x0,x1]

  y0=min(lat) & y1=max(lat)
  y=findgen(100)/99 * (y1-y0) + y0
  yrange=[y0,y1]

  set_plot,'ps'
  figname=figdir+'new_topo'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=5,ysize=3,/inches,$
    /helvetica

  position=[0.1,0.1,0.9,0.9]
  plot,x,y,/nodata,position=position,$
    xstyle=9,xmargin=0,ystyle=9,$
    xminor=1,yminor=1,$
    xticklen=xticklen,yticklen=yticklen,$
    xrange=xrange,yrange=yrange,$
    charsize=0.8

  loadct,0,/silent

  ;LAND-SEA BOUNDARY
    contour,iland,lon,lat,/follow,/overplot,/irregular,$
      levels=0.5,c_colors=0,c_charsize=0.7,c_thick=4,c_labels=0

  ;TOPOGRAPHY
    loadct,0,/silent
    cint=100
    clevs=[0.1,1,10,findgen(50)*cint+cint]
    contour,topog,lon,lat,/follow,/overplot,/irregular,$
      levels=clevs,c_colors=0,c_charsize=0.7;,c_labels=replicate(1,50)

  device,/close

  convert_png,figname,res=200,/remove_eps

print,'Done!'

end
