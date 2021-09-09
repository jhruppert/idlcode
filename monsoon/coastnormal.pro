; 
; Calculate the coast-normal wind speed and mask, average, and filter it.
;
; Called by run_imerg_windhist.pro
;
; James Ruppert
; 8/6/21
; 
pro coastnormal, icross, coast_thresh, u, v, eralon, eralat, uindex=uindex, $ ; Uindex is the main output
  unorm=unorm, xcross=xcross, ycross=ycross ; Optional returns

config_dir,dirs=dirs

dims=size(u,/dimen)
nxera=dims[0]
nyera=dims[1]
nd=dims[2]

;----DIRECTORIES--------------------

  metfil_dir=dirs.wkdir+'wrf_wps/wrfv4/output_myanmar/JJAS_2013-2017/era5/'
  coast_dist_fil=dirs.home+'idl/code/misc/dist2coast.nc'

;=====BEGIN READING=========================================================

  coast_dist = read_coastal_dist(eralon, eralat)
  lsmask = read_lsmask(eralon, eralat)

;=====CREATE INDEX FOR REGRESSION/BINNING=========================================================

  imerg_cross_settings, icross, xcross=xcross, ycross=ycross, bounds_ind=bounds_ind

  ;CHECK FOR REVERSE LATITUDINAL SECTION
    if ycross[1] lt ycross[0] then irev=1 else irev=0

  dx=xcross[1]-xcross[0]
  dy=ycross[1]-ycross[0]
  th = atan( dy / dx ); * 180./!pi

;  if icross eq 1 then 
  unorm = u*cos(th) + v*sin(th)
  ;if icross eq 2 then unorm2 = u*cos(th) + v*sin(th)
  if irev then unorm*=-1

;endfor

  ;SUBSET TO BOUNDING BOX
    ix=where((eralon ge bounds_ind[0]) and (eralon le bounds_ind[2]),bnx)
    iy=where((eralat ge bounds_ind[1]) and (eralat le bounds_ind[3]),bny)
    unorm2=unorm[ix,*,*] & unorm2=unorm2[*,iy,*]
    coast_dist2=coast_dist[ix,*] & coast_dist2=coast_dist2[*,iy]
    lsmask2=lsmask[ix,*] & lsmask2=lsmask2[*,iy]

  ;ADD TIME INDICES TO THESE TO SIMPLIFY WHERE FUNCTIONS
    coast_dist2=rebin(coast_dist2,[bnx,bny,nd])
    lsmask2=rebin(lsmask2,[bnx,bny,nd])

  ;IMPOSE RULES
    isea=where((coast_dist2 le coast_thresh) and (lsmask2 le 0.1),complement=inot)
    unorm2[inot]=!values.f_nan

  uindex=mean(mean(unorm2,dimension=1,/nan,/double),dimension=1,/nan,/double)

;save,uindex,filename='uindex'

;----PLOT WIND INDEX--------------------

iplot_ind=0
if iplot_ind then begin

  var_str='pw'
  setmax=20 & setmin='0.'
  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax
  figname=figdir+'test_wind_index_cross'+strtrim(icross,2)
  figspecs=create_struct(figspecs,'figname',figname)

  varm=max(abs(unorm2[*,*,jjas]),dimension=3,/nan)
  varm[where(finite(varm))]=20
  wrf_myanmar_map_plot, dirs, varm, eralon[ix], eralat[iy], figspecs
exit
endif


end
