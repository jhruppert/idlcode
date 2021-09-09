; 
; Diurnal composite Hovmoller diagrams from Myanmar/Meiyu IMERG data.
;
; James Ruppert
; 10/13/20
; 
pro run_imerg_dcomp_hov

message,'Update HOV prep to match run_imerg_windhist_dhov.pro'

config_dir,dirs=dirs

;NEED WRF SETTINGS JUST FOR MET_MET
;EXPERIMENT SETTINGS
  expname='myanmar'
  cases=['ctl'];,'morr_mynn']

;idir='9km'
idir='2dom'
if idir eq '2dom' then begin
  domtag='d01'
;  domtag='d02'
endif

;Set this if uneven n-output files
  nfils_set=241
casedir=dirs.scdir+expname+'/WRF/experiment/'+idir+'/'
dirs.figdir+=expname+'/'
config_wrfexp, casedir=casedir,cases=cases,dirs=dirs,$
  dims=dims, vars=vars, nfils_set=nfils_set, domtag=domtag;, /verbose

;OB DIRECTORIES
  maindir=dirs.scdir+'myanmar/'
  imergfil=maindir+'imerg/data/jjas_2013-2017_diurncomp_3B-HHR.MS.MRG.3IMERG.V06B.nc4';'imerg/20150620-20150629.nc4'
  npd_imerg=48
  ;era_fil=maindir+'era5/jjas_2013-2017/ERA5-JJAS_2013-2017_mean-pl.nc';'era5/ERA5-20150623-20150630-uv.nc'
  ;USING MET_EM FILES FOR DIURNAL COMPOSITE ERA5
  era_dir=dirs.wkdir+'wrf_wps/wrfv4/output_myanmar/JJAS_2013-2017/era5/'
  npd_era=24

;Local SOLAR time conversion
local=6;round(mean(dims.lon)/360.*24.) ; deg lon --> hours
print,'Adding +'+strtrim(local,2)+' for LT'

;----PLOT OPTIONS--------------------

if expname eq 'myanmar' then begin
  ;CROSS SECTIONS
  icross=1
  if icross eq 1 then begin
    ;SW-NE
    xcross=[77.8,97.4] ; lon,lat
    ycross=[12.,20.]   ; lon,lat
    width=3.0 ; width of cross section (degrees)
  endif else if icross eq 2 then begin
    ;NW-SE
    xcross=[83.,93.5] ; lon,lat
    ycross=[21.4,10]   ; lon,lat
    width=3.0 ; width of cross section (degrees)
  endif
endif else if expname eq 'meiyu' then begin
endif

ismooth=0 ; smooth output vars?

iwind=0 ; overlay wind?

;LEVEL SELECTION
  psel=925;700;800;700 700 is the winner of 600, 700, 800, 850, 925
  izlev=(where(dims.pres eq psel))[0]

;----READ OBS--------------------

  ;IMERG RAINFALL
    imerg=read_nc_var(imergfil,'precipitationCal') ; keep in mm/hr
    imerg=transpose(temporary(imerg))

  ;SHIFT TIME INDEX TO LAST POSITION
    imdims=size(imerg,/dim)
    im2=fltarr(imdims[1],imdims[2],imdims[0])
    for it=0,npd_imerg-1 do im2[*,*,it]=reform(imerg[it,*,*])
    imerg=im2 & im2=0

  ;IMERG LAT/LON
    lonim=read_nc_var(imergfil,'lon')
    latim=read_nc_var(imergfil,'lat')

;    npdim=48

  ;ERA5 WINDS

  eralon=dims.lon
  eralat=dims.lat

  if iwind then begin

    ;PLEVS
    ilevera=925;[925, 950, 975, 1000]
    spawn,'ls '+era_dir+'met_em.'+domtag+'*',era_fil
    era_fil=era_fil[0:23]
    p_era=reform(read_nc_var(era_fil[0],'PRES',count=[1,1,38,1],offset=[0,0,0,0])) ; Pa
    izlev_era=(where(p_era*1d-2 eq psel))[0]

    ;READ DIURNAL COMPOSITE
    u=fltarr(dims.nx,dims.ny,npd_era)
    v=u

    for it=0,npd_era-1 do begin
      fil=era_fil[it]
      count=[dims.nx+1,dims.ny,1,1] & offset=[0,0,izlev_era,0] ; x, y, p, t
      iu=reform(read_nc_var(fil,'UU',count=count,offset=offset))

    ;DESTAGGER
      ixin=indgen(dims.nx+1) & ixout=findgen(dims.nx)+0.5
      for iy=0,dims.ny-1 do u[*,iy,it]=interpol(reform(iu[*,iy]),ixin,ixout)
      count=[dims.nx,dims.ny+1,1,1] & offset=[0,0,izlev_era,0] ; x, y, p, t
      iv=reform(read_nc_var(fil,'VV',count=count,offset=offset))
      iyin=indgen(dims.ny+1) & iyout=findgen(dims.ny)+0.5
      for ix=0,dims.nx-1 do v[ix,*,it]=interpol(reform(iv[ix,*]),iyin,iyout)
    endfor

  endif

;READ TOPOGRAPHY AND LAND MASK FROM ERA5 MET_EM FILES
  spawn,'ls '+era_dir+'met_em.'+domtag+'*',era_fil
  fil=era_fil[0]
  topo=reform(read_nc_var(fil,'HGT_M'))*1e3 ; m --> km
  lsmask=reform(read_nc_var(fil,'LANDMASK'))

;----AVERAGE IN Y TO MAKE HOVMOLLER--------------------

;INTERPOLATE VARS ONTO HOV GRID
  cross_topo=cross_diag(topo,dims.lon,dims.lat,width,x_bounds=xcross,y_bounds=ycross,lonout=xlon,latout=xlat)
  topo=mean(cross_topo,dimension=2,/nan,/double) & cross_topo=0
  cross_lsmask=cross_diag(lsmask,dims.lon,dims.lat,width,x_bounds=xcross,y_bounds=ycross)
  lsmask=mean(cross_lsmask,dimension=2,/nan,/double) & cross_lsmask=0

;IMERG
  cross=cross_diag(imerg,lonim,latim,width,x_bounds=xcross,y_bounds=ycross,lonout=xlon_im,latout=xlat_im)
  imerg=mean(cross,dimension=2,/nan,/double) & cross=0

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
    xhov_era=fltarr(nxhov)
    for ix=0,nxhov-1 do $
      xhov_era[ix]=111.*sqrt( (xlon[ix]*cos(xlat[ix]*!dtor))^2 + xlat[ix]^2 ) ; Converts to km
    if icross eq 1 then xhov_era=reverse(xhov_era)
  ;ADJUST BY COASTLINE LOCATION
    if icross eq 1 then begin
      ixcheck=indgen(nxhov/2)+nxhov/2
      ixcoast=(where(lsmask[ixcheck] gt 0))[0]
    endif else if icross eq 2 then begin
      ixcheck=indgen(nxhov)
      ixcoast=(where(lsmask[ixcheck] eq 0))[0]
    endif
    xcoast=xhov_era[ixcheck[ixcoast]]
    xhov_era-=xcoast
  ;IMERG: DISTANCE FROM LAT/LON
    nxhov=n_elements(xlon_im)
    xhov_im=fltarr(nxhov)
    for ix=0,nxhov-1 do $
      xhov_im[ix]=111.*sqrt( (xlon_im[ix]*cos(xlat_im[ix]*!dtor))^2 + xlat_im[ix]^2 ) ; Converts to km
    if icross eq 1 then xhov_im=reverse(xhov_im)
    xhov_im-=xcoast

;----CREATE PLOTS--------------------

  var_str='RAINNC'
    setmin=0
    setmax=1.5

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figspecs.cbar_format='(f3.1)'
  figdir=dirs.figdir+'hov/'
  spawn,'mkdir '+figdir,tmp,tmpe
  figspecs=create_struct(figspecs,'figname',' ')
  titlesav=figspecs.title+' (IMERG)'

    figspecs.title=titlesav;+'; '+t_stamp+')'

    ;SMOOTH VARIABLES
    if ismooth then begin
      ndeg=0.5 ; output grid (degrees)
      nskip=round(ndeg/(dims.lon[1]-dims.lon[0]))
      ivar=gauss_smooth(temporary(ivar),replicate(nskip,2),/edge_truncate)
      iu=gauss_smooth(temporary(iu),replicate(nskip,2),/edge_truncate)
      iv=gauss_smooth(temporary(iv),replicate(nskip,2),/edge_truncate)
    endif

    figname=figdir+var_str+'_dcomp_imerg_'+domtag
    figspecs.figname=figname

;    ltim_imerg=findgen(npd_imerg)*24/npd_imerg+local
;    for it=0,npd_imerg-1 do ltim_imerg[it]-=24*(ltim_imerg[it] ge 24)

    ;REORDER TO PUT 0 LT AT FIRST INDEX
;    ltim_imerg=shift(temporary(ltim_imerg),local*npd_imerg/24)
    imerg=shift(temporary(imerg),[0,local*npd_imerg/24])

    if iwind then begin
      ltim_era=findgen(npd_era)*24/npd_era+local
      for it=0,npd_era-1 do ltim_era[it]-=24*(ltim_era[it] ge 24)
      ltim_era=shift(temporary(ltim_era),local*npd_era/24)
      u=shift(temporary(u),[0,local*npd_era/24])
      v=shift(temporary(v),[0,local*npd_era/24])
    endif

    stats,imerg

    if iwind then wind=create_struct('u',iu,'v',iv)

    wrf_monsoon_dcomp_hov, dirs, figspecs, imerg, indgen(npd_imerg), xhov_im, wind=wind, cvar=cvar

print,'DONE!!'
end
