; 
; Hovmoller diagrams from Myanmar/Meiyu IMERG data.
;
; This now also gives the option to plot diurnal comp for select time range
; See run_imerg_dcomp_hov.pro for diurnal composite that uses all-season composite file
;
; James Ruppert
; 11/20/20
; 
pro run_imerg_hov

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

dcomp=0 ; diurnal composite?
ismooth=0 ; smooth output vars?

;Not currently set up
  iwind=0 ; overlay wind?

;SELECT DATE RANGE
  yy_plot=[2013,2013]
  mm_plot=[6,7]
  dd_plot=[1,30] ; inclusive

expname='myanmar'
if expname eq 'myanmar' then begin
  ;CROSS SECTION
    xcross=[77.8,97.4] ; lon,lat
    ycross=[12.,20.]   ; lon,lat
    width=3.5 ; width of cross section (degrees)
endif else if expname eq 'meiyu' then begin
endif

;LEVEL SELECTION
;  psel=925;700;800;700 700 is the winner of 600, 700, 800, 850, 925
;  izlev=(where(dims.pres eq psel))[0]

;----SOME WRF STUFF TO GET TOPO--------------------

cases=['ctl']
idir='2dom'
if idir eq '2dom' then begin
  domtag='d01'
endif
casedir=dirs.scdir+expname+'/WRF/experiment/'+idir+'/'
dirs.figdir+=expname+'/'
config_wrfexp, casedir=casedir,cases=cases,dirs=dirs,$
  dims=dims, vars=vars, nfils_set=nfils_set, domtag=domtag;, /verbose

;----OB DIRECTORIES--------------------

  maindir=dirs.scdir+'myanmar/'
;  imergfil=maindir+'imerg/data/jjas_2013-2017_diurncomp_3B-HHR.MS.MRG.3IMERG.V06B.nc4';'imerg/20150620-20150629.nc4'
;  imergfil=maindir+'imerg/data/jjas_2013-2017_diurncomp_3B-HHR.MS.MRG.3IMERG.V06B_2.nc4';'imerg/20150620-20150629.nc4'
  npd_imerg=48
  ;era_fil=maindir+'era5/jjas_2013-2017/ERA5-JJAS_2013-2017_mean-pl.nc';'era5/ERA5-20150623-20150630-uv.nc'
  ;USING MET_EM FILES FOR DIURNAL COMPOSITE ERA5
  era_dir=dirs.wkdir+'wrf_wps/wrfv4/output_myanmar/JJAS_2013-2017/era5/'
;  npd_era=24

;Local SOLAR time conversion
local=6;round(mean(dims.lon)/360.*24.) ; deg lon --> hours
print,'Adding +'+strtrim(local,2)+' for LT'

;----READ OBS--------------------

  ;SELECTED TIME ARRAY
    time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
      final=julday(mm_plot[1],dd_plot[1],yy_plot[1],23,59,59),step_size=30,units='Minutes')
    nt=n_elements(time)
    nd=nt/npd_imerg

  ;USE FULL JJAS13-17 TIME ARRAY TO FIND SELECTION INDICES
;    tim_jjas=timegen(start=julday(6,1,2013,0,0,0),final=julday(9,30,2013,23,59,59),step_size=30,units='Minutes')
;    for iy=0,2017-2013-1 do $
;      tim_jjas = [ tim_jjas , timegen(start=julday(6,1,2013+iy+1,0,0,0),final=julday(9,30,2013+iy+1,23,59,59),step_size=30,units='Minutes') ]

  ;NEED TO READ DATA FROM HOURLY FILES THAT CONTAIN ENTIRE JJAS13-17 DATASET

  spawn,'ls '+maindir+'imerg/data/merge_*nc4',imfils

  ;DIMS OF HOURLY FILES
  nx=454 & ny=354

  rain=fltarr(nx,ny,nt)
  for ih=0,npd_imerg-1 do begin
    imtim=read_nc_var(imfils[ih],'time')
    imtim=julday(1,1,1970,0,0,temporary(imtim))
    ;READ ENTIRE SET OF DAYS
      t_read0=time[ih]
      it0=(where(abs(imtim-t_read0) lt 1e-9))[0]
;      caldat,imtim[it0],mm,dd,yy,hh,nn
;      print,mm,dd,yy,hh,nn
      count=[ny,nx,nd] & offset=[0,0,it0] ; y,x,t
      tmp=reform(read_nc_var(imfils[ih],'precipitationCal',count=count,offset=offset)) ; keep in mm/hr
      tmp2=transpose(tmp)
      for id=0,nd-1 do rain[*,*,id*npd_imerg+ih]=reform(tmp2[id,*,*])
  endfor

  if dcomp then begin
    rain2=fltarr(nx,ny,npd_imerg)
    for ih=0,npd_imerg-1 do $
      rain2[*,*,ih] = mean(rain[*,*,indgen(nd)*npd_imerg+ih],dimension=3,/nan,/double)
    rain=rain2 & rain2=0
  endif

  ;IMERG LAT/LON
    lonim=read_nc_var(imfils[0],'lon')
    latim=read_nc_var(imfils[0],'lat')

  ;ERA5 WINDS

;  eralon=dims.lon
;  eralat=dims.lat
;
;  if iwind then begin
;
;    ;PLEVS
;    ilevera=925;[925, 950, 975, 1000]
;    spawn,'ls '+era_dir+'met_em.'+domtag+'*',era_fil
;    era_fil=era_fil[0:23]
;    p_era=reform(read_nc_var(era_fil[0],'PRES',count=[1,1,38,1],offset=[0,0,0,0])) ; Pa
;    izlev_era=(where(p_era*1d-2 eq psel))[0]
;
;    ;READ DIURNAL COMPOSITE
;    u=fltarr(dims.nx,dims.ny,npd_era)
;    v=u
;
;    for it=0,npd_era-1 do begin
;      fil=era_fil[it]
;      count=[dims.nx+1,dims.ny,1,1] & offset=[0,0,izlev_era,0] ; x, y, p, t
;      iu=reform(read_nc_var(fil,'UU',count=count,offset=offset))
;
;    ;DESTAGGER
;      ixin=indgen(dims.nx+1) & ixout=findgen(dims.nx)+0.5
;      for iy=0,dims.ny-1 do u[*,iy,it]=interpol(reform(iu[*,iy]),ixin,ixout)
;      count=[dims.nx,dims.ny+1,1,1] & offset=[0,0,izlev_era,0] ; x, y, p, t
;      iv=reform(read_nc_var(fil,'VV',count=count,offset=offset))
;      iyin=indgen(dims.ny+1) & iyout=findgen(dims.ny)+0.5
;      for ix=0,dims.nx-1 do v[ix,*,it]=interpol(reform(iv[ix,*]),iyin,iyout)
;    endfor
;
;  endif

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
  cross=cross_diag(rain,lonim,latim,width,x_bounds=xcross,y_bounds=ycross,lonout=xlon_im,latout=xlat_im)
  rain=mean(cross,dimension=2,/nan,/double) & cross=0

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
  ;ADJUST BY COASTLINE LOCATION
    ixcheck=indgen(nxhov/2)+nxhov/2
    ixcoast=(where(lsmask[ixcheck] gt 0))[0]
    xhov_era=reverse(xhov_era)
    xcoast=xhov_era[ixcheck[ixcoast]]
    xhov_era-=xcoast
  ;IMERG: DISTANCE FROM LAT/LON
    nxhov=n_elements(xlon_im)
    xhov_im=fltarr(nxhov)
    for ix=0,nxhov-1 do $
      xhov_im[ix]=111.*sqrt( (xlon_im[ix]*cos(xlat_im[ix]*!dtor))^2 + xlat_im[ix]^2 ) ; Converts to km
    xhov_im=reverse(xhov_im)
    xhov_im-=xcoast

;----INDEX OF OFFSHORE RAIN--------------------

  ;BOUNDS FOR SELECTING ONSHORE V OFFSHORE RAINFALL
  ;USE COASTAL DISTANCE VARIABLE
    bounds=[-500,0,1000] ; km

  ;DAILY-MEAN RAINFALL
    days=indgen(nd)
    rain_dm=fltarr(nxhov,nd)
    ind=indgen(npd_imerg)
    for id=0,nd-1 do rain_dm[*,id]=mean(rain[*,id*npd_imerg+ind],dimension=2,/nan,/double)
    rain_dom_mn=mean(rain_dm,dimension=1,/nan,/double)

  ;STANDARDIZED RAIN INDEX
;    rain_dom_mn = (temporary(rain_dom_mn) - mean(rain_dom_mn,/nan,/double)) / stddev(rain_dom_mn)

  ;INDICES
    ind_off=where(xhov_im ge 0)
    ioffshore = mean(rain_dm[ind_off,*],dimension=1,/nan,/double); / rain_dom_mn
    ind_on=where(xhov_im lt 0)
    ionshore = mean(rain_dm[ind_on,*],dimension=1,/nan,/double); / rain_dom_mn
    ;STANDARDIZE
      ioffshore = (ioffshore - mean(ioffshore,/nan,/double)) / stddev(ioffshore,/nan,/double)
      ionshore = (ionshore - mean(ionshore,/nan,/double)) / stddev(ionshore,/nan,/double)

;----CREATE PLOTS--------------------

  var_str='RAINNC'
    setmin=0
    setmax=3.;1.5

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

    form2='(i2.2)'
    form4='(i4)'
    dat_str=string(mm_plot[0],format=form2)+string(dd_plot[0],format=form2)+strmid(strtrim(yy_plot[0],2),2,2)+'-'+$
            string(mm_plot[1],format=form2)+string(dd_plot[1],format=form2)+strmid(strtrim(yy_plot[1],2),2,2)

    if dcomp then begin

      figname=figdir+var_str+'_dcomp_imerg_'+dat_str

      ltim=findgen(npd_imerg)*24/npd_imerg+local
      for it=0,npd_imerg-1 do ltim[it]-=24*(ltim[it] ge 24)

      ;REORDER TO PUT 0 LT AT FIRST INDEX
      nshift = local * npd_imerg/24
      ltim=shift(temporary(ltim),nshift)
      rain=shift(temporary(rain),[0,nshift])

      if iwind then begin
        u=shift(temporary(u),[0,nshift])
        v=shift(temporary(v),[0,nshift])
      endif

    endif else $
      figname=figdir+var_str+'_imerg_'+dat_str

    figspecs.figname=figname
    print,figname

    if iwind then begin
      ltim_era=findgen(npd_era)*24/npd_era+local
      for it=0,npd_era-1 do ltim_era[it]-=24*(ltim_era[it] ge 24)
      ltim_era=shift(temporary(ltim_era),local*npd_era/24)
      u=shift(temporary(u),[0,local*npd_era/24])
      v=shift(temporary(v),[0,local*npd_era/24])
    endif

    if iwind then wind=create_struct('u',iu,'v',iv)

    keepopen=1

    if dcomp then $
      wrf_monsoon_dcomp_hov, dirs, figspecs, rain, ltim, xhov_im, wind=wind, cvar=cvar $
    else $
      wrf_monsoon_hov, dirs, figspecs, rain, time, npd_imerg, xhov_im, wind=wind, cvar=cvar, keepopen=keepopen

    if keyword_set(keepopen) then begin
      csize=0.3

;      axis,xaxis=1,xstyle=1,xminor=1,charsize=csize,xrange=[5,0],/save
;      oplot,rain_dom_mn,days,linestyle=0,thick=2,color=0

      axis,xaxis=1,xstyle=1,xminor=1,charsize=csize,xrange=[3,-3],/save
      plots,[0,0],!y.crange,linestyle=0,thick=1,color=0
      oplot,ioffshore,days,linestyle=2,thick=2,color=0
      oplot,ionshore,days,linestyle=1,thick=2,color=0

      device,/close
      convert_png,figspecs.figname,res=300,/remove_eps
    endif

print,'DONE!!'
end
