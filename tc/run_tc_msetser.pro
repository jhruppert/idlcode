; 
; TC-CENTERED DOMAIN SPECS
; 
pro tcloc_indices,dims,tcloc,nxread,nyread,keepx,keepy,keepx2,keepy2,nxread2,nyread2
;    keepx=where((dims.lon gt loc[0]-dx*0.5) and (dims.lon le loc[0]+dx*0.5),complement=nanx,nxread)
;    keepy=where((dims.lat lt loc[1]+dy*0.5) and (dims.lat ge loc[1]-dy*0.5),complement=nany,nyread)
  ;TC-CENTERED WITH CONSTANT NX,NY
    xdiff=abs(dims.lon-tcloc[0])
    locix=(where(xdiff eq min(xdiff)))[0]
    ydiff=abs(dims.lat-tcloc[1])
    lociy=(where(ydiff eq min(ydiff)))[0]
    keepx=indgen(nxread)+locix-round(0.5*nxread)
    keepy=indgen(nyread)+lociy-round(0.5*nyread)
    ;FOR IDENTIFYING OUT-OF-BOUNDS
      keepx2=where((keepx ge 0) and (keepx lt dims.nx),nxread2)
      keepy2=where((keepy ge 0) and (keepy lt dims.ny),nyread2)
end
; 
; 
; MSE budget time series from Cartesian TC output, TC-following
; 
; 
; James Ruppert
; 3/12/19
; 
pro run_tc_msetser

tcname='maria'
;tcname='haiyan'

if tcname eq 'maria' then begin
  tcyear='2017'
  hurdat=read_hurdat(tcname,tcyear)
endif else if tcname eq 'haiyan' then begin
  tcyear='2013'
  hurdat=read_jtwcdat(tcname) 
endif

;subdir='static_nest/'+tcname
subdir='redux/'+tcname
tc_sim_config, subdir, dirs=dirs, dims=dims, vars=vars;, /verbose
dirs.figdir+=tcname+'/'

do_smooth=1
  idaily=0 ; do daily-means? Only works with do_smooth=1

do_ddt=0 ; calculate DDT term?

;BOX SIZE FOR MSE CALCULATIONS
;  box_size=5.0;10. ; degrees
  box_size=10. ; Using 10 degrees for paper

;VORTEX TRACKING LEVEL
  psel=700 ; (hPa) level to locate vortex
  izsel=(where(dims.pres eq psel))[0]

;----TIME ARRAYS--------------------

  time=dims.time
  nt=dims.nt-1;n_elements(time)-1
  npd=dims.npd
  nhrs=1.*nt*npd/24.
  nd=(nhrs-(nhrs mod 24))/24.
  t_ind_sav=indgen(nt)

;----BEST TRACK SUBSET--------------------

if keyword_set(hurdat) then begin

  subset=where( (hurdat.jultim ge time[0] and hurdat.jultim le max(time)) and $
                (hurdat.lon ge dims.lon[0] and hurdat.lon le max(dims.lon)) and $
                (hurdat.lat ge dims.lat[0] and hurdat.lat le max(dims.lat)) ,nthd)
  hdtim=hurdat.jultim[subset]
  caldat,hdtim[0],mm,dd,yy,hh
  caldat,time[0],mm,ddt,yy,hht
  hdif=24*(dd-ddt)+hh-hht
  hdtim=indgen(nthd)*6+hdif
  hdwspd=hurdat.wspd[subset] * 0.51444 ; knots to m/s
  hdpsfc=hurdat.pres[subset]

endif

;----READ VARS--------------------

  g=9.81
  cp=1004.

  ;LAND MASK
    vtag='LANDMASK'
    lsfil=dirs.files_raw[0,2,0]
    mask=reform(read_nc_var(lsfil,'LANDMASK'))
    land_big=where(mask eq 1,nland_big)
    mask=0

hvar=fltarr(dirs.nc,nt)
hvar[*]=!values.f_nan
hpsefp=hvar
hplwp=hvar
hpswp=hvar
halfddt=hvar

for ic=0,dirs.nc-1 do begin
;for ic=0,0 do begin

  print,'CASE: ',dirs.cases[ic]

  t_ind=t_ind_sav
  if dirs.cases[ic] eq 'icrf_rst' then t_ind+=3*24+1 else $
  if ( strmatch(dirs.cases[ic],'*36h*') or strmatch(dirs.cases[ic],'icrf_*') or strmatch(dirs.cases[ic],'lwcr*') $
    or (dirs.cases[ic] eq 'lwswcrf') or (dirs.cases[ic] eq 'axisym') ) then t_ind+=36
  if strmatch(dirs.cases[ic],'*24h*') then t_ind+=24
  if strmatch(dirs.cases[ic],'*48h*') then t_ind+=48
  if strmatch(dirs.cases[ic],'*60h*') then t_ind+=60
  if strmatch(dirs.cases[ic],'*72h*') then t_ind+=72
  if strmatch(dirs.cases[ic],'*84h*') then t_ind+=84
  if strmatch(dirs.cases[ic],'*96h*') then t_ind+=96

  i_nt=nt-t_ind[0]
  t_ind=indgen(i_nt)+t_ind[0]
  it0=0
;  if dirs.cases[ic] eq 'ctl' then it0=24

  ;3D FILES
  ;T
    iv=where(strmatch(reform(dirs.files_post[ic,*]),'*/T.*'))
    ifil=dirs.files_post[ic,iv]
    tfil=ncdf_open(ifil,/nowrite)
  ;QV
    iv=where(strmatch(reform(dirs.files_post[ic,*]),'*QVAPOR*'))
    ifil=dirs.files_post[ic,iv]
    qfil=ncdf_open(ifil,/nowrite)
  ;QI
    iv=where(strmatch(reform(dirs.files_post[ic,*]),'*QICE*'))
    ifil=dirs.files_post[ic,iv]
    qifil=ncdf_open(ifil,/nowrite)
  ;QS
    iv=where(strmatch(reform(dirs.files_post[ic,*]),'*QSNOW*'))
    ifil=dirs.files_post[ic,iv]
    qsfil=ncdf_open(ifil,/nowrite)
  ;QG
    iv=where(strmatch(reform(dirs.files_post[ic,*]),'*QGRAUP*'))
    ifil=dirs.files_post[ic,iv]
    qgfil=ncdf_open(ifil,/nowrite)
  ;LW
    iv=where(strmatch(reform(dirs.files_post[ic,*]),'*RTHRATLW.*'))
    ifil=dirs.files_post[ic,iv]
    lwfil=ncdf_open(ifil,/nowrite)
  ;SW
    iv=where(strmatch(reform(dirs.files_post[ic,*]),'*RTHRATSW.*'))
    ifil=dirs.files_post[ic,iv]
    swfil=ncdf_open(ifil,/nowrite)

  ;VORTEX TRACKING
    ;READ ABSOLUTE VORTICITY
      iv=where(vars.vars eq 'AVOR')
      file=dirs.files_post[ic,iv]
      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,izsel,0] ; x,y,z,t
      avor=reform(read_nc_var(file,'AVOR',count=count,offset=offset))
    ;SMOOTH
      ixsmooth=round(111./3) ; 1-degree smoothing, run twice
      ismooth=[ixsmooth,ixsmooth,0]
      for i=1,2 do $
        avor=smooth(temporary(avor),ismooth,/edge_truncate,/nan)
    ;MASKING
      if nland_big gt 0 then avor[land_big]=!values.f_nan
      if tcname eq 'maria' then avor=wrf_maria_mask(temporary(avor),time[t_ind],hurdat,dims); else stop
      if tcname eq 'haiyan' then avor=wrf_haiyan_mask(temporary(avor),time[t_ind],hurdat,dims); else stop
    ;VORTEX TRACKING
      vloc=maria_vortex_locate(avor,dims);,/write)
      avor=0
    ;SMOOTH TRACKS
      for i=0,1 do $
        vloc=smooth(temporary(vloc),[0,3],/edge_truncate,/nan)

  ;LOAD FULL 2D VARS

  ;SEF
    count=[dims.nx,dims.ny,1,i_nt-it0] & offset=[0,0,0,it0] ; x,y,z,t
    iv=where(strmatch(reform(dirs.files_post[ic,*]),'*HFX*'))
    file=dirs.files_post[ic,iv]
    sef_all=reform(read_nc_var(file,'HFX',count=count,offset=offset)) ; W/m2
    iv=where(strmatch(reform(dirs.files_post[ic,*]),'*LH*'))
    file=dirs.files_post[ic,iv]
    sef_all+=reform(read_nc_var(file,'LH',count=count,offset=offset)) ; W/m2

  dx=1.*box_size
  dy=1.*box_size
  ;FOR CONSTANT NX,NY
    nxread=round(dx/(dims.lon[1]-dims.lon[0]))
    nyread=round(dy/(dims.lat[1]-dims.lat[0]))

  ;SAVE MSE FOR DDT
    hsav=fltarr(nxread,nyread,i_nt)
    hsav[*]=!values.f_nan

  ;FOR 3D VARS, ALL CALCULATIONS, LOOP OVER TIME
  for it=it0,i_nt-1 do begin
;  for it=96,i_nt-1 do begin
;  for it=44,44 do begin

  ;TC-CENTERED DOMAIN SPECS
    loc=reform(vloc[*,it])
    tcloc_indices,dims,loc,nxread,nyread,keepx,keepy,keepx2,keepy2,nxread2,nyread2

  ;LAND-SEA MASK
    lsfil=dirs.files_raw[0,2,0]
    count=[nxread2,nyread2,1] & offset=[keepx[keepx2[0]],keepy[keepy2[0]],0] ; x,y,z,t
    mask=reform(read_nc_var(lsfil,'LANDMASK',count=count,offset=offset))
    land=where(mask eq 1,nland,complement=sea)

  ;SEF
    sef=reform(sef_all[keepx[keepx2],*,it]) ; W/m2
    sef=sef[*,keepy[keepy2]]

  ;READ 3D VARS
    count=[nxread2,nyread2,dims.np,1] & offset=[keepx[keepx2[0]],keepy[keepy2[0]],0,it] ; x,y,z,t
  ;T
    ncdf_varget,tfil,1,tmpk,count=count,offset=offset ; K
  ;QV
;    qv=reform(read_nc_var(qfil,'QVAPOR',count=count,offset=offset)) ; kg/kg
    ncdf_varget,qfil,1,qv,count=count,offset=offset ; kg/kg
    ;CORRUPT VALUES
      inan=where(abs(qv) ge 1e5,nnan)
      if nnan gt 0 then qv[inan]=!values.f_nan
  ;QI
    ncdf_varget,qifil,1,qi,count=count,offset=offset ; kg/kg
  ;QS
    ncdf_varget,qsfil,1,qs,count=count,offset=offset ; kg/kg
  ;QG
    ncdf_varget,qgfil,1,qg,count=count,offset=offset ; kg/kg
  ;LW
    ncdf_varget,lwfil,1,lw,count=count,offset=offset ; K/s
    lw*=cp ; K/s --> J/kg/s
  ;SW
    ncdf_varget,swfil,1,sw,count=count,offset=offset ; K/s
    sw*=cp ; K/s --> J/kg/s

  ;ADD SNOW, GRAUPEL TO CLOUD ICE
  qi=temporary(qi)+qs+qg
  qs=0 & qg=0

  ;RHO
    tvirt = tmpk*(1.+0.61*qv)
    rho=tvirt
    for iz=0,dims.np-1 do rho[*,*,iz] = dims.pres[iz]*1e2 / ( 287. * tvirt[*,*,iz] )
    ;tvirt=0

  ;HEIGHT FROM HYDROSTATIC
    z = fltarr(nxread2,nyread2,dims.np) ; m
    for iz=1,dims.np-1 do $
      z[*,*,iz] = z[*,*,iz-1] + (dims.pres[iz-1]-dims.pres[iz])*1e2 / 9.81 / mean(reform(rho[*,*,iz-1:iz]),dimension=3,/double,/nan)
    ;rho=0

  ;REMOVE LAND
    if nland gt 0 then begin
      sef[land]=!values.f_nan
      maskz=fltarr(nxread2,nyread2,dims.np)
      for iz=0,dims.np-1 do maskz[*,*,iz]=mask
      landz=where(maskz eq 1)
      tmpk[landz]=!values.f_nan
      qv[landz]=!values.f_nan
      qi[landz]=!values.f_nan
      lw[landz]=!values.f_nan
      sw[landz]=!values.f_nan
      rho[landz]=!values.f_nan
      z[landz]=!values.f_nan
    endif

  ;MSE
    h = calc_mse_frozen(1d*qv,1d*qi,1d*tmpk,1d*z)
;    tmpk=0 & qv=0 & qi=0 & z=0

  ;VERTICALLY INTEGRATE
    dp=(dims.pres[2]-dims.pres[1])*(-1e2)
    iplev=where(dims.pres le 950)
    hint=total(h[*,*,iplev],3,/double)*dp/g
    h=hint ;& hint=0

  hsav[keepx2,keepy2,it]=h

;    lwint=total(lw[*,*,iplev],3,/double)*dp/g
    lwint=total(lw[*,*,*,*],3,/double)*dp/g
    lw=lwint ;& lwint=0
    swint=total(sw[*,*,*,*],3,/double)*dp/g
    sw=swint ;& lwint=0

  ;COSINE WEIGHTING FOR MEANS
    wgt=fltarr(nxread2,nyread2)
    ctmp=1d*cos(dims.lat[keepy[keepy2]]*!dtor)
    ctmp/=mean(ctmp,/double)
    for ix=0,nxread2-1 do wgt[ix,*]=ctmp

  ;VARIANCE TERMS
    hp   = h   *1d
    sefp = sef *1d
    lwp  = lw  *1d
    swp  = sw  *1d

    hp   -=mean(h*wgt,/nan,/double)
    sefp -=mean(sef*wgt,/nan,/double)
    lwp  -=mean(lw*wgt,/nan,/double)
    swp  -=mean(sw*wgt,/nan,/double)

    hv      = hp*hp
    ihpsefp = hp*sefp
    ihplwp  = hp*lwp
    ihpswp  = hp*swp

    hvar[ic,t_ind[it]]   = mean(hv*wgt,/nan,/double)
    hpsefp[ic,t_ind[it]] = mean(ihpsefp*wgt,/nan,/double)
    hplwp[ic,t_ind[it]]  = mean(ihplwp*wgt,/nan,/double)
    hpswp[ic,t_ind[it]]  = mean(ihpswp*wgt,/nan,/double)

  endfor ; it

  ;CLOSE FILES
  ncdf_close,tfil
  ncdf_close,qfil
  ncdf_close,qifil
  ncdf_close,lwfil
  ncdf_close,swfil

  ;CALCULATE DDT TERM

  if do_ddt then begin

    print,'Doing DDT'

    ;DDT
    ddt=dblarr(nxread,nyread,i_nt)
    ddt[*]=!values.f_nan
    for ix=0,nxread-1 do for iy=0,nyread-1 do begin
      ih=reform(hsav[ix,iy,*])
      ifin=where(finite(ih))
      ddt[ix,iy,ifin]=deriv(ih[ifin])
    endfor
    ddt*=1d/3600d ; per time step --> per sec

    ;REMOVE MEAN AND MULTIPLY BY H'
    for it=it0,i_nt-1 do begin
;    for it=60,60 do begin
      ;TC LOC
        loc=reform(vloc[*,it])
        tcloc_indices,dims,loc,nxread,nyread,keepx,keepy,keepx2,keepy2,nxread2,nyread2
      ;COSINE WEIGHTING FOR MEANS
        wgt=fltarr(nxread2,nyread2)
        ctmp=1d*cos(dims.lat[keepy[keepy2]]*!dtor)
        ctmp/=mean(ctmp,/double)
        for ix=0,nxread2-1 do wgt[ix,*]=ctmp
      iddt=reform(ddt[keepx2,keepy2,it])
      ih=reform(hsav[keepx2,keepy2,it])
      iddt-=mean(iddt*wgt,/double,/nan)
      ih-=mean(hp*wgt,/nan,/double)
      iddt*=ih
      halfddt[ic,t_ind[it]]=mean(iddt*wgt,/nan,/double)
    endfor

  endif

;fname=dirs.casedir[ic]+'post/msevar.nc'
;write_sing_ncvar,fname,hvar,'msevar',dimtag1='lon',dimtag2='lat',dimtag3='time'
;exit

  ;FIRST TIME STEP FOR H'LW' IS VERY DIFFERENT FROM THE REST (IS TAKEN FROM CTL)
  if dirs.cases[ic] ne 'ctl' and dirs.cases[ic] ne 'lwcrf' then begin
    hplwp[ic,t_ind[0]]=!values.f_nan
    hpswp[ic,t_ind[0]]=!values.f_nan
  endif

  ;FILL IN EACH CASE WITH CTL, THEN OVERWRITE LATER TIME STEPS
  if do_smooth then begin

    if idaily then begin

      npd=24
      nd=(dims.nt-1)/npd
      print,'DAILY-AVERAGING: ',(1.*dims.nt-1)/(1.*npd),' days'

      ;SHIFT INDEX TO NOON
        ;if tcname eq 'maria' then ishift=12 else ishift=0
      ishift=npd/2

      hvar_d=fltarr(dirs.nc,nd) & hvar_d[*]=!values.f_nan
      halfddt_d=hvar_d
      hpsefp_d=hvar_d
      hplwp_d=hvar_d
      hpswp_d=hvar_d

      for id=0,nd-1 do begin
        t_ind_av=indgen(npd)+id*npd
        hvar_d[ic,id]=mean(hvar[ic,t_ind_av],/nan)
        halfddt_d[ic,id]=mean(halfddt[ic,t_ind_av],/nan)
        hpsefp_d[ic,id]=mean(hpsefp[ic,t_ind_av],/nan)
        hplwp_d[ic,id]=mean(hplwp[ic,t_ind_av],/nan)
        hpswp_d[ic,id]=mean(hpswp[ic,t_ind_av],/nan)
      endfor

      hvar=hvar_d
      halfddt=halfddt_d
      hpsefp=hpsefp_d
      hplwp=hplwp_d
      hpswp=hpswp_d

    endif else begin

      if dirs.cases[ic] eq 'ctl' then $
        for icc=1,dirs.nc-1 do begin
          hvar[icc,t_ind]=hvar[0,t_ind]
          hpsefp[icc,t_ind]=hpsefp[0,t_ind]
;          hplwp[icc,t_ind]=hplwp[0,t_ind]
;          hpswp[icc,t_ind]=hpswp[0,t_ind]
        endfor
  
    ;SMOOTH TIME SERIES
      ;ismooth=6;12
      ;FOR GAUSSIAN FILTER
      ismooth=3
      ;for i=0,1 do begin
        hvar[ic,*]=gauss_smooth(reform(hvar[ic,*]),ismooth,/edge_truncate,/nan)
        halfddt[ic,*]=gauss_smooth(reform(halfddt[ic,*]),ismooth,/edge_truncate,/nan)
        hpsefp[ic,*]=gauss_smooth(reform(hpsefp[ic,*]),ismooth,/edge_truncate,/nan)
        hplwp[ic,*]=gauss_smooth(reform(hplwp[ic,*]),ismooth,/edge_truncate,/nan)
        hpswp[ic,*]=gauss_smooth(reform(hpswp[ic,*]),ismooth,/edge_truncate,/nan)
      ;endfor
  
    ;NOW REPLACE NANS FOR TIMES PRIOR TO TEST START
      if dirs.cases[ic] ne 'ctl' then begin
        hvar[ic,0:t_ind[0]-1]=!values.f_nan
        hpsefp[ic,0:t_ind[0]-1]=!values.f_nan
        hplwp[ic,0:t_ind[0]-1]=!values.f_nan
        hpswp[ic,0:t_ind[0]-1]=!values.f_nan
      endif

    endelse

  endif

endfor ; icase


;----CREATE PLOT--------------------


  figdir=dirs.figdir+'/mse_budget/'
  spawn,'mkdir '+figdir,tmp,tmpe

;stop

wrf_mse_vartser, tcname, hvar, halfddt, hpsefp, hplwp, hpswp, dirs.nc, nhrs, figdir, dirs.cases, idaily

print,'DONE!!'
end
