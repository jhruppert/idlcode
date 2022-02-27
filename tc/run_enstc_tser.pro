; 
; Time series of TC output from ensemble simulations

; James Ruppert
; 1/17/21
; 
pro run_enstc_tser

;tcname='maria'
tcname='haiyan'
case_str='ncrf';'ctl'

;if tcname eq 'maria' then begin
;  tcyear='2017'
;  hurdat=read_hurdat(tcname,tcyear)
;endif else if tcname eq 'haiyan' then begin
;  tcyear='2013'
;  hurdat=read_jtwcdat(tcname)
;endif

dom='d02'
tc_ens_config, tcname, case_str, dom, dirs=dirs, dims=dims, vars=vars;, /verbose
dirs.figdir+=tcname+'/'

type=1 ; 0 - pres,wnd; 1 - single var
  ; if type=1
  var_str='ri';'lh_rad';'w';'ike';'hfx';'lh'

iazim=0 ; use azimuthal wind?
do_smooth=1 ; only applies to type=1

;VORTEX TRACKING LEVEL
  psel=700 ; (hPa) level to locate vortex
  izsel=(where(dims.pres eq psel))[0]

;----OPTIONS FOR AZIM FILES (VTAN)--------------------

nrad=267 & naz=360
;if (tcname eq 'haiyan') or strmatch(subdir,'*redux*') then nrad=433
;WILL ONLY READ OUT TO RADIUS OF 800 KM
radius_range=[0,600] ; max is 1400 km xx800 km

;----TIME ARRAYS--------------------

  time=dims.time
  nt=dims.nt;n_elements(time)-1
;time=dims.time[24:nt-1]
;nt=n_elements(time)
  npd=dims.npd
  nhrs=1.*nt*npd/24.
  nd=(nhrs-(nhrs mod 24))/24.
  t_ind_sav=indgen(nt)
;t_ind_sav=indgen(nt)+24

;----BEST TRACK SUBSET--------------------

if keyword_set(hurdat) then begin

  subset=where( (hurdat.jultim ge time[0] and hurdat.jultim le max(time)) and $
                (hurdat.lon ge dims.lon[0] and hurdat.lon le max(dims.lon)) and $
                (hurdat.lat ge dims.lat[0] and hurdat.lat le max(dims.lat)) ,nthd)
  hdtim=hurdat.jultim[subset]
;  caldat,hdtim[0],mm,dd,yy,hh
  caldat,hdtim,mm,dd,yy,hh
;  caldat,time[0],mm,ddt,yy,hht
;  hdif=24.*(dd-ddt)+hh-hht
;  hdif=(dd-ddt)+hh-hht
;  hdtim=findgen(nthd)*6/24+dd+hh/24;hdif/24
  hdtim=dd+hh/24.
  hdwspd=hurdat.wspd[subset] * 0.51444 ; knots to m/s
  hdpsfc=hurdat.pres[subset]

endif

;----READ VARS--------------------


  ;LAND MASK
    vtag='LANDMASK'
;    file=dirs.files_raw[0,2,0]
    file=dirs.files_raw[0,2,0]
    mask=reform(read_nc_var(file,'LANDMASK'))
    land=where(mask eq 1,nland)
    mask=0

var1=fltarr(dirs.nc,nt)
var1[*]=!values.f_nan
var2=var1
var2car=var1

for ic=0,dirs.nc-1 do begin
;for ic=0,1 do begin

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
;print,t_ind
;continue
  if type eq 0 then begin

  ;MIN PRESSURE
    vtag='SLP';'PSFC'
;    iv=where(vars.vars eq vtag,count)
    iv=where(strmatch(reform(dirs.files_post[ic,*]),'*'+vtag+'*'),count)
    if count eq 0 then message,'File not found!'
    file=dirs.files_post[ic,iv]
    count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
    psfc=reform(read_nc_var(file,vtag,count=count,offset=offset))

    ;MASKING
;      if nland gt 0 then psfc[land]=!values.f_nan
      imask=fltarr(dims.nx,dims.ny,i_nt)
      for it=0,i_nt-1 do imask[*,*,it]=mask
      land=where(imask eq 1,nland)
      if nland gt 0 then psfc[land]=!values.f_nan

      if tcname eq 'maria' then psfc=wrf_maria_mask(temporary(psfc),time[t_ind],hurdat,dims)
      if tcname eq 'haiyan' then psfc=wrf_haiyan_mask(temporary(psfc),time[t_ind],hurdat,dims)

    psfc=min(min(temporary(psfc),dimension=1,/nan),dimension=1,/nan)
;    psfc*=1e-2

  ;SMOOTH TIME SERIES
    ismooth=3
    for i=0,1 do $
      psfc=smooth(temporary(psfc),ismooth,/edge_truncate,/nan)

    var1[ic,t_ind]=psfc

  if iazim then begin
  ;USE AZIMUTHAL WIND

    hr_fil=strtrim(t_ind[0],2)+'-'+strtrim(t_ind[i_nt-1],2)+'hr'
    file=dirs.casedir[ic]+'azim_U10_'+hr_fil+'.nc'

    radius=read_nc_var(file,'radius')
    azimuth=read_nc_var(file,'azmiuth')

    ilevwnd=0
    if ilevwnd eq 0 then begin
    ;10m
      ititle='Max avg(v!Dtan!N), 10 m' ; Title for wind calculation type
      file=dirs.casedir[ic]+'azim_U10_'+hr_fil+'.nc'
      u=reform(read_nc_var(file,'U10',count=[nrad,naz,1,i_nt],offset=[0,0,0,0]))
      file=dirs.casedir[ic]+'azim_V10_'+hr_fil+'.nc'
      v=reform(read_nc_var(file,'V10',count=[nrad,naz,1,i_nt],offset=[0,0,0,0]))
    endif else if ilevwnd eq 1 then begin
    ;PRESSURE LEVEL
      ititle='Max avg(v!Dtan!N), 1000 hPa' ; Title for wind calculation type
      psel=1000;950;700 ; (hPa)
      izsel=(where(dims.pres eq psel))[0]
      file=dirs.casedir[ic]+'azim_U_'+hr_fil+'.nc'
      u=reform(read_nc_var(file,'U',count=[nrad,naz,1,i_nt],offset=[0,0,izsel,0]))
      file=dirs.casedir[ic]+'azim_V_'+hr_fil+'.nc'
      v=reform(read_nc_var(file,'V',count=[nrad,naz,1,i_nt],offset=[0,0,izsel,0]))
    endif else if ilevwnd eq 2 then begin
    ;LOWEST MODEL LEVEL
      ititle='Max avg(v!Dtan!N), Lowest modlev' ; Title for wind calculation type
      file=dirs.casedir[ic]+'azim_UL_'+hr_fil+'.nc'
      u=reform(read_nc_var(file,'UL',count=[nrad,naz,1,i_nt],offset=[0,0,0,0]))
      file=dirs.casedir[ic]+'azim_VL_'+hr_fil+'.nc'
      v=reform(read_nc_var(file,'VL',count=[nrad,naz,1,i_nt],offset=[0,0,0,0]))
    endif

    wnd_azim=azim_wind_conv(u,v,azimuth) & u=0 & v=0
    wspd=wnd_azim.v_tan
    wnd_azim=0

  ;AZIMUTHALLY AVERAGE
    wspd=mean(temporary(wspd),dimension=2,/nan,/double)
  ;RADIUS SUBSET
    irad=where((radius ge radius_range[0]) and (radius le radius_range[1]) , nradtmp)
    wspd=wspd[irad,*]

    wspd=max(temporary(wspd),dimension=1,/nan)

  ;SAVE MAX CARTESIAN WIND
      file=dirs.casedir[ic]+'post/U10.nc'
      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
      u=reform(read_nc_var(file,'U10',count=count,offset=offset))
      file=dirs.casedir[ic]+'post/V10.nc'
      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
      v=reform(read_nc_var(file,'V10',count=count,offset=offset))
    wspdcar=sqrt(u^2+v^2)
    u=0 & v=0
    if nland gt 0 then wspdcar[land]=!values.f_nan
    if tcname eq 'maria' then wspdcar=wrf_maria_mask(temporary(wspdcar),time[t_ind],hurdat,dims)
    wspdcar=max(max(temporary(wspdcar),dimension=1,/nan),dimension=1,/nan)
    var2car[ic,t_ind]=wspdcar

  endif else begin
  ;USE REGULAR WIND

  ;MAX WIND
    ilevwnd=0
    if ilevwnd eq 0 then begin
      ititle='Max (u,v), 10 m' ; Title for wind calculation type
      file=dirs.casedir[ic]+'post/U10.nc'
      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
      u=reform(read_nc_var(file,'U10',count=count,offset=offset))
      file=dirs.casedir[ic]+'post/V10.nc'
      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
      v=reform(read_nc_var(file,'V10',count=count,offset=offset))
    endif else if ilevwnd eq 1 then begin
    ;PRESSURE LEVEL
      ititle='Max (u,v), 1000 hPa' ; Title for wind calculation type
      psel=1000;850;925 ; hPa
      iz=where(dims.pres eq psel)
      file=dirs.casedir[ic]+'post/U.nc'
      u=reform(read_nc_var(file,'U',count=count,offset=[0,0,iz,0]))
      file=dirs.casedir[ic]+'post/V.nc'
      v=reform(read_nc_var(file,'V',count=count,offset=[0,0,iz,0]))
    endif else if ilevwnd eq 2 then begin
    ;LOWEST MODEL LEVEL
      ititle='Max (u,v), Lowest modlev' ; Title for wind calculation type
      file=dirs.casedir[ic]+'post/UL.nc'
      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
      u=reform(read_nc_var(file,'UL',count=count,offset=offset))
      file=dirs.casedir[ic]+'post/VL.nc'
      v=reform(read_nc_var(file,'VL',count=count,offset=offset))
    endif

    wspd=sqrt(u^2+v^2)
    u=0 & v=0

    ;MASKING
      if nland gt 0 then wspd[land]=!values.f_nan
      if tcname eq 'maria' then wspd=wrf_maria_mask(temporary(wspd),time[t_ind],hurdat,dims)

    wspd=max(max(temporary(wspd),dimension=1,/nan),dimension=1,/nan)

  endelse

    print,'Max/min wind speed:'
    stats,wspd

  ;SMOOTH TIME SERIES
;    ismooth=3
;    for i=0,1 do $
;      wspd=smooth(temporary(wspd),ismooth,/edge_truncate,/nan)

    var2[ic,t_ind]=wspd

  endif else begin

;    if var_str eq 'ike' or var_str eq 'w' then begin
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
          if nland gt 0 then avor[land]=!values.f_nan
          if tcname eq 'maria' then avor=wrf_maria_mask(temporary(avor),time[t_ind],hurdat,dims); else stop
          if tcname eq 'haiyan' then avor=wrf_haiyan_mask(temporary(avor),time[t_ind],hurdat,dims); else stop
        ;VORTEX TRACKING
          vloc=maria_vortex_locate(avor,dims);,/write)
          avor=0
;    endif

    if var_str eq 'hfx' then begin
    ;TOTAL SURFACE FLUXES

      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
      iv=where(strmatch(reform(dirs.files_post[ic,*]),'*HFX*'),count)
      file=dirs.files_post[ic,iv]
      sef=reform(read_nc_var(file,'HFX',count=count,offset=offset))
      iv=where(strmatch(reform(dirs.files_post[ic,*]),'*LH*'),count)
      file=dirs.files_post[ic,iv]
      sef+=reform(read_nc_var(file,'LH',count=count,offset=offset))

      ;SPECS OF CALCULATION
        radius_thresh = 400 ; km

      ;SUBSET BY RADIUS
      sef_mean=fltarr(i_nt)
      for it=0,i_nt-1 do begin
        ivar = reform(sef[*,*,it])
        tcrad = radius_tc_ll(reform(vloc[*,it]),dims.lon,dims.lat)
        irad = where(tcrad le radius_thresh,count)
        ivar2 = ivar[irad]
        sef_mean[it] = mean(ivar2,/nan,/double) ; W/m2
      endfor

      var1[ic,t_ind]=sef_mean

    endif else if var_str eq 'ike' then begin
    ;INTEGRATED KINETIC ENERGY (IKE; POWELL AND REINHOLD 2007, BAMS)

      iv=where(strmatch(vars.vars,'U10',/fold_case))
      u=reform(read_nc_var(dirs.files_post[ic,iv],'U10'))
      iv=where(strmatch(vars.vars,'V10',/fold_case))
      v=reform(read_nc_var(dirs.files_post[ic,iv],'V10'))
      wsp=sqrt(u^2+v^2)
      u=0 & v=0

      ;SPECS OF CALCULATION
        radius_thresh = 800 ; km
        ws_thresh = 10.;18. ; TS = 18
        area = 1d*!pi*(radius_thresh*1e3)^2 ; m^2
        scale=1d-18
        area*=scale

      ;CALCULATE INTEGRATED KINETIC ENERGY
      ike=fltarr(i_nt)
      for it=0,i_nt-1 do begin

        ispd = reform(wsp[*,*,it])

        ;SUBSET FOR RADIUS THRESHOLD
          tcrad = radius_tc_ll(reform(vloc[*,it]),dims.lon,dims.lat)
          irad = where(tcrad le radius_thresh,count)
          ispd = ispd[irad]

        ;SUBSET FOR WIND SPEED THRESHOLD
          iwst = where(ispd ge ws_thresh,nwst)
          ispd = ispd[iwst]

        ike[it] = 0.5 * total(ispd^2,/nan,/double) * area
        ; 1/2 * mass * windspeed^2 (1 kg/m3  density and 1 m depth vertical layer implied)
        ; units of Joules (x some factor)

      endfor

      var1[ic,t_ind]=ike

    endif else if var_str eq 'w' then begin
    ;INTEGRATED VERTICAL MOTION xxx AT 500 HPA

      psel=500 ; (hPa)
      izsel=(where(dims.pres eq psel))[0]

      dp=deriv(dims.pres)*(-1e2)
      dpa=dp[0];fltarr(dims.nx,dims.ny,dims.np,i_nt)
;      for iz=0,dims,np-1 do dpa[*,*,iz,*]=dp[iz]

;      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,izsel,0] ; x,y,z,t
      count=[dims.nx,dims.ny,dims.np,i_nt] & offset=[0,0,0,0] ; x,y,z,t
      iv=where(strmatch(vars.vars,'W',/fold_case))
      w=reform(read_nc_var(dirs.files_post[ic,iv],'W',count=count,offset=offset))

      vmf=total(w*dpa,3,/nan,/double)/9.81 ; m/s --> kg/m/s

      ;SPECS OF CALCULATION
        radius_thresh = 300;200;300 ; km

      ;SUBSET BY RADIUS
      w_mean=fltarr(i_nt)
      for it=0,i_nt-1 do begin

;        ivar = reform(w[*,*,it])
        ivar = reform(vmf[*,*,it])

        tcrad = radius_tc_ll(reform(vloc[*,it]),dims.lon,dims.lat)
        irad = where(tcrad le radius_thresh,count)
        ivar2 = ivar[irad]

        w_mean[it] = mean(ivar2,/nan,/double) ; kg/m/s xxx m/s

      endfor

;      var1[ic,t_ind]=w_mean*1e2
      var1[ic,t_ind]=w_mean*1e-2

    endif else if var_str eq 'lh_rad' then begin
    ;LATENT HEAT IN VORTEX REGION

      iv=where(strmatch(vars.vars,'lh',/fold_case))
      file=dirs.files_post[ic,iv]
      ivar=reform(read_nc_var(file,'LH'))

      ;SPECS OF CALCULATION
        radius_thresh = 400 ; km

      ;CALCULATE INTEGRATED KINETIC ENERGY
      for it=0,i_nt-1 do begin

        ilh = reform(ivar[*,*,it])

        ;SUBSET FOR RADIUS THRESHOLD
          tcrad = radius_tc_ll(reform(vloc[*,it]),dims.lon,dims.lat)
          irad = where(tcrad le radius_thresh,count)
          ilh2 = ilh[irad]

        var1[ic,t_ind[it]] = mean(mean(ilh2,/nan,/double,dimension=1),/nan,/double,dimension=1)

      endfor

    endif else if var_str eq 'ri' then begin
    ;WIND SPEED INTENSIFICATION RATE
    ;DERIVATIVE IS TAKEN AFTER SMOOTHING WIND SPEED

      ititle=' '
      if iazim then begin
      ;USE AZIMUTHAL WIND
    
        hr_fil=strtrim(t_ind[0],2)+'-'+strtrim(t_ind[i_nt-1],2)+'hr'
        file=dirs.casedir[ic]+'azim_U10_'+hr_fil+'.nc'
    
        radius=read_nc_var(file,'radius')
        azimuth=read_nc_var(file,'azmiuth')
    
        ilevwnd=0
        if ilevwnd eq 0 then begin
        ;10m
          ititle='Max avg(v!Dtan!N), 10 m' ; Title for wind calculation type
          file=dirs.casedir[ic]+'azim_U10_'+hr_fil+'.nc'
          u=reform(read_nc_var(file,'U10',count=[nrad,naz,1,i_nt],offset=[0,0,0,0]))
          file=dirs.casedir[ic]+'azim_V10_'+hr_fil+'.nc'
          v=reform(read_nc_var(file,'V10',count=[nrad,naz,1,i_nt],offset=[0,0,0,0]))
        endif else if ilevwnd eq 1 then begin
        ;PRESSURE LEVEL
          ititle='Max avg(v!Dtan!N), 1000 hPa' ; Title for wind calculation type
          psel=1000;950;700 ; (hPa)
          izsel=(where(dims.pres eq psel))[0]
          file=dirs.casedir[ic]+'azim_U_'+hr_fil+'.nc'
          u=reform(read_nc_var(file,'U',count=[nrad,naz,1,i_nt],offset=[0,0,izsel,0]))
          file=dirs.casedir[ic]+'azim_V_'+hr_fil+'.nc'
          v=reform(read_nc_var(file,'V',count=[nrad,naz,1,i_nt],offset=[0,0,izsel,0]))
        endif else if ilevwnd eq 2 then begin
        ;LOWEST MODEL LEVEL
          ititle='Max avg(v!Dtan!N), Lowest modlev' ; Title for wind calculation type
          file=dirs.casedir[ic]+'azim_UL_'+hr_fil+'.nc'
          u=reform(read_nc_var(file,'UL',count=[nrad,naz,1,i_nt],offset=[0,0,0,0]))
          file=dirs.casedir[ic]+'azim_VL_'+hr_fil+'.nc'
          v=reform(read_nc_var(file,'VL',count=[nrad,naz,1,i_nt],offset=[0,0,0,0]))
        endif
    
        wnd_azim=azim_wind_conv(u,v,azimuth) & u=0 & v=0
        wspd=wnd_azim.v_tan
        wnd_azim=0
    
      ;AZIMUTHALLY AVERAGE
        wspd=mean(temporary(wspd),dimension=2,/nan,/double)
      ;RADIUS SUBSET
        irad=where((radius ge radius_range[0]) and (radius le radius_range[1]) , nradtmp)
        wspd=wspd[irad,*]
    
        wspd=max(temporary(wspd),dimension=1,/nan)
        var1[ic,t_ind]=wspd

      endif else begin
      ;MAX WIND
        count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
        file=dirs.casedir[ic]+'post/U10.nc'
        u=reform(read_nc_var(file,'U10',count=count,offset=offset))
        count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
        file=dirs.casedir[ic]+'post/V10.nc'
        v=reform(read_nc_var(file,'V10',count=count,offset=offset))
        wspd=sqrt(u^2+v^2)
        u=0 & v=0
    
        ;MASKING
          if nland gt 0 then wspd[land]=!values.f_nan
          if tcname eq 'maria' then wspd=wrf_maria_mask(temporary(wspd),time[t_ind],hurdat,dims)
    
        wspd=max(max(temporary(wspd),dimension=1,/nan),dimension=1,/nan)
        var1[ic,t_ind]=wspd;deriv(wspd)*24. ; m/s / day

      endelse

      if ~do_smooth then var1[ic,*]=deriv(wspd)*24. ; m/s / day

    endif else begin

      iv=where(strmatch(vars.vars,var_str,/fold_case))
      file=dirs.files_post[ic,iv]
      ivar=reform(read_nc_var(file,vars.vars[iv[0]]))
  ;    ivar=stddev(stddev(temporary(ivar),dimension=1,/nan),dimension=1,/nan)
      ;CALCULATE STDDEV
;      for it=0,i_nt-1 do begin
;        tmpv=reform(ivar[*,*,it])
;        var1[ic,t_ind[it]]=stddev(tmpv,/double,/nan)
;      endfor
      var1[ic,t_ind]=mean(mean(ivar,/nan,/double,dimension=1),/nan,/double,dimension=1)
  ;    var1[ic,t_ind]=iv

    endelse

  endelse

  ;FILL IN EACH CASE WITH CTL, THEN OVERWRITE LATER TIME STEPS
  if do_smooth then begin

    if dirs.cases[ic] eq 'ctl' then $
      for icc=1,dirs.nc-1 do begin
        var1[icc,t_ind]=var1[0,t_ind]
        var2[icc,t_ind]=var2[0,t_ind]
      endfor

  ;SMOOTH TIME SERIES
;    ismooth=[0,6];12]
    ;for i=0,1 do begin
;      var1[ic,*]=smooth(var1[ic,*],ismooth,/edge_truncate,/nan)
;      var2[ic,*]=smooth(var2[ic,*],ismooth,/edge_truncate,/nan)
    ;endfor
    if type eq 1 and var_str eq 'ri' then begin
      ;CALCULATE INTENSIFICATION RATE
      ;SMOOTH WINDS FIRST, THEN TAKE DERIVATIVE
    ismooth=6
      var1[ic,*]=gauss_smooth(reform(var1[ic,*]),ismooth,/edge_truncate,/nan)
      var1[ic,*]=deriv(var1[ic,*])*24. ; m/s / day
;      var1[ic,*]=gauss_smooth(reform(var1[ic,*]),3,/edge_truncate,/nan)
    endif else begin
    ismooth=3
      var1[ic,*]=gauss_smooth(reform(var1[ic,*]),ismooth,/edge_truncate,/nan)
      var2[ic,*]=gauss_smooth(reform(var2[ic,*]),ismooth,/edge_truncate,/nan)
    endelse

  ;NOW REPLACE NANS FOR TIMES PRIOR TO TEST START
    if dirs.cases[ic] ne 'ctl' then begin
      var1[ic,0:t_ind[0]-1]=!values.f_nan
      var2[ic,0:t_ind[0]-1]=!values.f_nan
    endif

;INTENSIFICATION RATE
;  inan=where(finite(reform(var2[ic,*])))
;  var2[ic,inan]=deriv(var2[ic,inan])*24. ; m/s / day
;  ;var2[ic,inan]=smooth(reform(var2[ic,inan]),13,/edge_truncate,/nan)
;  var2[ic,inan]=gauss_smooth(reform(var2[ic,inan]),6,/edge_truncate,/nan)

  endif

endfor ; icase


;SMOOTH TIME SERIES
;  ismooth=[0,3];24]
;  for i=0,1 do begin
;    var1=smooth(temporary(var1),ismooth,/edge_truncate,/nan)
;;    var2=smooth(temporary(var2),ismooth,/edge_truncate,/nan)
;  endfor


;----CREATE PLOT--------------------


  if type eq 0 then ftag='pwind' else ftag=var_str
  figname=dirs.figdir+'tser_'+ftag

  ;PLOT SPECS
    csize=0.8
    position=[0.14,0.18,0.89,0.89]
    xsize=4.2 & ysize=2
;    xtitle='Time [ hr ]'
    ytitle1='Min. Pressure [ hPa ]'
    iv2=1
    ytitle2='m s!U-1!N';'Max 10-m Wind [ m s!U-1!N ]'
    yrange=[890,1020]
    if tcname eq 'maria' then yrange[0]=910;920
    yrange2=[0,90];80]
    if tcname eq 'haiyan' then yrange2[1]=90
    title=ititle;strupcase(dirs.cases[ic])
    icplot=indgen(dirs.nc)

    if type eq 1 then begin
      if var_str eq 'lh' then begin
        ytitle1='Stddev(LH) [ W m!U-1!N ]'
        yrange=[0,250]
        ytitle2='Mean(LH) [ W m!U-2!N ]'
        yrange2=[0,250]
;        icplot=indgen(4)
      endif else if var_str eq 'hfx' then begin
        ytitle1='Mean(SEF) [ W m!U-2!N ]'
        yrange=[200,350]
        iv2=0
      endif else if var_str eq 'ike' then begin
        ytitle1=' [ 10!U18!N J ]'
        yrange=[0,15]
        if tcname eq 'haiyan' then yrange[1]=35
        iv2=0
      endif else if var_str eq 'w' then begin
;        ytitle1='[ cm/s ]'
;        yrange=[0,7]
        ytitle1='[ 10!U2!N kg/m/s ]'
        yrange=[0,7];15]
        if tcname eq 'haiyan' then yrange[1]=35
        iv2=0
      endif else if var_str eq 'lh_rad' then begin
        ytitle1='[ W m!U-2!N ]'
        yrange=[125,300]
;        if tcname eq 'haiyan' then yrange[1]=35
        iv2=0
      endif else if var_str eq 'ri' then begin
        ytitle1='m s!U-1!N day!U-1!N'
        yrange=[-10,40]
        iv2=0
      endif
    endif

  ncplot=n_elements(icplot)

  ;AXES
    x=findgen(nhrs)/npd
;xrange=[15.5,20.5]
    y=indgen(10)
    if ~keyword_set(yrange) then $
      yrange=[min(y),max(y)]

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  if tcname eq 'maria' then begin
    xticks=5
    xtickv=indgen(xticks+1)+15
;xticks=4
;xtickv=indgen(xticks+1)+16
    x+=14.5
    ;xtickname=strtrim(indgen(xticks+1)+15,2)
    xtitle='Day in Sept [ UTC ]'
    xrange=[14.5,20.5]
  endif else if tcname eq 'haiyan' then begin
    xticks=7
    xtickv=indgen(xticks+1)+1
;xticks=6
;xtickv=indgen(xticks+1)+2
    x+=1
    ;xtickname=strtrim(indgen(xticks+1)+1,2)
    xtitle='Day in Nov [ UTC ]'
    xrange=[1,8]
  endif

  if ~keyword_set(xrange) then $
    xrange=[min(x),max(x)]

  ;REPLACE PRESSURE WITH JUST WIND FOR PAPER FIG
  if type eq 0 then begin
    yrange=yrange2
    ytitle1=ytitle2
  endif

  plot,x,y,/nodata,position=position,$
    xstyle=9,ystyle=9,$
    xrange=xrange,yrange=yrange,yminor=2,$
    xticks=xticks,xtickv=xtickv,xticklen=0.034,xminor=2,$
    xtitle=xtitle,ytitle=ytitle1,$
    charsize=csize,$
    title=title

;SAFFIR-SIMPSON THRESHOLDS
  ls=1
  thick=1
  col=180
;  axis,yaxis=1,ystyle=5,yrange=yrange2,/noerase,/save
  ;TS
;    plots,!x.crange,replicate(18,2),linestyle=ls,color=col,thick=thick,/data
  ;CAT 1-5
    lim=[33,43,50,58,70]
;    lim=lim[[0,3,4]]
    nlim=n_elements(lim)
;    for i=0,nlim-1 do plots,!x.crange,replicate(lim[i],2),linestyle=ls,color=col,thick=thick,/data
;  axis,yaxis=1,ystyle=5,yrange=yrange,/noerase,/save

;DAY LINES
;  for id=1,nd-1 do $
;    plots,replicate(24*id,2),!y.crange,linestyle=0,thick=1,color=col

;VARIABLES

  ;HURDAT
    hdcol=160
;    if keyword_set(hdtim) then $
;      if type eq 0 then oplot,hdtim,hdpsfc,linestyle=1,thick=2.5,color=hdcol

  ;HURDAT INTENSIFICATION RATE
  if type eq 1 and var_str eq 'ri' then begin
    ;RAPID INTENSIFICATION THRESHOLD
    plots,!x.crange,replicate(15.65,2),linestyle=2,thick=1
;    nhdnew=max(hdtim)-min(hdtim)+1
    nhdnew=(max(hdtim)-min(hdtim))*24+1
    hdnew=findgen(nhdnew)/24+min(hdtim)
;    hdnew=interpol(hdtim,hdtim,hdnew)
    wspd=interpol(hdwspd,hdtim,hdnew)
    ismooth=6
;    wspd=smooth(temporary(wspd),ismooth,/edge_truncate,/nan)
    ddt=deriv(wspd)*24. ; m/s / day
    ddt=gauss_smooth(ddt,6,/edge_truncate)
    oplot,hdnew,ddt,linestyle=1,thick=2.5,color=hdcol
  endif

  cols=(findgen(dirs.nc)+1)/dirs.nc*220
  cols[0]=[0]
  lstyle=replicate(0,dirs.nc);indgen(dirs.nc)
;  lstyle[[4,5]]=[1,2]
  lthick=reverse((findgen(dirs.nc)+3)*4.5/dirs.nc)
lthick*=0.6
;  lthick_ctl=2.4
;  lthick[0]=lthick_ctl

  loadct,3,/silent

  ;VAR1
  if type ne 0 then $ ; Suppress plotting of PRESSURE for paper fig
  for ic=0,dirs.nc-1 do $
    oplot,x,var1[ic,*],linestyle=lstyle[ic],thick=lthick[ic],color=cols[ic]

  ;MARKER FOR ONSET OF RI
  ionset=0
  if ionset then begin
    if tcname eq 'maria' then t_onset=18.12
    if tcname eq 'haiyan' then t_onset=5.5
    plots,replicate(t_onset,2),!y.crange,linestyle=2,thick=1
  endif

;  if type eq 1 then begin
;    ;ADD HASHES FOR START TIME
;    for ic=1,dirs.nc-1 do begin
;      ivar=reform(var1[ic,*])
;      x0=min(where(finite(ivar)))
;      iv0=reform(ivar[x0])
;      wd=0.01
;      datc=convert_coord([0.2,position[1]+wd],/normal,/to_data)
;      plots,replicate(x[x0],2),[iv0-datc[1],iv0+datc[1]],linestyle=lstyle[ic],thick=lthick[ic],color=cols[ic],/data
;    endfor
;  endif

if iv2 then begin

  ;VAR2
  if type eq 0 then yst=5 else yst=9
  axis,yaxis=1,ystyle=yst,ytitle=ytitle2,charsize=csize,yrange=yrange2,yminor=2,/save

  ;CAT 1 THRESHOLD
    x0=!x.crange[1]
    buff=0.03
;    plots,[x0+buff,x0-buff],replicate(lim[0],2),linestyle=0,thick=2.0,color=0,/data

  ;HURDAT
    loadct,0,/silent
    if keyword_set(hdtim) then begin
      if type eq 0 then oplot,hdtim,hdwspd,linestyle=1,thick=2.5,color=hdcol
      max=max(hdwspd,mloc,/nan)
      print,'Max, HURDAT:',max,x[mloc]
      loc=max(x)
      plots,loc,max,psym=1,thick=1.3,symsize=0.8,color=hdcol,/data
    endif

  loadct,3,/silent

  for ic=0,dirs.nc-1 do $
    oplot,x,var2[ic,*],linestyle=lstyle[ic],thick=lthick[ic],color=cols[ic]
;  if type eq 0 then $
;    for ic=4,5 do $
;      oplot,x,var2[ic,*],linestyle=lstyle[ic],thick=lthick,color=180

  ;ADD MARKERS FOR MAX CARTESIAN WINDS
  ;(ONLY IF USING VTAN FOR MAIN VAR)
  if iazim then begin
    for ic=0,dirs.nc-1 do begin
      max=max(reform(var2car[ic,*]),mloc,/nan)
      print,'Max, ',dirs.cases[ic],', ',max,x[mloc]
      loc=max(x);x[mloc]
      plots,loc,max,psym=1,thick=1.3,symsize=0.8,color=cols[ic],/data
    endfor
  endif

  ;ADD HASHES FOR START TIME
  for ic=1,dirs.nc-1 do begin
    ivar=reform(var2[ic,*])
    x0=min(where(finite(ivar)))
    iv0=reform(ivar[x0])
;    wd=3
;    plots,replicate(x0,2),[iv0-wd,iv0+wd],linestyle=lstyle[ic],thick=lthick[ic],color=cols[ic]
    wd=0.028
    datc=convert_coord([0.2,position[1]+wd],/normal,/to_data)
    plots,replicate(x[x0],2),[iv0-datc[1],iv0+datc[1]],linestyle=lstyle[ic],thick=lthick[ic],color=cols[ic],/data
  endfor

endif; else begin

;  ;ADD HASHES FOR START TIME
;  for ic=1,dirs.nc-1 do begin
;    ivar=reform(var1[ic,*])
;    x0=min(where(finite(ivar)))
;    iv0=reform(ivar[x0])
;    wd=0.028
;    datc=convert_coord([0.2,position[1]+wd],/normal,/to_data)
;    plots,replicate(x0,2),[iv0-datc[1],iv0+datc[1]],linestyle=lstyle[ic],thick=lthick[ic],color=cols[ic]
;  endfor

;endelse

  loadct,0,/silent

  ;LEGEND
  ileg=1
  if ileg then begin
  loadct,3,/silent
;  if ileg then begin
    csize_fac=0.7
    margin=0.1
    pspacing=2.0 ; length of lines
    spacing=0.8 ; between lines
    leg_str=strupcase(dirs.cases[icplot])
    leg_style=lstyle
    leg_thick=lthick
    leg_color=cols;replicate(0,ncplot);cols
    legend2,leg_str,linestyle=leg_style,thick=leg_thick,COLORS=leg_color,$
      charsize=csize*csize_fac,/top_legend,/left_legend,/normal,$
      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.16,0.75]
;      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.56,0.35]
;  loadct,0,/silent
  endif

  device,/close

  convert_png,figname,res=200;,/remove_eps


print,'DONE!!'
end
