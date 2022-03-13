; 
; Maps from ensemble TC output

; James Ruppert
; 12/20/21
; 
pro run_enstc_maps

;tcname='maria'
tcname='haiyan'
case_str='ctl';'ncrf';'ctl'

if tcname eq 'maria' then begin
;  tcyear='2017'
;  hurdat=read_hurdat(tcname,tcyear)
  nt_ctl=145
endif else if tcname eq 'haiyan' then begin
;  tcyear='2013'
;  hurdat=read_jtwcdat(tcname)
  nt_ctl=169
endif

dom='d02'
tc_ens_config, tcname, case_str, dom, dirs=dirs, dims=dims, vars=vars;, /verbose
dirs.figdir+=tcname+'/'

;----PLOT OPTIONS--------------------

iwind=0 ; Plot 10m wind vectors?
itcloc=0 ; Plot TC location?

ishear=0 ; Plot shear vector?

allvars=['wspd','RAINNC','rainrate','avor','slp','lh','th_e','olr','pw']
allvars=['olr'];,'pw','olr','pw']
;allvars=['rainrate']
allvars=['avor']
allvars=['strat'] ; varout = 1 if convective, = 2 if stratiform, = 3 other, = 0 if no rain
;allvars=['ngms']
;allvars=['wspd']
;allvars=['z']
nvsel=n_elements(allvars)

;LEVEL SELECTION
  psel=700;800;700 700 is the winner of 600, 700, 800, 850, 925
  izlev=(where(dims.pres eq psel))[0]

;TIME SELECTION
;  hr_sel=[0,200] ; entire simulation
;  hr_sel=[0,60];72]
;  hr_sel=[70,72]
;  hr_sel=[24,48]


;----TIME SPECS--------------------


;FULL TIME SERIES
  time=dims.time
  nt_full=dims.nt
  npd=dims.npd
  nhrs=1.*nt_full*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nt_ctl)


;----READ VARS--------------------


loc=fltarr(dirs.nens,2,nt_full)
loc[*]=!values.f_nan

;LAND MASK
  vtag='LANDMASK'
  file=dirs.files_raw[0,2,0]
  mask=reform(read_nc_var(file,'LANDMASK'))
;  land=where(mask eq 1,nland)
;  mask=0


;for ic=0,dirs.nens-1 do begin
;for ic=0,dirs.nens-1,2 do begin
for ic=0,0 do begin

  print,'Memb: ',ic+1

    i_nt=nt_full
    offset=0
    if strmatch(case_str,'ncrf') then offset=36
    it_test=indgen(i_nt)+offset

;IVAR LOOP
for ivar_sel=0,nvsel-1 do begin
;for ivar_sel=0,0 do begin

var_str=allvars[ivar_sel]
print,'VAR: ',var_str


  iwmax=0

  if var_str eq 'th_e' then begin
    ;T
      iv=where(vars.vars eq 'T2')
      file=dirs.files_post[ic,iv]
      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
      t=reform(read_nc_var(file,'',varid='0',count=count,offset=offset))
    ;q
      iv=where(vars.vars eq 'Q2')
      file=dirs.files_post[ic,iv]
      q=reform(read_nc_var(file,'',varid='0',count=count,offset=offset))
    ;PSFC
      iv=where(vars.vars eq 'PSFC')
      file=dirs.files_post[ic,iv]
      p=reform(read_nc_var(file,'',varid='0',count=count,offset=offset))
    ;THETA-E (ASSUME R_TOT = R_V; EQN SOURCE: BRYAN AND FRITSCH 2002, P.2921)
      var = theta_e(t,p,q,q)
  endif else if var_str eq 'wspd' then begin
    ;U
      iv=where(vars.vars eq 'U10')
      file=dirs.files_post[ic,iv]
      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
      u=reform(read_nc_var(file,'U10',count=count,offset=offset))
    ;V
      iv=where(vars.vars eq 'V10')
      file=dirs.files_post[ic,iv]
      v=reform(read_nc_var(file,'V10',count=count,offset=offset))
    var = sqrt(u^2+v^2)
    u=0 & v=0
    iwmax=0;1
    if iwmax then begin
      ;USE CTL FOR TIMES UP TO RESTART
      if dirs.cases[ic] ne 'ctl' then begin
        ic2=where(dirs.cases eq 'ctl',count)
        if count eq 0 then message,'Need CTL for MAXWSPD!!'
        iv=where(vars.vars eq 'U10')
        file=dirs.files_post[ic2,iv]
        count=[dims.nx,dims.ny,1,int_ctl] & offset=[0,0,0,0] ; x,y,z,t
        u=reform(read_nc_var(file,'U10',count=count,offset=offset))
        iv=where(vars.vars eq 'V10')
        file=dirs.files_post[ic2,iv]
        v=reform(read_nc_var(file,'V10',count=count,offset=offset))
        varctl=sqrt(u^2+v^2)
      ;NOW CONCATENATE
        var=[[[varctl]],[[temporary(var)]]]
      endif
      t_sub=where(time_hrs ge 0)
      ;t_sub=where(time_hrs ge 84)
      var=max(temporary(var[*,*,t_sub]),dimension=3)
    endif
  endif else if var_str eq 'tqc' then begin
    ;QC
;      iv=where(vars.vars eq 'QCLOUD')
;      file=dirs.files_post[ic,iv]
;      count=[dims.nx,dims.ny,dims.np,i_nt] & offset=[0,0,0,0] ; x,y,z,t
;;      qc=reform(read_nc_var(file,'QCLOUD',count=count,offset=offset))
;      qc=reform(read_nc_var(file,'QCLOUD'))
    ;QI
      iv=where(vars.vars eq 'QICE')
      file=dirs.files_post[ic,iv]
      count=[dims.nx,dims.ny,dims.np,i_nt] & offset=[0,0,0,0] ; x,y,z,t
      qi=reform(read_nc_var(file,'QICE',count=count,offset=offset))
;    var=total(qc+qi,3,/double)*2500./9.81
    var=total(qi,3,/double)*2500./9.81 ; dp/g
  endif else if var_str eq 'lw' then begin
    ;VERTICALLY INTEGRATED LONGWAVE RADIATIVE HEATING
      iv=where(vars.vars eq 'RTHRATLW')
      file=dirs.files_post[ic,iv]
      count=[dims.nx,dims.ny,dims.np,i_nt] & offset=[0,0,0,0] ; x,y,z,t
      lw=reform(read_nc_var(file,'RTHRATLW',count=count,offset=offset))*1004. ; K/s --> J/kg/s = W/kg
      var=total(lw,3,/double)*2500./9.81 ; W/kg --> W/m2
      i_rm_mean=1
      if i_rm_mean then begin
        for it=0,i_nt-1 do begin
          varm=mean(var[*,*,it],/nan,/double)
          var[*,*,it]-=varm
        endfor
      endif
  endif else if var_str eq 'ngms' then begin
    ;NORMALIZED GROSS MOIST STABILITY (SOBEL ET AL. 2014, JAS)

      iv=where(vars.vars eq 'W')
      file=dirs.files_post[ic,iv]
      count=[dims.nx,dims.ny,dims.np,i_nt] & offset=[0,0,0,it_test[0]] ; x,y,z,t
      w=reform(read_nc_var(file,'W',count=count,offset=offset)) ; m/s
      iv=where(vars.vars eq 'QVAPOR')
      file=dirs.files_post[ic,iv]
      count=[dims.nx,dims.ny,dims.np,i_nt] & offset=[0,0,0,it_test[0]] ; x,y,z,t
      qv=reform(read_nc_var(file,'QVAPOR',count=count,offset=offset)) ; kg/kg
      iv=where(vars.vars eq 'T')
      file=dirs.files_post[ic,iv]
      count=[dims.nx,dims.ny,dims.np,i_nt] & offset=[0,0,0,it_test[0]] ; x,y,z,t
      tmpk=reform(read_nc_var(file,'T',count=count,offset=offset)) ; K

      ;RHO
        tvirt = tmpk*(1.+0.61*qv)
        rho=tvirt
        for iz=0,dims.np-1 do rho[*,*,iz,*] = dims.pres[iz]*1e2 / ( 287. * tvirt[*,*,iz,*] )
        tvirt=0
      ;HEIGHT FROM HYDROSTATIC
        z = rho & z[*]=0.
        for iz=1,dims.np-1 do $
          z[*,*,iz,*] = z[*,*,iz-1,*] + (dims.pres[iz-1]-dims.pres[iz])*1e2 / 9.81 / mean(reform(rho[*,*,iz-1:iz,*]),dimension=3,/double,/nan)
        rho=0
      dse = 1004.*tmpk + 9.81*z ; J/kg
      mse = dse + qv*2.5e6      ; J/kg

      ;SMOOTH INPUT FIELDS
        ;ix_sm=round(1./(dims.lat[1]-dims.lat[0])) ; 1-degree
        ix_sm=9
        ismooth=[ix_sm,ix_sm,0,0]
        nsm=10;29
        for i=0,nsm-1 do begin
          w=smooth(temporary(w),ismooth,/edge_truncate,/nan)
          dse=smooth(temporary(dse),ismooth,/edge_truncate,/nan)
          mse=smooth(temporary(mse),ismooth,/edge_truncate,/nan)
        endfor

      ;VERTICAL ADVECTION
        ddpd = tmpk & ddpd[*]=0.
        ddpm = ddpd
        for ix=0,dims.nx-1 do $
          for iy=0,dims.ny-1 do $
            for it=0,i_nt-1 do begin
              ddpd[ix,iy,*,it]=deriv(dims.pres*1e2,reform(dse[ix,iy,*,it]))
              ddpm[ix,iy,*,it]=deriv(dims.pres*1e2,reform(mse[ix,iy,*,it]))
            endfor
        vadvd = w * ddpd
        vadvm = w * ddpm

      vadvd_int = total(vadvd*100e2,3,/double)/9.81
      vadvm_int = total(vadvm*100e2,3,/double)/9.81
      var = vadvm_int / vadvd_int ; Unitless

  endif else begin

    if var_str eq 'rainrate' or var_str eq 'strat' then iv_str=var_str $
    else iv_str=strupcase(var_str)

    iv=where(vars.vars eq iv_str)
    file=dirs.files_post[ic,iv]

    if var_str eq 'avor' then klev=izlev else klev=0

    count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,klev,0] ; x,y,z,t
    var=reform(read_nc_var(file,iv_str,count=count,offset=offset))

;vloc=tc_track(var,dims,time,hurdat,kernel=kernel);,/write,trackfil=dirs.ensdir[ic]

;    if var_str eq 'olr' then begin
;      var[*,*,0]=!values.f_nan
;      ix=indgen(3)
;      var[[ix,dims.nx-1-ix],*,*]=!values.f_nan
;      var[*,[ix,dims.ny-1-ix],*]=!values.f_nan
;    endif

  endelse

  ;RAIN RATE
  if var_str eq 'rainrate' then var*=1./24 ; mm/d --> mm/h

  ;VERTICAL SHEAR
  if ishear then begin

    pshr=200 ; (hPa) level
    izshr2=(where(dims.pres eq pshr))[0]
    pshr=850 ; (hPa) level
    izshr8=(where(dims.pres eq pshr))[0]

    iv=where(vars.vars eq 'U')
    file=dirs.files_post[ic,iv]
    count=[dims.nx,dims.ny,1,i_nt]
    u200=reform(read_nc_var(file,'U',count=count,offset=[0,0,izshr2,0]))
    u850=reform(read_nc_var(file,'U',count=count,offset=[0,0,izshr8,0]))

    iv=where(vars.vars eq 'V')
    file=dirs.files_post[ic,iv]
    count=[dims.nx,dims.ny,1,i_nt] 
    v200=reform(read_nc_var(file,'V',count=count,offset=[0,0,izshr2,0]))
    v850=reform(read_nc_var(file,'V',count=count,offset=[0,0,izshr8,0]))
;v850[where(abs(v850) ge 5e2)]=!values.f_nan

    shr=fltarr(2,i_nt)

    ;WITHIN CERTAIN RADIUS, RADIUS-NORMALIZED
    for it=0,i_nt-2 do begin

      ;TC LOCATION
        ivloc = [ vloc[0,it] , vloc[1,it] ]
        tcrad = radius_tc_ll(ivloc,dims.lon,dims.lat)
        irad = where(((tcrad ge 0) and (tcrad le 800)))
        tcrad=1d*reform(tcrad[irad])
        tcradtot=total(tcrad,/double,/nan)

      u2=total( ( (reform(u200[*,*,it]))[irad] * tcrad) ,/double,/nan) / tcradtot
      v2=total( ( (reform(v200[*,*,it]))[irad] * tcrad) ,/double,/nan) / tcradtot
      u8=total( ( (reform(u850[*,*,it]))[irad] * tcrad) ,/double,/nan) / tcradtot
      v8=total( ( (reform(v850[*,*,it]))[irad] * tcrad) ,/double,/nan) / tcradtot

      shr[0,it]=u2-u8
      shr[1,it]=v2-v8

    endfor

  endif

  ;SST
  if var_str eq 'SST' then begin
    for it=0,i_nt-1 do begin
      isst=reform(var[*,*,it]) - 273.15 ; K --> C
      isst[land]=!values.f_nan
      var[*,*,it]=isst
    endfor
  endif

  ;WIND
    if iwind then begin
    ;U
      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
      iv=where(vars.vars eq 'U10')
      file=dirs.files_post[ic,iv]
      u=reform(read_nc_var(file,'',varid='0',count=count,offset=offset))
    ;V
      iv=where(vars.vars eq 'V10')
      file=dirs.files_post[ic,iv]
      v=reform(read_nc_var(file,'',varid='0',count=count,offset=offset))
;      wind=create_struct('u',u,'v',v)
    endif


;----CREATE PLOTS--------------------

  if var_str eq 'th_e' then begin
    setmin=340
    setmax=370
  endif else if var_str eq 'wspd' then begin
    setmin=18
    setmax=58
  endif else if var_str eq 'slp' then begin
    setmin=1000
    setmax=1015
setmax=5
setmin=-5
  endif else if var_str eq 'pw' then begin
    setmin=30
    setmax=70
  endif else if var_str eq 'olr' then begin
    setmin=90
    setmax=290
  endif else if var_str eq 'lh' then begin
    setmin=100
    setmax=300
  endif else if var_str eq 'avor' then begin
    setmin=0
    setmax=20;100
  endif else if var_str eq 'rainrate' then begin
    setmin=.1
    setmax=40
  endif else if var_str eq 'RAINNC' then begin
    setmin=10
    setmax=50
  endif else if var_str eq 'SST' then begin
    setmin=26
    setmax=32
  endif else if var_str eq 'ngms' then begin
    setmin=-2
    setmax=2
  endif
print,var_str
  tc_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figdir=dirs.figdir+var_str+'/'
  spawn,'mkdir '+figdir,tmp,tmpe
  figspecs=create_struct(figspecs,'figname',' ')
;  figspecs.title+=' ('+strupcase(dirs.cases[ic]);+')'
  titlesav=figspecs.title+' ('+strupcase(dirs.memb[ic]);+')'

  ;LOOP THROUGH TIMES

  skip=12
  t1=70
 ; for it=24,120,skip do begin
;  for it=96,120,6 do begin
  for it=t1,t1 do begin

  if var_str eq 'rainrate' then itop-=1

    ;TIME INFO
;    t0=it;*npd
    it_sel = where(time_hrs[it_test] eq it,count)
    utc=time[it_sel]
    time_zone_conv=-4d/24
    ;caldat,utc+time_zone_conv,mm,dd,yy,hh,nn
    ;t_stamp=string(hh,format='(i2.2)')+' L'
    caldat,utc,mm,dd,yy,hh,nn
    t_stamp=string(mm,format='(i2.2)')+'/'+string(dd,format='(i2.2)')+' '+string(hh,format='(i2.2)')+'00 UTC'
    print,t_stamp

;    if count eq 0 then continue

;    tim_tag='d'+string(it,format='(i2.2)')
;    tim_tag='h'+string(time_hrs[it]+hr_sel[0],format='(i2.2)')
    tim_tag='h'+string(time_hrs[it_test[it_sel]],format='(i3.3)')
    print,tim_tag

    figspecs.title=titlesav+'; '+string(time_hrs[it_test[it_sel]],format='(i3.3)')+' hr)'
;    figspecs.title=t_stamp+', '+strtrim(round(wspd_sav[it_sel]),2)+' m/s'

    if iwmax then ivar=var $
    else ivar = reform(var[*,*,it_sel])

    if var_str eq 'rainrate' then begin
      it_int=indgen(24)+it_sel[0]-23
      ivar = total(var[*,*,it_int],3,/double,/nan)
      print,'RAINRATE: LATEST 24 HR'
      figspecs.cbar_tag='[ mm d!U-1!N ]'
    endif

    if iwind then begin
;      iu = mean(u[*,*,it_sel],dimension=3,/nan,/double)
;      iv = mean(v[*,*,it_sel],dimension=3,/nan,/double)
      iu = reform(u[*,*,it_sel])
      iv = reform(v[*,*,it_sel])
      wind=create_struct('u',iu,'v',iv)
      cvar=sqrt(iu^2+iv^2)
    endif

    ;SMOOTHING
    if var_str eq 'slp' then begin
      ;ix_sm=round(111./3/2) ; 0.5-degree
      ;ix_sm=round(1./(dims.lat[1]-dims.lat[0])/2) ; 0.5-degree
      ix_sm=9;round(1./(dims.lat[1]-dims.lat[0])) ; 1-degree
      ismooth=[ix_sm,ix_sm]
      ivar -= mean(var,dimension=3,/nan,/double)
      nsm=5
      for i=1,nsm do $
        ivar=smooth(temporary(ivar),ismooth,/edge_truncate,/nan)
    endif

    figname=figdir+tim_tag+'_'+var_str+'_'+case_str+'_'+dirs.memb[ic]
    if var_str eq 'avor' then figname+='_'+string(psel,format='(i3.3)')
    figspecs.figname=figname
    print,figname

    ;HURDAT LOCATION
;    tdiff=abs(hurdat.jultim - time[it_test[it_sel]])
;    tloc=(where((tdiff eq min(tdiff)) and (tdiff lt 5d/24),count))[0]
;    if count then loc_pmin=[hurdat.lon[tloc],hurdat.lat[tloc]] else loc_pmin=0

;    if total(strmatch(['rainrate','RAINNC','wspd'],var_str)) then loc[*]=!values.f_nan

    stats,ivar

    ;REMOVE EXCEEDING POINTS
;    if var_str eq 'olr' or var_str eq 'pw' then begin
;      iless=where(ivar lt min(figspecs.levels),count)
;      if count gt 0 then ivar[iless]=min(figspecs.levels)
;    endif

    if ishear then ishr=reform(shr[*,it_sel])

    wrf_enstc_map_plot, dirs, ivar, dims.lon, dims.lat, figspecs, wind=wind, loc_pmin=reform(loc[ic,*,it_test[it_sel]]), cvar=cvar, $
      itcloc=itcloc, shr=ishr
;    wrf_tc_map_plot, dirs, ivar, dims.lon, dims.lat, figspecs, wind=wind, loc_pmin=reform(loc[ic,*,it_test[it_sel]]), itcloc=itcloc, shr=ishr
;    wrf_tc_smoothcomp_map, dirs, ivar, ivar_sm, dims.lon, dims.lat, figspecs, wind=wind, loc_pmin=reform(loc[ic,*,it_test[it_sel]]), cvar=filler, itcloc=itcloc, shr=ishr

  endfor

endfor ; icase

;endfor ; itfil (plotting via raw files)

endfor ; ivar

print,'DONE!!'
end
