; 
; TIME SERIES OF STRATIFORM RAIN FRACTION FOR ENSEMBLE TC RUNS

; Uses Rosi's stratiform rain fraction post-processing
;   varout = 1 if convective, = 2 if stratiform, = 3 other, = 0 if no rain

; TC MOTION IS ALSO CALCULATED HERE

; James Ruppert
; 3/15/22
; 
pro run_enstc_strat

;tcname='maria'
tcname='haiyan'

if tcname eq 'maria' then begin
  tcyear='2017'
  hurdat=read_hurdat(tcname,tcyear)
  nt_ctl=145
endif else if tcname eq 'haiyan' then begin
  tcyear='2013'
  hurdat=read_jtwcdat(tcname)
  nt_ctl=169
endif

dom='d02'
case_str='ctl'
tc_ens_config, tcname, case_str, dom, dirs=dirs, dims=dims, vars=vars;, /verbose

;underlay=0 ; plot field as underlay?
;  undervar_str='SST'
;  dom='01'

;TIME ARRAYS
  time=dims.time
  time_ctl=time
  nt_full=dims.nt;-1
  npd=dims.npd
  nhrs=1.*nt_ctl*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nt_ctl)

;----BEST TRACK SUBSET--------------------

radius_thresh = 1e3 ; km radius threshold for averaging

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

  subset_wrf_ctl = where(time ge hurdat.jultim[0] and time le max(hurdat.jultim),ntall)

  ;INTERPOLATE BEST TRACK ONTO SAME TIME INTERVAL AS WRF
  hdtim=interpol(hurdat.jultim,ntall)
  hdlon=interpol(hurdat.lon,ntall)
  hdlat=interpol(hurdat.lat,ntall)

endif

;----READ VARS--------------------

iv_str='strat'

nc=7
;0 - m1, ctl
;1 - m1, ncrf
;2 - m2, ncrf
;3 - m3, ncrf
;4 - m4, ncrf
;5 - m5, ncrf
;6 - m5, crfon

;locvort=fltarr(dirs.nens,3,nt)
;locvort[*]=!values.f_nan
strat=fltarr(nt_full,nc,4)
strat[*]=!values.f_nan

;TRACK ALL EVENTS BASED ON MEMB_01 CTL

  ic=0

  ;TC TRACKING
;    spawn,'ls '+dirs.ensdir[ic]+'track*',trackfiles,err
;    ntrack=n_elements(trackfiles)
;    if ~keyword_seT(trackfiles) then continue
;
;  for ifile=0,ntrack-1 do begin
;
;    temp=read_nc_var(trackfiles[ifile],'locvort')
;    temp=reform(temporary(temp),[4,dims.nt,1])
;
;    if ~keyword_set(locvort) then $
;      locvort=temp $
;    else $
;      locvort=[[[locvort]],[[temp]]]
;
;  endfor

  file_name1='/scratch/06040/tg853394/wrfenkf/ensemble/haiyan/'
  file_name2='post/d02/strat.nc'

for ic=0,nc-1 do begin

  if ic eq 0 then begin
    case_str='ctl'
    imem=1
  endif else if ic ge 1 and ic le 5 then begin
    case_str='ncrf'
    imem=ic
  endif else if ic eq 6 then begin
    case_str='crfon'
    imem=5
  endif

  print,'Member: ',imem
  print,'Case: ',case_str

  file = file_name1+'memb_'+string(imem,format='(i2.2)')+'/'+case_str+'/'+file_name2
  tc_ens_config, tcname, case_str, dom, imember=imem, dirs=dirs, dims=dims, vars=vars;, /verbose

;  subset_wrf = where(dims.time ge hdtim[0] and dims.time le max(hdtim),nttest)
  nttest=dims.nt
;  if nttest eq 0 then message,"This test doesn't have any!"

  ;EXCLUDE EDGES
  ixskip=round( 100./(111.*(dims.lat[1]-dims.lat[0])) )
  nx=n_elements(dims.lon)
  ny=n_elements(dims.lat)
  nxskip=dims.nx-2*ixskip
  nyskip=dims.ny-2*ixskip
  ix = indgen(nxskip)+ixskip
  iy = indgen(nyskip)+ixskip
  ;LOWER CORNER
;  ix = where(( dims.lon ge 155. ),nxx)
;  iy = where(( dims.lat le 8. ),nyy)

;  count=[dims.nx,dims.ny,1,nttest] & offset=[0,0,0,0] ; x,y,z,t
  count=[nxskip,nyskip,1,nttest] & offset=[ix[0],iy[0],0,0] ; x,y,z,t
;  count=[nxx,nyy,1,nttest] & offset=[ix[0],iy[0],0,0] ; x,y,z,t
  var=reform(read_nc_var(file,'strat',count=count,offset=offset))

  offset_ctl = min(where(time_ctl ge dims.time[0]))
;  offset_hd  = min(where(hdtim ge dims.time[0]))

  for it=0,nttest-1 do begin

    ;SUBSET BY RADIUS
      ivar = reform(var[*,*,it])
;      tcloc=[hdlon[it + offset_hd],hdlat[it + offset_hd]]
;      tcrad = radius_tc_ll(tcloc,dims.lon,dims.lat)
;      irad = where(tcrad le radius_thresh,count_all)
      irad = where(finite(ivar),count_all)
      ivar2 = ivar;[irad]

    ; varout = 1 if convective, = 2 if stratiform, = 3 other, = 0 if no rain
    for ist=0,3 do begin
      icstr = where(ivar2 eq ist,count)
      if count ne 0 then $
        strat[it + offset_ctl,ic,ist] = 1d * count; / (1d*count_all)
    endfor

    strat[it + offset_ctl,ic,*] /= count_all

  endfor

endfor ; icase


;----CREATE PLOT--------------------

  dirs.figdir+=tcname+'/'
  figname=dirs.figdir+'stratfrac'
;figname=dirs.figdir+'stratfrac_subdm'
;figname=dirs.figdir+'tracks_memb'+string(ic+1,format='(i2.2)')


  ;PLOT SPECS
    csize=0.8
    position=[0.14,0.1,0.92,0.92]
    xsize=4.2 & ysize=3.4
    ctbl=39

  x=time_hrs;[subset_wrf_ctl]
  xrange=[min(x),max(x)]

  title=''
  xtitle='Hours'
  ytitle='Area Fraction'

;  leg_str=['Raining','Conv','Strat','Anvil','S/C']
  leg_str=['CTL-M1','NCRF-M1','NCRF-M2','NCRF-M3','NCRF-M4','NCRF-M5','CRFON-M5']

;SMOOTHING
for ism=0,1 do $
  strat = smooth(strat,[5,0,0],/edge_truncate,/nan)

for istrat=0,4 do begin

  ;varout = 1 if convective, = 2 if stratiform, = 3 other, = 0 if no rain

  if istrat eq 0 then begin
    ifigname = figname+'_rain'
    title = 'Raining (any)'
    yrange=[0,0.3]
    ivar = reform(strat[*,*,1] + strat[*,*,2] + strat[*,*,3])
  endif else if istrat eq 1 then begin
    ifigname = figname+'_conv'
    title = 'Convective'
    yrange=[0,.02]
    ivar = reform(strat[*,*,istrat])
  endif else if istrat eq 2 then begin
    ifigname = figname+'_strat'
    title = 'Stratiform'
    yrange=[0,.2]
    ivar = reform(strat[*,*,istrat])
;ivar += reform(strat[*,*,3])
  endif else if istrat eq 3 then begin
    ifigname = figname+'_other'
    title = 'Other (anvil)'
    yrange=[0,.1]
    ivar = reform(strat[*,*,istrat])
  endif else if istrat eq 4 then begin
    ifigname = figname+'_stratfrac'
    title = 'Stratiform / Convective'
    yrange=[5,20]
;yrange=[10,40]
    ivar = reform( (strat[*,*,2] + strat[*,*,3]) / strat[*,*,1])
  endif

  print,title
  stats,ivar

  set_plot,'ps'
  epsname=ifigname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  plot,x,indgen(10),/nodata,position=position,$
    xstyle=9,ystyle=9,$
    xrange=xrange,yrange=yrange,yminor=2,$
    xticks=xticks,xtickv=xtickv,$;xticklen=0.034,xminor=2,$
    xtitle=xtitle,ytitle=ytitle,$
    charsize=csize,$
    title=title

  loadct,ctbl,/silent

  cols=findgen(nc)*235/(nc-1);[60,120,250]
;  lstyle=replicate(0,nc)
  lstyle=replicate(0,nc);[0,1,2]
  lthick=replicate(3,nc)

  for ic=0,nc-1 do begin

    if ic eq 0 then begin
    endif else if ic ge 1 and ic le 5 then begin
    endif else if ic eq 6 then begin
    endif

    oplot,x,reform(ivar[*,ic]),linestyle=lstyle[ic],thick=lthick[ic],color=cols[ic]

  endfor

  ;LEGEND
  ileg=1
  if ileg then begin
    csize_fac=0.6;0.7
    margin=0.1
    pspacing=2.0 ; length of lines
    spacing=0.8 ; between lines
    ;leg_str=strupcase(dirs.cases[icplot])
    leg_style=lstyle;replicate(0,nstrat)
    leg_thick=lthick;replicate(2,ncplot)
    leg_color=cols
    legend2,leg_str,linestyle=leg_style,thick=leg_thick,COLORS=leg_color,$
      charsize=csize*csize_fac,/top_legend,/left_legend,/normal,$
      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.7,0.82]
  endif

  device,/close

  convert_png,ifigname,res=200,/remove_eps

endfor ; istrat

print,'DONE!!'
end
