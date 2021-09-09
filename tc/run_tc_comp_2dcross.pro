; 
; Generate PW- or OLR-binned 2D cross section from TC output
;
; Mean water vapor profile used for experiment 'FIXQV' generated here (average over hr_sel=[0,24*2]).
;
; James Ruppert
; 10/2/19
; 
pro run_tc_comp_2dcross

tcname='maria';'haiyan';'edouard';
;tcname='haiyan'
tcyear='2017';'2013'
;tcyear='2013'
if tcname eq 'haiyan' then $
  hurdat=read_jtwcdat(tcname) $
else $
  hurdat=read_hurdat(tcname,tcyear)

;subdir='moving_nest/'+tcname
subdir='static_nest/'+tcname
tc_sim_config, subdir, dirs=dirs, dims=dims, vars=vars;, /verbose
dirs.figdir+=tcname+'/'

model_pres_levels_nc=dirs.casedir[0]+'mean_pres_tim0_ctl.nc'

;----PLOT OPTIONS--------------------


binvar_str='pw';'rainrate';'pw';'olr';'pw';
var2_str='LH';'OLR';

icont=2 ; 0=stream, 1=cloud, 2=LH, 3=AVOR

;TIME SELECTION
;  hr_sel=[0,200] ; entire simulation
;  hr_sel=[0,24*2]
  hr_sel=[49,72]
;  hr_sel=[49,96]
  hr_tag=string(hr_sel[0],format='(i3.3)')+'-'+string(hr_sel[1],format='(i3.3)')+'hr'

;VORTEX TRACKING LEVEL
  psel=700 ; (hPa) level to locate vortex
  izsel=(where(dims.pres eq psel))[0]

;RADIUS THRESHOLD TO CONSIDER
  radius_thresh=800 ; km


;----TIME SPECS--------------------


;FULL TIME SERIES
  time=dims.time
  nt=dims.nt-1
  npd=dims.npd
  nhrs=1.*nt*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nhrs)

;TIME SELECTION
  t_ind=where((time_hrs ge hr_sel[0]) and (time_hrs le hr_sel[1]))
  t_ind=t_ind[where(t_ind le nt-1)]
  nt_sel=n_elements(t_ind)
  ;OVERWRITE THESE
    nhrs=1.*nt_sel*npd/24.
    nd=(1.*nhrs-(nhrs mod 24))/24.
    time_hrs=indgen(nt_sel)
    time=time[t_ind]
    nt_sav=nt_sel


;----READ VARS--------------------


;LAND MASK
  vtag='LANDMASK'
  file=dirs.files_raw[0,2,0]
  mask=reform(read_nc_var(file,'LANDMASK'))
  land=where(mask eq 1,nland)
  mask=0


nbin=40
bin_count=dblarr(dirs.nc,nbin)

bin_var2=fltarr(dirs.nc,nbin)
bin_var2[*]=!values.f_nan
bin_olrc=bin_var2

for ic=0,dirs.nc-1 do begin
;for ic=0,3 do begin
;for ic=0,1 do begin

  print,'CASE: ',strupcase(dirs.cases[ic])

  hr0=0
  if strmatch(dirs.cases[ic],'*36h*') then hr0=36
  if strmatch(dirs.cases[ic],'*24h*') then hr0=24
  if strmatch(dirs.cases[ic],'*48h*') then hr0=48
  if strmatch(dirs.cases[ic],'*60h*') then hr0=60
  if strmatch(dirs.cases[ic],'*72h*') then hr0=72
  if strmatch(dirs.cases[ic],'*84h*') then hr0=84

  ;TIME SELECTION
    if hr0 gt hr_sel[0] then begin
      print,'Not enough times in test, so skipping...'
      continue
    endif
    t_offset=max([0,hr0-hr_sel[0]])
    ut_offset=max([0,hr_sel[0]-hr0])
    nt_test=nt_sel-t_offset
    t_ind_test=indgen(nt_test)+ut_offset
;    hrs_test=time_hrs[t_ind_test+hr0]

  it0 = t_ind_test[0]

  ;VORTEX TRACKING

    ;READ ABSOLUTE VORTICITY
      iv=where(vars.vars eq 'AVOR')
      file=dirs.files_post[ic,iv]
      count=[dims.nx,dims.ny,1,nt_test] & offset=[0,0,izsel,it0] ; x,y,z,t
      avor=reform(read_nc_var(file,'AVOR',count=count,offset=offset))

    ;SMOOTH
      ixsmooth=round(111./3) ; 1-degree smoothing, run twice
      ismooth=[ixsmooth,ixsmooth,0]
      for i=1,2 do $
        avor=smooth(temporary(avor),ismooth,/edge_truncate,/nan)

    ;MASKING
      if nland gt 0 then avor[land]=!values.f_nan
      if tcname eq 'maria' then avor=wrf_maria_mask(temporary(avor),time,hurdat,dims)
      if tcname eq 'haiyan' then avor=wrf_haiyan_mask(temporary(avor),time,hurdat,dims)

    ;VORTEX TRACKING
      vloc=maria_vortex_locate(avor,dims);,/write)
      avor=0

  ;BIN VAR
    iv=where(strmatch(vars.vars,binvar_str,/fold_case))
;where(vars.vars eq strupcase(binvar_str))
    file=dirs.files_post[ic,iv]
    count=[dims.nx,dims.ny,1,nt_test] & offset=[0,0,0,it0] ; x,y,z,t
    if binvar_str ne 'rainrate' then str=strupcase(binvar_str) else str=binvar_str
    binvar=reform(read_nc_var(file,str,count=count,offset=offset))
    if binvar_str eq 'rainrate' then begin
      binvar*=1./24 ; mm/h
    endif
;    if binvar_str eq 'olr' then begin
;      binvar[*,*,0]=!values.f_nan
;      ix=indgen(3)
;      binvar[[ix,dims.nx-1-ix],*,*]=!values.f_nan
;      binvar[*,[ix,dims.ny-1-ix],*]=!values.f_nan
;    endif

  ;VAR2
    iv=where(vars.vars eq var2_str)
    file=dirs.files_post[ic,iv]
    var2=reform(read_nc_var(file,var2_str,count=count,offset=offset))
;    if var2_str eq 'OLR' then begin
      iv=where(vars.vars eq 'OLRC')
      file=dirs.files_post[ic,iv]
      olrc=reform(read_nc_var(file,'OLRC',count=count,offset=offset))
;    endif


;----SUBSET TO CERTAIN RADIUS--------------------


    for it=0,nt_test-1 do begin

      ;TC LOCATION
        ivloc = [ vloc[0,it] , vloc[1,it] ]
        tcrad = radius_tc_ll(ivloc,dims.lon,dims.lat)
        irad = where((tcrad le radius_thresh),count_rad,complement=nan)

      ibv=reform(binvar[*,*,it])
      ibv[nan]=!values.f_nan
      binvar[*,*,it]=ibv

      ivar2=reform(var2[*,*,it])
      ivar2[nan]=!values.f_nan
      var2[*,*,it]=ivar2

;      if var2_str eq 'OLR' then begin
        iolrc=reform(olrc[*,*,it])
        iolrc[nan]=!values.f_nan
        olrc[*,*,it]=iolrc
;      endif

    endfor


;----BINNING--------------------


  if binvar_str eq 'pw' then begin
    min=25
    max=80
    xtitle='PW [ mm ]'
  endif else if binvar_str eq 'olr' then begin
    min=70
    max=330
    xtitle='OLR [ W m!U-2!N ]'
  endif else if binvar_str eq 'rainrate' then begin
    min=0
    max=75
    xtitle='Rainfall [ mm h!U-1!N ]'
  endif else message,'add bin var here'

  delta=1.*(max-min)/nbin
  bins=findgen(nbin)*delta+min ; BIN VARIABLE

  for ibin=0,nbin-1 do begin
    x0 = min + delta*ibin
    ix=where((binvar ge x0) and (binvar lt x0+delta),count)
    bin_count[ic,ibin]=count
;    if count le 25 then begin
    if count lt 2 then begin
      bins[ibin]=!values.f_nan
      continue
    endif
    bin_var2[ic,ibin]=mean(var2[ix],/nan,/double)
    bin_olrc[ibin]=mean(olrc[ix],/nan,/double)
  endfor

  ;NORMALIZE PDF
;  bin_count[ic,*]/=total(reform(bin_count[ic,*]),/double)

  if binvar_str eq 'olr' then begin
    bins=reverse(temporary(bins))
    bin_var2=reverse(temporary(bin_var2))
    bin_olrc=reverse(temporary(bin_olrc))
  endif


endfor ; icase

;----CREATE PLOTS--------------------


  tc_figspecs, 'AVOR', figspecs

  figdir=dirs.figdir+'/comp_cross/'
  spawn,'mkdir '+figdir,tmp,tmpe
  figname=figdir+'2d_pdf'
  figspecs=create_struct(figspecs,'figname',figname)

  ;PLOT SPECS
    csize=0.85
    position=[0.14,0.22,0.85,0.97]
    xsize=4.8 & ysize=1.7
    x=bins

  set_plot,'ps'
  epsname=figspecs.figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  ;PDFs
    plot,x,indgen(10),/nodata,position=position,ylog=1,$
      xstyle=9,ystyle=9,xminor=4,$
      yticklen=0.015,xticklen=0.05,$
      xrange=xrange,yrange=[10,1e6],$;yrange=[-1,6],$;yrange=[0,0.1],$
      xtitle=xtitle,$;yminor=4,$
      ytitle='N',charsize=csize

    for ic=0,dirs.nc-1 do $
      oplot,x,bin_count[ic,*],linestyle=ic,thick=2,color=0
;    for ic=0,dirs.nc-1 do $
;      oplot,x,(bin_count[ic,*]/bin_count[0,*]-1),linestyle=ic,thick=2,color=0

  ;VAR2
  if var2_str eq 'OLR' then begin
    axis,yaxis=1,yticks=2,ytickv=[100,200,300],ystyle=1,yminor=4,charsize=csize,yrange=[60,310],ylog=0,/save,$
      yticklen=0.015,ytitle='OLR [ W m!U-2!N ]'
  endif else if var2_str eq 'LH' then begin
    yrange2=[0,300]
    ytitle2='LH [ W m!U-2!N ]'
    axis,yaxis=1,ystyle=1,charsize=csize,yrange=yrange2,ylog=0,/save,$
      yticklen=0.015,ytitle=ytitle2
  endif

    for ic=0,dirs.nc-1 do $
      oplot,x,bin_var2[ic,*],linestyle=ic,thick=3,color=0
    if var2_str eq 'OLR' then $
      for ic=0,dirs.nc-1 do $
        oplot,x,bin_olrc[ic,*],linestyle=ic,thick=1,color=0

  device,/close

  convert_png,figspecs.figname,/remove_eps,res=400

;endfor ; icase

print,'DONE!!'
end
