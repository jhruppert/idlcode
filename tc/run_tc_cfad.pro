; 
; Make PDF from Cartesian (instead of azim) TC output.
;
; James Ruppert
; 4/29/19
; 
pro run_tc_cfad

tcname='maria';'edouard';
tcyear='2017'
hurdat=read_hurdat(tcname,tcyear)

;subdir='moving_nest/'+tcname
subdir='static_nest/'+tcname
tc_sim_config, subdir, dirs=dirs, dims=dims, vars=vars;, /verbose
dirs.figdir+=tcname+'/'

;TIME SELECTION
;  hr_sel=[0,200] ; entire simulation
  hr_sel=[49,72]
;hr_sel+=24
  hr_tag=strtrim(hr_sel[0],2)+'-'+strtrim(hr_sel[1],2)+'hr'

;VORTEX TRACKING LEVEL
  psel=700 ; (hPa) level to locate vortex
  izsel=(where(dims.pres eq psel))[0]

;---SETTINGS---------------------

;RADIUS THRESHOLD FOR HISTOGRAM
;  radius_thresh=[0,100]
  radius_thresh=[0,300] ; km
;  radius_thresh=[100,300] ; km
;  radius_thresh=[0,500]
  radius_thresh=[650,800]
;  radius_thresh=[700,800]

;VERTICAL LEVELS AND RANGE
  ytickv=[1000,700,500,400,300,200,150];,100];,50]
;  ytickv=[1000,925,850,700,600,500];,400,300,200,150];,100];,50]

;FIGURE TYPE
  type='lr';'wprm';'rh';
;  type='wprm'
;  type='rh'
;  type='ths'
;  type='the'
;  type='th'

;----TIME SPECS--------------------

;FULL TIME SERIES
  time=dims.time
  nt_full=dims.nt-1
  npd=dims.npd
  nhrs=1.*nt_full*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nhrs)
  nt_sav=dims.nt-1

;TIME SELECTION
  t_ind=where((time_hrs ge hr_sel[0]) and (time_hrs le hr_sel[1]))
  t_ind=t_ind[where(t_ind le nt_full-1)]
  nt_sel=n_elements(t_ind)
  time=time[t_ind]

;----READ VARS--------------------


;LAND MASK
  vtag='LANDMASK'
  file=dirs.files_raw[0,2,0]
  mask=reform(read_nc_var(file,'LANDMASK'))
  land=where(mask eq 1,nland)
  mask=0

;NEED LARGE PRES ARRAY FOR RH
;  if var1_str eq 'RH' then begin
;    prx=fltarr(dims.nx,dims.ny,dims.np,nt_sel)
;    for iz=0,dims.np-1 do prx[*,*,iz,*]=dims.pres[iz]*1e2
;  endif

for ic=0,dirs.nc-1 do begin
;for ic=0,0 do begin

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

  ;VORTEX TRACKING
  
  it0 = t_ind_test[0]

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
      if tcname eq 'maria' then avor=wrf_maria_mask(temporary(avor),time,hurdat,dims) else stop
  
    ;VORTEX TRACKING
      vloc=maria_vortex_locate(avor,dims);,/write)
      avor=0

  ;MAIN VARIABLES

  if type eq 'wprm' then begin
  ;W
    count=[dims.nx,dims.ny,dims.np,nt_test] & offset=[0,0,0,it0] ; x,y,z,t
    iv=where(vars.vars eq 'W')
    file=dirs.files_post[ic,iv]
    var=reform(read_nc_var(file,'W',count=count,offset=offset)) * 1e2 ; --> cm/s
;    iv=where(vars.vars eq 'QCLOUD')
;    file=dirs.files_post[ic,iv]
;    qc=reform(read_nc_var(file,'QCLOUD',count=count,offset=offset))
;    iv=where(vars.vars eq 'QICE')
;    file=dirs.files_post[ic,iv]
;    qi=reform(read_nc_var(file,'QICE',count=count,offset=offset))
;    cloud=(qc+qi)*1e6 & qc=0 & qi=0
;    cl_thresh=5
;    w_thresh=0.;5
;    icloud=where((cloud ge cl_thresh) and (w ge w_thresh),ncloud,complement=comp)
;;    icloud=where((cloud ge cl_thresh),ncloud,complement=comp)
;;print,1d*ncloud/n_elements(comp)*100
;;    ifin=where(abs(w) ge 0.5,complement=comp)
;    w[comp]=!values.f_nan
  endif else if type eq 'lr' then begin
  ;LAPSE RATE
    count=[dims.nx,dims.ny,dims.np,nt_test] & offset=[0,0,0,it0] ; x,y,z,t
    iv=where(vars.vars eq 'T')
    file=dirs.files_post[ic,iv]
    tmpk=reform(read_nc_var(file,'T',count=count,offset=offset))
;    iv=where(vars.vars eq 'QVAPOR')
;    file=dirs.files_post[ic,iv]
;    qv=reform(read_nc_var(file,'QVAPOR',count=count,offset=offset))
    pra=tmpk
    for iz=0,dims.np-1 do pra[*,*,iz,*]=dims.pres[iz];*1e2
    th=theta(tmpk-273.15,pra);*(1.+0.61*qv)
    pra=0 ;& tmpk=0 & qv=0
;    iv=where(vars.vars eq 'QCLOUD')
;    file=dirs.files_post[ic,iv]
;    qc=reform(read_nc_var(file,'QCLOUD',count=count,offset=offset))
;    iv=where(vars.vars eq 'QICE')
;    file=dirs.files_post[ic,iv]
;    qi=reform(read_nc_var(file,'QICE',count=count,offset=offset))
;    cloud=(qc+qi)*1e6 & qc=0 & qi=0
;    cl_thresh=5
;    icloud=where((cloud ge cl_thresh),ncloud,complement=comp)
;    cloud=0
;;print,1d*ncloud/n_elements(comp)*100
;    tmpk[comp]=!values.f_nan
  endif else if type eq 'rh' then begin
  ;RELATIVE HUMIDITY
    count=[dims.nx,dims.ny,dims.np,nt_test] & offset=[0,0,0,it0] ; x,y,z,t
    iv=where(vars.vars eq 'T')
    file=dirs.files_post[ic,iv]
    tmpk=reform(read_nc_var(file,'T',count=count,offset=offset))
    iv=where(vars.vars eq 'QVAPOR')
    file=dirs.files_post[ic,iv]
    qv=reform(read_nc_var(file,'QVAPOR',count=count,offset=offset))
    pra=tmpk
    for iz=0,dims.np-1 do pra[*,*,iz,*]=dims.pres[iz]*1e2
    var=calc_relh(qv,tmpk,pra)
    pra=0 & tmpk=0 & qv=0
  endif else if total(strmatch(['th','the','ths'],type)) then begin
    count=[dims.nx,dims.ny,dims.np,nt_test] & offset=[0,0,0,it0] ; x,y,z,t
    iv=where(vars.vars eq 'T')
    file=dirs.files_post[ic,iv]
    tmpk=reform(read_nc_var(file,'T',count=count,offset=offset))
    pra=tmpk
    for iz=0,dims.np-1 do pra[*,*,iz,*]=dims.pres[iz]*1e2
  ;THETA
    if type eq 'th' then var=theta(tmpk-273.15,pra*1e-2)
    if type eq 'the' then begin
  ;THETA-E
      iv=where(vars.vars eq 'QVAPOR')
      file=dirs.files_post[ic,iv]
      qv=reform(read_nc_var(file,'QVAPOR',count=count,offset=offset))
      var=theta_e(tmpk,pra,qv,qv)
    endif else if type eq 'ths' then begin
  ;THETA-ES
      qvs=mixr_sat(tmpk,pra*1e-2)*1e-3 ; --> kg/kg
      var=theta_e(tmpk,pra,qvs,qvs)
    endif
    pra=0 & tmpk=0 & qvs=0
  endif


;----CREATE PDF--------------------


    if type eq 'wprm' then begin
    ;UNITS OF CM/S
;      histmin=w_thresh*1d
      histmax=50d;.125d*1d2
      histmin=-1.*histmax
      binsize=3d;0.2*1d2
      nbins=(histmax-histmin)/binsize + 1
;      nbins=18;100
    endif else if type eq 'lr' then begin
    ;UNITS OF mK/hPa
      histmin=0d
      histmax=102.5d
      binsize=5
      nbins=(histmax-histmin)/binsize + 1
;      nbins=30
    endif else if type eq 'rh' then begin
    ;UNITS OF %
      histmin=30d
      histmax=110d
      binsize=5
      nbins=(histmax-histmin)/binsize + 1
    endif else if type eq 'ths' then begin
    ;UNITS OF K
      histmin=320d
      histmax=385d
      binsize=2
      nbins=(histmax-histmin)/binsize + 1
    endif else if type eq 'the' then begin
    ;UNITS OF K
      histmin=320d
      histmax=365d
      binsize=1;3;2
      nbins=(histmax-histmin)/binsize + 1
    endif else if type eq 'th' then begin
    ;UNITS OF K
      histmin=290d
      histmax=360d
      binsize=2
      nbins=(histmax-histmin)/binsize + 1
    endif

    var_prof=fltarr(nt_test,dims.np) ; To save mean or stddev profile

    nbins=round(nbins)
;    nbins=histmax/binsize
;    binsize=(histmax-histmin)/nbins
    bins=indgen(nbins)*binsize+histmin
    pdf=dblarr(nt_test,nbins,dims.np)
    pdf[*]=!values.d_nan
;    var_mean=fltarr(nt_test)
;    varian=fltarr(nt_test)

    min_points=50 ; Leave as NAN if fewer than this many data points

    for it=0,nt_test-1 do begin

      ;TC LOCATION
        ivloc = [ vloc[0,it] , vloc[1,it] ]
        tcrad = radius_tc_ll(ivloc,dims.lon,dims.lat)
        irad = where(((tcrad ge radius_thresh[0]) and (tcrad le radius_thresh[1])),count_rad)
        itcrad_sav=tcrad[irad]
        totrads=total(itcrad_sav,/double)

      ;CALCULATE LAPSE RATE
        if type eq 'lr' then begin
          s=fltarr(count_rad,dims.np)
          s[*]=!values.f_nan
          for i=0,count_rad-1 do begin
            ix=array_indices([dims.nx,dims.ny],irad[i],/dimensions)
            itmp = reform( tmpk[ix[0],ix[1],*,it] )
            ith  = reform(   th[ix[0],ix[1],*,it] )
            dthdz = deriv(dims.pres*1e2,ith)
;dthdz=smooth(dthdz,3,/edge_truncate)
            s[i,*] = -1. * itmp / ith * dthdz * 1e3 * 1e2 ; K/Pa --> mK / hPa
;s[i,*] = itmp
        ;for ip=0,dims.np-1 do print,dims.pres[ip],' ',ivar[i,ip]
          endfor
;s=smooth(temporary(s),[0,3],/edge_truncate,/nan)
        endif

      ;HEIGHT LOOP
      for iz=0,dims.np-1 do begin
  
        if type eq 'lr' then $
          ivar = reform(s[*,iz]) $
        else $
          ivar = (reform(var[*,*,iz,it]))[irad]

        ;SAVE MEAN/STDDEV PROFILE
;        if type eq 'wprm' then $
;          var_prof[it,iz] = total(ivar*itcrad_sav,/double,/nan) / totrads $
;        else if type eq 'lr' then $
          var_prof[it,iz] = total(ivar*itcrad_sav,/double,/nan) / totrads

        fin=where(finite(ivar),nfin)
        if nfin lt min_points then continue
        ivar=reform(ivar[fin])

        itcrad=itcrad_sav[fin]
  
;        if type eq 'wprm' then begin
;          ;REMOVE MEAN CIRCULATION
;  ;        ivar -= mean(ivar,/double,/nan)
;          ;FILTER FOR UPDRAFTS
;            loc=where((ivar ge w_thresh),nupdraft,complement=comp)
;            ivar=reform(ivar[loc])
;          if nupdraft lt min_points then continue
;        endif; else if type eq 'lr' then begin
;        ;endif
  
        ihist = histogram( ivar, nbins=nbins, binsize=binsize, min=histmin ,/nan ,/l64, REVERSE_INDICES=ri)

        ;NORMALIZE BY RADIUS
        totrad=total(itcrad,/double)
        fact=dblarr(nbins)
        for i=0,nbins-1 do begin
          if ri[i] eq ri[i+1] then continue
          fact[i] = total( itcrad[ ri[ri[i]:ri[i+1]-1] ] , /double )
        endfor
        fact/=totrad
        ihist*=fact

        ;NORMALIZE BY AREA UNDER CURVE
        area = total(ihist,/nan,/double,/preserve_type)
        pdf[it,*,iz] = ihist*100d/area
  ;      var_mean[it] = mean( ivar[irad], /nan, /double)
  ;      varian[it] = variance( ivar[irad], /nan, /double)
  
      endfor ; height loop

    endfor ; time loop

    ;TIME AVERAGE
    pdf = mean(temporary(pdf),dimension=1,/nan,/double) 
    var_prof = mean(temporary(var_prof),dimension=1,/nan,/double)
;for iz=0,dims.np-1 do print,dims.pres[iz],' ',var_prof[iz]-273.15
;p_zeroc=interpol(dims.pres,var_prof-273.15,0)
;print,p_zeroc
;exit

;    ;SMOOTH
;      ismooth=1
;      nsmooth=3
;      if ismooth then begin
;        var_t1=smooth(temporary(var_t1),nsmooth,/edge_truncate,/nan)
;        var_t2=smooth(temporary(var_t2),nsmooth,/edge_truncate,/nan)
;        var_t3=smooth(temporary(var_t3),nsmooth,/edge_truncate,/nan)
;      endif

    pdf = smooth(temporary(pdf),[0,3],/edge_truncate,/nan)

;  if dirs.cases[ic] eq 'ctl' then begin
  if ic eq 0 then begin
    pdf_ctl=pdf
    var_prof_ctl=var_prof
  endif else begin
;    pdf = temporary(pdf) - pdf_ctl
    pdf = pdf_ctl - temporary(pdf)
  endelse


;----CREATE PLOTS--------------------


  figdir=dirs.figdir+'cfads/'
  spawn,'mkdir '+figdir,tmp,tmpe
  figname=figdir+type+'_pdf_prof_'+dirs.cases[ic]
  radtag=string(radius_thresh[0],format='(i3.3)')+'-'+string(radius_thresh[1],format='(i3.3)')+'km'
  figname+='_'+radtag

  ;SET UP SHADING

    cbar_tag='[ % ]'
    if type eq 'wprm' then begin
      title='Histogram of !8w!X'
      ;FOR CTL
        ctblc=11;0
        maxc=30;16
        minc=1
        irevc=0
        ncolsc=21
        setlevs=0
      ;SENS TESTS
        ctbl=70
        irev=1
        max=2;4;6;8
        min=-1*max
        ncols=21;15
      ndivs=4
      cbar_format='(i2)'
      xtitle='[ cm/s ]'
;      xrange=[histmin,histmax-binsize/2]
      xrange=round(histmin)
      xrange=[xrange,-1*xrange]
;      xticks=3
    endif else if type eq 'lr' then begin
      title='Histogram of !8S!X'
      ;FOR CTL
        ctblc=11
        maxc=50
        minc=1
        irevc=0
        ncolsc=21;15;31
        setlevs=0
      ;SENS TESTS
        ctbl=70
        irev=1
        max=10
        min=-1*max
        ncols=21;15
      ndivs=4
      cbar_format='(i3)'
      xtitle='[ mK/hPa ]'
      xrange=[histmin,100];histmax-binsize/2]
    endif else if type eq 'rh' then begin
      title='Histogram of RH'
      ;FOR CTL
        ctblc=11
        maxc=60
        minc=1
        irevc=0
        ncolsc=21;15;31
        setlevs=0
      ;SENS TESTS
        ctbl=70
        irev=1
        max=10
        min=-1*max
        ncols=21;15
      ndivs=4
      cbar_format='(i3)'
      xtitle='[ % ]'
      xrange=[histmin,105];histmax-binsize/2]
    endif else if type eq 'ths' then begin
      title='Histogram of TH-S'
      ;FOR CTL
        ctblc=11
        maxc=60
        minc=1
        irevc=0
        ncolsc=21;15;31
        setlevs=0
      ;SENS TESTS
        ctbl=70
        irev=1
        max=10
        min=-1*max
        ncols=21;15
      ndivs=4
      cbar_format='(i3)'
      xtitle='[ K ]'
      xrange=[histmin,380];histmax-binsize/2]
    endif else if type eq 'the' then begin
      title='Histogram of TH-E'
      ;FOR CTL
        ctblc=11
        maxc=60
        minc=1
        irevc=0
        ncolsc=21;15;31
        setlevs=0
      ;SENS TESTS
        ctbl=70
        irev=1
        max=10
        min=-1*max
        ncols=21;15
      ndivs=4
      cbar_format='(i3)'
      xtitle='[ K ]'
      xrange=[histmin,360];histmax-binsize/2]
    endif else if type eq 'th' then begin
      title='Histogram of TH'
      ;FOR CTL
        ctblc=11
        maxc=60
        minc=1
        irevc=0
        ncolsc=21;15;31
        setlevs=0
      ;SENS TESTS
        ctbl=70
        irev=1
        max=10
        min=-1*max
        ncols=21;15
      ndivs=4
      cbar_format='(i3)'
      xtitle='[ K ]'
      xrange=[histmin,380];histmax-binsize/2]
    endif

;    if dirs.cases[ic] eq 'ctl' then begin
    if ic eq 0 then begin
      ctbl=ctblc
      max=maxc
      min=minc
      irev=irevc
      ncols=ncolsc
      if type eq 'lr' then ndivs=5 else if type eq 'wprm' then ndivs=3 else ndivs=3
      setlevs=[string(min,format='(i1)'),string(max*(findgen(ndivs)+1)/ndivs,format='(i2)')] 
    endif

    colors=findgen(ncols)/(ncols-1)*255
    if irev then colors=reverse(colors)
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;  figspecs=create_struct('figname',figname,'title',title,$
;    'cbar_tag',cbar_tag,'cbar_format',cbar_format,'ctbl',ctbl,'ndivs',ndivs,'levels',levels,'colors',colors)

  ;PLOT SPECS
    csize=0.8
    position=[0.22,0.16,0.81,0.92]
    xsize=2.9 & ysize=2.7
    ytitle1='Pressure [ hPa ]'
    yrange=[max(dims.pres),100]

  ;AXES
    x=bins
    y=dims.pres;[ilev]
    if ~keyword_set(xrange) then $
      xrange=[min(x),max(x)]
    if ~keyword_set(yrange) then $
      yrange=[min(y),max(y)]

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  loadct,0,/silent

  yticklen=-0.03
  xticklen=-0.028

  yrange=[max(y),min(ytickv)]
  yticks=n_elements(ytickv)-1

  plot,x,y,/nodata,position=position,/ylog,$
    xstyle=9,xminor=4,xticks=xticks,ystyle=9,yminor=0,yticks=yticks,ytickv=ytickv,$
    yticklen=yticklen,xticklen=xticklen,$
    xrange=xrange,yrange=yrange,$
    xtitle=xtitle,ytitle=ytitle1,$
    charsize=csize,$
    title=title

;VARIABLES

  loadct,ctbl,/silent

loc=where(pdf lt min(levels),nl)
;if nl gt 0 and dirs.cases[ic] ne 'ctl' then pdf[loc]=min(levels)
if nl gt 0 and ic ne 0 then pdf[loc]=min(levels)

  ;PDF
    for i=0,1 do $
      contour,pdf,x,y,/cell_fill,/overplot,$
        levels=levels,c_colors=colors

  loadct,0,/silent

;  if dirs.cases[ic] eq 'ctl' then col0=255 else col0=0
  if ic eq 0 then col0=255 else col0=0

  ;ZERO LINE
  if type eq 'wprm' then $
    plots,[0,0],yrange,linestyle=0,thick=1,color=col0,/data

  ;0C LINE
    if strmatch(dirs.cases[ic],'*LW*',/fold_case) then p_zero=557. else p_zero=568. ; VALID FOR CTL
    plots,!x.crange,replicate(p_zero,2),linestyle=1,thick=1.6,/data

  ;MEAN
    oplot,var_prof_ctl,y,linestyle=0,thick=5,color=col0
;  if type eq 'wprm' then begin
;    oplot,x[it_test],var_t1,linestyle=0,thick=2,color=0
;    oplot,x[it_test],var_t2,linestyle=1,thick=2,color=0
;    oplot,x[it_test],var_t3,linestyle=2,thick=2,color=0
;  endif

  loadct,0,/silent

  ;OLR
;  axis,yaxis=1,ystyle=9,ytitle=ytitle2,charsize=csize,yrange=yrange2,yminor=2,yticklen=yticklen,/save
;  oplot,x,olr,linestyle=2,thick=2,color=0

  ;BOX AROUND PLOT
    for i=0,0 do begin
      plots,!x.crange,replicate(yrange[i],2),linestyle=0,color=0,thick=1,/data
      plots,replicate(!x.crange[i],2),yrange,linestyle=0,color=0,thick=1,/data
    endfor

  ;COLOR BAR
  icbar=1
  if icbar then begin
    loadct,ctbl,/silent
    cpos= [ position[2]+0.018 ,$
            position[1]+0.2 ,$
            position[2]+0.041 ,$
            position[3]-0.2 ]
    colorbar2, colors=colors, range=[min(levels),max(levels)],divisions=ndivs,$
      charsize=csize*0.9, position=cpos, /right, /vertical, title=cbar_tag, $
      annotatecolor='black',format=cbar_format,setlevels=setlevs
    loadct,0,/silent
  endif

  ;LEGEND
  ileg=0
  if ileg then begin
;  if ileg then begin
    csize_fac=0.7
    margin=0.1
    pspacing=2.0 ; length of lines
    spacing=0.8 ; between lines
    leg_str=strupcase(dirs.cases[icplot])
    leg_style=lstyle[icplot]
    leg_thick=replicate(lthick,ncplot)
    leg_color=replicate(0,ncplot)
    legend2,leg_str,linestyle=leg_style,thick=leg_thick,COLORS=leg_color,$
      charsize=csize*csize_fac,/top_legend,/left_legend,/normal,$
      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.16,0.75]
;      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.56,0.35]
  endif

  device,/close

  convert_png,figname,/remove_eps,res=400


endfor ; icase

print,'DONE!!'
end
