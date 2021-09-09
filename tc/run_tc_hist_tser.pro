; 
; Make PDF from Cartesian (instead of azim) TC output.
;
; James Ruppert
; 4/29/19
; 
pro run_tc_hist_tser

tcname='maria';'edouard';
tcyear='2017'
hurdat=read_hurdat(tcname,tcyear)

;subdir='moving_nest/'+tcname
subdir='static_nest/'+tcname
tc_sim_config, subdir, dirs=dirs, dims=dims, vars=vars;, /verbose
dirs.figdir+=tcname+'/'

;TIME SELECTION
  hr_sel=[0,200] ; entire simulation
;  hr_sel=[24,96]
  hr_tag=strtrim(hr_sel[0],2)+'-'+strtrim(hr_sel[1],2)+'hr'

;VORTEX TRACKING LEVEL
  psel=700 ; (hPa) level to locate vortex
  izsel=(where(dims.pres eq psel))[0]

;---SETTINGS---------------------

;RADIUS THRESHOLD FOR HISTOGRAM
  radius_thresh = 500 ; km

;FIGURE TYPE
  type='lh';'wprm';'pw';

;----TIME SPECS--------------------

;FULL TIME SERIES
  time=dims.time
  nt_full=dims.nt-1
  npd=dims.npd
  nhrs=1.*nt_full*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nhrs)

;TIME SELECTION
  t_ind=where((time_hrs ge hr_sel[0]) and (time_hrs le hr_sel[1]))
  t_ind=t_ind[where(t_ind le nt_full-1)]
  nt_sel=n_elements(t_ind)

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
;for ic=0,2 do begin

  print,'CASE: ',strupcase(dirs.cases[ic])

    i_nt=nt_full
    if strmatch(dirs.cases[ic],'*36h*') then i_nt-=36
    if strmatch(dirs.cases[ic],'*24h*') then i_nt-=24
    if strmatch(dirs.cases[ic],'*48h*') then i_nt-=48
    if strmatch(dirs.cases[ic],'*60h*') then i_nt-=60
    if strmatch(dirs.cases[ic],'*72h*') then i_nt-=72
    if strmatch(dirs.cases[ic],'*84h*') then i_nt-=84
    it_test=indgen(i_nt)+nt_full-i_nt

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
      if tcname eq 'maria' then avor=wrf_maria_mask(temporary(avor),time[it_test],hurdat,dims) else stop
  
    ;VORTEX TRACKING
      vloc=maria_vortex_locate(avor,dims);,/write)
      avor=0

  ;MAIN VARIABLES

  if type eq 'pw' then begin
  ;PW
    count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
    iv=where(vars.vars eq 'PW')
    file=dirs.files_post[ic,iv]
    var=reform(read_nc_var(file,'PW',count=count,offset=offset))
  endif else if type eq 'wprm' then begin
  ;W
    ilev=(where(dims.pres eq '500'))[0]
    count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,ilev,0] ; x,y,z,t
;    ilev=where(dims.pres ge 100,nlev)
;    count=[dims.nx,dims.ny,nlev,i_nt] & offset=[0,0,ilev[0],0] ; x,y,z,t
    iv=where(vars.vars eq 'W')
    file=dirs.files_post[ic,iv]
    var=reform(read_nc_var(file,'W',count=count,offset=offset))
    iv=where(vars.vars eq 'QCLOUD')
    file=dirs.files_post[ic,iv]
    qc=reform(read_nc_var(file,'QCLOUD',count=count,offset=offset))
    iv=where(vars.vars eq 'QICE')
    file=dirs.files_post[ic,iv]
    qi=reform(read_nc_var(file,'QICE',count=count,offset=offset))
    cloud=(qc+qi)*1e6 & qc=0 & qi=0
;    icloud=where((cloud ge 3) and (var ge 1),ncloud,complement=comp)
    icloud=where((cloud ge 3) and (abs(var) ge 1),ncloud,complement=comp)
;print,1d*ncloud/n_elements(comp)*100
    var[comp]=!values.f_nan
  endif else if type eq 'lh' then begin
  ;LH
    count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
    iv=where(vars.vars eq 'LH')
    file=dirs.files_post[ic,iv]
    var=reform(read_nc_var(file,'LH',count=count,offset=offset))
  endif


;----CREATE PDF--------------------


  if type eq 'pw' or type eq 'lh' then begin

    if type eq 'pw' then begin
      histmin=0d
      histmax=100d
      binsize=1.0d
      nbins=histmax/binsize
    endif else if type eq 'lh' then begin
      histmin=0d;100d
      histmax=600d
      binsize=25.0d
      nbins=histmax/binsize
    endif

  ;CREATE HISTOGRAM
;    nbins=histmax/binsize
    nbins=round(nbins)
    bins=indgen(nbins)*binsize;+binsize*0.5
    pdf=fltarr(i_nt,nbins)
    var_mean=fltarr(i_nt)
    varian=fltarr(i_nt)

    for it=0,i_nt-1 do begin

      ;TC LOCATION
        ivloc = [ vloc[0,t_ind[it]] , vloc[1,t_ind[it]] ]
        tcrad = radius_tc_ll(ivloc,dims.lon,dims.lat)
        irad = where(tcrad le radius_thresh,count)

      ivar = reform(var[*,*,it])

      pdf[it,*] = histogram( ivar[irad], nbins=nbins, binsize=binsize, min=histmin ) * 1d/count*100
      var_mean[it] = mean( ivar[irad], /nan, /double)
      varian[it] = variance( ivar[irad], /nan, /double)

    endfor

  endif else if type eq 'wprm' then begin

  ;CREATE HISTOGRAM
    histmax=50d & binsize=1.0d
    nbins=histmax/binsize
    bins=indgen(nbins)*binsize;-histmax
    pdf=fltarr(i_nt,nbins)
    var_mean=fltarr(i_nt)
    varian=fltarr(i_nt)

    var_t1=fltarr(i_nt)
    var_t2=fltarr(i_nt)
    var_t3=fltarr(i_nt)

    pdf[*]=!values.f_nan
    for it=0,i_nt-1 do begin

      ;TC LOCATION
        ivloc = [ vloc[0,t_ind[it]] , vloc[1,t_ind[it]] ]
        tcrad = radius_tc_ll(ivloc,dims.lon,dims.lat)
        irad = where(tcrad le radius_thresh,count_rad)

      ivar = reform(var[*,*,it])

;      pdf[it,*] = histogram( ivar[irad], nbins=nbins, binsize=binsize, min=0.0d ,/nan) * 1d/count*100
;      var_mean[it] = mean( ivar[irad], /nan, /double)
;      varian[it] = variance( ivar[irad], /nan, /double)

      fin=where(finite(ivar),ntot)

      loc=where(ivar ge  5.0,count,complement=comp)
;      inorm=count/ntot
      inorm=count/count_rad
      var_t1[it] = 100d*count/count_rad;ntot
      loc=where(ivar ge 10.0,count,complement=comp)
;      inorm=count/ntot
      inorm=count/count_rad
      var_t2[it] = 100d*count/count_rad;ntot
      loc=where(ivar ge 15.0,count,complement=comp)
;      inorm=count/ntot
      inorm=count/count_rad
      var_t3[it] = 100d*count/count_rad;ntot

    endfor

    ;SMOOTH
      ismooth=1
      nsmooth=3
      if ismooth then begin
        var_t1=smooth(temporary(var_t1),nsmooth,/edge_truncate,/nan)
        var_t2=smooth(temporary(var_t2),nsmooth,/edge_truncate,/nan)
        var_t3=smooth(temporary(var_t3),nsmooth,/edge_truncate,/nan)
      endif

  endif


;----CREATE PLOTS--------------------


  figdir=dirs.figdir
  spawn,'mkdir '+figdir,tmp,tmpe

  ;SET UP SHADING

    if type eq 'pw' then begin
      title='PDF of Column Water Vapor'
      ctbl=9
      irev=1
      max=20
      min=0
      ncols=40;15
      ndivs=4
      cbar_tag='[ % ]'
      cbar_format='(i2)'
      ytitle1='PW [ mm ]'
      yrange=[20,85];75]
      ylog=0
    endif else if type eq 'wprm' then begin
      title='PDF of Cloudy Updrafts'
      ctbl=9
      irev=1
      max=0.5
      min=0
      ncols=40;15
      ndivs=4
      cbar_tag='[ % ]'
      cbar_format='(i2)'
      ;ytitle1='W-variance [ m!U2!N s!U-2!N ]'
      ;ytitle1='W [ m s!U-1!N ]'
      ytitle1='Percent above thresh [ % ]'
      yrange=[0.001,10]
      ;yrange=[0,0.2]
      ylog=1
    endif else if type eq 'lh' then begin
      title='PDF of LH'
      ctbl=9
      irev=1
      max=50
      min=0
      ncols=21;15
      ndivs=4
      cbar_tag='[ % ]'
      cbar_format='(i2)'
      ytitle1='LH [ W/m!U2!N ]'
      yrange=[histmin,histmax]
      ylog=0
    endif

  figname=figdir+type+'_pdf_'+dirs.cases[ic]

    colors=findgen(ncols)/(ncols-1)*255
    if irev then colors=reverse(colors)
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;  figspecs=create_struct('figname',figname,'title',title,$
;    'cbar_tag',cbar_tag,'cbar_format',cbar_format,'ctbl',ctbl,'ndivs',ndivs,'levels',levels,'colors',colors)

  ;PLOT SPECS
    csize=0.8
    position=[0.14,0.20,0.89,0.89]
    xsize=4.2 & ysize=2
;    xtitle='Time [ hr ]'
    xtitle='Date in Sept [ UTC ]'
    ytitle2='OLR [ W m!U-2!N ]'
    yrange2=[60,310]

  ;AXES
    x=time_hrs
    y=bins
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

  yticklen=-0.014

  xticks=5
  xtickv=indgen(xticks+1)*24+12
  xtickname=strtrim(indgen(xticks+1)+15,2)

  plot,x,y,/nodata,position=position,ylog=ylog,$
    xstyle=9,ystyle=9,$
    yticklen=yticklen,xticklen=-0.033,$
    xticks=xticks,xtickv=xtickv,xtickname=xtickname,$
    xrange=xrange,yrange=yrange,yminor=2,$
    xtitle=xtitle,ytitle=ytitle1,$
    charsize=csize,$
    title=title

;VARIABLES

  loadct,ctbl,/silent

  ;PDF
    for i=0,1 do $
      contour,pdf,x[it_test],y,/cell_fill,/overplot,$
        levels=levels,c_colors=colors

  ;MEAN
  loadct,0,/silent
  if type eq 'pw' then begin
    oplot,x[it_test],var_mean,linestyle=0,thick=2,color=0
    oplot,x[it_test],varian,linestyle=1,thick=2,color=0
  endif else if type eq 'wprm' then begin
    oplot,x[it_test],var_t1,linestyle=0,thick=2,color=0
    oplot,x[it_test],var_t2,linestyle=1,thick=2,color=0
    oplot,x[it_test],var_t3,linestyle=2,thick=2,color=0
  endif

  loadct,0,/silent

  ;OLR
;  axis,yaxis=1,ystyle=9,ytitle=ytitle2,charsize=csize,yrange=yrange2,yminor=2,yticklen=yticklen,/save
;  oplot,x,olr,linestyle=2,thick=2,color=0

  ;BOX AROUND PLOT
    for i=0,0 do begin
      plots,!x.crange,replicate(!y.crange[i],2),linestyle=0,color=0,thick=1,/data
      plots,replicate(!x.crange[i],2),!y.crange,linestyle=0,color=0,thick=1,/data
    endfor

  ;COLOR BAR
  icbar=1
  if icbar then begin
    loadct,ctbl,/silent
    cpos= [ position[2]+0.018 ,$
            position[1]+0.1 ,$
            position[2]+0.031 ,$
            position[3]-0.1 ]
    colorbar2, colors=colors, range=[min(levels),max(levels)],divisions=ndivs,$
      charsize=csize*0.9, position=cpos, /right, /vertical, title=cbar_tag, $
      annotatecolor='black',format=cbar_format
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

  convert_png,figname,/remove_eps,res=200


endfor ; icase

print,'DONE!!'
end
