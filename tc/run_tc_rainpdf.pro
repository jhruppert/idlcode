; 
; Make PDF from Cartesian (instead of azim) TC output, as function of radius instead of time.
;
; James Ruppert
; 9/2/21
; 
pro run_tc_rainpdf

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

;TIME SELECTION
;  hr_sel=[0,200] ; entire simulation
  t0=37
  hr_sel=[t0,t0+1*24]

;VORTEX TRACKING LEVEL
  psel=700 ; (hPa) level to locate vortex
  izsel=(where(dims.pres eq psel))[0]

;RADIUS THRESHOLD FOR PW PDF
    radius_thresh = 1000 ; km

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
  nt_sav=n_elements(t_ind)

;SELECT TIME-BINNING
;  nt_bin=24;6;12;24
;  nt_bin=48
  nt_bin=nt_sav

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

;PW
;histmax=100d
;binsize=1d
;RAINFALL
histmax=1000d
;binsize=0.1d

  ;CREATE HISTOGRAM
;    nbins=histmax/binsize
;    bins=indgen(nbins)*binsize+binsize*0.5

nbins=22
bins=10^(0.5*(findgen(nbins)-nbins+8))

  pdf=fltarr(dirs.nc,nbins-1)

;for ic=0,dirs.nc-1 do begin
for ic=0,2 do begin
;for ic=0,0 do begin

  print,'CASE: ',strupcase(dirs.cases[ic])

  hr0=0
  if strmatch(dirs.cases[ic],'*36h*') then hr0=36
  if strmatch(dirs.cases[ic],'*48h*') then hr0=48
  if strmatch(dirs.cases[ic],'*60h*') then hr0=60
  if strmatch(dirs.cases[ic],'*72h*') then hr0=72
  if strmatch(dirs.cases[ic],'*84h*') then hr0=84
  if strmatch(dirs.cases[ic],'*96h*') then hr0=96
  if strmatch(dirs.cases[ic],'lwcrf*') then hr0=36
;    i_add=0
;    if strmatch(dirs.cases[ic],'*36h*') then i_add=36
;    if strmatch(dirs.cases[ic],'*48h*') then i_add=48
;    if strmatch(dirs.cases[ic],'*60h*') then i_add=60
;    if strmatch(dirs.cases[ic],'*72h*') then i_add=72
;    if strmatch(dirs.cases[ic],'*84h*') then i_add=84
;    if strmatch(dirs.cases[ic],'*96h*') then i_add=96
;    if strmatch(dirs.cases[ic],'lwcrf*') then i_add=36
;    i_nt = nt_full - i_add
;    it_test=indgen(i_nt)+nt_full-i_nt

  ;TIME SELECTION
    t_offset=max([0,hr0-hr_sel[0]])
    ut_offset=max([0,hr_sel[0]-hr0])
    nt_test=nt_sav-t_offset
    t_ind_test=indgen(nt_test)+ut_offset
    hrs_test=time_hrs[t_ind_test+hr0]
    hr_fil=strtrim(hrs_test[0],2)+'-'+strtrim(hrs_test[nt_test-1],2)+'hr'
;    t_ind_plot=where(hrs_test ge hr_plot[0] and hrs_test le hr_plot[1],nt_plot)
;    hrs_plot=hrs_test[t_ind_plot]
;    hr_tag_plot=string(hrs_plot[0],format='(i3.3)')+'-'+string(hrs_plot[nt_plot-1],format='(i3.3)')+'hr'

  ;VORTEX TRACKING

    ;READ ABSOLUTE VORTICITY
      iv=where(vars.vars eq 'AVOR')
      file=dirs.files_post[ic,iv]
      count=[dims.nx,dims.ny,1,nt_test] & offset=[0,0,izsel,t_ind_test[0]] ; x,y,z,t
      avor=reform(read_nc_var(file,'AVOR',count=count,offset=offset))

    ;SMOOTH
      ixsmooth=round(111./3) ; 1-degree smoothing, run twice
      ismooth=[ixsmooth,ixsmooth,0]
      for i=1,2 do $
        avor=smooth(temporary(avor),ismooth,/edge_truncate,/nan)

    ;MASKING
      if nland gt 0 then avor[land]=!values.f_nan
      if tcname eq 'maria' then avor=wrf_maria_mask(temporary(avor),time[t_ind],hurdat,dims) else stop

    ;VORTEX TRACKING
      vloc=maria_vortex_locate(avor,dims);,/write)
      avor=0

  ;MAIN VARIABLES

  count=[dims.nx,dims.ny,1,nt_test] & offset=[0,0,0,t_ind_test[0]] ; x,y,z,t

  ;PW
;    iv=where(vars.vars eq 'PW')
;    file=dirs.files_post[ic,iv]
;    pw=reform(read_nc_var(file,'PW',count=count,offset=offset))

  ;RAINFALL
    iv=where(vars.vars eq 'rainrate')
    file=dirs.files_post[ic,iv]
    rain=reform(read_nc_var(file,'rainrate',count=count,offset=offset))
;    rain*=1./24 ; mm/d --> mm/h

  ;AVERAGE
;    pw_mean=mean(pw,dimension=2,/nan,/double)


;----TIME SUBSET--------------------

;nt_pw = i_nt/nt_bin
;
;;STRING ARRAY OF TIMES
;  hrs=indgen(i_nt/nt_bin)*nt_bin+nt_bin + i_add
;  hrs=string(hrs,format='(i3.3)')

;----PDF--------------------

  ;RADIUS SPECS
;    radbinsize=25 ; km
;    radmax=800 + radbinsize*0.5 ; km
;    nradbin=round(radmax/radbinsize)
;    radbins=findgen(nradbin)*radbinsize+radbinsize*0.5

;BEGIN TIME LOOP

;  for it=0,nt_pw-1 do begin
;  for it=0,0 do begin
;
;    print,'Hour: ',hrs[it]
;
;    pdf=fltarr(nbins)
;
;    for ir=0,nradbin-1 do begin
;
;      minrad = radbins[ir] - radbinsize*0.5
;      maxrad = radbins[ir] + radbinsize*0.5

      var_sav=rain

      for it=0,nt_bin-1 do begin

        ;TC LOCATION
          ivloc = [ vloc[0,it] , vloc[1,it] ]
          tcrad = radius_tc_ll(ivloc,dims.lon,dims.lat)
;          irad = where((tcrad gt minrad) and (tcrad le maxrad),count)
          irad = where(tcrad le radius_thresh,complement=nan,nradius)

         ivar=reform(var_sav[*,*,it])
         ivar[nan]=!values.f_nan

         var_sav[*,*,it]=ivar
 
      endfor ; itbin

      varreform=reform(var_sav,[1l*dims.nx*dims.ny*nt_test])

;      ihist = histogram( varreform, nbins=nbins, binsize=binsize, min=0.0d )
      ihist=fltarr(nbins-1)
      for ibin=0,nbins-2 do begin
        loc=where(varreform ge bins[ibin] and varreform lt bins[ibin+1],npts)
        ihist[ibin]=npts
      endfor

      ;nradius=70000.d
      nnorm=nradius*nt_test
      ihist *= (1d/nnorm)*1e2

      pdf[ic,*] = ihist

;    endfor ; ir

endfor ; icase

bins=bins[0:nbins-2]

;----CREATE PLOTS--------------------


  figdir=dirs.figdir+'azim_2d/'
  spawn,'mkdir '+figdir,tmp,tmpe
;  figname=figdir+dirs.cases[ic]+'_rain_pdf'
  figname=figdir+'rain_pdf'
;  figname+='_2d_i'+strtrim(nt_bin,2)+'h_'+hrs[it]

  ;SET UP SHADING

;    if var1_str eq 'PW' then begin
      title='PDF of Rainfall'
      ctbl=9
      irev=1
      max=20
      min=0
      ncols=40;15
      ndivs=4
      cbar_tag='[ % ]'
      cbar_format='(i2)'
;    endif

;    colors=findgen(ncols)/(ncols-1)*255
;    if irev then colors=reverse(colors)
;    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;  figspecs=create_struct('figname',figname,'title',title,$
;    'cbar_tag',cbar_tag,'cbar_format',cbar_format,'ctbl',ctbl,'ndivs',ndivs,'levels',levels,'colors',colors)

  ;PLOT SPECS
    csize=0.8
    position=[0.14,0.19,0.89,0.89]
    xsize=4.2 & ysize=2
    xtitle='Rain Rate [ mm d!U-1!N ]'
;    xrange=[1e-1,300]
    ytitle1='Frequency [ % ]'
    yrange=[1e-1,1e1]
;yrange=[1e2,1e6];max(pdf)]
ylog=1
yrange=[-0.1,0.1];max(pdf)]
ylog=0

  ;AXES
    x=bins
    y=indgen(5)
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

  plot,x,y,/nodata,/xlog,ylog=ylog,position=position,$
    xstyle=9,ystyle=9,$
    yticklen=yticklen,xticklen=-0.033,$
    xrange=xrange,yrange=yrange,$;yminor=2,$
    xtitle=xtitle,ytitle=ytitle1,$
    charsize=csize,$
    title=title

;VARIABLES

  loadct,ctbl,/silent

  ;PDF
  loadct,0,/silent
;  for ic=0,dirs.nc-1 do $
;    oplot,x,reform(pdf[ic,*]),linestyle=ic,thick=2,color=0
  for ic=1,dirs.nc-1 do begin
    ivar = reform(pdf[ic,*]) / reform(pdf[0,*]) - 1
    oplot,x,ivar,linestyle=ic,thick=2,color=0
  endfor

  loadct,0,/silent

  ;BOX AROUND PLOT
;    for i=0,0 do begin
;      plots,!x.crange,replicate(!y.crange[i],2),linestyle=0,color=0,thick=1,/data
;      plots,replicate(!x.crange[i],2),!y.crange,linestyle=0,color=0,thick=1,/data
;    endfor

  ;COLOR BAR
  icbar=0
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
  ileg=1
  if ileg then begin
;  if ileg then begin
    csize_fac=0.7
    margin=0.1
    pspacing=2.0 ; length of lines
    spacing=0.8 ; between lines
    ncplot=dirs.nc
    icplot=indgen(ncplot)
    leg_str=strupcase(dirs.cases[icplot])
    leg_style=indgen(ncplot)
    leg_thick=replicate(2,ncplot)
    leg_color=replicate(0,ncplot)
    legend2,leg_str,linestyle=leg_style,thick=leg_thick,COLORS=leg_color,$
      charsize=csize*csize_fac,/top_legend,/left_legend,/normal,$
      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.6,0.75]
;      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.56,0.35]
  endif

  device,/close

  convert_png,figname,/remove_eps,res=200


;endfor ; it

;endfor ; icase

print,'DONE!!'
end
