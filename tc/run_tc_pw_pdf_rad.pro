; 
; Make PDF from Cartesian (instead of azim) TC output, as function of radius instead of time.
;
; James Ruppert
; 4/29/19
; 
pro run_tc_pw_pdf_rad

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

;VORTEX TRACKING LEVEL
  psel=700 ; (hPa) level to locate vortex
  izsel=(where(dims.pres eq psel))[0]


;RADIUS THRESHOLD FOR PW PDF
    radius_thresh = 500 ; km

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

;SELECT TIME-BINNING
  nt_bin=24;6;12;24

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

    i_add=0
    if strmatch(dirs.cases[ic],'*36h*') then i_add=36
    if strmatch(dirs.cases[ic],'*24h*') then i_add=24
    if strmatch(dirs.cases[ic],'*48h*') then i_add=48
    if strmatch(dirs.cases[ic],'*60h*') then i_add=60
    if strmatch(dirs.cases[ic],'*72h*') then i_add=72
    if strmatch(dirs.cases[ic],'*84h*') then i_add=84
    i_nt = nt_full - i_add
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

  count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t

  ;PW
    iv=where(vars.vars eq 'PW')
    file=dirs.files_post[ic,iv]
    pw=reform(read_nc_var(file,'PW',count=count,offset=offset))

  ;AVERAGE
;    pw_mean=mean(pw,dimension=2,/nan,/double)


;----TIME SUBSET--------------------

nt_pw = i_nt/nt_bin

;STRING ARRAY OF TIMES
  hrs=indgen(i_nt/nt_bin)*nt_bin+nt_bin + i_add
  hrs=string(hrs,format='(i3.3)')

;----PDF of PW--------------------


  ;CREATE HISTOGRAM
    histmax=100d & binsize=1d
    nbins=histmax/binsize
    pwbins=indgen(nbins)*binsize;+binsize*0.5

  ;RADIUS SPECS
    radbinsize=25 ; km
    radmax=800 + radbinsize*0.5 ; km
    nradbin=round(radmax/radbinsize)
    radbins=findgen(nradbin)*radbinsize+radbinsize*0.5

;BEGIN TIME LOOP

  for it=0,nt_pw-1 do begin

    print,'Hour: ',hrs[it]

    it_bin = it*nt_bin + indgen(nt_bin)

    pw_pdf=fltarr(nradbin,nbins)
    pw_mean=fltarr(nradbin)

    for ir=0,nradbin-1 do begin

      minrad = radbins[ir] - radbinsize*0.5
      maxrad = radbins[ir] + radbinsize*0.5

      for itbin=0,nt_bin-1 do begin
  
        ;TC LOCATION
          ivloc = [ vloc[0,it_bin[itbin]] , vloc[1,it_bin[itbin]] ]
          tcrad = radius_tc_ll(ivloc,dims.lon,dims.lat)
          irad = where((tcrad gt minrad) and (tcrad le maxrad),count)
 
        ipw=reform(pw[*,*,it_bin[itbin]])
 
        if itbin eq 0 then $
          pwsav = reform(ipw[irad],count) $
        else $
          pwsav = [ pwsav , reform(ipw[irad],count) ]

      endfor ; itbin

      npts = n_elements(pwsav)
      pw_pdf[ir,*] = histogram( pwsav, nbins=nbins, binsize=binsize, min=0.0d ) * (1d/npts)*1e2
      pw_mean[ir] = mean(pwsav,/nan,/double)

    endfor ; ir


;----CREATE PLOTS--------------------


  figdir=dirs.figdir+'azim_2d/'
  spawn,'mkdir '+figdir,tmp,tmpe
  figname=figdir+dirs.cases[ic]+'_pw_pdf'
  figname+='_2d_i'+strtrim(nt_bin,2)+'h_'+hrs[it]

  ;SET UP SHADING

;    if var1_str eq 'PW' then begin
      title='PDF of Column Water Vapor'
      ctbl=9
      irev=1
      max=20
      min=0
      ncols=40;15
      ndivs=4
      cbar_tag='[ % ]'
      cbar_format='(i2)'
;    endif

    colors=findgen(ncols)/(ncols-1)*255
    if irev then colors=reverse(colors)
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;  figspecs=create_struct('figname',figname,'title',title,$
;    'cbar_tag',cbar_tag,'cbar_format',cbar_format,'ctbl',ctbl,'ndivs',ndivs,'levels',levels,'colors',colors)

  ;PLOT SPECS
    csize=0.8
    position=[0.14,0.19,0.89,0.89]
    xsize=4.2 & ysize=2
    xtitle='Radius [ km ]'
    xrange=[800,0]
    ytitle1='PW [ mm ]'
    ytitle2='OLR [ W m!U-2!N ]'
    yrange=[20,75]
    yrange2=[60,310]

  ;AXES
    x=radbins
    y=pwbins
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

  plot,x,y,/nodata,position=position,$
    xstyle=9,ystyle=9,$
    yticklen=yticklen,xticklen=-0.033,$
    xrange=xrange,yrange=yrange,yminor=2,$
    xtitle=xtitle,ytitle=ytitle1,$
    charsize=csize,$
    title=title

;VARIABLES

  loadct,ctbl,/silent

  ;PW PDF
    for i=0,1 do $
      contour,pw_pdf,x,y,/cell_fill,/overplot,$
        levels=levels,c_colors=colors

  ;MEAN PW
  loadct,0,/silent
    oplot,x,pw_mean,linestyle=0,thick=2,color=0

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


endfor ; it

endfor ; icase

print,'DONE!!'
end
