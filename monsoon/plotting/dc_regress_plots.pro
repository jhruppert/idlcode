; 
; Plotting call routines for diurnal regression analysis using IMERG (Myanmar/Meiyu IMERG data)
; Called by run_imerg_coast_dc_regress.pro
;
; James Ruppert
; 1/4/21
; 
pro dc_regress_plots, dirs, rtag, dims, var_plot, npd_imerg, local, ltim_imerg, psel_era, $
  xcross, ycross, width, xhov_im, $
  a1_all, a1_coast, a1_offshore, a1_smd, a1_cdiff, $;a1_onshore, $
  do_wind, u_all=u_all, u_coast=u_coast, u_offshore=u_offshore, u_onshore=u_onshore, $
  lonim=lonim, latim=latim, $
  lonera=lonera, latera=latera

  print,'Beginning plotting...'

;----REGRESSION PLOTS--------------------

;iplot_reg=0
;if iplot_reg then begin

;for itest=0,0 do begin
for itest=0,4 do begin

  print,'  iTest ',itest

  if var_plot eq 'rain' then begin
    var_str='RAINNC'
    setmax=200 & setmin='0.'
;FOR ALL-TIME-MEAN
;setmax=30 & setmin='0.'
    cbform='(i3)'
    cbtag='Rain [ mm d!U-1!N ]'
    lon=lonim
    lat=latim
  endif else if var_plot eq 'pw' then begin
    var_str=var_plot
    setmax=15 & setmin=-15
    cbform='(i3)'
    cbtag='PW [ mm ]'
    lon=eralon
    lat=eralat
  endif

  if itest eq 0 then begin
    if var_plot eq 'rain' then $
      var_plt=a1_all $
    else if var_plot eq 'pw' then message,"Haven't set up PW yet!";$
;      var_plt=pw_all
    if do_wind then begin
      u_plt=u_all
      v_plt=v_all
    endif
    itag='All'
  endif else if itest eq 1 then begin
    if var_plot eq 'rain' then $
      var_plt=a1_coast ;$
;    else if var_plot eq 'pw' then $
;      var_plt=pw_coast
    if do_wind then begin
      u_plt=u_coast
      v_plt=v_coast
    endif
    itag='Coast'
  endif else if itest eq 2 then begin
    if var_plot eq 'rain' then $
      var_plt=a1_offshore ;$
;    else if var_plot eq 'pw' then $
;      var_plt=pw_offshore
    if do_wind then begin
      u_plt=u_offshore
      v_plt=v_offshore
    endif
    itag='Offshore'
;  endif else if itest eq 3 then begin
;    if var_plot eq 'rain' then $
;      var_plt=a1_onshore ;$
;;    else if var_plot eq 'pw' then $
;;      var_plt=pw_onshore
;    if do_wind then begin
;      u_plt=u_onshore
;      v_plt=v_onshore
;    endif
;    itag='Onshore'
  endif else if itest eq 3 then begin
    if var_plot eq 'rain' then $
      var_plt=a1_smd ;$
;    else if var_plot eq 'pw' then $
;      var_plt=pw_smd
    if do_wind then begin
      u_plt=u_smd
      v_plt=v_smd
    endif
    itag='SMD'
  endif else if itest eq 4 then begin
    if var_plot eq 'rain' then $
      var_plt=a1_cdiff ;$
;    else if var_plot eq 'pw' then $
;      var_plt=pw_cdiff
    if do_wind then begin
      u_plt=u_cdiff
      v_plt=v_cdiff
    endif
    itag='cdiff'
  endif

  ;REGRESSION MAPS

    myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, setndivs=4, set_cint=3
    figdir=dirs.figdir+rtag+'/'
    figspecs=create_struct(figspecs,'figname',' ')
    figspecs.cbar_format=cbform
    figspecs.cbar_tag=cbtag
  
    ntplot=1;4;48
    nskip=npd_imerg/ntplot
    for ihr=0,npd_imerg-1,nskip do begin
;for ihr=0,0 do begin

      print,'    Map ',ihr,' of ',npd_imerg-1
  
      it_str=string(ihr,format='(i2.2)')
      lt_str=string(ltim_imerg[ihr],format='(f4.1)')

      figspecs.title='Regression: "'+itag+'" ('+lt_str+' LT)'
      figspecs.figname=figdir+rtag+'_reg_'+var_plot+'_'+strlowcase(itag)+'_'+strtrim(psel_era,2)+'_'+it_str+'it'
  
      ivar=reform(var_plt[*,*,ihr])

      if do_wind then begin
        ihr_era=1.*ihr*npd_era/npd_imerg
        ;only update winds if at a new ERA time step
        if ihr_era eq round(ihr_era) then begin
          iu_plt=reform(u_plt[*,*,ihr_era])
          iv_plt=reform(v_plt[*,*,ihr_era])
        endif
        wind=create_struct('u',iu_plt,'v',iv_plt,'x',eralon,'y',eralat)
        wspd=sqrt(iu_plt^2+iv_plt^2);*10.
        cvar=create_struct('cvar',wspd,'x',eralon,'y',eralat)
      endif
  
      wrf_myanmar_map_plot, dirs, ivar, lon, lat, figspecs, wind=wind, cvar=cvar
  
    endfor

  ;DIURNAL HOVMOLLERS

    if var_plot eq 'rain' then begin
      var_str='RAINNC'
      setmax=8 & setmin='0.'
;FOR ALL-TIME-MEAN
setmax=1.5 & setmin='0.'
setmax=3. & setmin=-1.*setmax
      cbform='(i1)'
cbform='(f4.1)'
      cbtag='Rain [ mm h!U-1!N ]'
      lon=lonim
      lat=latim
;    endif else if var_plot eq 'pw' then begin
;      var_str=var_plot
;      setmax=15 & setmin=-15
;      cbform='(i3)'
;      cbtag='PW [ mm ]'
;      lon=eralon
;      lat=eralat
    endif
  
    myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, setndivs=3, set_cint=3
    figdir=dirs.figdir+rtag+'/'
    figspecs=create_struct(figspecs,'figname',' ')
    figspecs.cbar_format=cbform
    figspecs.cbar_tag=cbtag

    ;INTERPOLATE ONTO HOVMOLLER CROSS
    cross=cross_diag(var_plt,lonim,latim,width,x_bounds=xcross,y_bounds=ycross,lonout=xlon_im,latout=xlat_im)
    imerg=mean(cross,dimension=2,/nan,/double) & cross=0
    imerg*=1./24 ; mm/d --> mm/h

    figspecs.figname=figdir+'../hov/'+rtag+'_reg_'+var_plot+'_'+strlowcase(itag)+'_'+strtrim(psel_era,2)
;figspecs.figname=figdir+'../hov/imerg_stddev_'+var_plot+'_'+strtrim(psel_era,2)
;figspecs.figname=figdir+'../hov/imerg_dcomptest_'+var_plot
    figspecs.title='Regression: "'+itag+'"'

    ;REORDER TO PUT 0 LT AT FIRST INDEX
    imerg=shift(temporary(imerg),[0,1.*local*npd_imerg/24])
    ;ltim_imerg_shift=shift(ltim_imerg,local*npd_imerg/24)

  ;    if iwind then begin
  ;      ltim_era=findgen(npd_era)*24/npd_era+local
  ;      for it=0,npd_era-1 do ltim_era[it]-=24*(ltim_era[it] ge 24)
  ;      ltim_era=shift(temporary(ltim_era),local*npd_era/24)
  ;      u=shift(temporary(u),[0,local*npd_era/24])
  ;      v=shift(temporary(v),[0,local*npd_era/24])
  ;    endif
  
  ;    if iwind then wind=create_struct('u',iu,'v',iv)
  
      wrf_monsoon_dcomp_hov, dirs, figspecs, imerg, indgen(npd_imerg), xhov_im, wind=wind, cvar=cvar

endfor ; itest
;endif


print,'DONE PLOTTING!!'
end
