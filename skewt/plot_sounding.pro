;---------------------------------------------------------------------------------------------------------------
;  Subroutine to plot a sounding provided with the thermodynamic data.
;  
;  Called by dyn_diurnal_var.pro
;  
;  CALL FORMAT:
;  sounding_average, tempc, dewpc, TOP, FIGNAME, TITLE, FIGPATH=FIGPATH, tlp=tlp
;  
;  INPUT
;  TOP:       Upper bound (in hPA) on skewt.
;  FIGNAME:   Names of figure including directory path.
;  TITLE: Station location as it will be used in title.
;  
;  PRES2, etc.: these will be overplotted in gray curves.
;  
;  OPTIONS
;  PLOTPARCELT: Plot the dry-bulb parcel curve.
;  PLOTPARCELV: Plot the virtual parcel temperature curve.
;  
;  
;  James Ruppert (ruppert@atmos.colostate.edu)
;  10/4/2013
;---------------------------------------------------------------------------------------------------------------
pro plot_sounding, pres, tempc, dewpc, TOP, figname=FIGNAME, title=TITLE, tlp=tlp, plotparcelt=plotparcelt, $
  plotparcelv=plotparcelv,$
  wspd=wspd, wdir=wdir, hght=hght, mld=mld,$
  pres2=pres2,tempc2=tempc2,dewpc2=dewpc2

on_error,0

nan=!values.f_nan

;stop  
;--PLOT SOUNDING------------------------------------------------------------------------------------------------
  
  set_plot,'ps'
  device,filename=figname+'.eps',/encapsulated,/color,bits=8,xsize=5.8,ysize=5.6,/inches,$
    /helvetica
  ;xyouts,0,0,'!5'
  loadct,4,/silent
  
  if top eq 600 then begin
    tmin=15
    tmax=35
  endif else if top eq 400 then begin
    tmin=0
    tmax=40
  endif else begin
    tmin=-40
    tmax=40
  endelse
  
  skewt,[tmin,tmax],[1040.,top],title=title,position=[0.135,0.105,0.94,0.94];position=[0.135,0.105,0.83,0.94]
  
  ind=where(finite(pres) eq 1 and pres ge top,npts)
  loadct,4,/silent
  
  if ~keyword_set(wspd) then begin
    wspd=replicate(nan,(size(tempc))[1],(size(tempc))[2])
    wdir=wspd
  endif
  
  if ~keyword_set(hght) then begin
    hght=replicate(nan,(size(tempc))[1],(size(tempc))[2])
  endif
  
  ;plot_skewt,tempc[ind],dewpc[ind],wspd[ind],wdir[ind],pres[ind],hght=hght[ind],$
;  plot_skewt,tempc[ind],dewpc[ind],pres[ind],$
  plot_skewt,tempc,dewpc,pres,$
  ;plot_skewt,tempv,dewpc[ind],wspd[ind],wdir[ind],pres[ind],hght=hght[ind],$
    col_t=150,col_dewpt=40,thick=7
  
  nlevs2=n_elements(pres2)
  if nlevs2 gt 0 then begin
    ind2=where(finite(pres2) eq 1 and pres2 ge top,npts)
    wspd2=replicate(!values.f_nan,nlevs2) & wdir2=wspd2
    loadct,0,/silent
    ;plot_skewt,tempc2[ind2],dewpc2[ind2],wspd2[ind2],wdir2[ind],pres2[ind],$
;    plot_skewt,tempc2[ind2],dewpc2[ind2],pres2[ind],$
    plot_skewt,tempc2,dewpc2,pres2,$
      col_t=150,col_dewpt=150,thick=3
  endif
  
  ;CAPE, LFC, AND DCAPE
;    mixr=mixr_sat(dewpc,pres)
;    if n_elements(where(finite(mixr[10:30]))) lt 18 then begin ; CHECK LOWER TROP
;      cape='NAN' & lfc=nan & dcape='NAN'
;    endif else begin
;      loc1=where(finite(mixr) and finite(tempc),nloc1)
;      ;loc1=loc1[1:n_elements(loc1)-1]
;      mixr=mixr/1000.
;      cape=cape_sound(pres[loc1],tempc[loc1]+273.15,mixr[loc1],/surf,tlp=tlp,tlvp=tlvp,lfc=lfc)
;      ;parcelcalc,tempc[loc1],dewpc[loc1],pres[loc1],cape2,lfc=lfctmp,/virtual
;      if cape le 0 then begin
;        cape=0. & lfc=nan
;      endif
;      if finite(cape) then cape=strtrim(round(cape),2)+' J/kg' else cape='NAN'
;      ;PLOT PARCEL CURVES
;      nans1=replicate(nan,nloc1)
;      if finite(cape) and keyword_set(plotparcelt) then plot_skewt,reform(tlp[0,*])-273.15,nans1,nans1,nans1,pres[loc1],thick=2.5
;      if finite(cape) and keyword_set(plotparcelv) then plot_skewt,reform(tlvp[0,*])-273.15,nans1,nans1,nans1,pres[loc1],thick=2.5
;    endelse
  
  ;PUT INFO ON SKEWT
    
    ;xloc=0.145 & yloc=0.91 & dy=0.035
    xloc=0.93 & yloc=0.90
    charsize=1.4 & charthick=3.4
;    xyouts,xloc,yloc,'CAPE = '+cape,$
;      align=1,charsize=charsize,charthick=charthick,/normal
;    if finite(lfc) then begin
;      plots,[20,35],[lfc,lfc],linestyle=0,thick=2.5,color=0,/data
;      xyouts,35,lfc-10,'LFC',align=1,charsize=charsize*.7,charthick=charthick,/data
;    endif
    
    if keyword_set(mld) then begin
      
      plots,[20,35],[mld,mld],linestyle=0,thick=2.5,color=0,/data
      xyouts,35,mld-10,'MLD',align=1,charsize=charsize*.7,charthick=charthick,/data
      
    endif
    
  device,/close
  convert_png,figname,res=130,/remove_eps


;print,'  DONE!!'
;stop

end