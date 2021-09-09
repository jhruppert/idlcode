; 
; Create y-averaged cross sections from ICON monsoon simulations.
; 
; James Ruppert
; 16.12.17
; 
pro yavg_cross_driver

;SOME SETTINGS

imist=0

itimeheight=0 ; plot time-height? otherwise cross sections

;exp_all=['24h','gscp1','gscp2','gscp1_303k']
;nexp=n_elements(exp_all)

;for iexp=0,nexp-1 do begin
;for iexp=0,0 do begin

 exp = 'terra_t3'

;  exp=exp_all[iexp]
  exp_name = exp;'nh_rce_'+strtrim(exp,2)

  print,'Running ',strupcase(exp)

;time range to include for diurnal composite [ in d ]
  t_sel=[15,20]

;EXPERIMENT SPECS
  icon_rce_exp_presets,exp_name,jult_fil0=jult_fil0,nfil=nfil,ndays=ndays,$
    delta_file=delta_file,grid_fil=grid_fil,npday=npday,$
    moddir=moddir,datdir=datdir,figdir=figdir,imist=imist,$
    npd2=npd2

  figdir+=exp+'/'

;WHOLE DAYS ONLY
  if n_elements(t_sel) gt 1 then $
    t_sel=[t_sel[0],t_sel[1]-1d/npday]

  ;GET TIME ARRAY
    ;icon_time_multi,datfils_3d,time=time
    time=dindgen(npday*ndays)/npday

;Y-AVG DATA OUTPUT FILES
  spawn,'ls '+moddir+'*yavg.nc',datfils,err
  if keyword_set(err) then stop
  t_ind1=max(where(time le t_sel[0]))
  if n_elements(t_sel) eq 2 then begin
    t_ends=[ max(where(time le t_sel[0])) , max(where(time le t_sel[1])) ]
    t_ind1=indgen(t_ends[1]-t_ends[0]+1)+t_ends[0]
  endif
  ntpfil=delta_file*npday
  t_ind_multi,time,time[t_ind1],ntpfil,$
    t_ind=t_ind,fil_ind=fil_ind

;HEIGHT
  icon_file_multi,nfil,jult_fil0,delta_file,'3d',moddir,exp_name,$
    datfils=datfils_3d
  nz=75
  hghtw=reform(read_nc_var(datfils_3d[0],'z_ifc',count=[1,nz+1])) ; [ m ]
  hght=interpol(hghtw,findgen(nz+1),findgen(nz)+0.5)


;START VARIABLE LOOP


  ;PLOT VARIABLES
    vars=['pres','w','u','rad','sw','lw','tmp',$
    ;        0    1   2    3    4    5     6
          'q1_ddt','q1_uddx','q1_wddz','q1_total','q1_qc']
    ;        7         8         9         10       11
    nvar=n_elements(vars)


;for ivar=0,nvar-1 do begin
for ivar=1,1 do begin

  plot_var = vars[ivar]
  print,'Plot var: ',plot_var

  ifigdir = figdir

  if itimeheight then begin
    if plot_var eq 'pres' then message,'Pick 3D var for time-height!'
  endif else $
    if plot_var ne 'pres' then ifigdir+=plot_var+'/'

  plot_type='3d'

  if plot_var eq 'pres' then begin

    message,"Haven't done this one yet"

    plot_type='2d'

    pres_plot='theta';'flux';'pres';

  endif; else if plot_var eq 'lw' or plot_var eq 'sw' or plot_var eq 'rad' then begin

;  endif else if plot_var eq 'cld' then begin
;
;  endif

  ;TIME LABEL
    tim_str = string(t_sel[0],format='(f5.2)');
    if n_elements(t_sel) gt 1 then $
      tim_str+='-'+string((t_sel[1]-(npday-1d)/npday),format='(i2.2)')

;CREATE FIGURE

  if plot_var eq 'pres' then $
    figname=ifigdir+pres_plot+'_'+tim_str $
  else $
    figname=ifigdir+plot_var+'_'+tim_str

  yavg_call_plots,itimeheight,exp,datfils,fil_ind,t_ind,plot_var,hght,hghtw,$
    figname=figname,pres_plot=pres_plot

endfor ; ivar

;endfor ; iexp

print,'Done with all!!'

end
