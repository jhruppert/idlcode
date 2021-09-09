; 
; Scatter diagrams from TC output
;
; James Ruppert
; 3/19/19
; 
pro run_tc_scatter

tc_sim_config, dirs=dirs, dims=dims, vars=vars;, /verbose
dirs.figdir+='edouard/'


;----PLOT OPTIONS--------------------


var1_str='pw'
var2_str='olr'

;TIME SELECTION
  hr_sel=[0,24*2]
  hr_tag=strtrim(hr_sel[0],2)+'-'+strtrim(hr_sel[1],2)+'hr'


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


;----READ VARS--------------------


print,'VAR1: ',var1_str
print,'VAR2: ',var2_str


;for ic=0,dirs.nc-1 do begin
;for ic=0,3 do begin
for ic=0,0 do begin

  print,'CASE: ',strupcase(dirs.cases[ic])

  ;VAR1
    iv=where(vars.vars eq var1_str)
    file=dirs.files_post[ic,iv]
    count=[dims.nx,dims.ny,1,nt_sel] & offset=[0,0,0,t_ind[0]] ; x,y,z,t
    var1=reform(read_nc_var(file,'',varid='0',count=count,offset=offset))
    var1=mean(temporary(var1),dimension=3,/nan,/double)

  ;VAR2
    iv=where(vars.vars eq var2_str)
    file=dirs.files_post[ic,iv]
    var2=reform(read_nc_var(file,'',varid='0',count=count,offset=offset))
    var1=mean(temporary(var2),dimension=3,/nan,/double)


;----CREATE PLOT--------------------


  figname=dirs.figdir+'scatter_'+dirs.cases[ic]+'_'+var1_str+'_'+var2_str

  ;PLOT SPECS
    csize=0.8
    position=[0.14,0.18,0.89,0.89]
    xsize=4.2 & ysize=2
    xtitle=var1_str;'Time [ hr ]'
    ytitle=var2_str;'Min. Pressure [ hPa ]'
    title='';strupcase(dirs.cases[ic])

  ;AXES
    x=var1
    y=var2
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

  plot,x,y,/nodata,position=position,$
    xstyle=9,ystyle=9,$
    xrange=xrange,yrange=yrange,yminor=2,$
    xtitle=xtitle,ytitle=ytitle1,$
    charsize=csize,$
    title=title

  oplot,x,y,psym=0,thick=1,color=0

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
      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.16,0.7]
  endif

  device,/close

  convert_png,figname,/remove_eps,res=200


endfor ; icase

print,'DONE!!'
end
