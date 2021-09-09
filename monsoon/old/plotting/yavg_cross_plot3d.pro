; 
; Plot cross sections from y-averaged ICON 3D output.
; 
; James Ruppert
; 16.12.17
; 
pro yavg_cross_plot3d,figname,title,radius,hght,pcps,$
  levels,colors,col_table,$
  icbar,cbar_tag,icbar_format,$
  icont,cstyle,$
  clevels1,clevels2,$
  shs=shs,lhs=lhs,t_gs=t_gs,$
  v1s=v1s,v2s=v2s,psi=psi

t_int=1 ; skip this many (1=do all)

plot_stmf=1 ; plot streamfunction?
sub_day_mean=0
sub_dom_mean=0
sub_dom_mean_v2=0
ismooth=1
ismooth_v2=0

if sub_day_mean then $
  print,'Subtracting daily mean!'
if sub_dom_mean then $
  print,'Subtracting domain mean!'

;CREATE FIGURE

  ;PLOT SPECS
    ;csize=0.7
    csize=1.1
    xticklen=-0.025 & yticklen=-0.015
    xsize=3.0 & ysize=2.5
    ysize=3.2
    position=[0.14,0.09,0.89,0.94]

  ;AXES
    x=radius*1e-3 ; km
    ;xrange=[min(x),max(x)]
;    xrange=[0,max(x)];+0.5*dr]
    xrange=[min(x),max(x)];500]
    nx = n_elements(x)
    y=hght*1e-3 ; m --> km
    yrange=[0,17]
    xtitle='!8x!X [ km ]'
    ytitle='Height [ km ]'
    p_ytitle='P!D0!N [ mm h!U-1!N ]'
    p_yrange=[0,4]
    p_yticks=2

  ;MULTI PLOT
    yb1=0.5
    position1 = [ position[0] , yb1 , position[2] , position[3] ]
    gap=0.025
    total=yb1-position[1]
    position2 = [ position[0] , position[1]+2*0.33*total , position[2] , position1[1]-gap ]
    position3 = [ position[0] , position[1]+0.33*total , position[2] , position2[1]-gap ]
    position4 = [ position[0] , position[1] , position[2] , position3[1]-gap ]


nv1=v1s
nv2=v2s
npsi=psi

specs=size(v1s)
if specs[0] eq 2 then ntime=1 else ntime=specs[3]
npday=24/ntime
dtimes=string(indgen(ntime)*npday,format='(i2.2)')

if sub_day_mean then begin
  day_mean=mean(v1,dimension=3,/nan,/double)
  for it=0,ntime-1 do v1[*,*,it] -= day_mean
endif


;LOOP OVER TIME

for it=-1,ntime-1,t_int do begin
;for it=0,ntime-1,6 do begin
;for it=12,12 do begin

print,'time = ',strtrim(it,2)

if it eq -1 then begin
  v1=mean(nv1,dimension=3,/double,/nan)
  v2=mean(nv2,dimension=3,/double)
  ipsi=mean(npsi,dimension=3,/double)
  pcp=mean(pcps,dimension=2,/double)
  sh=mean(shs,dimension=2,/double)
  lh=mean(lhs,dimension=2,/double)
  t_g=mean(t_gs,dimension=2,/double)
  ftim='mean'
  title_str=title+' (time mean)'
;  title_str=title+'*'
endif else begin
  v1=reform(nv1[*,*,it])
  v2=reform(nv2[*,*,it])
  ipsi=reform(npsi[*,*,it])
  pcp=reform(pcps[*,it])
  sh=reform(shs[*,it])
  lh=reform(lhs[*,it])
  t_g=reform(t_gs[*,it])
  ftim=string(it,format='(i2.2)')
  title_str=title+' (T='+dtimes[it]+'L)'
endelse

;pcp *= 24. ; mm/h --> mm/d
lh*=-1.
sh*=-1.
stop
meanv1=mean(v1,dimension=1,/nan,/double)
meanv2=mean(v2,dimension=1,/nan,/double)
if sub_dom_mean then $
  for ix=0,nx-1 do v1[ix,*] -= meanv1
if sub_dom_mean_v2 then $
  for ix=0,nx-1 do v2[ix,*] -= meanv2

;SMOOTHING
if ismooth then begin
  xsmooth=5;9;3
  v1=smooth(v1,[xsmooth,0],/edge_truncate,/nan)
;  for i=0,2 do $
;    v1=gauss_smooth(v1,[5,0],/edge_truncate,/nan)
;  v2=smooth(v2,[xsmooth,0],/edge_truncate,/nan)
  ipsi=smooth(ipsi,[xsmooth,0],/edge_truncate,/nan)
  pcp=smooth(pcp,xsmooth,/edge_truncate,/nan)
;  ih_trop=smooth(ih_trop,xsmooth,/edge_truncate,/nan)
endif
if ismooth_v2 then begin
  xsmooth=5;9;3
  v2=smooth(v2,[xsmooth,0],/edge_truncate,/nan)
endif

fname=figname+'_'+ftim
cbar_format=icbar_format


  set_plot,'ps'
  epsname=fname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  !p.multi=[0,0,4]

  loadct,0,/silent

;MAIN SECTION

  plot,x,y,/nodata,position=position1,$
    xstyle=9,ystyle=9,$
    ylog=ylog,ytickv=ytickv,yticks=yticks,$
;    xticks=5,xtickv=[45,30,20,15,10,5],$
    xtick_get=xtickv,ytick_get=ytickv,$
    xminor=2,yminor=2,$
    xticklen=xticklen,yticklen=yticklen,$
    xtickname=replicate(' ',10),$
    xrange=xrange,$
    yrange=yrange,$;yticks=4,yminor=6,ytickname=strtrim([0,12,0,12,24],2),
    ytitle=ytitle,$
    charsize=csize,title=title_str,color=0

  loadct,col_table,/silent

  ;FILL SHADING
    for i=0,1 do $
      contour,v1,x,y,/cell_fill,/overplot,$
        levels=levels,c_colors=colors

  loadct,0,/silent

  ;OVERPLOT AXIS LINES
    for i=0,1 do begin
      plots,!x.crange,!y.crange[[i,i]]
      plots,!x.crange[[i,i]],!y.crange
    endfor

  cthick=1.7

  ;STREAMFUNCTION
    if plot_stmf then begin
      cint=2
      psi_levs=findgen(200)*cint+cint
      contour,ipsi,x,y,/follow,/overplot,$
        levels=psi_levs,c_linestyle=0,c_thick=cthick,c_labels=0
      contour,ipsi,x,y,/follow,/overplot,$
        levels=-1.*reverse(psi_levs),c_linestyle=2,c_thick=cthick,c_labels=0
    endif

  ;CONTOURS
    csize2 = csize*0.85
    if icont[0] then begin
;      contour,v1,x,y,/follow,/overplot,$
;        levels=clevels1,c_charsize=csize2,c_linestyle=cstyle[0],c_thick=cthick,c_labels=replicate(1,20)
      ;cint=2;0.05
      ;clevs=findgen(40)*cint+cint
;      loadct,11,/silent
      cthick=0.5
      contour,v1,x,y,/follow,/overplot,$
        levels=clevels1,c_charsize=csize2,c_linestyle=0,c_thick=cthick,$;c_labels=0,$
        c_colors=0
      contour,v1,x,y,/follow,/overplot,$
        levels=-1*reverse(clevels1),c_charsize=csize2,c_linestyle=1,c_thick=cthick,$;c_labels=0,$
        c_colors=0
    endif
  ;LIQUID AND ICE CLOUD or CLC
    if icont[1] then begin
;      clevels2=[0.1,0.5,5]
      loadct,11,/silent
      contour,v2,x,y,/follow,/overplot,$
        levels=clevels2,c_charsize=csize2,c_linestyle=cstyle[1],c_thick=cthick,$
        c_labels=replicate(1,10),$
        c_colors=210
;        c_colors=0
    endif

  ;COLOR BAR
  if icbar then begin
    loadct,col_table,/silent
    cpos= [ 0.906 ,$
            position1[1]+0.1 ,$
            0.92 ,$
            position1[3]-0.1 ]
    colorbar2, colors=colors, range=[min(levels),max(levels)],$; divisions=ndivs,$
      charsize=csize2*0.7, position=cpos, /right, /vertical, title=cbar_tag, $
      annotatecolor='black',format=cbar_format,DIVISIONS=6
    loadct,0,/silent
  endif


lthick=2.0

;RAINFALL

  plot,x,y,/nodata,position=position2,$
;    /ylog,$
    yticks=p_yticks,$
    xstyle=9,ystyle=9,$
    xtickname=replicate(' ',10),$
;    xticks=5,xtickv=[45,30,20,15,10,5],$
    xminor=2,$;yminor=2,yticks=2,$
    xticklen=xticklen*3.,yticklen=yticklen,$
    xrange=xrange,yrange=p_yrange,$;yticks=4,yminor=6,ytickname=strtrim([0,12,0,12,24],2),
    ytitle=p_ytitle,$
    charsize=csize

  oplot,x,pcp,linestyle=0,thick=lthick

;FLUXES

;    axis,yaxis=1,ystyle=1,charsize=csize,yticklen=yticklen,$
;      yrange=[-10,100],ytitle='Fluxes [ W m!U-2!N ]',/save,/data
  plot,x,y,/nodata,position=position3,$
    yticks=p_yticks,xstyle=9,ystyle=9,xminor=2,$
    xtickname=replicate(' ',10),$
    xticklen=xticklen*3.,yticklen=yticklen,$
    xrange=xrange,yrange=[-50,250],$
    ytitle='F!D0!N [ W m!U-2!N ]',$
    charsize=csize

    oplot,x,sh,linestyle=2,thick=lthick
    oplot,x,lh,linestyle=0,thick=lthick

;SURFACE TEMP

  plot,x,y,/nodata,position=position4,$
    yticks=p_yticks,xstyle=9,ystyle=9,xminor=2,$
    xticklen=xticklen*3.,yticklen=yticklen,$
    xrange=xrange,yrange=[290,340],$
    xtitle=xtitle,ytitle='T!D0!N [ K ]',$
    charsize=csize

  oplot,x,t_g,linestyle=0,thick=lthick

  ;OVERPLOT AXIS LINES
;  for i=0,1 do begin
;    plots,xrange,p_yrange[[i,i]]
;    plots,xrange[[i,i]],p_yrange
;  endfor


  device,/close

  !p.multi=0

  convert_png,fname,res=200,/remove_eps

endfor

print,'Done plotting!!'

end
