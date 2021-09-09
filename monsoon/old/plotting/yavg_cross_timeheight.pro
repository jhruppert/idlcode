; 
; Plot PW cross sections from ICON 3D output.
; 
; James Ruppert
; 20.10.16
; 
pro yavg_cross_timeheight,figname,title,hght,pcps,$
  levels,colors,col_table,$
  icbar,cbar_tag,icbar_format,$
  icont,cstyle,$
  clevels1,clevels2,$
  v1s=v1s,v2s=v2s,psi=psi

plot_stmf=1 ; plot streamfunction?
sub_day_mean=0

if sub_day_mean then $
  print,'Subtracting daily mean!'

;CREATE FIGURE

  ;PLOT SPECS
    csize=0.7
    xticklen=-0.025 & yticklen=-0.015
    xsize=3.0 & ysize=2.5
    position=[0.14,0.14,0.89,0.92]

  ;AXES
    x=indgen(49) ; km
    xrange=[0,49]
    y=hght*1e-3 ; m --> km
    yrange=[0,15]
    xtitle='Time [ h ]'
    ytitle='Height [ km ]'
    p_ytitle='P [ mm h!U-1!N ]'
    p_yrange=[0,5]
    p_yticks=2

  ;MULTI PLOT
    position1 = [ position[0] , 0.37 , position[2] , position[3] ]
    gap=0.04
    position2 = [ position[0] , position[1] , position[2] , position1[1]-gap ]


nv1=v1s
nv2=v2s
npsi=psi

specs=size(v1s)
if specs[0] eq 2 then ntime=specs[2]
npday=24/ntime
dtimes=string(indgen(ntime)*npday,format='(i2.2)')
dtimes=[dtimes,dtimes,dtimes[0]]

if sub_day_mean then begin
  stop
  day_mean=mean(v1,dimension=3,/nan,/double)
  for it=0,ntime-1 do v1[*,*,it] -= day_mean
endif


;if it eq -1 then begin
;  v1=mean(nv1,dimension=3,/double,/nan)
;  v2=mean(nv2,dimension=3,/double)
;  ipsi=mean(npsi,dimension=3,/double)
;  pcp=mean(pcps,dimension=2,/double)
;  ftim='mean'
;  title_str=title+' (time mean)'
;;  title_str=title+'*'
;endif else begin
  v1=[[nv1],[nv1],[nv1[*,0]]]
  v2=[[nv2],[nv2],[nv2[*,0]]]
  ipsi=[[npsi],[npsi],[npsi[*,0]]]
  pcp=[pcps,pcps,pcps[0]]
  title_str=title


fname=figname
cbar_format=icbar_format


  set_plot,'ps'
  epsname=fname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  !p.multi=[0,0,2]

  loadct,0,/silent

;MAIN SECTION

  plot,x,y,/nodata,position=position1,$
    xstyle=9,ystyle=9,$
    ylog=ylog,ytickv=ytickv,yticks=yticks,$
    xticks=8,$;xtickv=[45,30,20,15,10,5],$
    xtick_get=xtickv,ytick_get=ytickv,$
    xminor=2,yminor=2,$
    xticklen=xticklen,yticklen=yticklen,$
    xtickname=replicate(' ',9),$
    xrange=xrange,$
    yrange=yrange,$;yticks=4,yminor=6,ytickname=strtrim([0,12,0,12,24],2),
    ytitle=ytitle,$
    charsize=csize,title=title_str,color=0

  loadct,col_table,/silent

  ;FILL SHADING
    for i=0,1 do $
      contour,transpose(v1),x,y,/cell_fill,/overplot,$
        levels=levels,c_colors=colors

  loadct,0,/silent

  ;OVERPLOT AXIS LINES
    for i=0,1 do begin
      plots,!x.crange,!y.crange[[i,i]]
      plots,!x.crange[[i,i]],!y.crange
    endfor

  cthick=1.7

  csize2 = csize*0.85

  ;STREAMFUNCTION
    if plot_stmf then begin
      cint=4
      psi_levs=findgen(200)*cint+cint
      contour,transpose(ipsi),x,y,/follow,/overplot,$
        levels=psi_levs,c_linestyle=0,c_charsize=csize2,c_thick=cthick;,c_labels=0
      contour,transpose(ipsi),x,y,/follow,/overplot,$
        levels=-1.*reverse(psi_levs),c_charsize=csize2,c_linestyle=2,c_thick=cthick;,c_labels=0
    endif

  ;CONTOURS
    if icont[0] then begin
;      contour,v1,x,y,/follow,/overplot,$
;        levels=clevels1,c_charsize=csize2,c_linestyle=cstyle[0],c_thick=cthick,c_labels=replicate(1,20)
      ;cint=2;0.05
      ;clevs=findgen(40)*cint+cint
;      loadct,11,/silent
      cthick=0.5
      contour,transpose(v1),x,y,/follow,/overplot,$
        levels=clevels1,c_charsize=csize2,c_linestyle=0,c_thick=cthick,$;c_labels=0,$
        c_colors=0
      contour,transpose(v1),x,y,/follow,/overplot,$
        levels=-1*reverse(clevels1),c_charsize=csize2,c_linestyle=1,c_thick=cthick,$;c_labels=0,$
        c_colors=0
    endif
  ;LIQUID AND ICE CLOUD or CLC
    if icont[1] then begin
;      clevels2=[0.1,0.5,5]
      loadct,11,/silent
      contour,transpose(v2),x,y,/follow,/overplot,$
        levels=clevels2,c_charsize=csize2,c_linestyle=cstyle[1],c_thick=cthick,$
        c_labels=replicate(1,20),$
        c_colors=210
;        c_colors=0
    endif
;    if icont[2] then begin
;;      contour,v3,x,y,/follow,/overplot,$
;;        levels=clevels3,c_charsize=csize2,c_linestyle=cstyle[2],c_thick=cthick,c_labels=0,$
;;        c_colors=210
;      ;cthick=0.8
;      cint=0.5
;      clevs=findgen(10)*cint+cint
;      contour,transpose(v2),x,y,/follow,/overplot,$
;        levels=clevs,c_charsize=csize2,c_linestyle=0,c_thick=cthick,$;c_labels=0,$
;        c_colors=0
;      clevs=-1*reverse(clevs)
;      contour,transpose(v2),x,y,/follow,/overplot,$
;        levels=clevs,c_charsize=csize2,c_linestyle=1,c_thick=cthick,$;c_labels=0,$
;        c_colors=0
;    endif

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


;RAINFALL

  plot,x,y,/nodata,position=position2,$
;    /ylog,$
    yticks=p_yticks,$
    xstyle=9,ystyle=9,$
    xticks=8,xtickname=dtimes[indgen(9)*6],$;xtickv=[45,30,20,15,10,5],$
    xminor=2,$;yminor=2,yticks=2,$
    xticklen=xticklen*3.,yticklen=yticklen,$
    xrange=xrange,yrange=p_yrange,$;yticks=4,yminor=6,ytickname=strtrim([0,12,0,12,24],2),
    xtitle=xtitle,ytitle=p_ytitle,$
    charsize=csize

  oplot,x,pcp,linestyle=0,thick=2.5

  ;OVERPLOT AXIS LINES
  for i=0,1 do begin
    plots,xrange,p_yrange[[i,i]]
    plots,xrange[[i,i]],p_yrange
  endfor


  device,/close

  !p.multi=0

  convert_png,fname,res=200,/remove_eps


print,'Done plotting!!'

end
