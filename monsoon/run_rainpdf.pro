; 
; Generate raifnall rate PDF.

; James Ruppert
; 8/14/20
; 
pro run_rainpdf

config_dir,dirs=dirs

;EXPERIMENT SETTINGS
  expname='myanmar'
  cases=['ctl','icloud0']

idir='9km'
;idir='2dom'
if idir eq '2dom' then domtag='d02'

;Set this if uneven n-output files
;  nfils_set=240-24
casedir=dirs.scdir+expname+'/WRF/experiment/'+idir+'/'
dirs.figdir+=expname+'/'
config_wrfexp, casedir=casedir,cases=cases,dirs=dirs,$
  dims=dims, vars=vars, nfils_set=nfils_set, domtag=domtag;, /verbose

;----TIME SPECS--------------------

;FULL TIME SERIES
  time=dims.time
  nt_full=dims.nt
  npd=dims.npd
  nhrs=1.*nt_full*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nhrs)

;----PDF SETTINGS--------------------

var_str='RAINNC'

;REGION SUBSET
  lonsel=[86.,101.5]
  latsel=[8.,25.]
  ixread=where((dims.lon ge lonsel[0]) and (dims.lon le lonsel[1]),nx_read)
  iyread=where((dims.lat ge latsel[0]) and (dims.lat le latsel[1]),ny_read)

;TIME SUBSET
  itsel=1
  if itsel then $
    t_ind=indgen(nt_full-npd)+npd $
  else $
    t_ind=indgen(nt_full)
  nt_read=n_elements(t_ind)

;HISTOGRAM SETTINGS
  histmin=2d
  histmax=60d
  binsize=2
  nbins=(histmax-histmin)/binsize + 1

;----READ VARS--------------------

rain_pdf=fltarr(dirs.nc,nbins)

for ic=0,dirs.nc-1 do begin
;for ic=0,0 do begin

  print,'CASE: ',strupcase(dirs.cases[ic])

    i_nt=nt_full
    if ( strmatch(dirs.cases[ic],'*36h*') or strmatch(dirs.cases[ic],'icrf_*') or strmatch(dirs.cases[ic],'lwcr*') $
      or (dirs.cases[ic] eq 'lwswcrf') or (dirs.cases[ic] eq 'axisym') ) then i_nt-=36
    if strmatch(dirs.cases[ic],'*24h*') then i_nt-=24
    if strmatch(dirs.cases[ic],'*48h*') then i_nt-=48
    it_test=indgen(i_nt)+nt_full-i_nt

    ;ACCUMULATED RAIN OR RAIN RATE
    iv_str=var_str
    iv=where(vars.vars eq iv_str)
    file=dirs.files_post[ic,iv]

    count=[nx_read,ny_read,1,nt_read] & offset=[ixread[0],iyread[0],0,t_ind[0]] ; x,y,z,t
    var=reform(read_nc_var(file,iv_str,count=count,offset=offset)) ; mm/d

    ;GET RAIN RATE USING CENTERED DIFFERENCE
      rain=fltarr(nx_read,ny_read,nt_read)
      rain[*]=!values.f_nan
      it0=0
      if t_ind[0] ne 0 then it0=1
      for it=it0,nt_read-2 do begin
        ind0=it-1
        ind1=it+1
        rain[*,*,it] = 0.5*(var[*,*,ind1] - var[*,*,ind0]) ; mm/time step (=mm/hr)
      endfor

  rain=reform(temporary(rain),1l*nx_read*ny_read*nt_read)

  ihist = histogram( var, nbins=nbins, binsize=binsize, min=histmin ,/nan ,/l64)
;  ihist /= total(ihist,/nan,/double)
  ihist=double(ihist)

  ihist=smooth(ihist,3,/edge_truncate,/nan)

  nmin=10 ; Minimum N to consider
  loc_setnan=where(ihist lt nmin)
  ihist[loc_setnan[0]:nbins-1]=!values.d_nan

  rain_pdf[ic,*] = ihist

endfor ; icase

;----CREATE PLOTS--------------------

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figdir=dirs.figdir
  figspecs=create_struct(figspecs,'figname',' ')
  figspecs.title='Rainfall Rate PDF'
  figname=figdir+'rain_pdf'
  figspecs.figname=figname

  ;PLOT SPECS
    csize=0.5;6
    position=[0.14,0.22,0.86,0.86]
    xsize=2.0 & ysize=1.0
    ;if do_sea then 
    xtitle='!8R!X [ mm hr!U-1!N ]'; else xtitle=''
    ytitle='N'
    ytitle2='Factors';'(!8p!Di!N!X-1)/!8p!X!DCTL!N'

  ;AXES
    x=indgen(nbins)*binsize+histmin
    y=indgen(10)
    yrange=[1e4,1e6]
;    xrange=[min(x)+0.4,max(x)-0.4]
    if ~keyword_set(xrange) then $
      xrange=[min(x),max(x)]
xrange[0]=0
    if ~keyword_set(yrange) then $
      yrange=[max(y),min(y)]

  set_plot,'ps'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=xsize,ysize=ysize,/inches,$
    /helvetica

  iytitle=ytitle

  ixtitle=xtitle
  ititle='';'Rainfall Rate PDF'

;  xtickv=indgen(5)*6
;  xtickname=string(xtickv,format='(i2.2)')
;  xticks=n_elements(xtickname)-1

  loadct,0,/silent

  plot,x,y,/nodata,position=position,/ylog,$
    xstyle=9,ystyle=9,$;xticks=xticks,xtickv=xtickv,xtickname=xtickname,$
;    yticklen=-0.016,
    xticklen=0.03,$
    xticks=3,xtickv=[0,20,40,60],xminor=4,$
    xrange=xrange,yrange=yrange,$;[0,10000],$
    xtitle=ixtitle,ytitle=iytitle,$
    charsize=csize,title=ititle

  lthick=replicate(3.,dirs.nc)
  lstyle=indgen(dirs.nc)
  lcolor=findgen(4)/3*250
  lcolor=reverse(lcolor)

  for ic=0,dirs.nc-1 do $
;    oplot,x,reform(rain_pdf[ic,*]),linestyle=lstyle[ic],thick=lthick[ic]*0.55,color=lcolor[ic]
    oplot,x,reform(rain_pdf[ic,*]),linestyle=([0,1])[ic],thick=1.2,color=0

  ;LEGEND
  ileg=0
  if ileg then begin
  loadct,3,/silent
;  if ileg then begin
    csize_fac=1.;0.7
    margin=0.1
    pspacing=1.8 ; length of lines
    spacing=0.5 ; between lines
    leg_str=strupcase(dirs.cases)
    leg_style=lstyle
    leg_thick=lthick
    leg_color=lcolor
    legend2,leg_str,linestyle=leg_style,thick=leg_thick,COLORS=leg_color,$
      charsize=csize*csize_fac,/top_legend,/left_legend,/normal,$
      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.16,0.85]
;      margin=margin,pspacing=pspacing,spacing=spacing,box=0,position=[0.56,0.35]
;  loadct,0,/silent
  endif

  device,/close

  convert_png,figname,res=400;,/remove_eps

print,'DONE!!'
end
