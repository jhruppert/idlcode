; 
; Maps of diurnal composite Myanmar WRF output side-by-side with IMERG
;
; Now using MET_EM BC files for ERA5 winds.
;
; James Ruppert
; 8/14/20
; 
pro myanmar_wrf_imerg_dcomp_maps

config_dir,dirs=dirs

;EXPERIMENT SETTINGS
  expname='myanmar'
  cases=['ctl'];,'morr_mynn']

;idir='9km'
idir='2dom'
if idir eq '2dom' then begin
  domtag='d01'
;  domtag='d02'
endif

;Set this if uneven n-output files
  nfils_set=241
casedir=dirs.scdir+expname+'/WRF/experiment/'+idir+'/'
dirs.figdir+=expname+'/'
config_wrfexp, casedir=casedir,cases=cases,dirs=dirs,$
  dims=dims, vars=vars, nfils_set=nfils_set, domtag=domtag;, /verbose

;OB DIRECTORIES
maindir=dirs.scdir+'myanmar/'
imergfil=maindir+'imerg/data/jjas_2013-2017_diurncomp_3B-HHR.MS.MRG.3IMERG.V06B.nc4';'imerg/20150620-20150629.nc4'
npd_imerg=48
;era_fil=maindir+'era5/jjas_2013-2017/ERA5-JJAS_2013-2017_mean-pl.nc';'era5/ERA5-20150623-20150630-uv.nc'
;USING MET_EM FILES FOR DIURNAL COMPOSITE ERA5
era_dir=dirs.wkdir+'wrf_wps/wrfv4/output_myanmar/JJAS_2013-2017/era5/'
npd_era=24

;Local SOLAR time conversion
local=round(mean(dims.lon)/360.*24.) ; deg lon --> hours
print,'Adding +'+strtrim(local,2)+' for LT'

;----PLOT OPTIONS--------------------

nday_avg=5 ; will composite over last 'nday_avg' days for WRF

iwind=1 ; Plot 10m wind vectors?

allvars=['RAINNC']
nvsel=n_elements(allvars)

;LEVEL SELECTION
  psel=925;700;800;700 700 is the winner of 600, 700, 800, 850, 925
  izlev=(where(dims.pres eq psel))[0]

;----TIME SPECS--------------------

;FULL TIME SERIES
  time=dims.time
  nt_full=dims.nt
  npd=dims.npd
  nhrs=1.*nt_full*npd/24.
  nd=(1.*nhrs-(nhrs mod 24))/24.
  time_hrs=indgen(nhrs)

;----READ OBS--------------------

  ;IMERG RAINFALL
    imerg=read_nc_var(imergfil,'precipitationCal'); Keep in mm/hr ; * 24 ; mm/hr --> mm/d
    imerg=transpose(temporary(imerg))

;    timeim=read_nc_var(imergfil,'time')
;    timeim=julday(1,1,1970,0,0,temporary(timeim))
;    ;FIND OFFSET (IMERG TIME SERIES STARTS EARLIER)
;      tdiff=abs(timeim-time[0])
;      it_offset=(where(tdiff eq min(tdiff)))[0]

;    specs=size(rain,/dimensions)
;    ntim=specs[0] & nxim=specs[1] & nyim=specs[2]

    ;IMERG LAT/LON
    lonim=read_nc_var(imergfil,'lon')
    latim=read_nc_var(imergfil,'lat')

;    npdim=48

  ;ERA5 WINDS

  eralon=dims.lon
  eralat=dims.lat

  if iwind then begin

    ;PLEVS
    ilevera=925;[925, 950, 975, 1000]
    spawn,'ls '+era_dir+'met_em.'+domtag+'*',era_fil
    era_fil=era_fil[0:23]
    p_era=reform(read_nc_var(era_fil[0],'PRES',count=[1,1,38,1],offset=[0,0,0,0])) ; Pa
    izlev_era=(where(p_era*1d-2 eq psel))[0]

    ;READ DIURNAL COMPOSITE
    uera=fltarr(dims.nx,dims.ny,npd_era)
    vera=uera

    for it=0,npd_era-1 do begin
      fil=era_fil[it]
      count=[dims.nx+1,dims.ny,1,1] & offset=[0,0,izlev_era,0] ; x, y, p, t
      iu=reform(read_nc_var(fil,'UU',count=count,offset=offset))

    ;DESTAGGER
      ixin=indgen(dims.nx+1) & ixout=findgen(dims.nx)+0.5
      for iy=0,dims.ny-1 do uera[*,iy,it]=interpol(reform(iu[*,iy]),ixin,ixout)
      count=[dims.nx,dims.ny+1,1,1] & offset=[0,0,izlev_era,0] ; x, y, p, t
      iv=reform(read_nc_var(fil,'VV',count=count,offset=offset))
      iyin=indgen(dims.ny+1) & iyout=findgen(dims.ny)+0.5
      for ix=0,dims.nx-1 do vera[ix,*,it]=interpol(reform(iv[ix,*]),iyin,iyout)
    endfor

  endif

;----WRF OUTPUT--------------------

;for ic=0,dirs.nc-1 do begin
for ic=0,0 do begin

  print,'CASE: ',strupcase(dirs.cases[ic])

    i_nt=nt_full
    if ( strmatch(dirs.cases[ic],'*36h*') or strmatch(dirs.cases[ic],'icrf_*') or strmatch(dirs.cases[ic],'lwcr*') $
      or (dirs.cases[ic] eq 'lwswcrf') or (dirs.cases[ic] eq 'axisym') ) then i_nt-=36
    if strmatch(dirs.cases[ic],'*24h*') then i_nt-=24
    if strmatch(dirs.cases[ic],'*48h*') then i_nt-=48
    it_test=indgen(i_nt)+nt_full-i_nt

nd=i_nt/npd

var_str=allvars[0]
print,'VAR: ',var_str

;  if var_str eq 'rainrate' or var_str eq 'RAINNC' then begin
    ;ACCUMULATED RAIN OR RAIN RATE
    iv_str=var_str
    iv=where(vars.vars eq iv_str)
    file=dirs.files_post[ic,iv]

    count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,0,0] ; x,y,z,t
    rain=reform(read_nc_var(file,iv_str,count=count,offset=offset)) ; mm or mm/d

    rr=rain
    rr[*]=!values.f_nan
    for it=1,i_nt-2 do rr[*,*,it]=0.5*(rain[*,*,it+1]-rain[*,*,it-1]) ; mm/h
    rain=rr & rr=0

  ;WIND
    if iwind then begin
    ;U
      count=[dims.nx,dims.ny,1,i_nt] & offset=[0,0,izlev,0] ; x,y,z,t
;      iv=where(vars.vars eq 'U10')
      iv=where(vars.vars eq 'U')
      file=dirs.files_post[ic,iv]
;      u=reform(read_nc_var(file,'',varid='0',count=count,offset=offset))
      u=reform(read_nc_var(file,'U',count=count,offset=offset))
    ;V
;      iv=where(vars.vars eq 'V10')
      iv=where(vars.vars eq 'V')
      file=dirs.files_post[ic,iv]
;      v=reform(read_nc_var(file,'',varid='0',count=count,offset=offset))
      v=reform(read_nc_var(file,'V',count=count,offset=offset))
;      wind=create_struct('u',u,'v',v)
    endif

;----DIURNAL COMPOSITE--------------------

  raindc=fltarr(dims.nx,dims.ny,npd)
  udc=raindc
  vdc=raindc
  
  nd=i_nt/npd
  for it=0,npd-1 do begin
    t_ind=(indgen(nday_avg)+nd-nday_avg)*npd+it
    raindc[*,*,it]=mean(rain[*,*,t_ind],dimension=3,/double,/nan)
    udc[*,*,it]=mean(u[*,*,t_ind],dimension=3,/double,/nan)
    vdc[*,*,it]=mean(v[*,*,t_ind],dimension=3,/double,/nan)
  endfor

;----TIME LOOP--------------------

iu=0 & iv=0
iuera=0 & ivera=0

hr_int=1;2
for it=0,npd-1,hr_int do begin
;for it=0,0 do begin

  irain = reform(raindc[*,*,it]) ; mm/h

  ;WRF WIND
    if iwind then begin
      iu=reform(udc[*,*,it])
      iv=reform(vdc[*,*,it])
    endif

  ;IMERG RAINFALL
  it_im=it*npd_imerg/npd
  rainim=reform(imerg[it_im,*,*])

  ;ERA5 WIND
    if iwind then begin
      it_era=it*npd_era/npd
      iuera=reform(uera[*,*,it_era])
      ivera=reform(vera[*,*,it_era])
    endif

;----CREATE PLOTS--------------------

;  if var_str eq 'RAINNC' then begin
    setmin=0
    setmax=5;36;50;200
;  endif

  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=5
  figdir=dirs.figdir+dirs.cases[ic]+'/'+var_str+'/dcomp/'
  spawn,'mkdir '+figdir,tmp,tmpe
  figspecs=create_struct(figspecs,'figname',' ')
;  figspecs.title+=' ('+strupcase(dirs.cases[ic]);+')'
  titlesav=figspecs.title+' ('+strupcase(dirs.cases[ic]);+')'

  ;TIME INFO
    ltim=it+local
    ltim-=24*(ltim ge 24)
    t_stamp=string(ltim,format='(i2.2)')+' L'
;    t_stamp_fig=string(it,format='(i2.2)')+'h'
    t_stamp_fig='jjas2013-2017_'+string(it,format='(i2.2)')+'h'

    print,t_stamp

    figspecs.title=titlesav+'; '+t_stamp+')'

    figname=figdir+t_stamp_fig+'_'+var_str+'_'+dirs.cases[ic];+'_dayavg'
    if var_str eq 'avor' then figname+='_'+string(psel,format='(i3.3)')
    figspecs.figname=figname

    s_wrf=create_struct('lon',dims.lon,'lat',dims.lat,'var',irain,'u',iu,'v',iv)
    s_obs=create_struct('lon',lonim,'lat',latim,'var',rainim,'lonwind',eralon,'latwind',eralat,'u',iuera,'v',ivera)

    plot_myanmar_obcomp_maps, dirs, figname, s_wrf, s_obs, figspecs

endfor ; icase

;endfor ; ivar

endfor ; iday

print,'DONE!!'
end
