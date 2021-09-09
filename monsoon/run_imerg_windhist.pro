; 
; Create composite plots as a function of filtered coast-normal wind speed.
;
; James Ruppert
; 7/22/21
; 
pro run_imerg_windhist

config_dir,dirs=dirs

;----PLOT OPTIONS--------------------

coast_thresh=200 ; km

plot_type='binned'
;plot_type='regression'

var_plot='rain';
;var_plot='pw';'rain';

icross=1 ; Which cross section for unorm calculation?

;ERAi SETTINGS
;LEVEL SELECTION
  psel_era=850;500;700;925

;BOUNDS FOR READ-IN
;  bounds=[78,6,104,29]
  bounds=[70,2,108,32]

;SELECT DATE RANGE
;  yy_plot=[2013,2017]
  yy_plot=[2000,2020]
yy_plot[0]=2001
  mm_plot=[1,12]
;  mm_plot=[6,12]
  dd_plot=[1,31] ; inclusive
;dat_str='2013-2017'
dat_str='2000-2020'

  ;DATE STRING
;    form2='(i2.2)'
;    form4='(i4)'
;    dat_str=string(mm_plot[0],format=form2)+string(dd_plot[0],format=form2)+strmid(strtrim(yy_plot[0],2),2,2)+'-'+$
;            string(mm_plot[1],format=form2)+string(dd_plot[1],format=form2)+strmid(strtrim(yy_plot[1],2),2,2)

;----DIRECTORIES--------------------

  figdir=dirs.figdir+'myanmar/imerg/wind_histogram/'
  ifigdir=figdir+'ind_unorm/'+dat_str+'/'

  maindir=dirs.wkdir
  im_fil=dirs.wkdir+'imerg/imerg/data/imerg_3B-HHR.MS.MRG.3IMERG.V06_daily_2000-2021_latlonsubset.nc4'
;  npd_imerg=48

  era_dir=maindir+'era5/'
  era_fil=era_dir+'ERA5-20000101-20201231-pl_dayavg.nc'
  era_pw=era_dir+'ERA5-20000101-20201231-pw_dayavg.nc'
;  npd_era=24

;Local SOLAR time conversion
;local=6;round(mean(dims.lon)/360.*24.) ; deg lon --> hours
;print,'Adding +'+strtrim(local,2)+' for LT'

;----TIME ARRAY FOR DAILY AVERAGE TIME SERIES--------------------

  ;SELECTED TIME ARRAY
    time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
      final=julday(mm_plot[1],dd_plot[1],yy_plot[1],23,59,59),step_size=1,units='Days');30,units='Minutes')
    nt=n_elements(time)
    nd=nt;/npd_imerg
    nyr=yy_plot[1]-yy_plot[0]+1

  ;SAVE JJAS INDICES
    caldat,time,mm,dd,yy
    jjas=where((mm ge 6) and (mm le 9))

;=====BEGIN READING=========================================================

;----READ RAIN--------------------

  if var_plot eq 'rain' then begin
    var=read_nc_imerg(time,im_fil,lon=lon,lat=lat,bounds=bounds) ; already in mm/d
    nx=n_elements(lon)
    ny=n_elements(lat)
  endif

;----READ ERA--------------------

  u=read_nc_era5(time,era_fil,'var131',plev=psel_era,lon=eralon,lat=eralat,bounds=bounds)
  v=read_nc_era5(time,era_fil,'var132',plev=psel_era,lon=eralon,lat=eralat,bounds=bounds)
  nxera=n_elements(eralon)
  nyera=n_elements(eralat)

  ;PW
  if var_plot eq 'pw' then begin;icalc_pw=1 else icalc_pw=0
    var=read_nc_era5(time,era_pw,'var137',lon=lon,lat=lat,bounds=bounds)
    nx=nxera
    ny=nyera
  endif

;=====CREATE INDEX FOR REGRESSION/BINNING=========================================================

  coastnormal, icross, coast_thresh, u, v, eralon, eralat, uindex=uindex

;----FILTER--------------------

  ;UNORM INDEX
  filter_monsoon, uindex, 'fft', var_bw=u_bw, var_intra=u_intra
;  filter_monsoon, uindex, 'rmean', var_bw=u_bw, var_intra=u_intra

  ;ALL VARS
;  filter_monsoon, var, 'rmean', var_bw=var_bw, var_intra=var_intra
;  filter_monsoon, u, 'rmean', var_bw=ua_bw, var_intra=ua_intra
;  filter_monsoon, v, 'rmean', var_bw=va_bw, var_intra=va_intra
;  filter_monsoon, unorm, 'rmean', var_bw=unorm_bw, var_intra=unorm_intra

;----PLOT WIND INDEX--------------------

iplot_ind=0
if iplot_ind then begin

  var_str='pw'
  setmax=20 & setmin='0.'
  myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax
  figname=figdir+'test_wind_index_cross'+strtrim(icross,2)
  figspecs=create_struct(figspecs,'figname',figname)

  varm=max(abs(unorm2[*,*,jjas]),dimension=3,/nan)
  varm[where(finite(varm))]=20
  wrf_myanmar_map_plot, dirs, varm, eralon[ix], eralat[iy], figspecs
exit
endif


;=====GENERATE COMPOSITES=========================================================

;----KEEP JJAS ONLY--------------------

  nd=n_elements(jjas)

  uindex=uindex[jjas]
  u_bw=u_bw[jjas]
  u_intra=u_intra[jjas]

;  unorm=unorm[*,*,jjas]

  var=var[*,*,jjas]
  u=u[*,*,jjas]
  v=v[*,*,jjas]

;----PLOT MEAN--------------------

iplot_mn=0
if iplot_mn then begin

  var_mn=mean(var,dimension=3,/nan,/double)
  umn=mean(u,dimension=3,/nan,/double)
  vmn=mean(v,dimension=3,/nan,/double)

;  cvarm=mean(unorm,dimension=3,/nan,/double)
  cvar=create_struct('cvar',cvarm,'x',eralon,'y',eralat)

  figname=ifigdir+'test_cross'+strtrim(icross,2)

  ;OVERLAY CROSS SECTION
  plt_cross=[xcross,ycross]

  plot_mean_map_myanmar, dirs, var_mn, lon, lat, var_plot, figname, $
    u=umn, v=vmn, eralon=eralon, eralat=eralat, cvar=cvar, cross=plt_cross
exit
endif

;----STANDARDIZE--------------------

  ;UNORM INDEX

    std_bw=stddev(u_bw,/nan,/double)
    std_intra=stddev(u_intra,/nan,/double)

    print,'Standard Deviations:'
    print,'BW:',std_bw
    print,'Intra:',std_intra
  
    u_bw    = (u_bw -    mean(u_bw,/nan,/double))    / std_bw
    u_intra = (u_intra - mean(u_intra,/nan,/double)) / std_intra
;set_plot,'x' 
;caldat,time[jjas],mm,dd,yy
;iy=where(yy ge 2013 and yy le 2015)
;plot,u_bw[iy]
;plot,u_intra[iy]
;stop
  
    ;N-CASES where |irain| > 1 sigma
    ibw = where(abs(u_bw) ge 1,nbw)
    iintra = where(abs(u_intra) ge 1,nintra)

  ;PRINT STATS
  
    print,'N where |index| >= 1:'
    print,'BW:',nbw
    print,'Intra:',nintra

  ;MAIN VARIABLES

i_std=1
if i_std then begin

    std_var = stddev(var,dimension=3,/nan,/double)
    std_var = rebin(temporary(std_var),[nx,ny,nd])

    std_u = stddev(u,dimension=3,/nan,/double)
    std_u = rebin(temporary(std_u),[nxera,nyera,nd])
    std_v = stddev(v,dimension=3,/nan,/double)
    std_v = rebin(temporary(std_v),[nxera,nyera,nd])
;    std_unorm = stddev(unorm,dimension=3,/nan,/double)
;    std_unorm = rebin(temporary(std_unorm),[nxera,nyera,nd])

    varm = mean(var,dimension=3,/nan,/double)
    varm = rebin(varm,nx,ny,nd)
    var = (var - varm)/std_var

    um = mean(u,dimension=3,/nan,/double)
    um = rebin(um,nxera,nyera,nd)
    u = (u - um)/std_u
    vm = mean(v,dimension=3,/nan,/double)
    vm = rebin(vm,nxera,nyera,nd)
    v = (v - vm)/std_v
;    unm = mean(unorm,dimension=3,/nan,/double)
;    unm = rebin(unm,nxera,nyera,nd)
;    unorm = (unorm - unm)/std_unorm

endif

;----CALCULATE COMPOSITES AND PLOT--------------------

if i_std then begin
  var_str='pw'
  setmax=0.5
;setmax=1.
  setmin=-1.*setmax
  cbform='(f4.1)'
endif else begin
  if var_plot eq 'rain' then begin
    var_str='RAINNC'
    setmax=40;24
    setmin='0.'
  endif else begin
    var_str=var_plot
    setmax=70
    setmin=30
  endelse
  cbform='(i2)'
endelse

myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=0.2
figspecs=create_struct(figspecs,'figname',' ')
figspecs.cbar_format=cbform
;figspecs.ndivs-=1
figspecs.title=''

if i_std and var_plot eq 'rain' then begin
  figspecs.cbar_tag=' '
  figspecs.col_table=71
  figspecs.colors=reverse(figspecs.colors)
endif

for iband=1,2 do begin ; Biweekly / Intraseasonal
;for iband=1,1 do begin ; Biweekly / Intraseasonal

  print,'iband: ',iband

  if iband eq 1 then begin
    bandtag='bw'
    uband=u_bw
    stdd=std_bw
  endif else if iband eq 2 then begin
    bandtag='intra'
    uband=u_intra
    stdd=std_intra
  endif

  if plot_type eq 'binned' then begin

    ;BINNING APPROACH

;      var_str='pw'
;;      setmax=0.5
;      setmax=1.
;      setmin=-1.*setmax
;      cbform='(f4.1)'
;
;      myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=0.2
;      figspecs=create_struct(figspecs,'figname',' ')
;      figspecs.cbar_format=cbform
;      ;figspecs.ndivs-=1
;      figspecs.title=''
;
;      figspecs.cbar_tag=' '
;      if var_plot eq 'rain' then begin
;        figspecs.col_table=71
;        figspecs.colors=reverse(figspecs.colors)
;      endif

    ;FOR NON-STANDARDIZED
;      var_str='RAINNC'
;      setmax=25 & setmin='0.'
;setmax=40
;      cbform='(i2)'
;      myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=2
;      figspecs=create_struct(figspecs,'figname',' ')
;      figspecs.cbar_format=cbform
;      figspecs.ndivs-=1
;      figspecs.title=''

    ;COMPOSITING THRESHOLDS
      bins = ['-1','-0.5','0.5','1']
bins = ['-1','-0.75','-0.5','-0.25','0','0.25','0.5','0.75','1']
bins = ['-1','-0.5','-0.25','0.25','0.5','1']
      nbin = n_elements(bins)

;    for ibin=0,nbin-2 do begin
    for ibin=0,nbin-1 do begin
  ;  for ibin=nbin-2,nbin-2 do begin

      ;AVERAGE VARIABLES OVER INDEX BINS
  
      print,'  ibin: ',ibin
  
      if ibin eq nbin-1 then begin
      ;COMPOSITE OVER ALL DATES
        it_sel=indgen(nd)
        figspecs.title='All Dates'
        figname=ifigdir+var_plot+'_'+bandtag+'_all_cross'+strtrim(icross,2)
        figspecs.figname=figname
      endif else begin
        bin_txt = [ bins[ibin] , bins[ibin+1] ]
        bin_p = float(bin_txt) * stdd
  ;      bin_n = -1*bin_p

        it_sel = where((uband ge bin_p[0]) and (uband le bin_p[1]),np_p)
  ;      it_n = where((uband le bin_n[0]) and (uband ge bin_n[1]),np_n)
        print,'Count-p:',np_p
  ;      print,'Count-n:',np_n 

        bintag=strtrim(ibin+1,2)

        figspecs.title=bin_txt[0]+' to '+bin_txt[1]+' sigma'
        figname=ifigdir+var_plot+'_'+bandtag+'_'+bintag+'_cross'+strtrim(icross,2)
        figspecs.figname=figname
      endelse

    ;MAP IT
    
;      for isign=0,1 do begin
;    
;        print,'    isign: ',isign
;  
;        if isign eq 0 then begin
;          it_sel=it_p 
;          s_tag='p'
;        endif else begin
;          it_sel=it_n
;          s_tag='n'
;        endelse
  
        var_mn=mean(var[*,*,it_sel],dimension=3,/nan,/double)
        umn=mean(u[*,*,it_sel],dimension=3,/nan,/double)
        vmn=mean(v[*,*,it_sel],dimension=3,/nan,/double)
;  scale=10
;  umn*=scale
;  vmn*=scale
;        unorm_mn=mean(unorm[*,*,it_sel],dimension=3,/nan,/double)
   
        wind=create_struct('u',umn,'v',vmn,'x',eralon,'y',eralat)
;        cvar=create_struct('cvar',unorm_mn,'x',eralon,'y',eralat)
   
        wrf_myanmar_map_plot, dirs, var_mn, lon, lat, figspecs, wind=wind, cvar=cvar, /noscalewind
    
;      endfor ; isign
  
    endfor ; ibin

  endif else if plot_type eq 'regression' then begin

    ;REGRESSION APPROACH

      var_str='pw'
      setmax=0.5
      ;setmax=1.
      setmin=-1.*setmax
      cbform='(f4.1)'
    
      myan_figspecs, var_str, figspecs, setmin=setmin, setmax=setmax, set_cint=0.2
      figspecs=create_struct(figspecs,'figname',' ')
      figspecs.cbar_format=cbform
      ;figspecs.ndivs-=1
      figspecs.title=''
      
      figspecs.cbar_tag=' '
      if var_plot eq 'rain' then begin
        figspecs.col_table=71
        figspecs.colors=reverse(figspecs.colors)
      endif

      dims=[nx,ny,nd]

      var_reg = regress_3d(uband,var,dims)

      dims=[nxera,nyera,nd]

      u_reg = regress_3d(uband,u,dims)
      v_reg = regress_3d(uband,v,dims)
;scale=10
;u_reg*=scale
;v_reg*=scale
;      unorm_reg = regress_3d(uband,unorm,dims)

      wind=create_struct('u',u_reg,'v',v_reg,'x',eralon,'y',eralat)
;      cvar=create_struct('cvar',unorm_reg,'x',eralon,'y',eralat)

      figname=ifigdir+var_plot+'_'+bandtag+'_reg_cross'+strtrim(icross,2)
      figspecs.figname=figname

      wrf_myanmar_map_plot, dirs, var_reg, lon, lat, figspecs, wind=wind, cvar=cvar, /noscalewind

  endif

endfor ; iband

print,'Done!!'
end
