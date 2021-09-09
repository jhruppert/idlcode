; 
; Calculate and write out PW from hourly ERA5 data for JJAS 2013-2017.
;
; James Ruppert
; 12/21/20
; 
pro run_calc_pw_era5

config_dir,dirs=dirs

;USE FULL JJAS 2013-2017
  yy_plot=[2013,2017]
  mm_plot=[6,9]
  dd_plot=[1,30] ; inclusive

  ;DATE STRING
    form2='(i2.2)'
    form4='(i4)'
    dat_str=string(mm_plot[0],format=form2)+string(dd_plot[0],format=form2)+strmid(strtrim(yy_plot[0],2),2,2)+'-'+$
            string(mm_plot[1],format=form2)+string(dd_plot[1],format=form2)+strmid(strtrim(yy_plot[1],2),2,2)


;----OB DIRECTORIES--------------------

  maindir=dirs.scdir+'myanmar/'
  era_dir=maindir+'era5/jjas_2013-2017/'
  npd_era=24


;----ONE TIME ARRAY FOR DAILY AVERAGE TIME SERIES--------------------

  ;ACCOMMODATE MULTIPLE YEARS
    if yy_plot[0] eq yy_plot[1] then $
      time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
        final=julday(mm_plot[1],dd_plot[1],yy_plot[1],23,59,59),step_size=1,units='hours') $
    else begin
      time=timegen(start=julday(mm_plot[0],dd_plot[0],yy_plot[0],0,0,0),$
        final=julday(mm_plot[1],dd_plot[1],yy_plot[0],23,59,59),step_size=1,units='hours')
      for iyy=yy_plot[0]+1,yy_plot[1] do begin
        itime=timegen(start=julday(mm_plot[0],dd_plot[0],iyy,0,0,0),$
          final=julday(mm_plot[1],dd_plot[1],iyy,23,59,59),step_size=1,units='hours')
        time=[time,itime]
      endfor
    endelse
    nt=n_elements(time)
    nd=nt/npd_era


;----READ ERA--------------------

  ;DIMENSIONS
  era_fil=era_dir+'ERA5-2013-06-01-2013-10-01-pl.nc'
  eralon=read_nc_var(era_fil,'lon')
  eralat=read_nc_var(era_fil,'lat')
  nxera=n_elements(eralon)
  nyera=n_elements(eralat)

  ;PLEVS
  p_era=reform(read_nc_var(era_fil,'plev'))*1d-2 ; Pa --> hPa
  nzera=n_elements(p_era)

  pw=fltarr(nxera,nyera,nt)
  it_year=0
  for year=yy_plot[0],yy_plot[1] do begin
;  for year=yy_plot[0],yy_plot[0] do begin

    print,'Running:',year

    era_fil=era_dir+'ERA5-'+strtrim(year,2)+'-06-01-'+strtrim(year,2)+'-10-01-pl.nc'
    era_sfil=era_dir+'ERA5-'+strtrim(year,2)+'-06-01-'+strtrim(year,2)+'-10-01-sl.nc'

    itime=timegen(start=julday(mm_plot[0],dd_plot[0],year,0,0,0),$
        final=julday(mm_plot[1],dd_plot[1],yy_plot[0],23,59,59),step_size=1,units='hours')
    i_nt=n_elements(itime)

    count=[nxera,nyera,nzera,i_nt] & offset=[0,0,0,0] ; x,y,z,t
    count2d=[nxera,nyera,i_nt] & offset2d=[0,0,0]

    ;CALC PW
      qv=reform(read_nc_var(era_fil,'var133',count=count,offset=offset)) ; kg/kg
      pr=fltarr(nxera,nyera,nzera,i_nt) & for iz=0,nzera-1 do pr[*,*,iz,*]=p_era[iz]
      pr*=1d2 ; --> Pa
      ;REVERSE VERTICAL
        qv=reverse(temporary(qv),3)
        pr=reverse(temporary(pr),3)
  
      ;FILL IN SURFACE VALUES
        prsfc=reform(read_nc_var(era_sfil,'var134',count=count2d,offset=offset2d)) ; Pa
        pr[*,*,0,*]=prsfc
  
      ipsel=indgen(nzera-1)+1
      for it=0,i_nt-1 do begin
      for ix=0,nxera-1 do begin
      for iy=0,nyera-1 do begin
        ip=where(pr[ix,iy,ipsel,it] lt pr[ix,iy,0,it])
        ip=[0,ip+1]
        dp=deriv(reform(pr[ix,iy,ip,it]))*(-1.)
        pw[ix,iy,it+it_year]=total(reform(qv[ix,iy,ip,it])*dp,/double)
      endfor
      endfor
      endfor

      it_year+=i_nt

  endfor ; iy

    pw/=9.81


;----WRITE OUT--------------------

fname=era_dir+'ERA5-JJAS_2013-2017-pw.nc'
write_sing_ncvar,fname,pw,'pw',dimtag1='lon',dimtag2='lat',dimtag3='time'


print,'DONE!!'
end
