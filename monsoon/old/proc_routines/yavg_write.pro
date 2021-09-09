; 
; Convert variables from irregular ICON grid) onto a 2D grid
; by averaging across y-points.
; 
; Also contains yavg_write_ncdf.
; 
; James Ruppert
; 16.12.17
; 
pro yavg_write

;SOME SETTINGS

  imist=1
  exp = 'terra_t3'

  t_sel=[0,30]

  newfil_tag = '_yavg'


;MODEL SPECS

  icon_rce_exp_presets,exp,jult_fil0=jult_fil0,nfil=nfil,ndays=ndays,$
    delta_file=delta_file,grid_fil=grid_fil,npday=npday,$
    moddir=moddir,datdir=datdir,figdir=figdir,imist=imist
  ntpfil=npday*delta_file

;  nfil+=1

  ;WHOLE DAYS ONLY
    if n_elements(t_sel) gt 1 then $
      t_sel=[t_sel[0],t_sel[1]-1d/npday]

  ;Z-LEVELS
    nz=75

  ;GET TIME ARRAY
    time=dindgen(npday*ndays)/npday
    nt=n_elements(time)

  ;DIMENSIONS
    dims=icon_domain_info(moddir+grid_fil,dx=dx)
    nx=n_elements(dims.x)
    xmax0=max(dims.x) & xmin0=min(dims.x) & lenx=xmax0-xmin0
    ymax0=max(dims.y) & ymin0=min(dims.y) & leny=ymax0-ymin0

  ;MODEL OUTPUT FILES
    icon_file_multi,nfil,jult_fil0,delta_file,'3d',moddir,exp,$
      datfils=datfils_3d
    icon_file_multi,nfil,jult_fil0,delta_file,'2d',moddir,exp,$
      datfils=datfils
    nfil=n_elements(datfils)

  ;SELECTED TIME(S) TO CALCULATE

    nt_plot = n_elements(t_sel) ; time range or single time?
    t_ind_calc = max(where(time le t_sel[0]))
    if nt_plot eq 2 then begin
      nt_calc = max(where(time le t_sel[1])) - t_ind_calc + 1
      ;t_ind_calc = [ t_ind_calc , max(where(time le t_sel[1])) ]
      t_ind_calc = indgen(nt_calc) + t_ind_calc
    endif else nt_calc=1

    t_ind_multi,time,time[t_ind_calc],ntpfil,$
      t_ind=t_ind,fil_ind=fil_ind
    uniq_fil = uniq(fil_ind)
    nfil_sel = n_elements(uniq_fil)

    t0_ind = [0,uniq_fil[0:nfil_sel-2]+1 ] ; file-block indices for output arrays

  ;Z-LEVEL SELECTION FOR ADDITIONAL 3D VARS
;    zsel=29 ; ~12 km
;    zsel=49 ; ~5 km
;    zsel=63 ; ~1500 m
    nz_read=nz
    if keyword_set(zsel) then z_offset=zsel else z_offset=nz-nz_read

  ;HEIGHT
    count=[1,nz_read+1] & offset=[0,z_offset]
    hghtw = read_nc_var(datfils_3d[0],'z_ifc',count=count,offset=offset)
    hghtw = reform(hghtw)
    count=[1,nz_read] & offset=[0,z_offset]
    hght = read_nc_var(datfils_3d[0],'z_mc',count=count,offset=offset)
    hght = reform(hght)


;PROCESS DATA

  ;NEW ARRAY SPECS
    nrad=round(lenx/(dx*2))
    ;Double the x-spacing to smooth slightly
    rad = findgen(nrad)*(2*dx)
    hbin = dx
    rad += hbin

;nfil = n_elements(ifil_sel)
;for ifil=0,nfil-1 do begin
for isel=0,nfil_sel-1 do begin

  ;READ-FILE INFO
;    ifil = ifil_sel[isel]
    ifil = fil_ind[ uniq_fil[isel] ]
    nt_fil = n_elements(where(fil_ind eq ifil))
    fil = datfils[ifil]
    fil3d = datfils_3d[ifil]

  ;WRITE-FILE INFO
    fil_tag = strsplit(fil,'/',/extract)
    fil_tag = fil_tag[n_elements(fil_tag)-1]
    fil_tag = (strsplit(fil_tag,'.',/extract))[0]
    outfil = moddir+fil_tag+newfil_tag+'.nc'

  ;Y-AVERAGED ARRAYS

    pres_pol = fltarr(nrad,nt_fil)
    pcp_pol = pres_pol
    ;    wsp_pol = pres_pol
    lh_pol = pres_pol
    sh_pol = pres_pol
    pw_pol = pres_pol
    t_g_pol = pres_pol
    w_so_pol = fltarr(nrad,8,nt_fil)
    t_so_pol = fltarr(nrad,9,nt_fil)
    pres_pol=0

    u_pol = fltarr(nrad,nz_read,nt_fil)
    v_pol = u_pol
    w_pol = fltarr(nrad,nz_read+1,nt_fil)
    lw_pol = u_pol
    sw_pol = u_pol
    ex_pol = u_pol
    tmpk_pol = u_pol
    qv_pol = u_pol
    qi_pol = u_pol
    qc_pol = u_pol

  it0 = t0_ind[isel]


;  for it=0,ntpfil-1 do begin
  for it=0,nt_fil-1 do begin

    tim = time[ t_ind[it] + ifil*ntpfil ]

    print,'TIME: ',tim

    ;READ VARIABLES

      t_read = t_ind[it0+it]

      ;2D VARS
        count=[nx,1] & offset=[0,t_read]
;        pcp = read_nc_var(fil,'tot_prec',count=count,offset=offset) ; mm
;        count=[nx,1] & offset=[0,t_read]
        lh = read_nc_var(fil,'lhfl_s',count=count,offset=offset) ; W/m2
        sh = read_nc_var(fil,'shfl_s',count=count,offset=offset) ; W/m2
        pw = read_nc_var(fil,'tqv_dia',count=count,offset=offset) ; mm
        t_g = read_nc_var(fil,'t_g',count=count,offset=offset) ; K

        count=[nx,8,1] & offset=[0,0,t_read]
        w_so = read_nc_var(fil,'w_so',count=count,offset=offset) ; kg/m2
        count=[nx,9,1] & offset=[0,0,t_read]
        t_so = read_nc_var(fil,'t_so',count=count,offset=offset) ; K

        ;PRECIP RATE VIA CENTERED DIFFERENCE
          tim_pcp = [ tim-1d/npday , tim+1d/npday ]
          t_ind_multi,time,tim_pcp,ntpfil,$
            t_ind=tind_pcp,fil_ind=fil_pcp
          pcp_m1 = read_nc_var(datfils[fil_pcp[0]],'tot_prec',$
            count=[nx,1],offset=[0,tind_pcp[0]])
          pcp_p1 = read_nc_var(datfils[fil_pcp[1]],'tot_prec',$
            count=[nx,1],offset=[0,tind_pcp[1]])
          pcp = 0.5*(pcp_p1 - pcp_m1) ; mm / time step
          pcp_p1=0 & pcp_m1=0

      ;3D VARS

        count=[nx,nz_read,1] & offset=[0,z_offset,t_read]
        u = read_nc_var(fil3d,'u',count=count,offset=offset)
        v = read_nc_var(fil3d,'v',count=count,offset=offset)

        lw = read_nc_var(fil3d,'ddt_temp_radlw',count=count,offset=offset)
        sw = read_nc_var(fil3d,'ddt_temp_radsw',count=count,offset=offset)
        ex = read_nc_var(fil3d,'exner',count=count,offset=offset)
        tmpk = read_nc_var(fil3d,'temp',count=count,offset=offset)
        qv = read_nc_var(fil3d,'tot_qv_dia',count=count,offset=offset)
        qi = read_nc_var(fil3d,'tot_qi_dia',count=count,offset=offset)
        qc = read_nc_var(fil3d,'tot_qc_dia',count=count,offset=offset)

        count=[nx,nz_read+1,1] & offset=[0,z_offset,t_read]
        w = read_nc_var(fil3d,'w',count=count,offset=offset)

    ;MODIFY X,Y TO PLACE P-MIN INTO CENTER
      xtmp = dims.x - xmin0

    ;AVERAGE OVER Y

    for irad=0,nrad-1 do begin
;    for irad=0,100 do begin

      loc_ave = where( (xtmp ge (rad[irad]-hbin)) and $
                       (xtmp lt (rad[irad]+hbin)) , nl )

      if nl le 2 then begin

        pcp_pol[irad,it] = !values.f_nan
        lh_pol[irad,it] = !values.f_nan
        sh_pol[irad,it] = !values.f_nan
        pw_pol[irad,it] = !values.f_nan
        t_g_pol[irad,it] = !values.f_nan

        w_so_pol[irad,*,it] = !values.f_nan
        t_so_pol[irad,*,it] = !values.f_nan

        u_pol[irad,*,it] = !values.f_nan
        v_pol[irad,*,it] = !values.f_nan
        lw_pol[irad,*,it] = !values.f_nan
        sw_pol[irad,*,it] = !values.f_nan
        ex_pol[irad,*,it] = !values.f_nan
        tmpk_pol[irad,*,it] = !values.f_nan
        qv_pol[irad,*,it] = !values.f_nan
        qc_pol[irad,*,it] = !values.f_nan
        qi_pol[irad,*,it] = !values.f_nan

        continue

      endif

      pcp_pol[irad,it] = mean(pcp[loc_ave],/double)
      lh_pol[irad,it] = mean(lh[loc_ave],/double)
      sh_pol[irad,it] = mean(sh[loc_ave],/double)
      pw_pol[irad,it] = mean(pw[loc_ave],/double)
      t_g_pol[irad,it] = mean(t_g[loc_ave],/double)

      w_so_pol[irad,*,it] = mean(w_so[loc_ave,*],dimension=1,/double)
      t_so_pol[irad,*,it] = mean(t_so[loc_ave,*],dimension=1,/double)

;      for iz=0,nz_read-1 do begin
        u_pol[irad,*,it] = mean(u[loc_ave,*],dimension=1,/double)
        v_pol[irad,*,it] = mean(v[loc_ave,*],dimension=1,/double)
        lw_pol[irad,*,it] = mean(lw[loc_ave,*],dimension=1,/double)
        sw_pol[irad,*,it] = mean(sw[loc_ave,*],dimension=1,/double)
        ex_pol[irad,*,it] = mean(ex[loc_ave,*],dimension=1,/double)
        tmpk_pol[irad,*,it] = mean(tmpk[loc_ave,*],dimension=1,/double)
        qv_pol[irad,*,it] = mean(qv[loc_ave,*],dimension=1,/double)
        qc_pol[irad,*,it] = mean(qc[loc_ave,*],dimension=1,/double)
        qi_pol[irad,*,it] = mean(qi[loc_ave,*],dimension=1,/double)
;      endfor

;      for iz=0,nz_read+1-1 do $
        w_pol[irad,*,it] = mean(w[loc_ave,*],dimension=1,/double)
;      if keyword_set(var1_tag) then $
;        for iz=0,nz_read+nzw-1 do $
;          var1_pol[irad,iz,it0+it] = mean(var1[loc_ave,iz],dimension=1,/double)
;
;      if keyword_set(var2_tag) then $
;        for iz=0,nz_read-1 do $
;          var2_pol[irad,iz,it0+it] = mean(var2[loc_ave,iz],dimension=1,/double)
;
;      if keyword_set(var3_tag) then $
;        for iz=0,nz_read-1 do $
;          var3_pol[irad,iz,it0+it] = mean(var3[loc_ave,iz],dimension=1,/double)
;
;      if keyword_set(var4_tag) then $
;        for iz=0,nz_read-1 do $
;          var4_pol[irad,iz,it0+it] = mean(var4[loc_ave,iz],dimension=1,/double)

    endfor ; irad

  endfor ; it

;  struct_vars = { pres_pol:pres_pol, pcp_pol:pcp_pol, $;wsp_pol:wsp_pol, $
;                  lh_pol:lh_pol, sh_pol:sh_pol, $
;                  u_pol:u_pol, v_pol:v_pol, $
;                  var1_pol:var1_pol, var2_pol:var2_pol, var3_pol:var3_pol, $
;                  var4_pol:var4_pol }


  ;WRITE OUT TO NETCDF FILE

  yavg_write_ncdf,outfil,rad,$
    pcp_pol, lh_pol, sh_pol, pw_pol, t_g_pol, w_so_pol, t_so_pol,$
    u_pol, v_pol, w_pol, lw_pol, sw_pol, ex_pol, tmpk_pol,$
    qv_pol, qi_pol, qc_pol


endfor ; isel (ifil)


end

; 
; Write out converted data to a netCDF file.
; 
pro yavg_write_ncdf,outfil,x,$
    pcp_pol, lh_pol, sh_pol, pw_pol, t_g_pol, w_so_pol, t_so_pol,$
    u_pol, v_pol, w_pol, lw_pol, sw_pol, ex_pol, tmpk_pol,$
    qv_pol, qi_pol, qc_pol

print,'Writing out!'

;GET DIMENSIONS
  specs=size(u_pol,/dimensions)
  nx = specs[0]
  nz = specs[1]
  nt = specs[2]
  nd = 8
  nd2 = 9


filename=outfil
rm_str='rm '+filename
spawn,rm_str,rm_out,rm_err
fid=ncdf_create(filename,/clobber,/netcdf4_format)
;fid=ncdf_create(filename,/netcdf4_format)

;DEFINE DIMENSIONS

;thid=ncdf_dimdef(fid,'theta',nth)
xid=ncdf_dimdef(fid,'x',nx)
zid=ncdf_dimdef(fid,'z',nz)
zid2=ncdf_dimdef(fid,'z2',nz+1)
tid=ncdf_dimdef(fid,'nt',nt)
did=ncdf_dimdef(fid,'nd',nd)
did2=ncdf_dimdef(fid,'nd2',nd2)

;DEFINE VARIABLES

varid_x=ncdf_vardef(fid,'x_yavg',xid,/float)

;2D VARS

varid_pcp=ncdf_vardef(fid,'pcp_yavg',[xid,tid],/float)
varid_lh=ncdf_vardef(fid,'lh_yavg',[xid,tid],/float)
varid_sh=ncdf_vardef(fid,'sh_yavg',[xid,tid],/float)
varid_pw=ncdf_vardef(fid,'pw_yavg',[xid,tid],/float)
varid_t_g=ncdf_vardef(fid,'t_g_yavg',[xid,tid],/float)

varid_w_so=ncdf_vardef(fid,'w_so_yavg',[xid,did,tid],/float)
varid_t_so=ncdf_vardef(fid,'t_so_yavg',[xid,did2,tid],/float)

;3D VARS

varid_u=ncdf_vardef(fid,'u_yavg',[xid,zid,tid],/float)
varid_v=ncdf_vardef(fid,'v_yavg',[xid,zid,tid],/float)
varid_w=ncdf_vardef(fid,'w_yavg',[xid,zid2,tid],/float)
varid_lw=ncdf_vardef(fid,'lw_yavg',[xid,zid,tid],/float)
varid_sw=ncdf_vardef(fid,'sw_yavg',[xid,zid,tid],/float)
varid_ex=ncdf_vardef(fid,'ex_yavg',[xid,zid,tid],/float)
varid_tmp=ncdf_vardef(fid,'tmpk_yavg',[xid,zid,tid],/float)
varid_qv=ncdf_vardef(fid,'qv_yavg',[xid,zid,tid],/float)
varid_qi=ncdf_vardef(fid,'qi_yavg',[xid,zid,tid],/float)
varid_qc=ncdf_vardef(fid,'qc_yavg',[xid,zid,tid],/float)

;ADD ATTRIBUTES

ncdf_attput,fid,varid_x,'def','new x-dimension'
ncdf_attput,fid,varid_x,'units','m'

;2D VARS

ncdf_attput,fid,varid_pcp,'def','accum. precipitation'
ncdf_attput,fid,varid_pcp,'units','mm'
ncdf_attput,fid,varid_lh,'def','latent heat flux'
ncdf_attput,fid,varid_lh,'units','W/m2'
ncdf_attput,fid,varid_sh,'def','sensible heat flux'
ncdf_attput,fid,varid_sh,'units','W/m2'
ncdf_attput,fid,varid_pw,'def','precip water'
ncdf_attput,fid,varid_pw,'units','mm'
ncdf_attput,fid,varid_t_g,'def','surface temperature'
ncdf_attput,fid,varid_t_g,'units','K'

ncdf_attput,fid,varid_w_so,'def','soil water content'
ncdf_attput,fid,varid_w_so,'units','kg/m2'
ncdf_attput,fid,varid_t_so,'def','soil temperature'
ncdf_attput,fid,varid_t_so,'units','K'

;3D VARS

ncdf_attput,fid,varid_u,'def','u-wind'
ncdf_attput,fid,varid_u,'units','m/s'
ncdf_attput,fid,varid_v,'def','v-wind'
ncdf_attput,fid,varid_v,'units','m/s'
ncdf_attput,fid,varid_w,'def','vertical motion (on mass levels)'
ncdf_attput,fid,varid_w,'units','m/s'
ncdf_attput,fid,varid_lw,'def','Longwave heat source'
ncdf_attput,fid,varid_lw,'units','K/d'
ncdf_attput,fid,varid_sw,'def','Shortwave heat source'
ncdf_attput,fid,varid_sw,'units','K/d'
ncdf_attput,fid,varid_ex,'def','Exner function'
ncdf_attput,fid,varid_ex,'units','unitless'
ncdf_attput,fid,varid_tmp,'def','temperature'
ncdf_attput,fid,varid_tmp,'units','K'
ncdf_attput,fid,varid_qv,'def','water vapor mixing ratio'
ncdf_attput,fid,varid_qv,'units','kg/kg'
ncdf_attput,fid,varid_qc,'def','Specific cloud water content'
ncdf_attput,fid,varid_qc,'units','kg/kg'
ncdf_attput,fid,varid_qi,'def','Specific cloud ice content'
ncdf_attput,fid,varid_qi,'units','kg/kg'


;END DEFINE MODE

ncdf_control,fid,/endef

;INSERT VARIABLES

ncdf_varput,fid,varid_x,x

;2D VARS

ncdf_varput,fid,varid_pcp,pcp_pol
ncdf_varput,fid,varid_lh,lh_pol
ncdf_varput,fid,varid_sh,sh_pol
ncdf_varput,fid,varid_pw,pw_pol
ncdf_varput,fid,varid_t_g,t_g_pol

ncdf_varput,fid,varid_w_so,w_so_pol
ncdf_varput,fid,varid_t_so,t_so_pol

;3D VARS

ncdf_varput,fid,varid_u,u_pol
ncdf_varput,fid,varid_v,v_pol
ncdf_varput,fid,varid_w,w_pol
ncdf_varput,fid,varid_lw,lw_pol
ncdf_varput,fid,varid_sw,sw_pol
ncdf_varput,fid,varid_ex,ex_pol
ncdf_varput,fid,varid_tmp,tmpk_pol
ncdf_varput,fid,varid_qv,qv_pol
ncdf_varput,fid,varid_qi,qi_pol
ncdf_varput,fid,varid_qc,qc_pol

;DONE

ncdf_close,fid


print,'Done writing out!'


end
