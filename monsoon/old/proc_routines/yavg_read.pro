;
; Read converted ICON data written to netCDF file by yavg_write_ncdf
;
; James Ruppert
; 16.12.17
;
pro yavg_read,azim_fils,fil_ind,t_ind,struct_vars ; return vars


  nt = n_elements(t_ind)

  uniq_fil = uniq(fil_ind)
  nfil_sel = n_elements(uniq_fil)

  t0_ind = [0,uniq_fil[0:nfil_sel-2]+1 ] ; file-block indices for output arrays

  ;GET DIM SPECS
  fid=ncdf_open(azim_fils[0],/nowrite)
  dimid = ncdf_dimid(fid,'x')
  ncdf_diminq,fid,dimid,name,nrad
  dimid = ncdf_dimid(fid,'z')
  ncdf_diminq,fid,dimid,name,nz
  dimid = ncdf_dimid(fid,'nd')
  ncdf_diminq,fid,dimid,name,nd
  dimid = ncdf_dimid(fid,'nd2')
  ncdf_diminq,fid,dimid,name,nd2
  ncdf_close,fid

  x_yavg = read_nc_var(azim_fils[fil_ind[ uniq_fil[0] ]],'x_yavg') ; [ m ]

  pcp_yavg=fltarr(nrad,nt)
  lh_yavg=fltarr(nrad,nt)
  sh_yavg=fltarr(nrad,nt)
  pw_yavg=fltarr(nrad,nt)
  t_g_yavg=fltarr(nrad,nt)

  w_so_yavg=fltarr(nrad,nd,nt)
  t_so_yavg=fltarr(nrad,nd2,nt)

  u_yavg=fltarr(nrad,nz,nt)
  v_yavg=fltarr(nrad,nz,nt)
  w_yavg=fltarr(nrad,nz+1,nt)
  lw_yavg=fltarr(nrad,nz,nt)
  sw_yavg=fltarr(nrad,nz,nt)
  ex_yavg=fltarr(nrad,nz,nt)
  tmpk_yavg=fltarr(nrad,nz,nt)
  qv_yavg=fltarr(nrad,nz,nt)
  qi_yavg=fltarr(nrad,nz,nt)
  qc_yavg=fltarr(nrad,nz,nt)

for isel=0,nfil_sel-1 do begin

  ;READ-FILE INFO
    ifil = fil_ind[ uniq_fil[isel] ]
    fil = azim_fils[ifil]
    nt_fil = n_elements(where(fil_ind eq ifil))

  t_ind_get = indgen(nt_fil)

  count=[nrad,nt_fil] & offset=[0,t_ind_get[0]]

  pcp_yavg[*,t_ind_get] = read_nc_var(fil,'pcp_yavg',count=count,offset=offset) ; [ mm / time step ]
  lh_yavg[*,t_ind_get] = read_nc_var(fil,'lh_yavg',count=count,offset=offset) ; [ W/m2 ]
  sh_yavg[*,t_ind_get] = read_nc_var(fil,'sh_yavg',count=count,offset=offset) ; [ W/m2 ]
  pw_yavg[*,t_ind_get] = read_nc_var(fil,'pw_yavg',count=count,offset=offset) ; [ mm ]
  t_g_yavg[*,t_ind_get] = read_nc_var(fil,'t_g_yavg',count=count,offset=offset) ; [ K ]

  count=[nrad,nd,nt_fil] & offset=[0,0,t_ind_get[0]]
  w_so_yavg[*,*,t_ind_get] = read_nc_var(fil,'w_so_yavg',count=count,offset=offset) ; [ kg/m2 ]
  count=[nrad,nd2,nt_fil] & offset=[0,0,t_ind_get[0]]
  t_so_yavg[*,*,t_ind_get] = read_nc_var(fil,'t_so_yavg',count=count,offset=offset) ; [ K ]


  count=[nrad,nz,nt_fil] & offset=[0,0,t_ind_get[0]]

  u_yavg[*,*,t_ind_get] = read_nc_var(fil,'u_yavg',count=count,offset=offset) ; [ m/s ]
  v_yavg[*,*,t_ind_get] = read_nc_var(fil,'v_yavg',count=count,offset=offset) ; [ m/s ]
  lw_yavg[*,*,t_ind_get] = read_nc_var(fil,'lw_yavg',count=count,offset=offset) ; [ K/s ]
  sw_yavg[*,*,t_ind_get] = read_nc_var(fil,'sw_yavg',count=count,offset=offset) ; [ K/s ]
  ex_yavg[*,*,t_ind_get] = read_nc_var(fil,'ex_yavg',count=count,offset=offset) ; [ Km2/kg/s ]
  tmpk_yavg[*,*,t_ind_get] = read_nc_var(fil,'tmpk_yavg',count=count,offset=offset) ; [ K ]
  qv_yavg[*,*,t_ind_get] = read_nc_var(fil,'qv_yavg',count=count,offset=offset) ; [ kg/kg ]
  qi_yavg[*,*,t_ind_get] = read_nc_var(fil,'qi_yavg',count=count,offset=offset) ; [ kg/kg ]
  qc_yavg[*,*,t_ind_get] = read_nc_var(fil,'qc_yavg',count=count,offset=offset) ; [ kg/kg ]


  count=[nrad,nz+1,nt_fil] & offset=[0,0,t_ind_get[0]]

  w_yavg[*,*,t_ind_get] = read_nc_var(fil,'w_yavg',count=count,offset=offset) ; [ m/s ]

endfor ; isel


  struct_vars = { x_yavg:x_yavg, pcp_yavg:pcp_yavg, $
                  lh_yavg:lh_yavg, sh_yavg:sh_yavg, pw_yavg:pw_yavg,$
                  t_g_yavg:t_g_yavg, w_so_yavg:w_so_yavg, t_so_yavg:t_so_yavg,$
                  u_yavg:u_yavg, v_yavg:v_yavg, w_yavg:w_yavg, $
                  lw_yavg:lw_yavg, sw_yavg:sw_yavg, $
                  ex_yavg:ex_yavg, tmpk_yavg:tmpk_yavg, qv_yavg:qv_yavg, $
                  qi_yavg:qi_yavg, qc_yavg:qc_yavg }


end
