pro readall_wrf_varlon,filelist

nfil=n_elements(filelist)

for ifil=0,nfil-1 do begin
  temp=' '
  lon=read_nc_var(filelist[ifil],'XLONG',count=[10,10,1],offset=[0,0,0])
;  stats,lon
endfor

end
;------------------------------
pro readall_post_varlon,filelist

nfil=n_elements(filelist)

for ifil=0,nfil-1 do begin
  temp=' '
lon=read_nc_var(filelist[ifil],'XLONG',count=[10,10,1],offset=[0,0,0])
help,lon
if ifil eq 5 then exit
endfor

end
;------------------------------
; 
; Read a single variable to be output from TC model output.

; James Ruppert
; 3/12/19
; 
pro run_tc_readall

for itc=0,1 do begin

if itc eq 0 then $
  tcname='maria' $
else $
  tcname='haiyan'

if tcname eq 'maria' then begin
  tcyear='2017'
  hurdat=read_hurdat(tcname,tcyear)
endif else if tcname eq 'haiyan' then begin
  tcyear='2013'
  hurdat=read_jtwcdat(tcname)
endif

subdir='redux/'+tcname
config_dir, dirs=dirs
upperdir=dirs.scdir+'tc/output/'+subdir+'/'
spawn,'ls '+upperdir,cases
nc=n_elements(cases)

for ic=0,nc-1 do begin
  casedir=upperdir+cases[ic]
  spawn,'echo '+casedir

  spawn,'ls '+casedir+'/wrf*',filelist,err
  if keyword_set(err) then continue
  readall_wrf_varlon,filelist
continue
  spawn,'ls '+casedir+'/post',filelist,err
  if keyword_set(err) then continue

exit
endfor
exit

endfor ; itc

print,'DONE!!'
end
