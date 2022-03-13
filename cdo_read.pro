pro cdo_read


;Directories

maindir='/scratch/06040/tg853394/'
dir = [$
       'era5/',$
       'imerg/',$
       'tc/output/redux/']

dir=maindir+dir

print,dir

ndir=n_elements(dir)

nwrite=0
for idir=0,ndir-1 do begin

adir=dir[idir]

  if strmatch(adir,'*tc/output*') then dotc=1 else dotc=0
print,dotc

  spawn,'ls '+adir,files,err
  print,files
  
    nfil=n_elements(files)
    for ifil=0,nfil-1 do begin
    ;  NCDF_LIST, files[ifil], /DIMENSIONS, /GATT, OUT=textout, /QUIET, /VARIABLES
      fid = ncdf_open(adir+files[ifil],/nowrite)
      sout = ncdf_inquire(fid)
      ncdf_close,fid
      openw,1,'var_list_'+string(nwrite,format='(i3.3)')
        printf,1,sout
      close,1
      nwrite+=1
    endfor

;cdostr='cdo sinfov '+erdir+'/ERA5* > cdo_txt1.out 2>&1'

;cdostr='cdo sinfov '+imdir+'/*'
;cdostr='cdo sinfov '+imdir+'/3B-HHR.MS.MRG.3IMERG.2000-2020_JJAS_diurncomp.V06B.nc4'
;cdostr='/home1/06040/tg853394/builds/cdo/bin/cdo --help'
;print,cdostr
;spawn,cdostr,out,err
;print,out
;print
;print,err

;cdo sinfov ${tcdir}/*/*/azim* > cdo_azim.out 2>&1
;cdo sinfov ${tcdir}/*/*/wrf* > cdo_wrf.out 2>&1
;cdo sinfov ${tcdir}/*/*/post/* > cdo_post.out 2>&1

endfor

end
