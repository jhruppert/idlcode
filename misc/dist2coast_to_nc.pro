pro dist2coast_to_nc

ascii='dist2coast.txt'

npts=40500000
print,npts
cdist=fltarr(3,npts)

tmp=' '
openr,1,ascii
  for ip=0,npts-1 do begin
    readf,1,tmp
    cdist[*,ip]=float(strsplit(tmp,/extract))
  endfor
close,1

ncfil='dist2coast.nc'
write_sing_ncvar,ncfil,cdist,'coast_dist';,dim1=3,dim2=npts,dimtag1='x-y-distkm',dimtag2='npts'

end
