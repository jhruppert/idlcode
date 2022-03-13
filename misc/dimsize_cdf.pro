; 
; Returns the size of a dimension, given it's name (as a string), in a netCDF file.
; 
; DIMSTR can be either a single string name or an array of names.
; If the latter, an array of dimsizes will be returned.
; 
; James Ruppert
; 28.10.16
; 

function dimsize_cdf , fname, dimstr
  
  nstr = size(dimstr,/dimensions)
  
  fid = ncdf_open(fname,/nowrite)
  
  dimid = ncdf_dimid(fid,dimstr[0])
  ncdf_diminq, fid, dimid, name, dimsize
  
  for istr=1,nstr[0]-1 do begin
    dimid = ncdf_dimid(fid,dimstr[istr])
    ncdf_diminq, fid, dimid, name, idimsize
    dimsize = [dimsize,idimsize]
  endfor
  
  ncdf_close,fid
  
  return,dimsize
  
end
