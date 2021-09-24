;
; Convert an eps file to a png file.
;
; Input:
;
;   Filename  =  complete filename up to and excluding extension (i.e., without the ".eps").
;   Res       =  (optional) Integer value for output png file resolution. Res = 100 if no value provided.
;
; Options:
;
;   remove_eps = set this to remove the corresponding eps file.
;
; 3/7/14
; James Ruppert
; ruppert@atmos.colostate.edu

pro convert_png , filename , res = res, remove_eps=remove_eps
  
  if ~keyword_set(res) then res=100
  res=round(res)
  
  convert='/usr/bin/gs -q -dBATCH -dNOPROMPT -dNOPAUSE -dEPSCrop -dGraphicsAlphaBits=4 -dTextAlphaBits=4 -sDEVICE=png16m '+$
    '-r'+strtrim(res,2)+' -sOutputFile='+filename+'.png '+filename+'.eps'
;  convert='convert -density '+strtrim(res,2)+' '+filename+'.eps '+filename+'.png'
  spawn,convert,out
  
  if keyword_set(remove_eps) then $
    spawn,'rm '+filename+'.eps'
  
end
