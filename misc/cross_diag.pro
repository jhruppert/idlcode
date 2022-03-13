;
; Remap an input array from 2D to a 1D diagonal cross section, using the entire
; x- and y-span of the array.
;
; Method: take an average of the diagonal of a given width; trim off ends.
;
; Input array should have x,y dimensions first, followed by any
; additional dimensions, which will be looped over.
;
; input:      (input) input array (input = input[x,y[,d1,d2])]).
; XXX iang:       (input) 0: LL-->UR, 1: UL-->LR
; x,y:        (input) x,y dimensions
; buffer:     (input) averaging distance in perpendicular dimension (in units of x,y).
; x,y_bounds: (input) bounds of cross section
; lonout:     (optional return) lon of cross section (center of transect)
; latout:     (optional return) lat " ".
;
; Returned array is of the form array[xnew[,d1,d2]], which is averaged across the perpendicular dimension.
;
; James Ruppert
; 2/5/19
function cross_diag, input, x, y, buffer, x_bounds=x_bounds, y_bounds=y_bounds, lonout=xlon0, latout=xlat0

;INPUT SPECS
  dims=size(input,/dimensions)
  nd=n_elements(dims)
;  x0=indgen(dims[0])
;  y0=indgen(dims[1])
  nx=dims[0]
  ny=dims[1]
  extdim = (nd gt 2)
  if n_elements(x) ne nx then message,'Check dimensions!'
  if n_elements(y) ne ny then message,'Check dimensions!'

;TRIANGULATION FOR REMAPPING
;  x2=fltarr(nx,ny) & y2=x2
;  for iy=0,ny-1 do x2[*,iy]=x
;  for ix=0,nx-1 do y2[ix,*]=y
;  triangulate,x2,y2,tri
  x2d=fltarr(nx,ny)
  y2d=x2d
  for ix=0,nx-1 do y2d[ix,*]=y
  for iy=0,ny-1 do x2d[*,iy]=x
  triangulate,x2d,y2d,tri

;CHECK FOR REVERSE LATITUDINAL SECTION
  if y_bounds[1] lt y_bounds[0] then irev=1 else irev=0

;SPECS FOR CROSS SECTION
;  nc=round(sqrt(1.+nx^2+ny^2)) ; PRESERVES APPROXIMATE RESOLUTION OF INPUT FIELD
  nxc=n_elements(where(x ge x_bounds[0] and x le x_bounds[1]))
  if irev then $
    nyc=n_elements(where(y ge y_bounds[1] and y le y_bounds[0])) $
  else $
    nyc=n_elements(where(y ge y_bounds[0] and y le y_bounds[1]))
  nc=round(sqrt(1.+nxc^2+nyc^2)) ; PRESERVES APPROXIMATE RESOLUTION OF INPUT FIELD
;  xc=interpol(x,nc)
;  yc=interpol(y,nc)
  ;REVERSE Y IF NECESSARY
;  if iang eq 1 then yc=reverse(temporary(yc))
;  isign=1
;  if iang eq 1 then isign=-1

;DETERMINE ANGLE (DEG) FROM INPUT COORDINATES
;  dx=max(x)-min(x)
;  dy=max(y)-min(y)
  dx=x_bounds[1]-x_bounds[0]
  dy=y_bounds[1]-y_bounds[0]
  if irev then dy*=-1.
  th = atan( dy / dx ) * 180./!pi

;  ;FOR TRIMMING EXCESS
;    th_trim = th
;    if th gt 45 then th_trim = 90 - th

;HALF-LENFTH OF BUFFER
  hbuff=0.5*buffer
  nhbuff = round(hbuff/(x[1]-x[0]))
  ibuff = hbuff*findgen(nhbuff-1)/(nhbuff-2) ; fraction of half buffer length
  ibuff = [-1*reverse(ibuff),0,ibuff]
  nbuff = n_elements(ibuff)
;  ntrim = round(1.*nhbuff/tan(th_trim*!pi/180.)) ; points trimmed at each end of cross section

;TRIM ENDS OF DIMENSIONAL ARRAYS
;  nc-=2*ntrim
;  xc=xc[ntrim:ntrim+nc-1]
;  yc=yc[ntrim:ntrim+nc-1]
;  x_bounds=xc[[0,nc-1]]
;  y_bounds=yc[[0,nc-1]]

;NEW CROSS SECTION LAT/LON COORDINATES
  xlon0=findgen(nc)*(x_bounds[1]-x_bounds[0])/(nc-1)+min(x_bounds)
  y_diff = y_bounds[1]-y_bounds[0]
;  if irev then y_diff *= -1.
  xlat0=findgen(nc)*y_diff/(nc-1)+min(y_bounds)
  if irev then xlat0=reverse(-1.*findgen(nc)*y_diff/(nc-1)+min(y_bounds))
  xlon=fltarr(nc,nbuff)
  xlat=xlon
  for ic=0,nc-1 do begin
      xlon[ic,*] = xlon0[ic] + ibuff * sin(th*!dtor)
      xlat[ic,*] = xlat0[ic] + ibuff * cos(th*!dtor)
  endfor

if extdim then $
  cross=fltarr([nc,nbuff,dims[2:nd-1]]) $
else $
  cross=fltarr(nc,nbuff)


;RUN REMAP

if extdim then begin

  if nd eq 3 then begin

;    for ic=0,nc-1 do begin
;      xband = xc[ic] +         ibuff * sin(th*!dtor)
;      yband = yc[ic] - isign * ibuff * cos(th*!dtor)
      for id1=0,dims[2]-1 do begin
;        tmp=griddata(x2,y2,reform(input[*,*,id1]),/linear,triangles=tri,xout=xband,yout=yband,/grid,missing=!values.f_nan)
;        cross[ic,id1]=mean(tmp[indgen(nbuff),indgen(nbuff)],/double,/nan)
        tmp=griddata(x2d,y2d,reform(input[*,*,id1]),/linear,triangles=tri,xout=xlon,yout=xlat,missing=!values.f_nan)
        cross[*,*,id1]=reform(tmp,nc,nbuff) & tmp=0
      endfor
;    endfor

  endif else if nd eq 4 then begin

;    for ic=0,nc-1 do begin
;      xband = xc[ic] +         ibuff * sin(th*!dtor)
;      yband = yc[ic] - isign * ibuff * cos(th*!dtor)
      for id1=0,dims[2]-1 do begin
        for id2=0,dims[3]-1 do begin
;          tmp=griddata(x2,y2,reform(input[*,*,id1,id2]),/linear,triangles=tri,xout=xband,yout=yband,/grid,missing=!values.f_nan)
;          cross[ic,id1,id2]=mean(tmp[indgen(nbuff),indgen(nbuff)],/double,/nan)
          tmp=griddata(x2d,y2d,reform(input[*,*,id1,id2]),/linear,triangles=tri,xout=xlon,yout=xlat,missing=!values.f_nan)
          cross[*,*,id1,id2]=reform(tmp,nc,nbuff) & tmp=0
        endfor
      endfor
;    endfor

  endif else message,'Check dimensions. Might have an extra?'

endif else begin

;  for ic=0,nc-1 do begin
;     xband = xc[ic] +         ibuff * sin(th*!dtor)
;     yband = yc[ic] - isign * ibuff * cos(th*!dtor)
;    xints = interpol(indgen(nx),x,xband)
;    yints = interpol(indgen(ny),y,yband)
;    tmp=griddata(x2,y2,input,/linear,triangles=tri,xout=xband,yout=yband,/grid,missing=!values.f_nan)
;    cross[ic]=mean(tmp[indgen(nbuff),indgen(nbuff)],/double,/nan)
    tmp=griddata(x2d,y2d,input,/linear,triangles=tri,xout=xlon,yout=xlat,missing=!values.f_nan)
    cross=reform(tmp,nc,nbuff) & tmp=0
;  endfor

endelse

return,cross

end
