;
; Create a one or two-column array of positions for multi-panel plotting.
;
; Returns array in arrangement top-bottom, left-right
;
; James Ruppert
; 11 Nov 2011
;
function position_multi, position_in, ncols, nrows, xbuff=xbuff, ybuff=ybuff

; position_in = position for entire plotting window
; ncols = n columns
; nrows = n rows

iposition=fltarr(ncols,nrows,4)

  ;Set buffer between panels
    if ~keyword_set(xbuff) then $
      xbuff=0.03
    if ~keyword_set(ybuff) then $
      ybuff=0.03

  ;Panel dimensions
    dx = ( (position_in[2]-position_in[0]) - xbuff*(ncols-1)) / ncols
    dy = ( (position_in[3]-position_in[1]) - ybuff*(nrows-1) ) / nrows

  for icol=0,ncols-1 do begin
    for irow=0,nrows-1 do begin
      iposition[icol,irow,0]=position_in[0] + (dx + xbuff) * icol
      iposition[icol,irow,1]=position_in[3] - (dy + ybuff) * irow - dy
      iposition[icol,irow,2]=iposition[icol,irow,0] + dx
      iposition[icol,irow,3]=iposition[icol,irow,1] + dy
    endfor
  endfor

return, iposition

end
