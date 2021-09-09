; 
; Produce variable details for WRF TC simulatoins.
; 
pro gfs_figspecs, var_str, figspecs, setmax=setmax, setmin=setmin, setscale=setscale, tag2=tag2, setndivs=setndivs,$
  idcomp=idcomp, set_cint=set_cint

scale=1.
if ~keyword_set(idcomp) then idcomp=0

;VAR = PMSL
  if var_str eq 'PMSL' then begin

    cbar_format='(i4)'
    cbar_tag='PMSL [ hPa ]'
    title='PMSL'

    scale=1e-2

    col_table=11
    icbar=1

    ndivs=6
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=50;15
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    colors=reverse(colors)

    max=1030.
    if keyword_set(setmax) then max=setmax
    min=980
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = Rel Vorticity
  endif else if var_str eq 'VOR' then begin

    cbar_format='(f5.1)'
    cbar_tag='Rel. Vor. [ 10!U-5!N s!U-1!N ]'
    title='VOR'

    scale=1e5

    col_table=71
    icbar=1

    ndivs=6
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=50;15
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    colors=reverse(colors)

    max=10.
    if keyword_set(setmax) then max=setmax
    min=-1.*max
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = PW
  endif; else if var_str eq 'pw' then begin

  ;CONTOUR VARIABLE
    col_table_c=11
    colors_c=200;200
    cint=5;3;0.5
    if keyword_set(set_cint) then cint=set_cint
    clevs=(findgen(50)+1)*cint

;  clevs=2^(indgen(9))
;;  cint=2;0.5
;;  clevs=(findgen(50)+1)*cint
;  clevs=[1,2,indgen(10)*5+5];5,10,20,30]
;  cint=5;0.5
;  clevs=(findgen(50)+1)*cint
  clevs2=clevs

  cint=0.5
  clevs3=(findgen(50)+1)*cint


  ;VARIABLE 2
;THESE WILL BE OVERWRITTEN BY HOV CALLING PROGRAM
    col_table2=70
    ncols=11
    maxv=10
    levels2=findgen(ncols)/(ncols-1)*2*maxv-maxv
    colors2=findgen(ncols)/(ncols-1)*255
    colors2=reverse(colors2)
    ndivs2=4
    cbar_tag2='U [ m s!U-1!N ]'
    cbar_format2='(i3)'


  figspecs={cbar_format:cbar_format,cbar_tag:cbar_tag,title:title,$
    col_table:col_table,icbar:icbar,colors:colors,levels:levels,levels2:levels2,colors2:colors2,$
    clevs:clevs,clevs2:clevs2,col_table2:col_table2,ndivs2:ndivs2,cbar_tag2:cbar_tag2,cbar_format2:cbar_format2,$
    clevs3:clevs3,ndivs:ndivs,$
    scale:scale,col_table_c:col_table_c,colors_c:colors_c}

end
