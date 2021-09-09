; 
; Calculate the diurnally resolved regressions of a given input variable yvar(x,y,hr,day),
;   given 4 daily-mean input rainfall indices.
;
; James Ruppert
; 12/31/20
; 
pro calc_dc_regress, yvar, $
  irain_all, irain_coast, irain_offshore, irain_onshore, irain_smd, irain_cdiff, $ ; input daily-mean rainfall indices
  a1_all=a1_all, a1_coast=a1_coast, a1_offshore=a1_offshore, a1_smd=a1_smd, a1_cdiff=a1_cdiff
;, a1_onshore=a1_onshore ; return regression coefficients

  ;DIMENSIONS
    dims=size(yvar,/dimen) ; [ x, y, hour, day ]
    nx=dims[0]
    ny=dims[1]
    npd=dims[2]
    nd=dims[3]
;    nt=dims[2]
;    nd=nt/npd

  ;REGRESS RAIN

  ;TIME-MEAN AS A FUNCTION OF HOUR
;  ymean=mean(yvar,dimension=4,/nan,/double) ; (sluggish job)
  ;REMOVE DAILY MEAN
  ymean=mean(yvar,dimension=3,/nan,/double) ; (sluggish job) - returns ymean(x,y,day)
  dims=[nx,ny,nd,npd]
  yp_all = yvar - transpose( rebin(ymean,dims), [0,1,3,2] )

  ;prep x-terms once (fast job)
  dims=[nd,nx,ny]
  xp_all      = transpose(rebin(irain_all,     dims),[1,2,0])
  xp_coast    = transpose(rebin(irain_coast,   dims),[1,2,0])
  xp_offshore = transpose(rebin(irain_offshore,dims),[1,2,0])
;  xp_onshore  = transpose(rebin(irain_onshore, dims),[1,2,0])
  xp_smd      = transpose(rebin(irain_smd,     dims),[1,2,0])
  xp_cdiff    = transpose(rebin(irain_cdiff,   dims),[1,2,0])

  ;regressed vars
  a1_all=fltarr(nx,ny,npd)
  a1_coast=a1_all
  a1_offshore=a1_all
;  a1_onshore=a1_all
  a1_smd=a1_all
  a1_cdiff=a1_all

  ;for indexing by hour
;  indhr=lindgen(nd)*npd

  for ih=0,npd-1 do begin
;for ih=0,5 do begin

    print,'  ih =',ih,' of',npd-1

;    yp=reform(yvar[*,*,ih,*])
;    ymean=mean(yp,dimension=3,/nan,/double)
;    dims=[nx,ny,nd]
;    yp = temporary(yp) - rebin(reform(ymean[*,*,ih]),dims)
    ;Standardize
;    yp = temporary(yp) / rebin(reform(stdev[*,*,ih]),dims)

    ;FOR YP WITH DAILY-MEAN REMOVED
    yp=reform(yp_all[*,*,ih,*])

    a1_all[*,*,ih]      = total((xp_all      * yp),3,/nan,/double)/nd
    a1_coast[*,*,ih]    = total((xp_coast    * yp),3,/nan,/double)/nd
    a1_offshore[*,*,ih] = total((xp_offshore * yp),3,/nan,/double)/nd
;    a1_onshore[*,*,ih]  = total((xp_onshore  * yp),3,/nan,/double)/nd
    a1_smd[*,*,ih]      = total((xp_smd      * yp),3,/nan,/double)/nd
    a1_cdiff[*,*,ih]    = total((xp_cdiff    * yp),3,/nan,/double)/nd

  endfor


end
