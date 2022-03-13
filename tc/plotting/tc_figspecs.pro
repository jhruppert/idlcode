; 
; Produce variable details for WRF TC simulatoins.
; 
pro tc_figspecs, var_str, figspecs, setmax=setmax, setmin=setmin, setscale=setscale, tag2=tag2, setndivs=setndivs,$
  idcomp=idcomp, set_cint=set_cint

scale=1.
if ~keyword_set(idcomp) then idcomp=0

;VAR = THETA-E
  if var_str eq 'th_e' then begin

    cbar_format='(i3)'
    cbar_tag='EPT [ K ]'
    title='Equiv Pot Temp'

    scale=1.

    col_table=66
    icbar=1

    ndivs=6
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=15;5
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
;    colors=reverse(colors)

    max=380.
    if keyword_set(setmax) then max=setmax
    min=320
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = PW
  endif else if var_str eq 'pw' then begin

    cbar_format='(i2)'
    cbar_tag='[ mm ]'
    title='PW'

    scale=1.

    col_table=66;69
    icbar=1

    ndivs=4
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=15;71;15
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
;    colors=reverse(colors)

    max=80
    if keyword_set(setmax) then max=setmax
    min=30
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = OLR
  endif else if var_str eq 'olr' then begin

    cbar_format='(i3)'
    cbar_tag='[ W m!U-2!N ]'
    title='OLR'

    scale=1.

    col_table=76 ; custom bar
    icbar=1

    ndivs=4
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=71;255;9;15;5
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    colors=reverse(colors)

    max=320
    if keyword_set(setmax) then max=setmax
    min=80
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = lw
  endif else if var_str eq 'lw' then begin

    cbar_format='(i4)'
    cbar_tag='[ W m!U-2!N ]'
    title='<LW>'
    icbar=1

    scale=1.

    col_table=71
    irev=1

    ndivs=4
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=71
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    if irev then colors=reverse(colors)

    max=150
    if keyword_set(setmax) then max=setmax
    min=-1.*max;-300
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = strat
  endif else if var_str eq 'strat' then begin

    cbar_format='(i4)'
    cbar_tag='[ 1-C, 2-S, 3-o ]'
    title='Morph'
    icbar=1

    scale=1.

    col_table=11

    ndivs=2

    ncols=3
    colors=[254,254/2,254./3]

    levels=indgen(ncols)+1

;VAR = TQC
  endif else if var_str eq 'tqc' then begin

    cbar_format='(f3.1)'
    cbar_tag='[ mm ]'
    title='TQC'
    icbar=1

    scale=1.

    col_table=77;78 ; 77-light (based on 11), 78-dark (rainbow)
    irev=0

    ndivs=1
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=255;71;41
    colors=findgen(ncols)/(ncols-1)*255
;    colors=findgen(ncols)/(ncols-1)*255/2+255/2
    if irev then colors=reverse(colors)

    max=0.6
    if keyword_set(setmax) then max=setmax

    if col_table eq 77 then $
      min=1e-3 $ ; for cbar 77 (light)
    else $
      min=0. ; for cbar 78 (dark)

    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = NGMS
  endif else if var_str eq 'ngms' then begin

    cbar_format='(f4.1)'
    cbar_tag='[ - ]'
    title='NGMS'
    icbar=1

    scale=1.

    col_table=71
    irev=0

    ndivs=1
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=15
    colors=findgen(ncols)/(ncols-1)*255
    if irev then colors=reverse(colors)

    max=2
    if keyword_set(setmax) then max=setmax
    min=-1.*max
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = MADV
  endif else if var_str eq 'madv' then begin

    cbar_format='(i3)'
    cbar_tag='[ g/kg / hr ]'
    title='MADV'

    scale=1e3*3600

    col_table=66;69
    icbar=1

    ndivs=4
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=15;71;15
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
;    colors=reverse(colors)

    max=20;80
    if keyword_set(setmax) then max=setmax
    min=-1.*max
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = AVOR
  endif else if var_str eq 'AVOR' or var_str eq 'avor' then begin

    cbar_format='(i3)'
    cbar_tag='AVOR [ 10!U-5!N s!U-1!N ]'
    title='AVOR'

    scale=1.

    col_table=16;13;66;16
    icbar=1

    ndivs=6
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=21;15;5
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
;    colors=reverse(colors)

    max=60;20
    if keyword_set(setmax) then max=setmax
    min=0;-1.*max;20
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = THETA-E
  endif else if var_str eq 'theta_e' then begin

    cbar_format='(i3)'
    cbar_tag='Theta-e [ K ]'
    title='theta_e'

    scale=1.

    col_table=66;16
    icbar=1

    ndivs=6
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=15;5
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
;    colors=reverse(colors)

    max=360
    if keyword_set(setmax) then max=setmax
    min=320
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = SST
  endif else if var_str eq 'SST' then begin

    cbar_format='(i2)'
    cbar_tag='SST [ C ]'
    title='SST'

    scale=1.

    col_table=39
    icbar=1

    ndivs=6
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=15
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
;    colors=reverse(colors)

    max=20
    if keyword_set(setmax) then max=setmax
    min=-20
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = H_DIABATIC
  endif else if var_str eq 'H_DIABATIC' then begin

    cbar_format='(f5.1)'
cbar_format='(i3)'
    cbar_tag='[ mm hr!U-2!N ]'
;    cbar_tag='[ K d!U-1!N ]'
    title='HDIA'

    scale=1.;3600.*24

    col_table=70
    icbar=1

    ndivs=5
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=51;15;5
;ncols=45
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    colors=reverse(colors)

    max=5.
    if keyword_set(setmax) then max=setmax
    min=0.;-1.*max
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = Radial wind
  endif else if var_str eq 'u_rad' then begin

    cbar_format='(i3)'
    cbar_tag='[ m s!U-1!N ]'
    title='URad'

    scale=1.

    col_table=69
    icbar=1
    irev=0

    ndivs=4
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=13;15;5
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    if irev then colors=reverse(colors)

    max=6.
    if keyword_set(setmax) then max=setmax
    min=-1.*max
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = Tangential wind
  endif else if var_str eq 'v_tan' then begin

    cbar_format='(i3)'
    cbar_tag='[ m s!U-1!N ]'
    title='Vtan'
    icbar=1

    scale=1.

    col_table=71;75;11;71;70
    irev=1

    ndivs=6
;ndivs=8
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=15;25;15;5
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    if irev then colors=reverse(colors)

    max=30.
max=15.
    if keyword_set(setmax) then max=setmax
    min=-1.*max
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = SEF surface enthalpy flux
  endif else if var_str eq 'SEF' then begin

    cbar_format='(i4)'
    cbar_tag='[ W m!U-2!N ]'
    title='SEF'

    scale=1.

    col_table=75;11;71;70
    irev=0
    icbar=1

    ndivs=6
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=31;15;5
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    if irev then colors=reverse(colors)

    max=30.
    if keyword_set(setmax) then max=setmax
    min=0
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = LAPSE RATE
  endif else if var_str eq 'LR' then begin

    cbar_format='(i3)'
    cbar_tag='LR [ 10!U-2!N K hPa!U-1!N ]'
    title='LR'

    scale=1.

    col_table=70
    icbar=1

    ndivs=6
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=15;5
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    colors=reverse(colors)

    max=10.
    if keyword_set(setmax) then max=setmax
    min=4;-1.*max
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = LW or SW
  endif else if strmatch(var_str,'*RTHRA*') then begin

;    cbar_format='(f4.1)'
;cbar_format='(f5.1)'
    cbar_format='(i4)'
    cbar_tag='LW [ K d!U-1!N ]'
;    cbar_tag='[ K d!U-1!N ]'
    title='Radheat'

    scale=3600d*24

    col_table=70
    icbar=1

    ndivs=4;6
    if keyword_set(setndivs) then ndivs=setndivs

;FOR AZIM HOVMOLLER OF LW
;  cbar_tag='[ W m!U-2!N ]'
;  scale=1.
;  ndivs=4
;  cbar_format='(i4)'

    ncols=15;21;15;5
;ncols=45
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    colors=reverse(colors)

    max=5
max=3
    if keyword_set(setmax) then max=setmax
    min=-1.*max;5
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = RH
  endif else if var_str eq 'RH' then begin

    cbar_format='(i3)'
    cbar_tag='RH [ % ]'
    title='RH'

    scale=1.

    col_table=66
    icbar=1

    ndivs=7
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=15;5
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
;    colors=reverse(colors)

    max=100
    if keyword_set(setmax) then max=setmax
    min=30
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = PSFC
  endif else if var_str eq 'slp' or var_str eq 'SLP' then begin

    cbar_format='(i4)'
    cbar_tag='SLP [ hPa ]'
    title='SLP'

    scale=1.;e-2

    col_table=66
    icbar=1

    ndivs=6
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=15;5
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    colors=reverse(colors)

    max=120
    if keyword_set(setmax) then max=setmax
    min=0
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = Z - Geopotential height
  endif else if var_str eq 'z' or var_str eq 'Z' then begin

    cbar_format='(i4)'
    cbar_tag='Height [ m ]'
    title='Height'

    scale=1.;e-2

    col_table=66
    icbar=1

    ndivs=6
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=15;5
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
;    colors=reverse(colors)

    max=120
    if keyword_set(setmax) then max=setmax
    min=0
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = Latent Heat
  endif else if var_str eq 'lh' then begin

    cbar_format='(i3)'
    cbar_tag='LH [ W d!U-2!N ]'
    title='LH'

    scale=1.

    col_table=11
    icbar=1

    ndivs=4
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=15;5
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
;    colors=reverse(colors)

    max=5
    if keyword_set(setmax) then max=setmax
    min=0
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = RAIN
  endif else if var_str eq 'rainrate' or var_str eq 'trmm_rain' then begin

    cbar_format='(i3)'
    cbar_tag='[ mm hr!U-1!N ]'
    title='Rainfall Rate'

    scale=1.; keep mm/d     /24 ; mm/d --> mm/h
    if var_str eq 'trmm_rain' then scale=24. ; --> mm/d     1. ; already in mm/h

    col_table=75;0;70
    icbar=1

    ndivs=9
    if keyword_set(setndivs) then ndivs=setndivs

    levels=[10,25,50,75,100,150,300]
    ncols=n_elements(levels)
    ndivs=ncols-1

;    ncols=10;15;5
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
;    colors=reverse(colors)

    colors[0:1]=[5,30]
;    colors[1]=20

    max=70.
    if keyword_set(setmax) then max=setmax
    min=10;0.5
    if keyword_set(setmin) then min=setmin
;    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = RAINNC
  endif else if var_str eq 'RAINNC' then begin

    cbar_format='(i3)'
    cbar_tag='Accum. Rain [ cm ]'
    title='Accumulated Rain'

    scale=0.1 ; mm --> cm

    col_table=75;0;70
    icbar=1

    ndivs=4
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=15;5
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
;    colors=reverse(colors)

    max=500.
    if keyword_set(setmax) then max=setmax
    min=0
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = CRH - column rel hum
  endif else if var_str eq 'CWV' or var_str eq 'CRH' then begin

    cbar_format='(i3)'
    cbar_tag='[ % ]'
    title='CRH'

    scale=1e2

    col_table=66
    icbar=1
    irev=0

    ndivs=4
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=21;15;21;5
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    if irev then colors=reverse(colors)

    max=70.
    if keyword_set(setmax) then max=setmax
    min=10;0.5
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = TEMP
  endif else if var_str eq 'T'  then begin

    cbar_format='(f4.1)'
    cbar_tag='[ K ]'
    title="T'"
;    title="Thv'"
    icbar=1

    scale=1.

    col_table=75;11
    irev=0

    ndivs=4
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=25;15
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    if irev then colors=reverse(colors)

    max=1.
    if keyword_set(setmax) then max=setmax
    min=-1.*max
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = U-WIND
  endif else if var_str eq 'U' or var_str eq 'V' then begin

    cbar_format='(f4.1)'
    cbar_tag='[ m s!U-1!N ]'
    title=var_str+"'"

    scale=1.
    if keyword_Set(setscale) then scale=setscale

    if scale eq 10 then cbar_tag='[ 10!U-1!N m s!U-1!N ]'

    col_table=70;11
    icbar=1

    ndivs=4
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=21
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    colors=reverse(colors)

    max=10.
    if keyword_set(setmax) then max=setmax
    min=-1.*max
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = W-WIND
  endif else if var_str eq 'W' then begin

    cbar_format='(f5.1)'
    cbar_tag='[ cm s!U-1!N ]'
;cbar_tag='[ 10!U2!N kg m!U-1!N s!U-1!N ]'
;    cbar_tag='[ hPa h!U-1!N ]'
    title='W'
;title='rhoW'
    icbar=1

    scale=1;e-2;1e2
    if keyword_Set(setscale) then scale=setscale

    col_table=71
    irev=1

    ndivs=4
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=15;13;9;15
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    if irev then colors=reverse(colors)

    max=12.
    if keyword_set(setmax) then max=setmax
    min=-1.*max
    if keyword_set(setmin) then min=setmin
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

    ;EXPONENTIAL
      cbar_tag='[ kg m!U-1!N s!U-1!N ]'
      cbar_format='(i5)'
      ;levels=10^(findgen(6)-2)
      levels=[10,50,100,500,1e3,5*1e3]
      levels=[-1.*reverse(levels),0,levels]
      ncols=n_elements(levels)
      colors=findgen(ncols)/(ncols-1)*255;/2+255/2
      if irev then colors=reverse(colors)

;VAR = MFZ (Radial mass flux)
  endif else if var_str eq 'mfz' then begin

    cbar_format='(i3)'
    cbar_tag='[ 10!U5!N kg m!U-1!N s!U-1!N ]'
    title='Radial Mass Flux'

    scale=1.e-5
    if keyword_Set(setscale) then scale=setscale

    col_table=70
    icbar=1

    ndivs=6
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=19;15
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    colors=reverse(colors)

    max=15.
    if keyword_set(setmax) then max=setmax
    min=-1.*max
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = MFR (Vertical mass flux)
  endif else if var_str eq 'mfr' then begin

    cbar_format='(i3)'
    cbar_tag='[ 10!U3!N kg m!U-1!N s!U-1!N ]'
    title='Vertical Mass Flux'

    scale=1.e-3
    if keyword_Set(setscale) then scale=setscale

    col_table=70
    icbar=1

    ndivs=4
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=19;15
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
    colors=reverse(colors)

    max=10.
    if keyword_set(setmax) then max=setmax
    min=-1.*max
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = QVAPOR
  endif else if var_str eq 'QVAPOR' then begin

    cbar_format='(f4.1)'
    cbar_tag='[ % ]'
    title=var_str+"'"

    scale=1.
    if keyword_Set(setscale) then scale=setscale

    col_table=66
    icbar=1

    ndivs=4
    if keyword_set(setndivs) then ndivs=setndivs

    ncols=19;15
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2

    max=10.
    if keyword_set(setmax) then max=setmax
    min=0
    levels=findgen(ncols)/(ncols-1)*(max-min)+min

;VAR = 10m WIND
  endif else if var_str eq 'U10' or var_str eq 'V10' or var_str eq 'wspd' then begin

;    cbar_format='(f4.1)'
    cbar_format='(i2)'
    cbar_tag='[ m/s ]'
;    title=var_str+"'"
    title='10m Wind Speed'
    icbar=1

    scale=1.
    if keyword_Set(setscale) then scale=setscale

    if scale eq 10 then cbar_tag='[ 10!U-1!N m s!U-1!N ]'

    col_table=75;10;75;11;70
    irev=0

;    ndivs=4
    if keyword_set(setndivs) then ndivs=setndivs

;    levels=[18,33,43,50,58,70]
    max=60;58
    min=18
    cint=3.
    nlevs=(max-min)/cint+1
    levels=indgen(nlevs)*cint+min
;    nlevs=255
;    levels=findgen(nlevs)/(nlevs-1)*(max-min)+min
;    ndivs=nlevs-1
    ndivs=3

    ncols=n_elements(levels);21
    colors=findgen(ncols)/(ncols-1)*255;/2+255/2
;colors=reverse(colors)
;colors[0:round(255./2)]=reverse(colors[0:round(255./2)])
;    colors=findgen(ncols)/(ncols-1)*230

    if irev then colors=reverse(colors)

;    max=10.
;    if keyword_set(setmax) then max=setmax
;    min=-1.*max
;    if keyword_set(setmin) then min=setmin
;    levels=findgen(ncols)/(ncols-1)*(max-min)+min

  endif

  ;CONTOUR VARIABLE
    col_table_c=11
    colors_c=200;200
    cint=5;3;0.5
    if keyword_set(set_cint) then cint=set_cint
    clevs=(findgen(50)+1)*cint
clevs=[15,30]

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
