;========================================================================
;  ROUTINE TO CALCULATE PARCEL TEMPERATURE, CAPE, CIN, AND LFC 
;  
;  Rob Seigel
;
PRO PARCELCALC, t, td, p, cape, tp=tp, cin=cin, lcl=lcl, lfc=lfc, el=el, $
                virtual=vir

  ON_Error,2

  ;  Find number of data levels
  szd     = size(t)
  nlevels = szd[1]

  ;begin calcs
  cape=0.
  cin=0.
  lcl=0.
  lfc=0.
  wp=mixr_sat(td,p)
  tp=t
  tv=((tp+273.15)*(1.+0.61*(wp/1000.)))-273.15
  tvp=tv
  wsp=MIXR_SAT(tp,p)
  tpot=(t[0]+273.15)*((1000./p[0])^(287.05/1004.))
  lcl=0.
  for k=1,nlevels-1 do begin
     tp[k]=(tpot*((p[k]/1000.)^(287.05/1004.)))-273.15
     wp[k]=wp[k-1]
     tvp[k]=((tp[k]+273.15)*(1.+0.61*(wp[k]/1000.)))-273.15
     wsp[k]=MIXR_SAT(tp[k],p[k])
     if (wsp[k] LT 1.e-30) then wsp[k]=0.
     if (wp[k] ge wsp[k]) then begin
        if (lcl eq 0.) then begin
           ;linearly interpolate between layers to find LCL
           lowvap=wp[k-1]-wsp[k-1]
           highvap=wp[k]-wsp[k]
           vapInterp=[lowvap,highvap]
           pInterp=[p[k-1],p[k]]
           lclInterp=INTERPOL(pInterp,vapInterp,0)
           lcl=lclInterp
           lclInd=k

           ; interpolate env and parcel vars to LCL
           ; for more accurate parcel calcs
           p[k]=lcl
           tv[k]=INTERPOL([tv[k-1],tv[k]],vapInterp,0)
           t[k]=INTERPOL([t[k-1],t[k]],vapInterp,0)
           tp[k]=INTERPOL([tp[k-1],tp[k]],vapInterp,0)
           td[k]=INTERPOL([td[k-1],td[k]],vapInterp,0)

           saT=OS(tp[k],lcl)
        endif
        tp[k]=TSA(saT,p[k])-273.15
        wp[k]=MIXR_SAT(tp[k],p[k])
        wsp[k]=MIXR_SAT(tp[k],p[k])
        tvp[k]=((tp[k]+273.15)*(1.+0.61*(wsp[k]/1000.))-273.15)
        tvdiff=tvp[k]-tv[k]
        tdiff=(tp[k]+273.15)-(t[k]+273.15)

        ;calculate cape and/or cin
        if (keyword_set(vir)) then begin
           ;VIRTUAL ADJUSTMENT
           if (tvdiff ge 0. and lcl gt 0) then begin ;CAPE
              if (cape eq 0) then lfc=p[k]
              cape=cape+287.*tvdiff*alog(p[k-1]/p[k])
              
              el=p[k]
           endif
           if (tdiff lt 0. and (cape eq 0. or p[k] gt 500)) then begin ;CIN
              cin=cin+287.*tvdiff*alog(p[k-1]/p[k])
           endif
        endif else begin
           ;NON-VIRTUAL CALCS
           if (tdiff ge 0.) then begin ;CAPE
              if (cape eq 0) then lfc=p[k]
              cape=cape+287.*tdiff*alog(p[k-1]/p[k])
              
              el=p[k]
           endif
           if (tdiff lt 0. and (cape eq 0. or p[k] gt 500)) then begin ;CIN
              cin=cin+287.*tdiff*alog(p[k-1]/p[k])
           endif
        endelse
     endif
  endfor

  ;reset temp profile and parcel profile to virtual if asked
  if (keyword_set(vir)) then begin
     t=tv
     tp=tvp
  endif

  ;if CIN is less than 1 from interpolation errors, then make it 0
  if (abs(cin) lt 1) then cin=0

  ;return calculations as string
  cape=strmid(strtrim(string(cape),1),0,5)
  cin=strmid(strtrim(string(cin),1),0,5)
  lcl=strmid(strtrim(string(lcl),1),0,5)
  lfc=strmid(strtrim(string(lfc),1),0,5)
  if n_elements(el) ne 0 then el=strmid(strtrim(string(el),1),0,5) else el='nan'

END