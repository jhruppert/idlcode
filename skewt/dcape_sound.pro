;+
; NAME:
; DCAPE_SOUND
; PURPOSE:
; CAlculate dowdraft CAPE (DCAPE) for a vertical sounding of
;     pressure, temperature and mixing ratio.
;
; INPUT PARAMETERS:
;   P: pressures (hPa) at the N levels
;   T: Temperatures (K) at the N levels
;   R: Water vapor mixing ratio R (kg H2O per kg of dry air) at the N levels
;   Note: Relation between R and specific humidity S: R=S/(1-S)
;   
; OPTIONS:
;  MIN:  Set keyword to calculate DCAPE from level of minimum wet bulb
;  DLEV: Set to a vertical index to calculated DCAPE from.
;
; OUTPUT PARAMETERS:
;  Funtion returns DCAPE from a specified level (or from the level of minimum
;  wet bulb temperature; see below).
;  DCAPEP: DCAPE of each air parcel
;
; EXAMPLE:
;   PRES=[ 1000., 850., 700., 500.,  300.,  200.,  100.] ; pressures
; TEMP=[ 30.,  20.,  15., 0.,  -40.,  -60.,  -70.] ; temp (Celsius)
; DEWPT=[  25.,  18.,   14., -20.,  -60., -90., -110.] ; dew points (Celsius) 
; R=MIXR_SAT(dewpt,pres)/1000.
; or    R=[0.0107849,0.00911408,0.00731244,0.00156797,4.04480e-05,7.08816e-07,2.71149e-08]
; print,cape(pres,temp+273.15,R)
; HISTORY:
;   Based on the CAPE_SOUND function of Kerry Emanuel and Dominik Brunner. Written
;   by James Ruppert. ruppert@atmos.colostat.edu
;   4/25/12
;-

;USING FUNCTION ESA FROM CAPE_SOUND (NEED TO COMPILE THAT FIRST)

FUNCTION DCAPE_SOUND,P,T,R,DCAPEP=DCAPEP,TLP=TLP,TLVP=TLVP,TVPDIF=TVPDIF,BUF=BUF,layer=layer

;
;BUF: A BUFFER (IN HPA) INDICATING THE DEPTH (P+BUF TO P-BUF) TO AVERAGE
;     THE DOWNDRAFT PARCEL WETBULB
;

ON_ERROR,3

; check for input
IF n_params() LT 3 THEN BEGIN
   message,'Usage: Cape_sound,p,T,R'
   return,-1
ENDIF

N=n_elements(p) ; number of vertical levels
IF (n_elements(T) NE N) OR (n_elements(R) NE N) THEN BEGIN
   print,'Error: arrays p,T,R must have same size'
   return,-1
ENDIF

IF N EQ 0 THEN return,-1

if keyword_set(buf) then begin
  dz=mean(deriv(p[5:15]),/nan)
  ibuf=fix(buf)/fix(-1*dz)
  if ibuf eq 0 then message,'BUF is not divisible by DZ'
endif

itop=N-1

TLP=fltarr(N,N)
TLVP=fltarr(N,N)
TVPDIF=fltarr(N,N)
DCAPEP=fltarr(N)
NAP=fltarr(N)
PAP=fltarr(N)

;
;   ***   ASSIGN VALUES OF THERMODYNAMIC CONSTANTS     ***
;
CPD=1005.7  ; specific heat dry air
RV=461.5  ; gas constant water vapor = R*/M_H2O
RD=287.04 ; gas constant dry air = R*/M_dry
EPS=RD/RV

;
;  *** Water vapor pressure EV and saturation vapor pressure ES ***
;
;EV=R*P/(EPS+R)  ; vapor pressure
ES=ESA(T-273.15)  ; saturation vapor pressure
RS=EPS*ES/(P-ES)

;
;   ***   Find theta-wetbulb for downdraft   ***
;

if keyword_set(buf) then begin
  
  layer=[itop-2*ibuf,itop]
  thetaw=fltarr(layer[1]-layer[0]+1)
  for i=layer[0],layer[1] do begin
    if ~finite(p[i]) then begin
      thetaw[i]=!values.f_nan
      continue
    endif
    RH=R[I]/RS[I]  ; relative humidity
    RH=(RH LT 1.0)*RH+(RH GE 1.0)
    CHI=T[I]/(1669.0-122.0*RH-T[I])
    PLCL=(RH GT 0.0)*P[i]*RH^CHI+(RH LE 0.)
    TLCL=T[I]*(PLCL/P[I])^(RD/CPD)
    thetaw[i-layer[0]]=os(tlcl,plcl)
  endfor
  thetaw=mean(thetaw,/nan)
  
endif else begin
  
  i=itop
  RH=R[I]/RS[I]  ; relative humidity
  RH=(RH LT 1.0)*RH+(RH GE 1.0)
  CHI=T[I]/(1669.0-122.0*RH-T[I])
  PLCL=(RH GT 0.0)*P[i]*RH^CHI+(RH LE 0.)
  TLCL=T[I]*(PLCL/P[I])^(RD/CPD)
  thetaw=os(tlcl,plcl)
  
endelse

;   ***  Starting downdraft parcel loop   ***

FOR I=0,0 DO BEGIN
  
  FOR J=I,N-1 DO BEGIN
     
     TLP[I,J]=_TSA(THETAW,P[J])
     TLVP[I,J]=TLP[I,J]*(1.+RS[J]/EPS)/(1.+RS[J])
     TVPDIF[I,J]=TLVP[I,J]-T[J]*(1.+R[J]/EPS)/(1.+R[J])
     ;TVPDIF[I,J]=TLP[I,J]-T[J]
     
  ENDFOR
  
  FOR J=(I+1),itop DO BEGIN
      TVM=0.5*(TVPDIF[I,J]+TVPDIF[I,J-1])
      PM=0.5*(P[J]+P[J-1])
      IF TVM LE 0.0 THEN $
         NAP[I]=NAP[I]-RD*TVM*(P[J-1]-P[J])/PM $
      ELSE BEGIN
         PAP[I]=PAP[I]+RD*TVM*(P[J-1]-P[J])/PM
      ENDELSE
  ENDFOR
  DCAPEP[I]=NAP[I];PAP[I]-
  
ENDFOR  ; loop over air parcel origins

return,DCAPEP[0]

END