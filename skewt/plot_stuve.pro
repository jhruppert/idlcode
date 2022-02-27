;+
; NAME:
;	plot_stuve
;
; PURPOSE:
;    Stuve (pseudo-adiabat) Atmospheric Thermodynamic Diagram
;    The axes are temperature and -pressure^(R/Cp)
;
;    A temperature and dew point sounding may optionally be plotted
;    on the diagram.
;   
; CATEGORY:
;
; CALLING SEQUENCE:
;
;   plot_stuve,pres=pres,Temp=Temp,td=td,col_temp=col_temp,$
;	col_td=col_td,wspeed=wspeed,wangle=wangle,title=title,$
;	linestyle=linestyle,/over
;
; EXAMPLE:
;
; INPUTS:
;
; OPTIONAL INPUT PARAMETERS:
;  pres		Fltarr(nlev): pressures (hPa) of the sounding at nlev different levels
;  temp		Fltarr(nlev): sounding temperatures (deg C)
;  td	        Fltarr(nlev): sounding dew point temperatures (deg C)
;  col_temp	Byte: color index for temperture line
;  col_td	Byte: color index for dew point temperature line
;  linestyle	Integer: linestyle index for sounding lines (default=0=solid line)
;  wspeed	Fltarr(nlev): wind speed at nlev levels (m/s)
;  wangle	Fltarr(nlev): wind angle at nlev levels (0 deg = wind
;			from the north, clockwise)
;  title	String: Plot title (default is 'Stuve')
;
; KEYWORD INPUT PARAMETERS:
;  over		Set keyword to plot sounding over existing graph
;
; OUTPUTS:
;
; COMMON BLOCKS:
;
; SIDE EFFECTS:
;
; RESTRICTIONS:
;
; PROCEDURE:
;	
; MODIFICATION HISTORY:
;
; Program in IDL      October 1997
; Frank Evans  University of Colorado    evans@nit.colorado.edu
; Dominik Brunner, KNMI, Feb 2000, brunner@atmos.umnw.ethz.ch
;	converted to a subroutine.
;	The sounding can now be passed as paramters.
;-
PRO plot_stuve,pres=pres,Temp=Temp,td=td,col_temp=col_temp,$
	col_td=col_td,wspeed=wspeed,wangle=wangle,title=title,$
	linestyle=linestyle,over=over

; number of levels of vertical sounding
nlev=n_elements(pres)
; wind data supplied?
windok=(n_elements(wspeed) EQ nlev) AND (n_elements(wangle) EQ nlev) $
	AND (nlev GT 0)

IF n_elements(title) EQ 1 THEN plottitle=title ELSE plottitle='Stuve Diagram'

; Set character sizes for title and labels
titlesize=1.3
labsize=0.8
legendsize=1.0

;  Initialize constants that control where lines drawn and labeled
;if (atmospart eq 1) then begin
Tlo=-80.  &  Thi=40.
Plo=100.  &  Phi=1050.
thetastartlab=200.  &  thetaendlab=450.  &  thetasteplab=10.
Pstartlab=1000.  &  Pendlab=Plo       &  Psteplab=-100.
;endif
;if (atmospart eq 2) then begin	; lower troposphere only
;Tlo=-40.  &  Thi=40.
;Plo=500.  &  Phi=1050.
;thetastartlab=240.  &  thetaendlab=370.  &  thetasteplab=10.
;  Pstartlab=1000.  &  Pendlab=Plo       &  Psteplab=-100.
;endif

Pylo=Plo^0.286   &  Pyhi=Phi^0.286
Tstart=Tlo       &  Tend=Thi          &  Tstep=10.
Tstartlab=Tlo    &  Tendlab=Thi       &  Tsteplab=Tstep
Pstart=Phi       &  Pend=Plo          &  Pstep=-50.
thetastart=190.  &  thetaend=450.     &  thetastep=10.
thetaestart=190. &  thetaeend=400.    &  thetaestep=10.

;  Set up plot area
IF windok THEN xhi=0.89 ELSE xhi=0.95
xlo=0.15  &  ylo=0.15 &  yhi=0.95

; overlay sounding over exising graph?
IF keyword_set(over) THEN goto,soundingonly

plot, [0], [0], position=[xlo,ylo,xhi,yhi], /nodata, ticklen=0.0, $
    xrange=[Tlo,Thi], yrange=[Pyhi,Pylo], /xstyle, /ystyle, $
    xticks=1, xtickname=[' ',' '], yticks=1, ytickname=['    ','    '], $
    xtitle='Temperature (C)', ytitle='Pressure (hPa)',thick=1.5
XYOUTS,(xlo+xhi)/2.,yhi+0.02,plottitle, charsize=titlesize,$
	alignment=0.5,/normal

;  Make the temperature lines and label them
for T = Tstart, Tend, Tstep do begin
    plots, [T,T], [Pylo,Pyhi],  thick=1.5
endfor
Pylabel = Pyhi + (Pyhi-Pylo)*0.025
for T = Tstartlab, Tendlab, Tsteplab do begin
    xyouts, T,Pylabel, alignment=0.5, orient=0, /noclip, $
    charsize=1.5*labsize, charthick=1.5, string(T,format='(I3)')
endfor

;  Make the pressure lines and label them
for P = Pstart, Pend, Pstep do begin
    Py=P^0.286
    plots, [Tlo,Thi], [Py,Py], thick=1.5
endfor
for P = Pstartlab, Pendlab, Psteplab do begin
    Py = P^0.286 + (Pyhi-Pylo)*0.005
    xyouts, Tlo, Py, alignment=1.1, orient=0, /noclip,  $
    charsize=1.5*labsize, charthick=1.5, string(P,format='(I4)')
endfor

; Make the theta curves - dry adiabats
; Use theta = T*(1000/p)^0.286
thetaslopefac=-26.0*(yhi-ylo)*(Thi-Tlo)*1000^0.286 / (20.*(xhi-xlo)*(Pyhi-Pylo))
Tarr=Tlo+(Thi-Tlo)*findgen(100)/99.0
for theta = thetastart, thetaend, thetastep do begin
    Parr = 1000.* ((Tarr+273.15)/theta)^(1./0.286)
    Py = Parr^0.286
    plots, Tarr, Py, noclip=0, linestyle=0
endfor
for theta = thetastartlab, thetaendlab, thetasteplab do begin
    T1 = theta*(Plo/1000.)^0.286 - 273.15
    P1 = 1000.* ( (Tlo+273.15)/theta )^(1./0.286)
    if P1 lt Plo then T=T1 else T=Tlo
    P = 1000.* ( (T+273.15)/theta )^(1./0.286)
    Py = P^0.286
    angle = (180./!PI)*atan(thetaslopefac/theta)
    xyouts, T, Py, alignment=-0.1, orient=angle,  $
    charsize=1.2*labsize, charthick=1.5, string(theta,format='(I3)')
endfor

;  Make the saturated specific humidity curves
;    q_s = 0.622*e_s/p
qs=[30,20,15,10,7,5,3,2,1.5,1.0,0.7,0.5,0.3,0.2,0.1]
qslab=['30','20','15','10','7','5','3','2','1.5','1.0','0.7','0.5','0.3','0.2','0.1']
Pylabel = Pyhi - (Pyhi-Pylo)*0.007
t0=273.15
Tarr=-80+findgen(120)
for i = 0, n_elements(qs)-1 do begin
    Tk = 273.15 + Tarr
    es = exp( -6763.6/Tk - 4.9283*alog(Tk) + 54.2190 )
    Pyarr = ( 0.622*es/(0.001*qs(i)) )^0.286
    plots, Tarr, Pyarr, noclip=0, linestyle=1
    es = 0.001*qs(i)*Phi/0.622
    T = T0/(1.-T0*alog(es/6.11)/5417) - 273.15
    if T gt Tlo+10 then $
      xyouts, T, Pylabel, alignment=0.5, orient=0,  $
    charsize=1.0*labsize, charthick=1.5, qslab(i)
endfor
xyouts, Tlo, Pylabel, alignment=-0.2, charsize=1.2*labsize, charthick=1.5, $
    'q!Ds!N (g/kg)'

;  Make the saturated adiabat curves
;    Use dtheta = -L*theta/(Cp*T) dqs and integrate from T=-60 C
nstep=100
Tarr=fltarr(nstep)  &  Pyarr=fltarr(nstep)
for theta = thetaestart, thetaeend, thetaestep do begin
   thetap = theta
   T = -60
   Tk = T + 273.15
   P = 1000.*(Tk/thetap)^(1./0.286)
   es = exp( -6763.6/Tk - 4.9283*alog(Tk) + 54.2190 )
   qs0 = 0.622*es/P
   Tarr(0) = T
   Pyarr(0) = P^0.286
   for i = 1, nstep-1 do begin
     T = T + 1.0
     Tk = T + 273.15
     es = exp( -6763.6/Tk - 4.9283*alog(Tk) + 54.2190 )
     qs = 0.622*es/P
     thetap = thetap - (2.5E6/1004)* (thetap/Tk) * (qs-qs0)
     P = 1000.*(Tk/thetap)^(1./0.286)
     Tarr(i) = T
     Pyarr(i) = P^0.286
     qs0 = qs
   endfor
   plots, Tarr, Pyarr, noclip=0, linestyle=2
endfor

plots, [0.25,0.30], [0.050,0.050], /normal, linestyle=0
xyouts, 0.31, 0.045, alignment=0.0,  '!5dry adiabats !9q!5', $
	charsize=1.2*legendsize,charthick=1.5, /normal
plots, [0.55,0.60], [0.050,0.050], /normal, linestyle=2
xyouts, 0.61, 0.045, alignment=0.0,  '!5saturated adiabats', $
	charsize=1.2*legendsize,charthick=1.5, /normal
plots, [0.55,0.60], [0.030,0.030], /normal, linestyle=1
xyouts, 0.61, 0.025, alignment=0.0,  '!5sat. specific humidity q!Ds!N', $
	charsize=1.2*legendsize,charthick=1.5, /normal

soundingonly:

; plot the optional sounding
IF nlev GT 0 THEN BEGIN
  IF n_elements(linestyle) NE 1 THEN linestyle=0. ; solid line is default
  Pyarr=pres^0.286
  IF n_elements(temp) EQ nlev THEN BEGIN
     ;   Plot the temperature profile
     IF n_elements(col_temp) NE 0 THEN $
        plots, temp, Pyarr, noclip=0, linestyle=linestyle, thick=3.,color=col_temp $
     ELSE plots, temp, Pyarr, noclip=0, linestyle=linestyle, thick=3.
  ENDIF
  IF n_elements(td) EQ nlev THEN BEGIN
     ;   Plot the dewpoint temperature profile
     IF n_elements(col_td) NE 0 THEN $
        plots, td, Pyarr, noclip=0, linestyle=linestyle, thick=3.,color=col_td $
     ELSE plots, td, Pyarr, noclip=0, linestyle=linestyle, thick=3.
  ENDIF
ENDIF

; plot optional wind barbs
IF windok THEN BEGIN
   pyarr=pres^0.286
   p=convert_coord(Thi,pyarr,/to_normal)
   index=WHERE(p[1,*] LE yhi)	; only wind barbs within plot range
   FOR i=0,n_elements(index)-1 DO wind_barb,wspeed[index[i]],size=1.7,$
	wangle[index[i]],p[0,index[i]],p[1,index[i]],/nocircle
ENDIF


END
