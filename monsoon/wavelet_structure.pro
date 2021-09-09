; 
; Use wavelet software based on WAVETEST from the following package to
; generate a structure containing all pertinent fields.
;
; See "http://paos.colorado.edu/research/wavelets/"
; Written January 1998 by C. Torrence
;
; James Ruppert
; 6/12/21
; 
function wavelet_structure,x,var,dx,pad,s0,dj,j1,mother,iavg,$ ; <-- inputs
  iavg2=iavg2 ; <-- inputs 
              ;(see WAVETEST and other wavelet routines for more info)

  nx=n_elements(x)

  recon_var=var

  ;estimate lag-1 autocorrelation, for red-noise significance tests
  lag1 = (A_CORRELATE(var,1) + SQRT(A_CORRELATE(var,2)))/2.

  ;Wavelet transform:

    wave = WAVELET(recon_var,dx,PERIOD=period,SCALE=scale,S0=s0, $
      PAD=1,COI=coi,DJ=dj,J=j1,MOTHER=mother,/RECON)

    power = (ABS(wave))^2  ; compute wavelet power spectrum
    global_ws = TOTAL(power,1)/nx   ; global wavelet spectrum (GWS)
    J = N_ELEMENTS(scale) - 1

  ;Significance levels:
    signif = WAVE_SIGNIF(var,dx,scale,0, $
      LAG1=lag1,SIGLVL=0.95,MOTHER=mother)
    signif = REBIN(TRANSPOSE(signif),nx,J+1)  ; expand signif --> (J+1)x(N) array
    signif = power/signif   ; where ratio > 1, power is significant
  
  ;GWS significance levels:
    dof = nx - scale   ; the -scale corrects for padding at edges
    global_signif = WAVE_SIGNIF(var,dx,scale,1, $
      LAG1=lag1,DOF=dof,MOTHER=mother,CDELTA=Cdelta,PSI0=psi0)
  
  ;check total variance (Parseval's theorem) [Eqn(14)]
    scale_avg = REBIN(TRANSPOSE(scale),nx,J+1)  ; expand scale-->(J+1)x(N) array
    power_norm = power/scale_avg
    variance = (MOMENT(var))[1]
    recon_variance = dj*dx/(Cdelta*nx)*TOTAL(power_norm)  ; [Eqn(14)]

;     IF (N_ELEMENTS(recon_var) GT 1) THEN BEGIN
;       recon_variance = (MOMENT(recon_var))[1]
;   ; RMS of Reconstruction [Eqn(11)]
;       rms_error = SQRT(TOTAL((var - recon_var)^2)/nx)
;       PRINT
;       PRINT,'        ******** RECONSTRUCTION ********'
;       PRINT,'original variance =',variance,' sigma^2'
;       PRINT,'reconstructed var =',FLOAT(recon_variance),' sigma^2'
;       PRINT,'Ratio = ',recon_variance/variance
;       PRINT,'root-mean-square error of reconstructed var = ',rms_error,' sigma'
;       PRINT
;       IF (mother EQ 'DOG') THEN BEGIN
;         PRINT,'Note: for better reconstruction with the DOG, you need'
;         PRINT,'      to use a very small s0.'
;       ENDIF
;       PRINT
;     ENDIF

  ;Scale-average over a set range of periods
    iavgin=iavg
    avg = WHERE((scale GE iavg[0]) AND (scale LT iavg[1]))
;    scale_avg = dj*dx/Cdelta*TOTAL(power_norm[*,avg],2)  ; [Eqn(24)]
scale_avg = dj*dx/Cdelta*TOTAL(power[*,avg],2)  ; [Eqn(24)]
    scaleavg_signif = WAVE_SIGNIF(var,dx,scale,2, $
      LAG1=lag1,SIGLVL=0.95,DOF=iavgin,MOTHER=mother)

;----RETURN STRUCTURE

  wavelet={wave:wave,recon_var:recon_var,period:period,scale:scale,power:power,global_ws:global_ws,signif:signif,$
           coi:coi,global_signif:global_signif,scale_avg:scale_avg,scaleavg_signif:scaleavg_signif}

  ;Scale-average 2
  if keyword_set(iavg2) then begin

    iavgin=iavg2
    avg = WHERE((scale GE iavg2[0]) AND (scale LT iavg2[1]))
    scale_avg2 = dj*dx/Cdelta*TOTAL(power_norm[*,avg],2)  ; [Eqn(24)]
;scale_avg2 = dj*dx/Cdelta*TOTAL(power[*,avg],2)  ; [Eqn(24)]
    scaleavg_signif2 = WAVE_SIGNIF(var,dx,scale,2, $
      LAG1=lag1,SIGLVL=0.95,DOF=iavgin,MOTHER=mother)

  wavelet=create_struct('scale_avg2',scale_avg2,wavelet)

  endif

  return,wavelet                                    

end
