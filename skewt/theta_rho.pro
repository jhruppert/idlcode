function theta_rho,theta,qv,ql,qi

; FUNCTION TO CALCULATE DENSITY POTENTIAL TEMPERATURE,
; WHICH IS CONVENIENT FOR CALCULATING BUYANCY (BASED ON CM1)
; 
; INPUTS:
; 
;   THETA:   POTENTIAL TEMPERATURE (K)
;   QV:      VAPOR MIXING RATIO (KG/KG)
;   QL:      LIQUID WATER MIXING RATIO (I.E., QC+QR; KG/KG)
;   QI:      ICE MIXING RATIO (KG/KG)
; 
; RETURNS:
; 
;   DENSITY POTENTIAL TEMPERATURE (K)
; 
; James Ruppert, ruppert@atmos.colostate.edu
; 12/10/12
; 


eps = 287.045/461.53 ; Rd/Rv ~ 0.622

theta_r = theta * ( 1 + qv/eps ) / ( 1 + qv + ql + qi )

return,theta_r

end