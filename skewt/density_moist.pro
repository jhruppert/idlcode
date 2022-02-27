function density_moist,t,p,qv

; CALCULATE DENSITY OF MOIST AIR PROVIDED WITH THE FOLLOWING:
; **ASSUMING FIELDS ARE PROVIDED IN SI UNITS**
;   TEMPERATURE (K)
;   PRESSURE (Pa)
;   QV, water vapor mass mixing ratio (kg/kg)
; 
; Program assumes you have the suite of thermodynamic programs by Dominik Brunner used below,
;   which can be found at http://www.iac.ethz.ch/staff/dominik/idltools/idl_atmosphys.html
; 
; James Ruppert, 12/13/12
; ruppert.@atmos.colostate.edu

rd=287.04
cp=1004.0
eps=rd/cp

stop
;FROM SUE'S NOTES: p0 * pi ^ (cv/Rd) / (th_v * Rd)
rho2 = 1e5*(ex^(717./287)) / (thv*287.) ; [ kg/m3 ]

;t_virt=t_virtual( t , p*1.e-2 , esat( rh2tdew(t,mixr2rh(qv*1.e3,p*1.e-2,t)) ) )

virt_corr = (1. + qv/eps)/(1.+qv)
;virt_corr = (1. + 0.61*qv)

density = p / ( rd * t * virt_corr )
;density = p / ( rd * t_virt )

return,density

end