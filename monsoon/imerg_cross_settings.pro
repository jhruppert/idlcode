; 
; Provide the settings for cross sections.
;
; Called by coastnormal.pro
;
; James Ruppert
; 8/11/21
; 
pro imerg_cross_settings, icross, xcross=xcross, ycross=ycross, bounds_ind=bounds_ind

  if icross eq 1 then begin

  ;BOB

  ;SW-NE
;    xcross=[77.8,97.4] ; lon,lat
;    ycross=[12.,20.]   ; lon,lat
    xcross=[77.8,96.2] ; lon,lat
    ycross=[12.,21.5]   ; lon,lat

  ;SUBSET BOX
    bounds_ind=[89.,15.5,95.,22.5] ; Northern Myanmar coastline
;bounds_ind=[89.,10.,102,22.5] ; Northern Myanmar coastline

  endif else if icross eq 2 then begin

  ;WG

  ;SW-NE
    xcross=[77.8,96.2] -13. ; lon,lat
    ycross=[12.,21.5] - 2.5  ; lon,lat

  ;SUBSET BOX
    bounds_ind=[69.5,8.,77.5,21.] ; Western Ghats

  endif else if icross eq 3 then begin

  ;BOB

  ;NW-SE
    xcross=[83.,93.5] ; lon,lat
    ycross=[21.4,10]   ; lon,lat
      xcross=reverse(xcross)
      ycross=reverse(ycross)

  ;SUBSET BOX
    bounds_ind=[78.,14.,89.7,24.5] ; Northwestern BoB coastline

  endif


end
