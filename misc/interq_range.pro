; 
; PROGRAM TO CALCULATE THE INTERQUARTILE RANGE OF A DATASET
; 
; TAKEN FROM DAVID FANNING'S "BOXPLOT.PRO"
; 
; JAMES RUPPERT
; 26.1.2017

function interq_range, data

  ; Sort the data.
  sortedData = data[Sort(data)]
;  IF N_Elements(sortedData) MOD 2 EQ 0 THEN BEGIN
;    index = N_Elements(sortedData)/2
;    medianData = (sortedData[index-1] + sortedData[index]) / 2.0
;    lowerGroup = sortedData[0:index-1]
;    higherGroup = sortedData[index:N_Elements(data)-1]
;  ENDIF ELSE BEGIN ; The middle point belongs to both upper and lower quartiles.
    index = N_Elements(sortedData)/2
    medianData = sortedData[index]
    lowerGroup = sortedData[0:index]
    higherGroup = sortedData[index:N_Elements(data)-1]
;  ENDELSE
  
  ; Find the quartiles.
  quartile_25 = Median(lowerGroup, /EVEN)
  quartile_75 = Median(higherGroup, /EVEN)
  
  ; Calculate IQR
  iqr = quartile_75 - quartile_25
  
  return, iqr

end