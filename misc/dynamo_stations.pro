;
; Return a structure containing the DYNAMO sounding station info.
;
; James Ruppert
; 24.11.16
;
function dynamo_stations


;SOUNDING STATIONS
    ;stn_info=[[$
    ;         0            1              2          3          4            5
    stn = [ 'Gan' , 'R/V !8Revelle!X' , 'Male' , 'Colombo' , 'Marai' , 'Diego Garcia' ];,[$
    ;Lat
    lat = [ -0.69 , 0.0 , 4.19 , 6.91 , -8.0 , -7.31 ];,[$
    ;Lon
    lon = [ 73.15 , 80.5 , 73.53 , 79.87 , 80.5 , 72.43 ];]
    stn_info = {stn:stn,lat:lat,lon:lon}


return,stn_info

end
