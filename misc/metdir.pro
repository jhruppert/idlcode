;**********************************************************************************
;*  METDIR                                                                        *
;*  A program that accepts vector components and returns resultant direction      *
;*    with 0 degrees at north, the meteorological convention.                     *
;*                                                                                *
;*  INPUT:                                                                        *
;*    - u,v (1 or more element arrays of vector components)                       *
;*  OUTPUT:                                                                       *
;*    - dir (same size array as input arrays, in degrees, with 0 at north)        *
;*                                                                                *
;*  Written by Brian McNoldy (mcnoldy@atmos.colostate.edu)                        *
;*  October 2008                                                                  *
;*                                                                                *
;*  To get the corresponding resultant magnitude from the u and v components,     *
;*    use sqrt(u^2+v^2)... this program doesn't do that for you.                  *
;*  To get the resultant direction using the mathematically conventional method,  *
;*    use atan(v/u)... this program won't do that for you either.                 *
;**********************************************************************************

pro metdir, u, v, dir

  size=size(u)
  xx=long(size[1]);n_elements(u))
  yy=long(size[2]);n_elements(v))
  arg=fltarr(xx,yy)
  dir=fltarr(xx,yy)

  for i=0L,xx-1 do begin
    for j=0L,yy-1 do begin

      if (abs(u[i,j]) gt abs(v[i,j])) then begin
        cas=1
        dx=v[i,j]
        dy=u[i,j]
      endif else begin
        cas=0
        dx=u[i,j]
        dy=v[i,j]
      endelse

      arg[i,j]=abs(dx/dy)
      if (u[i,j] gt 0) then begin
        if (v[i,j] gt 0) then begin          ;upper right quadrant
          if (cas eq 0) then begin
            dir[i,j] = 180. + atan(arg[i,j])/!dtor
          endif else begin
            dir[i,j] = 270. - atan(arg[i,j])/!dtor
          endelse
        endif else begin                     ;lower right quadrant
          if (cas eq 0) then begin
            dir[i,j] = 360. - atan(arg[i,j])/!dtor 
          endif else begin
            dir[i,j] = 270. + atan(arg[i,j])/!dtor 
          endelse 
        endelse 
      endif else begin
        if (v[i,j] gt 0) then begin          ;upper left quadrant
          if (cas eq 0) then begin
            dir[i,j] = 180. - atan(arg[i,j])/!dtor 
          endif else begin
            dir[i,j] = 90. + atan(arg[i,j])/!dtor 
          endelse
        endif else begin                     ;lower left quadrant
          if (cas eq 0) then begin
            dir[i,j] = 0. + atan(arg[i,j])/!dtor 
          endif else begin
            dir[i,j] = 90. - atan(arg[i,j])/!dtor 
          endelse
        endelse
      endelse

    endfor
  endfor
end