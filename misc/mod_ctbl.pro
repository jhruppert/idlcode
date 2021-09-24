;
; Modify my own color table.
;
; James Ruppert
; 17 April 2019
pro mod_ctbl

  config_dir, dirs=dirs
  ctbl_file='/home1/06040/tg853394/idl/code/misc/ctbl_jhr_17Apr2019.tbl'
  ;this table was copied from the default IDL color table file.

  dopng=1 ; Generate png of new color bar?

  name='rainfall' ; index = 75
   ;^ used for wspd history
;  name='olr'  ; 76
;  name='satir' ; 77 and 78


print,'NEW CBAR: ',strupcase(name)


ncol=256


if name eq 'rainfall' then begin


;RAINFALL, ALSO USED FOR WSPD HISTORY MAPS

  index=75 ; output index

;;  ctbl_ref=11 ; USE A PRIOR-TABLE'S DEFAULT AS STARTING POINT
;  ctbl_ref=33;4;13 ; USE A PRIOR-TABLE'S DEFAULT AS STARTING POINT
;  loadct,ctbl_ref,rgb_table=rgb_ref
;  loadct,0,rgb_table=rgb_gray
;  rgb_table=rgb_ref
;
;  ;REPLACE LOW INDICES WITH WHITE-GRAY
;  top_mod=40;20;35
;;  gray_ind = round( 215-findgen(top_mod)*2.5 )
;  gray_ind = round( 200-findgen(top_mod)*2.5 )
;  rgb_table[ 0:top_mod-1 , * ] = rgb_gray[ gray_ind , * ]
;  ;OVERWRITE BOTTOM INDEX WITH WHITE
;rgb_table[ 0 , * ] = rgb_gray[ 255 , * ]
;
;  ;ADD SOME BLUE HUE (BY REDUCING RED)
;;  for i=0,top_mod-1 do $
;  for i=1,top_mod-1 do $
;    rgb_table[i,0] -= round( i*2.0 + 2 )

      ;SET ANCHORS AND BLEND

      rgb_table=fltarr(ncol,3)

      ;ANCHORS
        anchors = [$
          [255,255,255], $ ; white
        ;gray-blue
          [220,220,235], $ ;R,G,B
          [200,220,235], $ ;R,G,B
          [163,197,225], $ ;R,G,B
          [96,158,205], $ ;R,G,B
          [38,127,205], $ ;R,G,B
          [35,100,200], $ ;R,G,B
        ;greens
          [0,120,130], $ ;R,G,B
          [0,140,100], $ ;R,G,B
          [35,176,70], $ ;R,G,B
          [35,220,20], $ ;R,G,B
          [10,245,10], $ ;R,G,B
;          [], $ ;R,G,B
;          [255,255,255], $ ;R,G,B
;          [], $ ;R,G,B
;          [], $ ;R,G,B
        ;yellow
          [255,245,30], $ ;R,G,B
          [255,185,25], $ ;R,G,B
          [239,166,20], $ ;R,G,B
          [230,104,16], $ ;R,G,B
          [230,40,10], $ ;R,G,B
          [230,0,0] ] ; red
      ;ANCHOR LOCATIONS (INPUT AS FRACTIONS OF 1
        nanc=1.*(size(anchors,/dim))[1]
        anloc=round(1.*ncol*findgen(nanc)/(nanc-1))

      ;BLEND AND FILL COLOR TABLE
      nanc=n_elements(anloc)
      nloc=anloc[indgen(nanc-1)+1]-anloc[indgen(nanc-1)]
      for i=0,2 do begin
        for a=0,nanc-2 do begin
          ind = interpol([anchors[i,a:a+1]],nloc[a]+1*((a lt nanc-2)))
          cind=indgen(n_elements(ind))+anloc[a]
          rgb_table[cind,i] = round(ind)
        endfor
      endfor

endif else if name eq 'olr' then begin


;OLR - red-black

  index=76

  ctbl_ref=71
  loadct,ctbl_ref,rgb_table=rgb_ref
  rgb_table=rgb_ref

  ;REVERSE TABLE TO PUT RED AT TOP
  rgb_table = reverse(rgb_table,1)
  ;REVERSE GRAY-BLACK SECTION
  ixrev=indgen(ncol/2)
;  rgb_table[ixrev,*] = reverse(rgb_table[ixrev,*],1)

  ;SELECT RANGE FROM REFERENCE
    minsel=45
    maxsel=ncol-1
    rangesel=maxsel-minsel
    ix = indgen(ncol)

  ;STRETCH BOTTOM
    ;EXPONENTIAL FUNCTION SPECS
      pow=2.6
      xmin=1 & xmax=4
      xrange=xmax-xmin
    ;FUNCTION
      ix_new = (findgen(ncol)/(ncol-1) * xrange + xmin )^pow
      ix_new -= min(ix_new)
      ix_new = ix_new/max(ix_new) * rangesel + minsel
print,ix_new

  for i=0,2 do $
    rgb_table[*,i] = interpol(reform(rgb_table[*,i]),ix,ix_new)


endif else if name eq 'satir' then begin


;SATIR

  itype=3 ; 1 - single clr table
          ; 2 - Classic IR color table (NESDIS)
          ; 3 - Based on NASA Worldview for MODIS Tb
          ; 4 - 2nd attempt at NESDIS table (based on Jerry's code)

  if itype eq 1 then begin
  ;SINGLE EXPONENTIAL

  index=77

      ;NOW USING CTBL 11 WITH BOTTOM REPLACED BY LIGHT-GRAYS (TAKING FROM RAINFALL ROUTINE ABOVE)

      ctbl_ref=11 ; USE A PRIOR-TABLE'S DEFAULT AS STARTING POINT
      loadct,ctbl_ref,rgb_table=rgb_ref
      rgb_table=rgb_ref

      ;REPLACE LOW INDICES WITH WHITE-GRAY

          loadct,0,rgb_table=rgb_gray
        
          top_mod=35
          gray_ind = round( 215-findgen(top_mod)*2.5 )
          rgb_table[ 0:top_mod-1 , * ] = rgb_gray[ gray_ind , * ]
        
          ;ADD SOME BLUE HUE (BY REDUCING RED)
          for i=0,top_mod-1 do $
            rgb_table[i,0] -= round( i*2.0 + 2 )

      ;REVERSE TABLE TO PUT BLACK AT TOP
      rgb_table = reverse(rgb_table,1)
    
      ;SELECT RANGE FROM REFERENCE
        minsel=.2*255;0
        maxsel=255
        rangesel=maxsel-minsel
        ix = indgen(ncol)
    
      ;STRETCH BOTTOM
        ;EXPONENTIAL FUNCTION SPECS
          pow=3.2;2.6
          xmin=1 & xmax=4
          xrange=xmax-xmin
        ;FUNCTION
          ix_new = (findgen(ncol)/(ncol-1) * xrange + xmin )^pow
          ix_new -= min(ix_new)
          ix_new = ix_new/max(ix_new) * rangesel + minsel
    
      for i=0,2 do $
        rgb_table[*,i] = interpol(reform(rgb_table[*,i]),ix,ix_new)

      ;REVERSE TABLE AGAIN
      rgb_table = reverse(rgb_table,1)

  endif else if itype eq 2 then begin
  ;CLASSIC IR COLORMAP (NESDIS STYLE)

  index=78

      ;SET ANCHORS AND BLEND

      rgb_table=fltarr(ncol,3)

      ;ANCHORS
        anchors = [$
          [114,26,120], $ ;R,G,B
          [238,138,200], $ ;R,G,B
          [230,230,230], $ ;R,G,B
          [5,5,5], $ ;R,G,B
          [235,53,35], $ ;R,G,B
          [244,173,60], $ ;R,G,B
          [255,254,85], $ ;R,G,B
          [139,250,77], $ ;R,G,B
          [59,132,65], $ ;R,G,B
          [1,23,101], $ ;R,G,B
          [105,231,241], $ ;R,G,B
          [115,252,253], $ ;R,G,B
          [193,193,193], $ ;R,G,B
          [0,0,0] ] ; black
      ;ANCHOR LOCATIONS (INPUT AS FRACTIONS OF 1
        anloc=round(1.*ncol*[0,1,1.1,2,3,3.5,4,4.5,5.5,6,8,8.1,8.2,16]/16.)

      ;BLEND AND FILL COLOR TABLE
      nanc=n_elements(anloc)
      nloc=anloc[indgen(nanc-1)+1]-anloc[indgen(nanc-1)]
      for i=0,2 do begin
        for a=0,nanc-2 do begin
          ind = interpol([anchors[i,a:a+1]],nloc[a]+1*((a lt nanc-2)))
          cind=indgen(n_elements(ind))+anloc[a]
          rgb_table[cind,i] = round(ind)
        endfor
      endfor

  endif else if itype eq 3 then begin
  ;BASED ON NASA WORLDVIEW MODIS TB COLOR BAR

  index=79

      ;SET ANCHORS AND BLEND

      rgb_table=fltarr(ncol,3)

      ;ANCHORS
        anchors = [$
;          [2,2,22], $ ; near-black
;          [41,51,118], $ ;R,G,B
;          [83,74,140], $ ;R,G,B
;          [137,98,149], $ ;R,G,B
;          [170,119,142], $ ;R,G,B
;          [207,157,149], $ ;R,G,B
;          [242,222,195], $ ;R,G,B
;          [255,251,230], $ ;R,G,B
;          [255,255,255] ] ; white
;Revised
          [2,2,22], $ ; near-black
          [30,34,103], $ ;R,G,B
          [68,52,126], $ ;R,G,B
          [123,73,134], $ ;R,G,B
          [164,96,123], $ ;R,G,B
          [203,133,127], $ ;R,G,B
          [241,214,177], $ ;R,G,B
          [254,250,221], $ ;R,G,B
          [255,255,255] ] ; white
      ;ANCHOR LOCATIONS (INPUT AS FRACTIONS OF 1
        anloc=round(1.*ncol*findgen(9)/8)

      ;BLEND AND FILL COLOR TABLE
      nanc=n_elements(anloc)
      nloc=anloc[indgen(nanc-1)+1]-anloc[indgen(nanc-1)]
      for i=0,2 do begin
        for a=0,nanc-2 do begin
          ind = interpol([anchors[i,a:a+1]],nloc[a]+1*((a lt nanc-2)))
          cind=indgen(n_elements(ind))+anloc[a]
          rgb_table[cind,i] = round(ind)
        endfor
      endfor

  endif else if itype eq 4 then begin
  ;REDUX OF CLASSIC (NESDIS) STYLE IR TABLE

  index=80

      ;SET ANCHORS AND BLEND

      rgb_table=fltarr(ncol,3)

           red = [[0.0,  0.5, 0.5],$ ; R
                 [1./28, 1.0, 1.0],$
                 [1./14, 1.0, 0.8],$
                 [2./14, 0.0, 0.0],$
                 [3./14, 1.0, 1.0],$
                 [4./14, 1.0, 1.0],$
                 [5./14, 0.0, 0.0],$
                 [7./14, 0.0, 1.0],$
                 [1.0,  0.0, 0.0]]
         green = [[0.0,  0.0, 0.0],$ ; G
                 [1./28, 0.0, 0.0],$
                 [1./14, 0.5, 0.8],$
                 [2./14, 0.0, 0.0],$
                 [3./14, 0.0, 0.0],$
                 [4./14, 1.0, 1.0],$
                 [5./14, 1.0, 1.0],$
                 [6./14, 0.0, 0.0],$
                 [7./14, 1.0, 1.0],$
                 [1.0,  0.0, 0.0]]
          blue = [[0.0,  0.5, 0.5],$ ; B
                 [1./28, 1.0, 1.0],$
                 [1./14, 1.0, 0.8],$
                 [2./14, 0.0, 0.0],$
                 [5./14, 0.0, 0.0],$
                 [6./14, 1.0, 1.0],$
                 [7./14, 1.0, 1.0],$
                 [1.0,  0.0, 0.0]]

      red  *=ncol-1
      green*=ncol-1
      blue *=ncol-1

      red  =round(red)
      green=round(green)
      blue =round(blue)

      ;RED
      nanc=(size(red,/dim))[1] - 1
      for ia=0,nanc-1 do begin
        nind=red[0,ia+1]-red[0,ia]+1
        ind=indgen(nind)+red[0,ia]
        rgb_table[ind,0] = interpol( reform([ red[2,ia] , red[1,ia+1] ]), nind )
      endfor

      ;GREEN
      nanc=(size(green,/dim))[1] - 1
      for ia=0,nanc-1 do begin
        nind=green[0,ia+1]-green[0,ia]+1
        ind=indgen(nind)+green[0,ia]
        rgb_table[ind,1] = interpol( reform([ green[2,ia] , green[1,ia+1] ]), nind ) 
      endfor

      ;BLUE
      nanc=(size(blue,/dim))[1] - 1
      for ia=0,nanc-1 do begin
        nind=blue[0,ia+1]-blue[0,ia]+1
        ind=indgen(nind)+blue[0,ia]
        rgb_table[ind,2] = interpol( reform([ blue[2,ia] , blue[1,ia+1] ]), nind ) 
      endfor

  endif

endif

;WRITE OUT
  modifyct,index,name,reform(rgb_table[*,0]),reform(rgb_table[*,1]),reform(rgb_table[*,2]),file=ctbl_file

;WRITE PNG
if dopng then begin
  set_plot,'ps'
  figname=dirs.figdir+'newcbar'
  epsname=figname+'.eps'
  !p.font=0
  device,filename=epsname,/encapsulated,/color,bits=8,xsize=7,ysize=1.,/inches,/helvetica
  x=indgen(ncol)
  y=indgen(5)
  loadct,0,/silent,file=dirs.ctfil
  plot,x,y,/nodata,position=[0.02,0.21,0.98,0.93],charsize=0.8,xminor=0,yminor=0,xstyle=1,xrange=[0,ncol-1]
  var=intarr(ncol,5) & for i=0,4 do var[*,i]=x
  loadct,index,/silent,file=dirs.ctfil
  for i=0,1 do contour,var,x,y,/cell_fill,/overplot,c_colors=x,levels=x
  loadct,0,/silent
  xyouts,0.5,0.5,'CTABLE = '+strtrim(index,2),align=0.5,charsize=0.8,/normal
  device,/close
  convert_png,figname,res=200,/remove_eps
endif

end
