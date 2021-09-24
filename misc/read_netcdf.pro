FUNCTION read_netcdf, file,GMT=gmt,INFO=info
;file='/users/ifmnfs2a/spreen/tmp/gmt_orig.grd'

; +
; INPUT:
;   /GMT   netcdf file is a old (< Ver4.1) GMT grd file
;   /INFO  show some information
;-
; Sep 2005 new version,   Gunnar Spreen
;

;info=1
; ---- End user variables

id=NCDF_OPEN(file[0])

dinfo=NCDF_INQUIRE(id)

IF KEYWORD_SET(info) THEN BEGIN
    HELP,dinfo,/STRUCTURE
    FOR i=0,dinfo.nvars-1 DO BEGIN
        varinf=NCDF_VARINQ(id,i)
        HELP,varinf,/STRUCTURE
        PRINT,'DIM',varinf.dim
        FOR j=0,varinf.natts-1 DO BEGIN
            attname= NCDF_ATTNAME(id,i,j)
            PRINT,attname
        ENDFOR
    ENDFOR
	name = NCDF_ATTNAME(id,/GLOBAL,0)
	NCDF_ATTGET,id,/GLOBAL,name,description
	result = NCDF_ATTINQ(id,/GLOBAL, name)
	help,name,description,result,/STRUCTURE
	indescription = string(description)

	outinfo = strarr(1)
	outinfo(0) = indescription
	GOTO, output_time
ENDIF

; read data
IF KEYWORD_SET(GMT) THEN BEGIN
    NCDF_VARGET,id,'x_range',xr
    help,xr
    NCDF_VARGET,id,'y_range',yr
    help,yr
    NCDF_VARGET,id,'z_range',zr
    help,zr
    NCDF_VARGET,id,'spacing',spacing
    help,spacing
    NCDF_VARGET,id,'dimension',dim
    help,dim
    NCDF_VARGET,id,'z',z
    help,z
ENDIF ELSE BEGIN
    FOR i=0,dinfo.nvars-1 DO BEGIN
        varinf=NCDF_VARINQ(id,i)
        NCDF_VARGET,id,i,tdat
        natts=varinf.natts
        stest=0
        otest=0
        utest=0
        FOR ii=0,natts-1 DO BEGIN
            attname = NCDF_ATTNAME(id,i,ii)
            CASE attname OF
                'scale_factor': stest=1
                'add_offset': otest=1
                'units': utest=1
                ELSE:
            ENDCASE
        ENDFOR
        IF stest THEN BEGIN
            NCDF_ATTGET,id,i,'scale_factor',sfact
            tdat=sfact*TEMPORARY(tdat)
        ENDIF
        IF otest THEN BEGIN
            NCDF_ATTGET,id,i,'add_offset',offset
            tdat=TEMPORARY(tdat)+offset
        ENDIF
        IF utest THEN BEGIN
            NCDF_ATTGET,id,i,'units',unit
            tunit=STRING(unit)
        ENDIF

        tname=varinf.name
        tname = IDL_VALIDNAME(tname,/CONVERT_ALL)
        IF i EQ 0 THEN BEGIN
            IF utest THEN $
              data=CREATE_STRUCT(tname,tdat,tname+'_Unit',tunit) $
            ELSE $
              data=CREATE_STRUCT(tname,tdat)
        ENDIF ELSE BEGIN
            IF utest THEN $
              data=CREATE_STRUCT(TEMPORARY(data),tname,tdat,tname+'_Unit',tunit) $
            ELSE $
              data=CREATE_STRUCT(TEMPORARY(data),tname,tdat)
        ENDELSE
    ENDFOR
ENDELSE

NCDF_CLOSE,id
GOTO, output_data

output_time:
NCDF_CLOSE,id
RETURN,outinfo

GOTO, exit_routine

output_data:
RETURN,data

exit_routine:

END
