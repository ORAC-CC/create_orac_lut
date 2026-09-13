;+
 PRO LOADXYZ,Filename,a,b,c,X,Y,Z,nskip=nskip
;
;
; General purpose routine for loading x,y and z arrays of data  from a '.TABLE'
; type file.  On input a, b and c are assumed to point to the columns containing 
; the X, Y and Z data (counting from 0). nskip is the number of lines to skip over.
; Date        Author Comment
; 2  FEB 1995   RGG   Created
; 22 NOV 1995   RGG   Updated. Checks to see if the file exists. Allows comment 
;                      fields (lines that start with a '*' throughout the file.
;                      Dynamically assigns own unit number.
; 08 MAY 2009   AMS   Added ability to skip nskip lines at the start
; 23 SEP 2014   AJAS  Set the lines to the number of lines in the file, using
;                      FILE_LINES(). Replaced obsolete STR_SEP with STRSPLIT.
; 26 MAY 2015   AJAS  Explicitly check that file exists. Cleaned up error checking.
;
;-
    COMPILE_OPT IDL2
    ON_ERROR, 2

    IF ~FILE_TEST( STRING(filename) ) THEN MESSAGE,'File not found: '+STRING(filename)

    IF NOT KEYWORD_SET(nskip) THEN nskip=0

    ;; Number of lines in the file.
    nl = FILE_LINES( filename )
    IF nl EQ 0 THEN MESSAGE,'Empty file: '+filename 
    
    ;; Open the file, read it all into a string array, and then close the file.
    OPENR, Unit, Filename, /GET_LUN
    Lines = STRARR( nl )
    Readf, Unit, Lines
    Close, Unit
    Free_LUN, Unit

    ;; Remove commented lines (those starting with "*" or "#".
    Q = WHERE(Lines NE '' AND STRMID(Lines,0,1) NE '*' AND $
              STRMID(lines,0,1) NE '#', Count )

    IF (Count GT 0) THEN BEGIN
       Line = STRSPLIT(STRTRIM(STRCOMPRESS(Lines[Q[0]]),2),/EXTRACT)
       COLUMNS = N_Elements(Line)
       format = '("Invalid ",A1," column specification in '+Filename+'")'
       IF ((a LT 0) or (a GT COLUMNS)) THEN MESSAGE,STRING('X',FORMAT=format)
       IF ((B Lt 0) or (b GT COLUMNS)) THEN MESSAGE,STRING('Y',FORMAT=format)
       IF ((C Lt 0) or (C GT COLUMNS)) THEN MESSAGE,STRING('Z',FORMAT=format)
       X = DBLARR(Count-nskip)
       Y = DBLARR(Count-nskip)
       Z = DBLARR(Count-nskip)
       FOR I = 0L, Count - 1-nskip DO BEGIN
          Line = STRSPLIT(STRTRIM(STRCOMPRESS(Lines[Q[I+nskip]]),2),' ',/EXTRACT)
          X[I] = DOUBLE(Line[A])
          Y[I] = DOUBLE(Line[B])
          Z[I] = DOUBLE(Line[C])
       ENDFOR
    ENDIF
    RETURN


 END

