;+
 Pro LOADXY,Filename,a,b,X,Y,nskip=nskip,lines=lines
;
; General purpose routine for loading x,y and z arrays of data  from a '.TABLE'
; type file.  On input a, b and c are assumed to point to the columns 
; containing the X and Y data (counting from 0). nskip is the number of lines to skip over.
; 
; The lines keyword is depricated but still exists for compatibility.
;
; Date        Author Comment
; 2  FEB 1995   RGG   Created
; 22 NOV 1995   RGG   Updated
;                     Checks to see if the file exists
;                     Allows comment fields (lines that start with a '*' throughout 
;                     the file
;                     Dynamically assigns own unit number                
; 08 MAY 2009   AMS   Added ability to skip nskip lines at the start
; 27 JAN 2011   DMP   Added Lines keyword to allow larger file reads
; 06 AUG 2014   AJAS  Dynamically decide number of lines using FILE_LINES() internal fn.
; 19 MAR 2015   AJAS  Rewrote to just call loadxyz.pro.
;-

   LOADXYZ, Filename, a,b,0, X,Y,dont_care_about_this_variable, nskip=nskip

   RETURN
END
