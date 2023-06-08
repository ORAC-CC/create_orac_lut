; procedure load_lutstr
;
; Read an LUT grid definition file.
;
; File should contain the following lines, with values separated by spaces.
;  (1) Number of optical depth LUT values; then
;      'linear'             if the step size is even in linear space
;      'uneven_linear'      if the step size is generally even in linear space (and should be interpolated in this space)
;      'logarithmic'        if the step size is even in natural log space
;      'uneven_logarithmic' if the step size is generally even in natural log space (and should be interpolated in this space)
;  (2) Optical depth LUT values themselves
;  - Repeat of the above two lines, but for:
;     - effective radius
;     - solar zenith
;     - satellite zenith
;     - relative azimuth
; Comment lines can be included and must have "#" as the first character in the
; line. Comments must not be placed between the (1) and (2) lines.
;
; INPUT ARGUMENTS:
; file (string) Path to file.
;
; INPUT KEYWORDS:
; None
;
; OUTPUT ARGUMENTS:
; LUTstr (structure) Structure with grid information.
;
; HISTORY:
; 02/03/21  RGG Use maximumzenith angle to limit number of zenith angles to those observed.  
; 17/02/21  RGG Generalised to state spacing for all dimensions rather than use difference between 1st two elements.  
; 21/06/13, G Thomas: Original version.

Function lut_quadrature, N, spacing, values , max_val=max_val
  Case 1 of 
  Spacing Eq 'linear': if keyword_set(max_val) Then X = values(0)+(max_val-values(0))*findgen(N)/(N-1) else X = values(0)+(values(1)-values(0))*findgen(N)/(N-1)
  Spacing Eq 'logarithmic':if keyword_set(max_val) then X = 10^(alog10(values(0))+alog10(max_val/values(0))*findgen(N)/(N-1)) else X = 10^(alog10(values(0))+alog10(values(1)/values(0))*findgen(N)/(N-1))
  Spacing Eq 'uneven_linear' or Spacing Eq 'uneven_logarithmic': X = values
  Else:  message,"Invalid spacing descriptor in LUT definition. Must be one of: 'linear', 'uneven_linear', 'logarithmic' or 'uneven_logarithmic'"    
  EndCase
  Return, X
End

pro load_lutstr, file, max_sat_zenith, LUTstr

; Read in the lut description
  Lines = strarr(File_lines(file))
  openr, lun, file, /get_lun
    readf, lun, lines
  close, lun
  free_lun, lun 
; Remove commented lines (those starting with "*" or "#".
  Q = WHERE(Lines NE '' AND STRMID(Lines,0,1) NE '*' AND STRMID(lines,0,1) NE '#', Count )
  Lines  = Lines[Q]    


; First two data lines should be optical depth values
  Words = strsplit(lines(0),' ',/extract)
  Opd_N = fix(Words[0])
  opd_Spacing = strlowcase(Words[1])
  opd = lut_quadrature(Opd_N, opd_spacing,float(strsplit(lines(1),' ',/extract))) 

; Next two data lines should be effective radius values
  Words = strsplit(lines(2),' ',/extract)
  Efr_N = fix(Words[0])
  Efr_Spacing = strlowcase(Words[1])
  Efr = lut_quadrature(Efr_N, Efr_spacing, float(strsplit(lines(3),' ',/extract)))

; Next two data lines should be solar zenith angle values
  Words = strsplit(lines(4),' ',/extract)
  Soz_N = fix(Words[0])
  Soz_Spacing = strlowcase(Words[1])
  Soz = lut_quadrature(Soz_N,Soz_spacing,float(strsplit(lines(5),' ',/extract)) )

; Next two data lines should be effective radius values
  Words = strsplit(lines(6),' ',/extract)
  Saz_N = fix(Words[0])
  Saz_Spacing = strlowcase(Words[1])
  Saz = lut_quadrature(  Saz_N,Saz_spacing,float(strsplit(lines(7),' ',/extract)),max_val=max_sat_zenith) 

; Next two data lines should be effective radius values
  Words = strsplit(lines(8),' ',/extract)
  Raa_N = fix(Words[0])
  Raa_Spacing = strlowcase(Words[1])
  If (Raa_Spacing Ne 'linear') Then stop,'Relative azimuth spacing must be linear' 
  Raa = lut_quadrature(  Raa_N,Raa_spacing,float(strsplit(lines(9),' ',/extract)) )

;  Build the output structure
   LUTstr = { Opd_N      : Opd_N      , Efr_N:       Efr_N      , Soz_N       : Soz_N      , Saz_N      : Saz_N      , Raa_N      : Raa_N, $
              OPD_Spacing: OPD_Spacing, Efr_Spacing: Efr_Spacing, Soz_Spacing : Soz_Spacing, Saz_Spacing: Saz_Spacing, Raa_Spacing: Raa_Spacing,  $
              opd        : opd        , EFR:         EfR        , SOz         : Soz        , Saz        : Saz        , Raa        : Raa    }

end
