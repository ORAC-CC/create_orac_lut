; procedure read_presdat
;
; Read a pressure profile definition file.
;
; The pressure profile driver file simply contains three columns:
;    Height (in km, top to bottom of atmosphere)
;    Pressure(in hPa)
;    Temperature (in K)
; NB. This has not changed from the old LUT code, except that there now
;    be comment lines, beginning with "#",  at the start of the file.
;
; INPUT ARGUMENTS:
; file (string) Path to file.
;
; INPUT KEYWORDS:
; None
;
; OUTPUT ARGUMENTS:
; presstr (structure) Structure with meteorological profile information.
;
; HISTORY:
; 21/06/13, G Thomas: Original version.

pro read_presdat, file, presstr

   line = ''
   openr, lun, file, /get_lun
;  Comment lines start with "#"
   readf,lun, line
   while strmid(line,0,1) eq '#' do readf,lun, line
   row = float(strsplit(line,/extract))
   H = row[0]
   P = row[1]
   T = row[2]
   row = fltarr(3)
   while ~EOF(lun) do begin
      readf,lun, row
      H = [H,row[0]]
      P = [P,row[1]]
      T = [T,row[2]]
   endwhile
   free_lun, lun
   presstr = { NLevels: n_elements(H), $
               height: H, pressure: P, temperature: T }
end

; procedure load_atmstr
;
; Read a pressure profile definition file (rfm atm file).
;
; INPUT ARGUMENTS:
; file (string) Path to file.
; atmospheres (integer) Code of desired atmosphere 0 = old style read.
;
; INPUT KEYWORDS:
; None
;
; OUTPUT ARGUMENTS:
; presstr (structure) Structure with meteorological profile information.
;
; HISTORY:
; 21/06/13, G Thomas: Original version.

pro load_atmstr, file, atmospheres, atmstr

  If (Atmospheres Eq 'midsatm.dat') Then begin
    read_presdat, file, atmstr
    
  EndIf Else Begin
   line = ''
   OPENR, lun, file, /GET_LUN
;  Comment lines start with "!"
     READF,lun, line
     WHILE strmid(line,0,1) eq '!' DO READF,lun, line
     nlevels = (STRSPLIT(line,/EXTRACT)) (0)   
     H = FLTARR(nlevels)
     P = FLTARR(nlevels)
     T = FLTARR(nlevels)
     READF,lun, line
     READF,lun, H
     READF,lun, line
     READF,lun, P
     READF,lun, line
     READF,lun, T
   FREE_LUN, lun

   Q = WHERE (h LE 100, nlevels) ;  MODTRAN Code limited to 100 km
   atmstr = { nlevels: nlevels, height: REVERSE(h[q]), pressure: REVERSE(p[q]), temperature: REVERSE(t[q])} 
   EndElse
   end

