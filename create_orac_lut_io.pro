; This library of contains the procedures for reading the various driver files
; needed by create_orac_lut, as well as the code for writing the ORAC lookup
; tables themselves.
;
; HISTORY:
; 21/06/13, G Thomas: Original version.




; procedure read_instdat
;
; Read an intrument definition file.
;
; This file should contain the instrument name on the first non-comment line,
; followed by a line giving the number of spectral channels and then four
; columns: channel number; centre wavelength; 1=thermal 0=not-thermal; 1=solar
; 0=not-solar.
;
; INPUT ARGUMENTS:
; file (string) Path to file.
;
; INPUT KEYWORDS:
; None
;
; OUTPUT ARGUMENTS:
; inststr (structure) Structure with instrument information.
;
; HISTORY:
; 21/06/13, G Thomas: Original version.

pro read_instdat, file, inststr
   line = ''
   openr, lun, file, /get_lun
   readf, lun, line
;  Skip over comment lines
   while strmid(line,0,1) eq '#' do readf, lun, line
;  Read the inststrrument name
   name = strtrim(line,2)
;  Read the number of channels
   readf, lun, line
   NChan = fix(line)
;  Define the channel arrays
   ChanName = intarr(NChan)
   ChanWl   = fltarr(NChan)
   ChanEm   = intarr(NChan)
   ChanSl   = intarr(NChan)
;  Read the channel data
   dat = fltarr(4)
   for i=0,NChan-1 do begin
      readf,lun,dat
      ChanName[i] = fix(dat[0])
      ChanWl[i]   = dat[1]
      ChanEm[i]   = fix(dat[2])
      ChanSl[i]   = fix(dat[3])
   endfor
   free_lun, lun

   inststr = {name: name, NChan: NChan, $
              ChanNum: ChanName,        $
              ChanWl:  ChanWl,          $
              ChanEm:  ChanEm,          $
              ChanSol: ChanSl           }
end





; procedure read_srfdat
;
; Read a channel spectral response function file.
;
; Following optional comment lines (starting with '#') this is a file with two
; columns: wavelength in microns and a relative spectral response as a function
; of wavelength.
;
; INPUT ARGUMENTS:
; file (string) Path to file.
;
; INPUT KEYWORDS:
; None
;
; OUTPUT ARGUMENTS:
; srfstr (structure) Structure with spectral response function information.
;
; HISTORY:
; 14/04/18, G McGarragh: Original version.

pro read_srfdat, file, srfstr

   line = ''
   openr, lun, file, /get_lun
   readf, lun, line
;  while strmid(line,0,1) eq '#' do readf, lun, line
   row = fltarr(2)
   readf, lun, row
   wl  = row[0]
   srf = row[1]
   while ~EOF(lun) do begin
      readf, lun, row
      wl  = [wl, row[0]]
      srf = [srf,row[1]]
   endwhile
   free_lun, lun
   srfstr= { NWls: n_elements(wl), wl: wl, srf: srf }
end


; LUT output procedures follow...
;
; The following routines write the LUT files, with one procedure per file
; dimension (1d, 2d, 3d or 5d LUTs).
;
; INPUT ARGUMENTS:
; lutfile (string) Full path to the file to write to.
; Axis    (array) Either 1, 2, 3 or 5-D vectors defining the axes of the LUT
;                 (OPD, EfR, Solar ZA, Sat ZA, RAA).
; Data    (array) An array containing the first set of LUT values.
;
; INPUT KEYWORDS:
; Wl=value    The wavelength corresponding to the LUT.
; Data2=array A second LUT data array to write to the file.
; /logOPD     The optical depth values should be output as log10(OPD).
; /logEfR     The effective radius values should be output as log10(EfR).
;
; OUTPUT ARGUMENTS:
; None
;
; HISTORY:
; 21/06/13, G Thomas: Original version.
; XX/XX/15, G McGarragh: Add write_orac_lut_1d for 1-D output.

pro write_orac_lut_1d, lutfile, EfR, Data, Wl=Wl, logEfR=logEfR
   if keyword_set(logEfR) then outEfR = alog10(EfR) else outEfR = EfR

   NEfR = n_elements(EfR)
   DelEfR = outEfR[1] - outEfR[0]

   openw, lun, lutfile, /get_lun
   if n_elements(Wl) gt 0 then printf,lun, Wl
   lutformat = '(10E14.6)'
   printf, lun, NEfR, DelEfR
   printf, lun, outEfR, format=lutformat
   printf, lun, Data, format=lutformat
   free_lun, lun
end ; write_orac_lut_1d


pro write_orac_lut_2d, lutfile, OPD, EfR, Data, Wl=Wl, Data2=Data2, $
                       logAOD=logAOD, logEfR=logEfR
   if keyword_set(logAOD) then outOPD = alog10(OPD) else outOPD = OPD
   if keyword_set(logEfR) then outEfR = alog10(EfR) else outEfR = EfR

   NOPD = n_elements(OPD)
   DelOPD = outOPD[1] - outOPD[0]
   NEfR = n_elements(EfR)
   DelEfR = outEfR[1] - outEfR[0]

   openw, lun, lutfile, /get_lun
   if n_elements(Wl) gt 0 then printf,lun, Wl
   lutformat = '(10E14.6)'
   printf, lun, NOPD, DelOPD
   printf, lun, outOPD, format=lutformat
   printf, lun, NEfR, DelEfR
   printf, lun, outEfR, format=lutformat
   printf, lun, Data, format=lutformat
   if n_elements(Data2) gt 1 then $
      printf,lun, Data2, format=lutformat
   free_lun, lun
end ; write_orac_lut_2d


pro write_orac_lut_3d, lutfile, OPD, Sol, EfR, Data, Wl=Wl, $
                       Data2=Data2, logAOD=logAOD, logEfR=logEfR
   if keyword_set(logAOD) then outOPD = alog10(OPD) else outOPD = OPD
   if keyword_set(logEfR) then outEfR = alog10(EfR) else outEfR = EfR

   NOPD = n_elements(OPD)
   DelOPD = outOPD[1] - outOPD[0]
   NSol = n_elements(Sol)
   DelSol = Sol[1] - Sol[0]
   NEfR = n_elements(EfR)
   DelEfR = outEfR[1] - outEfR[0]

   openw, lun, lutfile, /get_lun
   if n_elements(Wl) gt 0 then printf,lun, Wl
   lutformat = '(10E14.6)'
   printf, lun, NOPD, DelOPD
   printf, lun, outOPD, format=lutformat
   printf, lun, NSol, DelSol
   printf, lun, Sol, format=lutformat
   printf, lun, NEfR, DelEfR
   printf, lun, outEfR, format=lutformat

   printf,lun, Data, format=lutformat
   if n_elements(Data2) gt 1 then begin
      sz = size(Data2)
      if sz[0] eq 3 then printf,lun, Data2, format=lutformat $
      else printf,lun, Data2, format=lutformat
   endif
   free_lun, lun
end ; write_orac_lut_3d


pro write_orac_lut_5d, lutfile, OPD, Sat, Sol, Azi, EfR, Data, Wl=Wl, $
                       Data2=Data2, logAOD=logAOD, logEfR=logEfR
   if keyword_set(logAOD) then outOPD = alog10(OPD) else outOPD = OPD
   if keyword_set(logEfR) then outEfR = alog10(EfR) else outEfR = EfR

   NOPD = n_elements(OPD)
   DelOPD = outOPD[1] - outOPD[0]
   NSat = n_elements(Sat)
   DelSat = Sat[1] - Sat[0]
   NSol = n_elements(Sol)
   DelSol = Sol[1] - Sol[0]
   NAzi = n_elements(Azi)
   DelAzi = Azi[1] - Azi[0]
   NEfR = n_elements(EfR)
   DelEfR = outEfR[1] - outEfR[0]

   openw, lun, lutfile, /get_lun
   if n_elements(Wl) gt 0 then printf,lun, Wl
   lutformat = '(10E14.6)'
   printf, lun, NOPD, DelOPD
   printf, lun, outOPD, format=lutformat
   printf, lun, NSat, DelSat
   printf, lun, Sat, format=lutformat
   printf, lun, NSol, DelSol
   printf, lun, Sol, format=lutformat
   printf, lun, NAzi, DelAzi
   printf, lun, Azi, format=lutformat
   printf, lun, NEfR, DelEfR
   printf, lun, outEfR, format=lutformat

   printf, lun, Data, format=lutformat
   if n_elements(Data2) gt 1 then begin
      sz = size(Data2)
      if sz[0] eq 3 then printf,lun, Data2, format=lutformat $
      else printf,lun, Data2, format=lutformat
   endif
   free_lun, lun
end ; write_orac_lut_5d
