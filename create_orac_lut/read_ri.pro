function read_ri, Filename
; Procedure to read in refractive index data file Filename (These
; files should have the filename extension .ri).
; The returned structure contains
; .wvl the vavelngth in um
; .wvn the wavenumbe in cm^-1
; n - the real part of the refractive index
; k - imaginary part of the refractive index (NOTE this value is positive)
; dn - 1 sigma uncertainty in n
; dk - 1 sigma uncertainty in k
; One of n,dn,k,dk needs to be present in the input file but not all
; values need to be.
;
;
; 06-12-09 RGG Created 
; 11-01-10 RGG Modified to include FORMAT Keyword

  Get_LUN, Unit
 
  If (STRMID(Filename, 2, /REVERSE_OFFSET) EQ '.gz') Then $
    OPENR, Unit, Filename, /COMPRESS $
  ELSE $
    OPENR, Unit, Filename

  Line =''

  Lines = 0 ; number of lines in file
  HL = 0    ; number of lines in header
  HH = 0    ; number of continuation lines in header

; First count the number of lines and lines starting with # and ##
; to calculate the number of Tags and size of header
  while Not EOF(Unit) do begin
    Readf,Unit,Line
    Lines = Lines + 1
    If (strmid(Line,0,1) Eq '#') Then HL = HL + 1
    If (strmid(Line,0,2) Eq '##') Then HH = HH + 1
  endwhile
 
  Tags = HL - HH ; Number of Tags
  Vals = Lines - HL 

  Tagname = strarr(Tags)
  Tagvlue = strarr(Tags)


; Create output structure

  RI = {DESCRIPTION: '', $
      DISTRIBUTEDBY: '', $
          SUBSTANCE: '', $
         SAMPLEFORM: '', $
        TEMPERATURE: '', $
      CONCENTRATION: '', $
          REFERENCE: '', $
                DOI: '', $
             SOURCE: '', $
            CONTACT: '', $
            COMMENT: '', $
             FORMAT: '', $ 
               Vals: Vals,$
               Wavn: Dblarr(Vals),$
               Wavl: Dblarr(Vals),$
                  n: Dblarr(Vals),$
                 dn: Dblarr(Vals),$
                  k: Dblarr(Vals),$
                 dk: Dblarr(Vals)}

  TN=TAG_NAMES(RI)

; Reset the file to the start
  Point_Lun, Unit,0

; Read Header
  For Count = 0, HL - 1 do begin
    Readf,Unit,Line
    If (strmid(Line,0,2) Eq '##') Then $
      RI.(Index) = RI.(Index) + ' '+ strtrim(strmid(Line,2,strlen(Line)-2),2)$
    Else Begin
      I = strpos(Line,'=')
      Tagname = strupcase(strtrim(strmid(Line,1,I-1),2))
      Tagvlue = strtrim(strmid(Line,I+1, strlen(Line)-I),2)
      Index = where(Tagname Eq TN, Match)
      If (Match Eq 0) Then $
        print,'Warning in read_ri:',Tagname,' tag unknown so ignored' $
      Else $
        RI.(Index) = Tagvlue
    EndElse  
  EndFor 
; Read Values
  FMT = strupcase(strsplit(Ri.format,Count=Count,/extract))
  Values = dblarr(Count,Vals)
  Readf,Unit,Values
  Close, Unit
  Free_LUN, Unit

  For I = 0, Count-1 do begin
    case FMT(i) of
    'WAVN': RI.Wavn = Values(I,*)
    'WAVL': RI.Wavl = Values(I,*)
       'N': RI.N = Values(I,*)
      'DN': RI.dN = Values(I,*)
       'K': RI.K = Values(I,*)
      'DK': RI.dK = Values(I,*)
    Endcase
   Endfor
; Convert wavl to wavn or vice versa
  If (RI.Wavl(0) Eq 0) then Ri.Wavl = 1D4/RI.Wavn
  If (RI.Wavn(0) Eq 0) then Ri.Wavn = 1D4/RI.Wavl
  Return, RI
End
