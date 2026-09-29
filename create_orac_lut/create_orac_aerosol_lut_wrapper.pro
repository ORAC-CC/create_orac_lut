; procedure create_orac_aerosol_lut_wrapper
;
; A wrapper procedure (plus two I/O functions) for calling the main procedure
; create_orac_lut. The create_orac_lut_wrapper procedure accepts no arguments,
; so is compatible with creating a save file for use with the IDL virtual
; machine or runtime mode.
;
; The procedure reads a text file containing the arguments needed for calling
; create_orac_lut, the path to which is read from the environment variable
; 'CREATE_ORAC_LUT_DRIVER'. This file can contain comment lines (beginning with
; a '#' character) and gives the compulsory arguments in order:
;    - input path
;    - instrument driver file name
;    - pressure profile driver file name
;    - scattering properties driver file name
;    - lut definition driver file name
;    - output path
;
; Optional parameters can also be included in this file, as they would appear in
; a call to create_orac_lut. Eg.
;    - ChannelID=['01','02','06','20','31','32']
;    - /no_rayleigh
;    - /no_screen
;    - version='20'
;
; ... see the header for function create_orac_lut for the complete list.
;
; See the header for create_orac_lut for a full list of optional parameters.
;
; INPUT ARGUMENTS:
; None
;
; INPUT KEYWORDS:
; Keyword inputs correspond directly with the arguments and keyword inputs of
; procedure create_orac_lut and overwrite values in the diver file. See the
; header for create_orac_lut for a list and description of the these inputs.

; OUTPUT ARGUMENTS:
; None
;
; HISTORY:
; 21/06/13, G Thomas: Original version.
; 18/10/13, G Thomas: Added driver keyword.
; 31/10/13, G Thomas: Added error handler to ensure IDL exits gracefully if code
;                     crashes.
; 05/11/13, G Thomas: Added keyword input for all parameters which can be passed
;                     to create_orac_lut, for easier use with the RAL vm_control
;                     script.
; XX/XX/16, G McGarragh: Add support for new create_orac_lut keywords to the
;    driver file reader and the overwriting keyword list.
; 02/09/16, G McGarragh: Add full stack trace with line numbering to the error
;    hander output. Before handling of the error did not indicate which file/
;    procedure/line it occurred at. Now the output is the same as that without
;    the error handler.
; 01/03/18, G McGarragh: Add parse_multi_line() (which calls parse_line()) so
;    string arrays can be put on multiple lines with the continuation operator
;    '$'.  This is useful for instruments with lots of channels like MODIS.
; 14/04/18, G McGarragh: Add optional inputs srf_quad and srfdat to support
;    integration over channel spectral response functions (SRFs).
; 20/03/20, G Thomas: Add optional input opt_prop_luts
; 17/02/21, RGG: Rewritten. Added gas keyword, Channle ID now an integer not string to be consistent with instrument definition file
;                     

function parse_string_array, word
; Remove leading and trailing spaces
  word = strtrim(word,2)
; remove array brackets
  word = strmid(word,1,strlen(word)-2)
  out = strsplit(word,',',/extract)
  return, out
end

pro create_orac_aerosol_lut_wrapper, driver        = driver, $
                             in_path       = k_in_path, $
                             instfile      = k_instfile, $
                             mmfile        = k_mmfile, $
                             lutfile       = k_lutfile, $
                             out_path      = k_out_path, $
                             atmospheres   = k_atmospheres, $
                             channelID     = k_ChannelID, $
                             force_n       = k_force_n, $
                             force_k       = k_force_k, $
                             mie           = k_mie, $
                             gas           = k_gas, $
                             no_rayleigh   = k_no_rayleigh, $
                             no_screen     = k_no_screen, $
                             n_theta       = k_n_theta, $
                             srf_quad      = k_srf_quad, $  
                             reuse_scat    = k_reuse_scat, $
                             scat_only     = k_scat_only, $
                             srfdat        = k_srfdat, $
                             tmatrix_path  = k_tmatrix_path, $
                             opt_prop_luts = k_opt_prop_luts, $
                             version       = k_version

; Use the driver keyword for the location of the driver file, if it exists,
; otherwise use the environment variable

  if ~keyword_set(driver) then driver = getenv('CREATE_ORAC_AEROSOL_LUT_DRIVER')
  
; create outpath in luts directory based on driver filename 
  words = strsplit(driver,/extract,'/')
  out_path = 'luts/'+(strsplit(words(-1),/extract,'.'))[0]
  
; Load the driver file, which contains all the arguments for the
; create_orac_lut function. Comments can be added after a "#" sign.
; Read in the instrument description
  Lines_Read = strarr(File_lines(driver))
  openr, lun, driver, /get_lun
    readf, lun, lines_read
  close, lun
  free_lun, lun  

; Remove leading and trailing spaces
  lines_read = strtrim(lines_read,2)
  
; Remove comment lines
  First_Character = strmid(lines_read,0,1)
  OK = where(First_CharaCTER ne '#', LineCount)
  lines_ok = lines_read(OK)

; remove trailing comments 
  lines_clean = strarr(LineCount)
  For I = 0,LineCount-1 do Lines_clean[i]= (strsplit(lines_ok[i],'#',/extract))[0]

; Merge continued lines  
  Last_Character = strmid(lines_ok,0,1,/Reverse_Offset)
  Unique = where(Last_CharaCTER Eq '$', Continued_Count)
  Number_of_Lines = LineCount - Continued_Count
  Lines = strarr(Number_of_Lines)
  
  J = 0
  For I = 0, Number_of_Lines -1 do begin
    Current_Line = Lines_clean(J)
    While (strmid(Current_Line,0,1,/Reverse_Offset) Eq '$') do begin
      J= J + 1
      Current_Line =  strmid(Current_Line,0, strlen(Current_Line)-1) + Lines_clean(J)  
    EndWhile
    Lines(I) = Current_line
    J= J + 1
  endfor
   
; Process the lines  
  in_path  = lines[0]
  instfile  = lines[1]
  mmfile   = lines[2]
  lutfile   = lines[3]
  atmospheres = lines[4] 

; Now check for keywords  
  For I = 5, Number_of_Lines-1 do begin
    words = strtrim(strsplit(lines(I),'=',/extract),2)
    first_word = strlowcase(words[0])
    case first_word of
      'channelid'    : channelID    = fix(parse_string_array(words[1]))
      'force_n'      : force_n      = float(words[1])
      'force_k'      : force_k      = float(words[1])          
      'mie'          : mie          = 1
      'gas'          : gas          = 1
      'no_rayleigh'  : no_rayleigh  = 1
      'no_screen'    : no_screen    = 1
      'srf_quad'     : srf_quad     = fix(words[1])
      'n_theta'      : n_theta      = fix(words[1])
      'reuse_scat'   : reuse_scat   = 1
      'scat_only'    : scat_only    = 1
      'srfdat'       : srfdat       = words[1]
      'tmatrix_path' : tmatrix_path = words[1]
      'opt_prop_luts': opt_prop_luts= words[1]
      'version'      : version      = fix(words[1])
      'null'         :
      else           : print,'Keyword ',first_word,' not recognised'
    endcase
  endfor

;  Now check if any _other_ parameters have been passed by keyword. These values
;  will override any driver file settings
   if n_elements(k_in_path)       gt 0 then in_path     = k_in_path
   if n_elements(k_instfile)       gt 0 then instfile     = k_instfile
   if n_elements(k_mmfile)        gt 0 then mmfile      = k_mmfile
   if n_elements(k_lutfile)        gt 0 then lutfile      = k_lutfile
   if n_elements(k_out_path)      gt 0 then out_path    = k_outpath
   if n_elements(k_atmospheres)   gt 0 then atmospheres = k_atmospheres

   if n_elements(k_ChannelID)     gt 0 then channelID     = k_ChannelID
   if n_elements(k_force_n)       gt 0 then force_n       = k_force_n
   if n_elements(k_force_k)       gt 0 then force_k       = k_force_k
   if n_elements(k_mie)           gt 0 then mie           = k_mie
   if n_elements(k_gas)           gt 0 then gas           = k_gas
   if n_elements(k_no_rayleigh)   gt 0 then no_rayleigh   = k_no_rayleigh
   if n_elements(k_no_screen)     gt 0 then no_screen     = k_no_screen
   if n_elements(k_srf_quad)      gt 0 then srf_quad      = k_srf_quad
   if n_elements(k_n_theta)       gt 0 then n_theta       = k_n_theta
   if n_elements(k_reuse_scat)    gt 0 then reuse_scat    = k_reuse_scat
   if n_elements(k_scat_only)     gt 0 then scat_only     = k_scat_only
   if n_elements(k_srfdat)        gt 0 then srfdat        = k_srfdat
   if n_elements(k_tmatrix_path)  gt 0 then tmatrix_path  = k_tmatrix_path
   if n_elements(k_opt_prop_luts) gt 0 then opt_prop_luts = k_opt_prop_luts
   if n_elements(k_version)       gt 0 then version       = k_version


;  Call the create_orac_lut function itself making sure there are no spaces on the directory or file names
   stat = create_orac_aerosol_lut(strtrim(in_path,2),    $
                          strtrim(instfile,2),    $
                          strtrim(mmfile,2),     $
                          strtrim(lutfile,2),     $
                          strtrim(out_path,2),   $
                          strtrim(atmospheres,2),$                   
                          channelID     = channelID,     $
                          force_n       = force_n,       $
                          force_k       = force_k,       $
                          mie           = mie,           $
                          gas           = gas,           $
                          no_rayleigh   = no_rayleigh,   $
                          no_screen     = no_screen,     $
                          n_theta       = n_theta,       $
                          srf_quad      = srf_quad,      $
                          reuse_scat    = reuse_scat,    $
                          scat_only     = scat_only,     $
                          tmatrix_path  = tmatrix_path,  $
                          opt_prop_luts = opt_prop_luts, $
                          version       = version,       $
                          driver        = driver)
                          
;  Check the output status of create_orac_lut
   if stat ne 0 then print,'create_orac_lut failed with code: ', strtrim(stat,2)

   skip_all:
end
