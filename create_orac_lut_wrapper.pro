; procedure create_orac_lut_wrapper
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
;    - channels=['01','02','06','20','31','32']
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
; 14/04/18, G McGarragh: Add optional inputs n_srf_points and srfdat to support
;    integration over channel spectral response functions (SRFs).

function parse_line, lun
   line = '#'
   while (strmid(line,0,1) eq '#') and ~eof(lun) do begin
      readf,lun,line
   endwhile
   if strmid(line,0,1) eq '#' then begin
      return, 'null'
   endif else begin
      return, (strsplit(line,'#',/extract))[0]
   endelse
end


function parse_multi_line, lun
   line = ''
   line2 = parse_line(lun)
   while (strpos(line2, '$') ge 0) and ~eof(lun) do begin
      line = line + (strsplit(line2,'$',/extract))[0]
      line2 = parse_line(lun)
   endwhile
   return, line + line2
end


function parse_string_array, instr
   pos = strsplit(instr,"'")
   if (n_elements(pos) eq 1) and (pos[0] eq 0) then $
      pos = strsplit(instr,'"')
   if (n_elements(pos) eq 1) and (pos[0] eq 0) then return, instr
;  The following line assumes that the final character of the input string is
;  either a ' or ".
   if n_elements(pos) eq 1 then begin
      return, [strmid(instr,pos[0],strlen(instr)-pos[0]-1)]
   endif else if n_elements(pos) eq 2 then begin
      return, [strmid(instr,pos[0],pos[1]-pos[0]-1)]
   endif

   n_els = (n_elements(pos)-1)/2
   out = strarr(n_els)
   j = 0
   for i=1,n_elements(pos)-2,2 do begin
      out[j] = strmid(instr,pos[i],pos[i+1]-pos[i]-1)
      j = j+1
   endfor
   return, out
end


pro create_orac_lut_wrapper, driver=driver, in_path=k_in_path, $
                             instdat=k_instdat, miedat=k_miedat, $
                             lutdat=k_lutdat, presdat=k_presdat, $
                             out_path=k_out_path, channels=k_channels, $
                             force_n=k_force_n, force_k=k_force_k, $
                             gasdat=k_gasdat, mie=k_mie, $
                             no_rayleigh=k_no_rayleigh, $
                             no_screen=k_no_screen, n_theta=k_n_theta, $
                             n_srf_points=k_n_srf_points, $
                             reuse_scat=k_reuse_scat, scat_only=k_scat_only, $
                             srfdat=k_srfdat, tmatrix_path=k_tmatrix_path, $
                             version=k_version

;  Set up an error handler
   catch, error_index
   if error_index ne 0 then begin
      help, /last_message, output=error_message
      print, error_message

;     print, '*** create_orac_lut encountered an error. Index: ', $
;            strtrim(error_index,2)
;     print, '*** Error message: ', !ERROR_STATE.MSG
      print, '*** Aborting run'
      goto, skip_all
   endif

;  Use the driver keyword for the location of the driver file, if it exists,
;  otherwise use the environment variable
   if ~keyword_set(driver) then $
      driver = getenv('CREATE_ORAC_LUT_DRIVER')

;  Load the driver file, which contains all the arguments for the
;  create_orac_lut function. Comments can be added after a "#" sign.
   if n_elements(driver) gt 0 then begin
      openr,lun,driver,/get_lun
      in_path  = parse_line(lun)
      instdat  = parse_line(lun)
      miedat   = parse_line(lun)
      lutdat   = parse_line(lun)
      presdat  = parse_line(lun)
      out_path = parse_line(lun)
;     Now check for keywords
      while ~eof(lun) do begin
         line = parse_multi_line(lun)
         words = strsplit(line,'=',/extract)
         kw = strsplit(words[0],'/',/extract)
         case strlowcase(kw[0]) of
            'channels'     : channels     = parse_string_array(words[1])
            'force_n'      : force_n      = parse_string_array(words[1])
            'force_k'      : force_k      = parse_string_array(words[1])
            'gasdat'       : gasdat       = parse_string_array(words[1])
            'mie'          : mie          = 1
            'no_rayleigh'  : no_rayleigh  = 1
            'no_screen'    : no_screen    = 1
            'n_srf_points' : n_srf_points = uint(words[1])
            'n_theta'      : n_theta      = uint(words[1])
            'reuse_scat'   : reuse_scat   = 1
            'scat_only'    : scat_only    = 1
            'srfdat'       : srfdat       = parse_string_array(words[1])
            'tmatrix_path' : tmatrix_path = parse_string_array(words[1])
            'version'      : version      = parse_string_array(words[1])
            'null'         :
            else           : print,'Keyword ',kw[0],' not recognised'
         endcase
      endwhile
      free_lun, lun
   endif

;  Now check if any _other_ parameters have been passed by keyword. These values
;  will override any driver file settings
   if n_elements(k_in_path)      gt 0 then in_path=k_in_path
   if n_elements(k_instdat)      gt 0 then instdat=k_instdat
   if n_elements(k_miedat)       gt 0 then miedat=k_miedat
   if n_elements(k_lutdat)       gt 0 then lutdat=k_lutdat
   if n_elements(k_presdat)      gt 0 then presdat=k_presdat
   if n_elements(k_out_path)     gt 0 then out_path=k_outpath

   if n_elements(k_channels)     gt 0 then channels=k_channels
   if n_elements(k_force_n)      gt 0 then force_n=k_force_n
   if n_elements(k_force_k)      gt 0 then force_k=k_force_k
   if n_elements(k_gasdat)       gt 0 then gasdat=k_gasdat
   if n_elements(k_mie)          gt 0 then mie=k_mie
   if n_elements(k_no_rayleigh)  gt 0 then no_rayleigh=k_no_rayleigh
   if n_elements(k_no_screen)    gt 0 then no_screen=k_no_screen
   if n_elements(k_n_srf_points) gt 0 then n_srf_points=k_n_srf_points
   if n_elements(k_n_theta)      gt 0 then n_theta=k_n_theta
   if n_elements(k_reuse_scat)   gt 0 then reuse_scat=k_reuse_scat
   if n_elements(k_scat_only)    gt 0 then scat_only=k_scat_only
   if n_elements(k_srfdat)       gt 0 then srfdat=k_srfdat
   if n_elements(k_tmatrix_path) gt 0 then tmatrix_path=k_tmatrix_path
   if n_elements(k_version)      gt 0 then version=k_version

;  Call the create_orac_lut function itself
   stat = create_orac_lut(in_path, instdat, miedat, lutdat, presdat, $
                          out_path, channels=channels, force_n=force_n, $
                          force_k=force_k, gasdat=gasdat, mie=mie, $
                          no_rayleigh=no_rayleigh, no_screen=no_screen, $
                          n_theta=n_theta, n_srf_points=n_srf_points, $
                          reuse_scat=reuse_scat, scat_only=scat_only, $
                          srfdat=srfdat, tmatrix_path=tmatrix_path, $
                          version=version, driver=driver)

;  Check the output status of create_orac_lut
   if stat ne 0 then print,'create_orac_lut failed with code: ', strtrim(stat,2)

   skip_all:
end
