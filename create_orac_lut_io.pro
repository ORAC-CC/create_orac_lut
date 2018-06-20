; This library of contains the procedures for reading the various driver files
; needed by create_orac_lut, as well as the code for writing the ORAC lookup
; tables themselves.
;
; HISTORY:
; 21/06/13, G Thomas: Original version.


; procedure read_lutdat
;
; Read an LUT grid definition file.
;
; File should contain the following lines, with values separated by spaces.
;  (1) Number of aerosol optical depth LUT values; 1 for logged output, 0 for
;      linear
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
; 21/06/13, G Thomas: Original version.

pro read_lutdat, file, LUTstr
   line = ''
   openr, lun, file, /get_lun
   readf, lun, line
   for i=0,4 do begin
;  Skip over any comment lines
      while strmid(line,0,1) eq '#' do readf, lun, line
;     First two data lines should be aerosol optical depth values
      dat1 = strsplit(line,' ',/extract)
      Ndat = fix(dat1[0])
      if dat1[1] ne '#' then Ldat = fix(dat1[1])
      dat = fltarr(Ndat)
      readf, lun, dat
;     Now store the data in the correct output values
      case i of
         0: begin
            NAOD = Ndat
            LAOD = Ldat
            AOD  = dat
         end
         1: begin
            NEfR = Ndat
            LEfR = Ldat
            EfR  = dat
         end
         2: begin
            NSat = Ndat
            Sat  = dat
         end
         3: begin
            NSol = Ndat
            Sol  = dat
         end
         4: begin
            NAzi = Ndat
            Azi  = dat
         end
      endcase
      if ~eof(lun) then readf,lun, line
   endfor
   free_lun, lun

;  Build the output structure
   LUTstr = { NAOD: NAOD, NEFR: NEfR,             $
              NSOL: NSol, NSAT: NSat, NAZI: NAzi, $
              LAOD: LAOD, LEfR: LEfR,             $
              AOD: AOD,   EFR: EfR,               $
              SOL: Sol,   Sat: Sat,   Azi: Azi    }
end


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


; procedure read_miedat
;
; Read a microphysical definition file.
;
; These are the most complex of the driver files and contain the microphysical
; properties of the aerosol/cloud. Values in the file are denoted with ident-
; ifying strings at the start of each line (or block), so the order of entries
; in the file is not important.
;
; The labels are:
;
; lutname: The "long name" of the LUT (eg. A05_Maritime_Clean).
;
; outname: The "short name" of the LUT - this is used as the identifying string
;    in the output filenames of the LUT files themselves (eg. A05).
;
; lutdescription: A description string of the LUT. Allows a human readable
;    description of the particular LUT.
;
; profile: The relative AOD profile with height that will be interpolated onto
;    the pressure grid:
;    4         : Number of layers
;    4.5   0.0 : Layer height in km and relative AOD (AOD units are arbitrary.)
;    3.5   0.0 : "
;    2.5  15.0 : "
;    1.5  85.0 : "
;
; component: Denotes the start of a description of a new component. Each
;    component can be described by an OPAC style data file, Baran or Baum ice
;    crystal scattering properties or in-line in the driver file. Following
;    this line may be a number of data groups described later.
;
;    To use an OPAC component, the word "opac" should appear after "component".
;       followed by the name of the OPAC data file.  In this case the "mixing
;       ratio" and "scattering code" (and "aspect ratios", if T-matrix is being
;       used) data groups must be specified.  Note: Only the OPAC files for
;       aerosol are supported.  OPAC Files for a modified-gamma distribution of
;       cloud particles may be supported in the future if there is demand.
;
;    To use Baran ice crystal scattering properties the word "baran" should
;       appear after "component" followed by the path to the directory with the
;       optical properties. In this case no data groups need to be specified.
;
;    To use Baum ice crystal scattering properties the word "baum" should
;       appear after "component" followed by the full path to the NetCDF file:
;       GeneralHabitMixture_SeverelyRough_AllWavelengths_FullPhaseMatrix.nc,
;       and, optionally, if the Baum instrument specific scattering properties
;       are to be used, followed by the full path to the instrument specific
;       NetCDF file containing the optical properties for GeneralHabitMixture/
;       SeverelyRough for each channel.  In this case no data groups need to be
;       specified.
;
;    For user defined components, the word "user" should follow "component", as
;       well as an (optional) component name. In this case the "mixing ratio",
;       "scattering code", "size", "refractive index" (and "aspect ratios", if
;       T-matrix is being used) data groups must be specified.
;
;    The possible data groups denoted by a label beginning with "* " are:
;
;       * mixing ratio: A single value on the following line provides the
;            number-mixing ratio of the component in the class as a whole.
;
;       * scattering code: A value of "Mie" or "tmatrix" on the following line
;            defines which scattering code should be used for this component.
;
;       * size: The name of the size distribution on the following line followed
;            by a line of one or more size distribution parameters separated by
;            spaces. For 'log_normal' two values are required: the mode radius
;            and the log-normal spread (defined such that sigma(ln r) = ln S).
;            For 'modified_gamma' (according to Hansen and Travis 1974) the
;            values are: 'a', 'b', minimum radius and maximum radius.
;
;       * refractive index: In a format similar to the height profile, the
;            following lines should give the number of RI values, followed by
;            three columns giving the wavelength in microns, the real, and the
;            imaginary refractive index values. Again, these will be
;            interpolated to the require output wavelengths.
;
;       * aspect ratio: Again, similarly to the height profile, the number of
;            aspect ratio values, followed by two columns giving the aspect
;            ratios and their relative numbers. This is only needed if T-matrix
;            scattering is being used.
;
; end: Denotes the end of the driver file.
;
; INPUT ARGUMENTS:
; file (string) Path to file.
;
; INPUT KEYWORDS:
; None
;
; OUTPUT ARGUMENTS:
; scatstr (structure) Structure with aerosol/cloud microphysical information.
;
; HISTORY:
; 21/06/13, G Thomas: Original version.
; XX/XX/15, G McGarragh: Add a line to the 'size' data group to indicate the
;    name of the size distribution and generalize the input size distribution
;    parameters as a vector with length depending on the distribution type. See
;    the documentation for details.
; XX/XX/15, G McGarragh: The word "opac" is now required after "component" to
;    indicate that the next token will be an OPAC data file.
; XX/XX/15, G McGarragh: Add support for ice crystal components from either the
;    Bryan Baum or Anthony Baran datasets. See the documentation for details.

pro read_miedat, file, scat

   line = ''
;  Define some variables used in reading the file
   NComp       = 0
   NWl         = 0
   NAR         = 0
   NAlt        = 0
   outname     = ''
   lutname     = ''
   description = ''
   code        = ''
   MRat        = [1.0] & MRat1      = 0.0
   distname    = ['']  & distname1  = ''
   comptype    = ['']  & comptype1  = ''
   compname    = ['']  & compname1  = ''
   compname2   = ['']  & compname21 = ''
   Rm          = [0.0] & Rm1        = 0.0
   S           = [0.0] & S1         = 0.0

;  Define a dummy output structure which we populate below
   scat = { lutname:     '', $
            outname:     '', $
            description: '', $
            NComp:        0  }

   openr,lun, file, /get_lun

;  Comment lines start with "#", data descriptors start with "*"
   readf,lun, line
   while strmid(line,0,1) eq '#' do readf,lun, line

;  Extract the first word from the non-comment lines
   while strlowcase(strtrim(line,2)) ne 'end' do begin
      chunks = strsplit(line,' ',/extract)
      case strlowcase(chunks[0]) of
         'component': begin
            NComp = NComp+1 ; Update number of components
            code = ''
            has_ri = 0      ; Flag for the existence of refractive index
            has_as = 0      ; Flag for the existence of aspect ratio

            if n_elements(chunks) eq 1 then begin $
               message,'Scattering file components must specifiy the component' + $
                       'type as the second field'
            endif
            comptype1  = chunks[1]

            case strlowcase(comptype1) of
               'user': begin
                  if n_elements(chunks) lt 3 then begin $
                     message,"Scattering file component type 'user' " + $
                             "requires a component name as its first+ field."
                  endif
                  compname1=strjoin(chunks[2:*],' ',/single)
               end
               'opac': begin
                  if n_elements(chunks) ne 3 then begin $
                     message,"Scattering file component type 'opac' " + $
                             "requires an OPAC OptDat file as its first field."
                  endif
                  compname1=strjoin(chunks[2:2],' ',/single)
               end
               'baran': begin
                  if n_elements(chunks) ne 3 then begin $
                     message,"Scattering file component type 'baran' " + $
                             "requires the path to the directory with the " + $
                             "optical properties as its first field."
                  endif
                  compname1=strjoin(chunks[2:2],' ',/single)
               end
               'baum': begin
                  if n_elements(chunks) ne 3 and n_elements(chunks) ne 4 then begin $
                     message,"Scattering file component type 'baum' " + $
                             "requires the path to the generic optical " + $
                             "properties file as field one optionally " + $
                             "followed by the path to the imager specific " + $
                             "file as field 2."
                  endif
                  compname1=strjoin(chunks[2:2],' ',/single)
                  compname21 = ''
                  if n_elements(chunks) eq 4 then begin $
                     compname21=strjoin(chunks[3:3],' ',/single)
                  endif
               end
               else: begin
                  message,'Invalid component type: ' + comptype1
               end
            endcase ; End of component data case statement

;           If the component ID isn't "user" we assume it is an OPAC/GADS
;           component and load the data from the database.
            if strlowcase(comptype1) eq 'opac' then begin
               has_ri = 1
               rd_optdat_p,compname1,distname1,rmd,rm1,rml,rmh,s1,wl,Be,Bs,Ba, $
                           g,w,n,k,ang,phase
               distname1 = 'log_normal'
            endif
            readf,lun, line

;           Now we inspect each data descriptor and read the data
;           accordingly.
            while strmid(line,0,1) eq '*' do begin
               case strtrim(strlowcase(strmid(line,1)),2) of
                  'size': begin
                     readf,lun, distname1
                     readf,lun, Rm1, S1
                     distname1 = strtrim(distname1,2)
                  end
                  'scattering code': begin
                     readf,lun, Code
                     Code = strtrim(Code,2)
                  end
                  'mixing ratio': readf,lun, MRat1
                  'refractive index': begin
                     has_ri = 1
                     readf,lun, NWl
                     wl = fltarr(NWl)
                     n  = fltarr(NWl)
                     k  = fltarr(NWl)
                     row = fltarr(3)
                     for i=0,NWl-1 do begin
                        readf,lun, row
                        wl[i] = row[0]
                        n[i] = row[1]
                        k[i] = row[2]
                     endfor
                  end
                  'aspect ratio': begin
                     has_as = 1
                     readf,lun, NAR
                     eps  = fltarr(NAR)
                     neps = fltarr(NAR)
                     row = fltarr(2)
                     for i=0,NAR-1 do begin
                        readf,lun, row
                        eps[i]  = row[0]
                        neps[i] = row[1]
                     endfor
                  end
               endcase ; End of component 'user' data case statement
               if ~eof(lun) then readf,lun,line else break
            endwhile ; End of loop for reading each component

;           The size parameters and mixing ratios of each component are stored
;           in simple vectors, which makes them easy to reference in the main
;           routine.
            if NComp eq 1 then begin
               MRat[0]     = MRat1
               distname [0] = distname1
               comptype [0] = comptype1
               compname [0] = compname1
               compname2[0] = compname21
               Rm[0]       = Rm1
               S[0]        = S1
            endif else begin
               MRat      = [MRat,     MRat1]
               distname  = [distname,  distname1]
               comptype  = [comptype,  comptype1]
               compname  = [compname,  compname1]
               compname2 = [compname2, compname21]
               Rm        = [Rm,        Rm1]
               S         = [S,         S1]
            endelse

;           The refractive index, scattering code label and (optionally)
;           asymmetry data for each component are stored in their own sub-
;           structure. (This is because the refractive index and asymmetry
;           could have different numbers of elements for different components).
            tmp = {code : code}
            if has_ri then tmp = create_struct(tmp, 'WL', WL, $
                                                    'CM', complex(n,k))
            if has_as then tmp = create_struct(tmp, 'eps', eps, 'neps', neps)
            scat = create_struct(scat, 'Comp'+strtrim(NComp,2), tmp)
         end ; End of component case
         'profile': begin
            readf,lun,NAlt
            H   = fltarr(NAlt)
            REx = fltarr(NAlt)
            row = fltarr(2)
            for i=0,NAlt-1 do begin
               readf,lun, row
               H[i] = row[0]
               REx[i] = row[1]
            endfor
            scat = create_struct(scat, 'NLayer', NAlt, 'height', H, 'RExt', REx)
            if ~eof(lun) then readf,lun, line else break
         end ; End of profile case
         'outname': begin
            scat.outname = strjoin(chunks[1:*],'_',/single)
            if ~eof(lun) then readf,lun, line else break
         end
         'lutname': begin
            scat.lutname = strjoin(chunks[1:*],'_',/single)
            if ~eof(lun) then readf,lun, line else break
         end
         'lutdescription': begin
            scat.description = strjoin(chunks[1:*],' ',/single)
            if ~eof(lun) then readf,lun, line else break
         end
;        If we don't recognise the current line, break out of the case
;        statement.
         else: begin
            message,/info, 'Warning: unknown label line found and skipped: ' + $
                           line
            readf,lun, line
            break
         end
      endcase	; End of class case statement
   endwhile	; End of main while loop

;  We've found the end of the data file
   free_lun, lun

;  Overwrite the dummy NComp value in the output structure with the actual
;  value
   scat.NComp = NComp

;  Add the component mixing ratios, size distribution name, and the required
;  size distribution parameters to the output structure
   scat = create_struct(scat, 'MRat', MRat, 'distname', distname, $
                        'comptype', comptype, 'compname', compname, $
                        'compname2', compname2, 'Rm', Rm, 'S', S)

;  If the 'outname' wasn't defined, set to the be same as 'lutname'
   if scat.outname eq '' then begin
      if scat.lutname eq '' then $
         message, 'Either "outname" or "lutname" must be defined'
      scat.outname = scat.lutname
   endif else if lutname eq '' then scat.lutname = scat.outname
end


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
               H: H, P: P, T: T }
end


; procedure read_gasdat
;
; Read a gas optical depth profile file.
;
; The gas optical depth files can start with comment lines beginning with a "#"
; character. This is then followed by two possible labels:
;    '* instrument': The following line contains the instrument name for this
;                    file.
;    '* channel':    The following line contains the channel number  for this
;                    file.
; These should be followed by two columns giving the height and gas optical
; depth for each layer (as defined by the pressure/height data file). Note that
; the height grid for pressure and gas optical depth have to match.
;
; INPUT ARGUMENTS:
; file (string) Path to file.
;
; INPUT KEYWORDS:
; None
;
; OUTPUT ARGUMENTS:
; gasstr (structure) Structure with gas optical depth profile information.
;
; HISTORY:
; 21/06/13, G Thomas: Original version.

pro read_gasdat, files, gasstr

   line = ''
   ngasdat = n_elements(files)
   gasstr = {Nchan: ngasdat}
   for f=0,ngasdat-1 do begin
      openr,lun,files[f],/get_lun
      readf,lun, line
      while strmid(line,0,1) eq '#' do readf,lun, line
      while strmid(line,0,1) eq '*' do begin
         case strlowcase(strtrim(strmid(line,1),2)) of
            'instrument': begin
               readf,lun, line
               inst = strtrim(line,2)
            end
            'channel': begin
               readf,lun, line
               Chan = fix(line)
            end
         endcase
         readf,lun, line
      endwhile
      row = float(strsplit(line,' ',/extract))
      H   = row[0]
      OPD = row[1]
      row = fltarr(2)
      while ~EOF(lun) do begin
         readf,lun, row
         H   = [H,row[0]]
         OPD = [OPD,row[1]]
      endwhile
      free_lun, lun

      if f eq 0 then gasstr = create_struct(gasstr, 'H', H) $
      else if ~array_equal(H, gasstr.H) then $
         message, 'Height profile of file '+files[f]+' is different from ' + $
                  files[0]
      ChanName = 'C'+string(Chan,format='(i02)')
      Chanstr = {Chan: Chan, Inst: Inst, OPD: OPD}
      gasstr = create_struct(gasstr, ChanName, Chanstr)
   endfor
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
