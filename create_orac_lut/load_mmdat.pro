
; procedure read_mmdat
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
; substance: The substance being modelled (eg. liquid-water, water-ice).
; shortname: short microphysical code number starting
;            a - aerosol
;            c - cloud
;            v - vocanic ash  
; description: A description string of the microphysical model.
;
; profile: The relative optical depth  profile with height that will be interpolated onto
;    the pressure grid:
;    4         : Number of layers
;    4.5   0.0 : Layer height in km and relative optical depth (optical depth units are arbitrary.)
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
;       * refractive index: It is followed by a file of type ri that is issumed to live in the input_files/ri directory.
;            The values will be interpolated to the require output wavelengths.
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
; mmstr (structure) Structure with aerosol/cloud microphysical information.
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

pro load_mmdat, file, mmstr

   line = ''
;  Define some variables used in reading the file
   NComp       = 0
   NWl         = 0
   NAR         = 0
   NAlt        = 0
   substance     = ''
   shortname     = ''
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
   mmstr = {substance:   '', $
            shortname:   '', $
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
               message,'Scattering file components must specify the component' + $
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
                     rifilename = ' '
                     readf,lun, rifilename
                     ri = read_ri('input_files/ri/'+ rifilename)
                     nwl = ri.Vals
                     wl=ri.Wavl
                     n=ri.n  
                     k=-ri.k
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
            mmstr = create_struct(mmstr, 'Comp'+strtrim(NComp,2), tmp)
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
            mmstr = create_struct(mmstr, 'NLayer', NAlt, 'height', H, 'RExt', REx)
            if ~eof(lun) then readf,lun, line else break
         end ; End of profile case
         'substance': begin
            mmstr.substance = strjoin(chunks[1:*],'_',/single)
            if ~eof(lun) then readf,lun, line else break
         end
         'description': begin
            mmstr.description = strjoin(chunks[1:*],' ',/single)
            if ~eof(lun) then readf,lun, line else break
         end
         'shortname': begin
            mmstr.shortname = strjoin(chunks[1:*],' ',/single)
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
   mmstr.NComp = NComp

;  Add the component mixing ratios, size distribution name, and the required
;  size distribution parameters to the output structure
   mmstr = create_struct(mmstr, 'MRat', MRat, 'distname', distname, $
                        'comptype', comptype, 'compname', compname, $
                        'compname2', compname2, 'Rm', Rm, 'S', S)

end
