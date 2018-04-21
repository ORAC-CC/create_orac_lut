; function create_orac_lut
;
; A routine for generating ORAC look-up tables (LUTs) from aerosol/cloud optical
; properties. The code does the following steps:
; * Reads in a series of driver text files
; * Performs the necessary scattering calculations to produce scattering
;   properties for the given aerosol/cloud class
; * Calls DISORT to perform the radiative transfer calculations
; * Writes the LUT "static application data" (SAD) files needed by ORAC
;
; INPUT ARGUMENTS:
; driver_path (string) The path in which the code can find the driver files
; instdat     (string) The name of the instrument data driver file to use
; miedat      (string) The name of the aerosol/cloud optical/microphysical
;                      properties driver file to use
; lutdat      (string) The name of the LUT dimension and vertex driver file to
;                      use
; presdat     (string) The name of the atmospheric pressure profile driver file
;                      to use
; out_path    (string) The directory in which to save the generated LUTs
;
; INPUT KEYWORDS:
; channels=strarr      Provide a list of channel numbers to produce LUTs for.
;                      These must be a subset of the channel numbers listed in
;                      the instdat file. (By default the code will produce LUTs
;                      for all the channels defined in instdat).
; force_n=strarr       Allows the user to set a fixed value for the real part
;                      of the refractive index for 0.55um (element 0) and for
;                      each channel (elements 1 to n channels). As an array of
;                      strings if the value for a particular channel is ' '
;                      then the value is not changed for that channel.
; force_k=strarr       Allows the user to set a fixed value for the imaginary
;                      part of the refractive index for 0.55um (element 0) and
;                      for each channel (elements 1 to n channels). As an array
;                      of strings if the value for a particular channel is ' '
;                      then the value is not changed for that channel.
; gasdat=strarr        The base of optional files containing the optical depth
;                      due to trace gas absorption for all/some instrument
;                      channels. File matching the path: driver_path/gasdat*
;                      will be read.
; /mie                 Force Mie calculation, overriding specification in
;                      miedat file.
; /no_rayleigh         Do not include Rayleigh scattering.
; /no_screen           By default the code uses string(13b) to prevent new-
;                      lines at the end of print-to-screen commands (for
;                      aesthetic reasons): these make a real mess if output is
;                      piped to a file. Setting this keyword suppresses these
;                      control characters
; n_theta=integer      Set the number of angles at which to calculate the phase
;                      function moments. Defaults to 1000.
; n_srf_points=integer Set the number of points for the integration across
;                      spectral response function.
; /reuse_scat          Reuse the scattering computations from a previous run
;                      stored in the file scatfile.sav created with the IDL
;                      save procedure. This file is saved whenever /reuse_scat
;                      is not set. And may only be used by a run with the same
;                      microphysical configuration and set of channels.
; /scat_only           Only compute the scattering properties skipping over
;                      the radiative transfer calculations. In this case the
;                      file scatfile.sav and the Bext, BextRat, w and g LUTS
;                      will be output whereas the reflection, transmission and
;                      emission LUTS will not.
; srfdat=strarr        The files containing the spectral response functions for
;                      each channel.  Requires option n_srf_points to be set.
; tmatrix_path=string  Path to the Dubovik T-Matrix LUT base directory.
;                      Required only when T-Matrix calculations are to to be
;                      made.
; version=string       Set a version string to place in the output LUT file
;                      names (defaults to nothing).
;
; RETURN VALUE:
; At the moment, if the code finishes successfully, the value "0" will always be
; returned. Error/warning codes are a future option.
;
; OUTPUT ARGUMENTS:
; None
;
; HISTORY:
; 05/02/13, G Thomas: Original version.
; 08/02/13, G Thomas: Added channels keyword. Added code to make use of the
;                     log flag for AOD and EfR in lutdat files.
; 19/02/13, G Thomas: Debugging completed. Added check on solar component for
;                     calculation of direct-beam RT.
; 01/05/13, G Thomas: All scattering parameters calculated at the 550nm
;                     reference wavelength are now stored.
;                     Scattering parameters output in the form needed by RAL
;                     LUT code.
; 07/05/13, G Thomas: Added mie keyword.
; 10/05/13, G Thomas: Changed emissivity calculation so that it is done by
;                     DISORT rather than estimated from SSA.
; 17/05/13, G Thomas: Added the force_k keyword.
; 30/05/13, G Thomas: Added the no_screen keyword.
; 26/06/13, G Thomas: Altered the wavelength threshold for the force_k keyword
;                     to 2.0 microns, for producing LUTs for the SMASH project.
; XX/XX/15, G McGarragh: Add support for a gamma distribution as a single mode.
;    Useful for liquid water cloud.
; XX/XX/15, G McGarragh: Add support for a single mode of ice crystals according
;    Braun Baum et al.
; XX/XX/15, G McGarragh: Add support for a single mode of ice crystals according
;    Anthony Baran et al.
; XX/XX/15, G McGarragh: Change Bext output to be that actual Bext and not the
;    ratio with that of the reference wavelength.  Put the actual ratio output
;    into the BextRat LUT.  Output single scattering albedo into the 'w' LUT,
;    asymmetry parameter into the 'g' LUT and the average volume per particle
;    into the 'Vavg' LUT.
; XX/XX/15, G McGarragh: Add the /scat_only keyword to skip the time consuming
;    RT calculations if the scattering output is all that is desired.
; XX/XX/15, G McGarragh: Add the /reuse_scat keyword to reuse the scattering
;    calculations from a previous run.  See the documentation for details.
; XX/XX/15, G McGarragh: The force_k keyword is now a vector where the first
;    element is for the reference wavelength and the other elements are for the
;    channels.  Also added a corresponding force_n keyword.
; 08/06/16, G McGarragh: Add the /no_rayleigh keyword for the option to not
;    include Rayleigh scattering.  Useful for building LUTs for multilayer
;    algorithms.
; 31/08/16, G McGarragh: Add the tmatrix_path keyword as the path to the Dubovik
;    T-Matrix LUTs was fixed before.
; 12/10/16, G McGarragh: The phase functions interpolated from the Baran and
;    Baum ice crystal scattering properties must be normalized.
; 20/10/16, G McGarragh: Write the svn version and termination timestamp to
;    separate files.
; 21/10/16, G McGarragh: Copy the driver files to the output directory using the
;    LUT output base name.
; 14/04/18, G McGarragh: Add support for integration over channel spectral
;    response functions (SRFs).


; Include the libraries of procedures for reading and writing driver files and
; Look-Up Tables; setting up and running the scattering code; and calling the
; DISORT DLM.
; NOTE: These are assumed to be in the same directory as the main function.
@create_orac_lut_io.pro
@create_orac_lut_scat.pro
@create_orac_lut_disort.pro
@read_baran.pro
@read_baum.pro

; Begin the main LUT generation function
function create_orac_lut, driver_path, instdat, miedat, lutdat, presdat, $
                          out_path, channels=channels, force_n=force_n, $
                          force_k=force_k, gasdat=gasdat, mie=mie, $
                          no_rayleigh=no_rayleigh, no_screen=no_screen, $
                          n_theta=n_theta, n_srf_points=n_srf_points, $
                          reuse_scat=reuse_scat, scat_only=scat_only, $
                          srfdat=srfdat, tmatrix_path=tmatrix_path, $
                          version=version, driver=driver

;  -----------------------------------------------------------------------------
;  Test input and output files and directories
;  -----------------------------------------------------------------------------

   print,'Reading input files...'

;  Check for the existence of the driver file path and the required files
   ok = file_test(driver_path, /directory, /read, /write)
   if not ok then message, 'Driver_path not found, or is read/write protected'

;  Instrument definition file
   ok = file_test(instdat, /read)
   if not ok then message, 'instdat file not readable: '+instdat

;  Microphysical definition file
   ok = file_test(miedat, /read)
   if not ok then message, 'miedat file not readable: '+miedat

;  LUT grid definition file
   ok = file_test(lutdat, /read)
   if not ok then message, 'lutdat file not readable: '+lutdat

;  Pressure profile definition file
   ok = file_test(presdat, /read)
   if not ok then message, 'presdat file not readable: '+presdat

;  Gas optical depth profile files, if they exist
   ngasdat = n_elements(gasdat)
   if ngasdat gt 0 then begin
      gasdats = driver_path+'/'+gasdat
      for i=0,ngasdat-1 do begin
         ok = file_test(gasdats[i], /read)
         if not ok then message, 'gasdat file not readable: '+gasdats[i]
      endfor
   endif

;  Channel spectral response functions, if they exist
   nsrfdat = n_elements(srfdat)
   if nsrfdat gt 0 then begin
      srfdats = srfdat
      for i=0,nsrfdat-1 do begin
         ok = file_test(srfdats[i], /read)
         if not ok then message, 'srfdat file not readable: '+srfdats[i]
      endfor
   endif

;  Finally, check that the output directory exists
   ok = file_test(out_path, /directory, /read, /write)
   if not ok then message, 'Out_path not found, or is read/write protected'


;  -----------------------------------------------------------------------------
;  Read in files
;  -----------------------------------------------------------------------------

;  **** Read the LUT parameters file
   read_lutdat, lutdat, lutstr

;  **** Read the Instrument parameters file
   read_instdat, instdat, inststr

;  If the channels keyword has been set, then alter the inststr data
;  accordingly, keeping only those channels listed by the keyword.
   if n_elements(channels) gt 0 then begin
      match = replicate(0, inststr.NChan)
      for i=0,n_elements(channels)-1 do begin
         matchch = where(inststr.ChanNum eq channels[i])
         if matchch[0] lt 0 then $
            message,/info, 'Warning: channel '+strtrim(channels[i])+ $
                           ' not found in '+instdat $
         else match[matchch[0]] = 1
      endfor
      matchi = where(match,matchn)
      inststr = create_struct('name',    inststr.name, $
                              'NChan',   matchn, $
                              'ChanNum', inststr.ChanNum[matchi], $
                              'ChanWl',  inststr.ChanWl[matchi],  $
                              'ChanEm',  inststr.ChanEm[matchi],  $
                              'ChanSol', inststr.ChanSol[matchi])
   endif

;  Now set:
;  NSRF:      Number of SRF integration points
;  ChanSRFWl: 2-d array (NSRF, matchn) of the wavelengths for each integration
;             point for each channel
;  ChanSRF:   2-d array (NSRF, matchn) of the spectral response function values
;             for each integration point for each channel

;  For the default, when no SRF files are provided, one integration point
;  (monochomatic) is set.  This will produce results exactly the same as before
;  this code got SRF support.
   if nsrfdat eq 0 then begin
      NSRF = 1
      ChanSRFWl = reform(inststr.ChanWl, NSRF, matchn)
      ChanSRF   = replicate(1, 1, matchn)
;  Now in the case the SRF files are provided but we must check that other
;  required information is provided.
   endif else begin
      if not n_elements(n_srf_points) then begin
         message, 'Must specify n_srf_points if spectral response files are ' + $
                  'provided.'
      endif
      if n_srf_points mod 2 eq 0 then begin
         message, 'n_srf_points must be an odd number in order to resolve the ' + $
                  'central wavelength.'
      endif

      NSRF = n_srf_points
      ChanSRFWl = fltarr(NSRF, matchn)
      ChanSRF   = fltarr(NSRF, matchn)
      for i=0,matchn-1 do begin
         read_srfdat, srfdats[i], srfstr
         ChanWl1 = srfstr.wl[0]
         ChanWl2 = srfstr.wl[n_elements(srfstr.wl) - 1]
         ChanSRFWl[*,i] = ChanWl1 + findgen(NSRF) * $
                          (ChanWl2 - ChanWl1) / (NSRF - 1)
         ChanSRF  [*,i] = interpol(srfstr.srf, srfstr.wl, ChanSRFWl[*,i])
      endfor
   endelse

;  **** Read the scattering parameters file
   read_miedat, miedat, scatstr

;  **** Read the pressure profile file
   read_presdat, presdat, presstr

;  **** Read the GAS OPD files
   if ngasdat gt 0 then begin
      read_gasdat, gasdats, gasstr
;     Check that the Gas OPD values are on the same grid as the pressure
;     profile
      if ~array_equal(presstr.H, gasstr.H) then $
         message, 'The pressure and gas OPD profiles must be on the same ' + $
                  'altitude grid'
   endif
   print,'LUT calculation for class ',scatstr.lutname,' for the ', $
         inststr.name,' instrument'
   print,'Dimensions are:'
   print,'     Class components ', strtrim(scatstr.NComp,2)
   print,'  Instrument channels ', strtrim(inststr.NChan,2)
   print,'        Optical depth ', strtrim(lutstr.NAOD,2)
   print,'     Effective radius ', strtrim(lutstr.NEfR,2)
   print,'         Solar zenith ', strtrim(lutstr.NSol,2)
   print,'    Instrument zenith ', strtrim(lutstr.NSat,2)
   print,'     Relative azimuth ', strtrim(lutstr.NAzi,2)


;  -----------------------------------------------------------------------------
;  Some miscellaneous setup.
;  -----------------------------------------------------------------------------

;  If the no_screen keyword has been set, we suppress the control character
;  used to prevent a new-line for some print statements (see below). This is
;  useful if the output is being piped into a file, for example (where the
;  control character just messes up the formatting).
   if keyword_set(no_screen) then begin
      newlinechar = ' '
      newlinestr  = ')'
   endif else begin
      newlinechar = string(13b)
      newlinestr  = ',$)'
   endelse

;  If we have a version number, incorporate it into the output filenames
   if n_elements(version) gt 0 then verstrng = '_v'+version $
   else verstrng = ''


;  Generate the output file name pattern
   lutbase = out_path+'/'+strupcase(inststr.name)+'_'+scatstr.outname

;  Check the size of the channel string needed...
   Chfmt = '(i0)'
;  if max(inststr.ChanNum) ge 100 then Chfmt = '(i03)' $
;  else Chfmt = '(i02)'


;  -----------------------------------------------------------------------------
;  Output svn version and copy the driver files to the output directory using
;  the LUT output base name.
;  -----------------------------------------------------------------------------
   spawn, '; svnversion ' + file_dirname((routine_info('create_orac_lut', $
          /function, /source)).path) + ' > ' + lutbase + '_svn_version.txt'

   file_copy, driver,  lutbase + '.driver', /overwrite
   file_copy, instdat, lutbase + '.inst',   /overwrite
   file_copy, miedat,  lutbase + '.mie',    /overwrite
   file_copy, lutdat,  lutbase + '.lut',    /overwrite
   file_copy, presdat, lutbase + '.pres',   /overwrite
   for i=0,ngasdat-1 do begin
      ChStrng = 'Ch'+string(inststr.ChanNum[i], format=Chfmt)
      file_copy, gasdats[i], lutbase + '_' + ChStrng + '.gas', /overwrite
   endfor


;  -----------------------------------------------------------------------------
;  Interpolate the aerosol profile layers onto the atmos. pressure and gas OPD
;  layers, and the aerosol refractive index onto the channel wavelengths.
;  -----------------------------------------------------------------------------

;  Firstly, we have to define the height of the layers, which lie between each
;  pressure level...
   NLayers = presstr.NLevels -1
   HLayers = (presstr.H[0:NLayers-1] + presstr.H[1:NLayers]) / 2.0
   AerRelTau = interpol(scatstr.RExt, scatstr.Height, HLayers)
   AerRelTau = AerRelTau/total(AerRelTau)

;  Interpolate the components refractive index values onto the 0.55 micron
;  reference wavelength and the instrument channels.
   AerM550 = complexarr(scatstr.NComp)
   AerM = complexarr(NSRF, inststr.NChan, scatstr.NComp)

;  Which member of the scatstr structure is the first component with refractive
;  index/asymmetry information substructure?
   scatoffset = where(tag_names(scatstr) eq 'COMP1')
   scatoffset = scatoffset[0]
   if scatoffset eq -1 then message, "Can't find COMP1 substructure in scatstr"

   if scatstr.(scatoffset).code ne '' then begin
      for c=0,scatstr.NComp-1 do begin
         cc = scatoffset + c
         AerM550[c]  = interpol(scatstr.(cc).Cm, scatstr.(cc).wl, 0.55)
         AerM[*,*,c] = interpol(scatstr.(cc).Cm, scatstr.(cc).wl, ChanSRFWl)
      endfor

;     If the force_n keyword has been specified, replace real RI values with
;     those given in force_n.
      if n_elements(force_n) gt 0 then begin
         if force_n[0] ne ' ' then begin
            for c=0,scatstr.NComp-1 do begin
               AerM550[c] = AerM550[c] - complex(float(AerM550[c]), 0.0) + $
                            complex(float(force_n[0]), 0.0)
            endfor
         endif

         for i=0,inststr.NChan-1 do begin
            if force_n[i+1] ne ' ' then begin
               for c=0,scatstr.NComp-1 do begin
                  AerM[*,i,c] = AerM[*,i,c] - complex(float(AerM[*,i,c]), 0.0) + $
                                complex(float(force_n[i+1]), 0.0)
               endfor
            endif
         endfor
      endif

;     If the force_k keyword has been specified, replace imaginary RI values
;     with those given in force_k.
      if n_elements(force_k) gt 0 then begin
         if force_k[0] ne ' ' then begin
            for c=0,scatstr.NComp-1 do begin
               AerM550[c] = AerM550[c] - complex(0.0,imaginary(AerM550[c])) + $
                            complex(0.0,float(force_k[0]))
            endfor
         endif

         for i=0,inststr.NChan-1 do begin
            if force_k[i+1] ne ' ' then begin
               for c=0,scatstr.NComp-1 do begin
                  AerM[*,i,c] = AerM[*,i,c] - complex(0.0,imaginary(AerM[*,i,c])) + $
                                complex(0.0,float(force_k[i+1]))
               endfor
            endif
         endfor
      endif
   endif


;  -----------------------------------------------------------------------------
;  Determine scattering properties either through calling a scattering code for
;  Mie or t-matrix or loading 'baum' or 'buran' or ice crystal properties, or by
;  reloading properties saved from a previous run.
;  -----------------------------------------------------------------------------

   if keyword_set(reuse_scat) then begin

;     **** Read the scattering parameters for the class as a whole for
;          reuse
      restore, out_path+'/scatfile.sav'

   endif else begin

;     **** Generate the range of component mixing-ratios/mode radii required to
;          provide the required effective radii

      if scatstr.comptype[0] eq 'opac' or scatstr.comptype[0] eq 'user' then begin
         if scatstr.distname[0] eq 'log_normal' then begin
;           create_range returns (NComp, NEfR) arrays of Mixing ratio and Rm
            create_range, scatstr.MRat, scatstr.Rm, scatstr.S, lutstr.EfR, $
                          lut_MRat, lut_Rm
         endif else if scatstr.distname[0] eq 'modified_gamma' then begin
            lut_MRat = Dblarr(N_elements(scatstr.MRat),N_elements(lutstr.EfR))
            lut_Rm   = Dblarr(N_elements(scatstr.MRat),N_elements(lutstr.EfR))

            lut_MRat[0,*] = scatstr.MRat
            lut_Rm  [0,*] = lutstr.EfR
         endif
      endif else begin
         lut_MRat = Dblarr(N_elements(scatstr.MRat),N_elements(lutstr.EfR))

         lut_MRat[0,*] = scatstr.MRat
      endelse

;     **** Generate the quadrature points for the scattering phase function

;     Check if the NMom keyword has been set, if it hasn't we use the default
;     value of 1000.
      if n_elements(n_theta) eq 0 then begin
         NMom = 1000
;        x = 2. * !pi * 240. / .47;
;        NMom = fix(2 * (x + 4.05 * x^(1./3.) + 8))
      endif else begin
         NMom = n_theta
      endelse

;     The quadrature procedure gives us our phase function angles
      quadrature, 'g', NMom, Abscissas, Weights
;     Note that QV = cos(scattering_angle)
      QV0 =  1.0 ; theta = 0
      QV1 = -1.0 ; theta = 180
      QV = ((QV1-QV0)*Abscissas + (QV0+QV1)) / 2d0
      PTheta = acos(QV)

;     **** Call the scattering code for the required range of components and
;          mode radii

;     Define the arrays which hold the scattering parameters
      Vavg_c    = fltarr(scatstr.NComp, lutstr.NEfR) ; Average volume per particle

      ; At the reference wavelength
      Bext550_c = fltarr(scatstr.NComp, lutstr.NEfR) ; Extinction coefficient
      w550_c    = fltarr(scatstr.NComp, lutstr.NEfR) ; Single scattering albedo
      g550_c    = fltarr(scatstr.NComp, lutstr.NEfR) ; Asymmetry parameter
      Phs550_c  = fltarr(NMom, scatstr.NComp, lutstr.NEfR) ; Phase function

      ; At the individual channels
      Bext_c    = fltarr(NSRF, inststr.NChan, scatstr.NComp, lutstr.NEfR)
      w_c       = fltarr(NSRF, inststr.NChan, scatstr.NComp, lutstr.NEfR)
      g_c       = fltarr(NSRF, inststr.NChan, scatstr.NComp, lutstr.NEfR)
      Phs_c     = fltarr(NMom, NSRF, inststr.NChan, scatstr.NComp, lutstr.NEfR)

      for c=0,scatstr.NComp-1 do begin
         if scatstr.comptype[0] eq 'opac' or scatstr.comptype[0] eq 'user' then begin
            for r=0,lutstr.NEfR-1 do begin
               print,' Performing scattering calculations for component ' + $
                     scatstr.compname[c] + ', EfR: ', lutstr.EfR[r], newlinechar, $
                     format='(A,f8.4,A'+newlinestr

;              Note that the mode radii of each component don't usually change
;              from one effective radius to the next (it's the mixing ratio that
;              changes). Here we check if the component has changed size and
;              only call the scattering code if it has.
               calculated = 0
               if r gt 0 then begin
                  if (lut_Rm[c,r] eq lut_Rm[c,r-1]) and $
                     (Bext_c[0,0,c,r-1] ne 0) then calculated = 1
               endif
               if calculated then begin
                  Vavg_c[c,r]     = Vavg_c[c,r-1]
                  Bext550_c[c,r]  = Bext550_c[c,r-1]
                  w550_c[c,r]     = w550_c[c,r-1]
                  g550_c[c,r]     = g550_c[c,r-1]
                  Phs550_c[*,c,r] = Phs550_c[*,c,r-1]

                  Bext_c[*,*,c,r]  = Bext_c[*,*,c,r-1]
                  w_c[*,*,c,r]     = w_c[*,*,c,r-1]
                  g_c[*,*,c,r]     = g_c[*,*,c,r-1]
                  Phs_c[*,*,*,c,r] = Phs_c[*,*,*,c,r-1]
               endif else begin ; This is a new mode radius, do the calculation
                  if lut_MRat[c,r] gt 0 then begin
                     cc = scatoffset + c

;                    If the Mie keyword has been set, we use Mie scattering for
;                    all components, regardless of the driver settings.
                     if keyword_set(mie) then scode = 'mie' $
                     else scode = scatstr.(cc).code

;                    Set up the values of eps and neps, which only exist in the
;                    scatstr structure if tmatrix scattering is to be used.
                     if strlowcase(scode) eq 'tmatrix' then begin
                        epsvals = scatstr.(cc).eps
                        nepsvals = scatstr.(cc).neps
                     endif else begin
                        epsvals = 0
                        nepsvals = 0
                     endelse

;                    Calculate Bext at 550 nm (the reference wavelength) and
;                    Vavg (average volume per particle).
                     Vavg1 = 0.
                     create_bwgp, scatstr.distname[c], lut_Rm[c,r], $
                                  scatstr.S[c], AerM550[c], 0.55, QV, $
                                  Bext1, w1, g1, Phs1, scode=scode, $
                                  tmatrix_path=tmatrix_path, eps=epsvals, $
                                  neps=nepsvals, Vavg=Vavg1
                     Vavg_c[c,r]     = Vavg1
                     Bext550_c[c,r]  = Bext1
                     w550_c[c,r]     = w1
                     g550_c[c,r]     = g1
                     Phs550_c[*,c,r] = Phs1

;                    Calculate Bext, w (single scatter albedo), g (asymmetry
;                    parameter) and Phs (phase function) for each instrument
;                    channel.
                     create_bwgp, scatstr.distname[c], lut_Rm[c,r], $
                                  scatstr.S[c], AerM[*,*,c], ChanSRFWl, $
                                  QV, Bext1, w1, g1, Phs1, scode=scode, $
                                  tmatrix_path=tmatrix_path, eps=epsvals, $
                                  neps=nepsvals
                     Bext_c[*,*,c,r]  = reform(Bext1, NSRF, inststr.NChan)
                     w_c[*,*,c,r]     = reform(w1,    NSRF, inststr.NChan)
                     g_c[*,*,c,r]     = reform(g1,    NSRF, inststr.NChan)
                     Phs_c[*,*,*,c,r] = reform(Phs1,  NMom, NSRF, inststr.NChan)
                  endif
               endelse
            endfor

         endif else if scatstr.comptype[c] eq 'baran' then begin
            init_baran, scatstr.compname[c], baran

            read_baran, baran, 0.55, lutstr.EfR, Bext1, w1, g1, $
                        PTheta * 180. / !pi, Phs1
            Bext550_c[c,*]  = Bext1
            w550_c[c,*]     = w1
            g550_c[c,*]     = g1
            Phs550_c[*,c,*] = Phs1

            ; The Baran optical property data set for thermal channels does not
            ; have g in which case it sets g to zero but g is used below *only*
            ; to distinguish particle layers from layers without particles so we
            ; set g to a nonzero value if it is zero.
            g550_c[where(g_c eq 0)] = -999.

            ; Normalize
            for r=0,lutstr.NEfR-1 do begin
               Phs550_c[*, c, r] /= total(Phs550_c[*, c, r] * weights) / 2.
            endfor

            read_baran, baran, ChanSRFWl, lutstr.EfR, Bext1, w1, g1, $
                        PTheta * 180. / !pi, Phs1
            Bext_c[*,*,c,*]  = reform(Bext1, NSRF, inststr.NChan, lutstr.NEfR)
            w_c[*,*,c,*]     = reform(w1,    NSRF, inststr.NChan, lutstr.NEfR)
            g_c[*,*,c,*]     = reform(g1,    NSRF, inststr.NChan, lutstr.NEfR)
            Phs_c[*,*,*,c,*] = reform(Phs1,  NMom, NSRF, inststr.NChan, lutstr.NEfR)

            g_c[where(g_c eq 0)] = -999.

            ; Normalize
            for r=0,lutstr.NEfR-1 do begin
               for l=0,inststr.NChan-1 do begin
                  for m=0,NSRF-1 do begin
                     Phs_c[*, m, l, c, r] /= total(Phs_c[*, m, l, c, r] * weights) / 2.
                  endfor
               endfor
            endfor

         endif else if scatstr.comptype[c] eq 'baum' then begin
            read_baum_lambda, scatstr.compname[c], 0.55, lutstr.EfR, Bext1, $
                              w1, g1, PTheta * 180. / !pi, Phs1
            Bext550_c[c,*]  = Bext1
            w550_c[c,*]     = w1
            g550_c[c,*]     = g1
            Phs550_c[*,c,*] = Phs1

            ; Normalize
            for r=0,lutstr.NEfR-1 do begin
               Phs550_c[*, c, r] /= total(Phs550_c[*, c, r] * weights) / 2.
            endfor

            ; Choose between spectral or instrument/channel specific properties
            if scatstr.compname2[c] eq '' then begin
               read_baum_lambda, scatstr.compname[c], ChanSRFWl, lutstr.EfR, $
                                 Bext1, w1, g1, PTheta * 180. / !pi, Phs1
            endif else begin
               read_baum_channel, scatstr.compname2[c], fix(channels), lutstr.EfR, $
                                  Bext1, w1, g1, PTheta * 180. / !pi, Phs1
            endelse
            Bext_c[*,*,c,*]  = reform(Bext1, NSRF, inststr.NChan, lutstr.NEfR)
            w_c[*,*,c,*]     = reform(w1,    NSRF, inststr.NChan, lutstr.NEfR)
            g_c[*,*,c,*]     = reform(g1,    NSRF, inststr.NChan, lutstr.NEfR)
            Phs_c[*,*,*,c,*] = reform(Phs1,  NMom, NSRF, inststr.NChan, lutstr.NEfR)

            ; Normalize
            for r=0,lutstr.NEfR-1 do begin
               for l=0,inststr.NChan-1 do begin
                  for m=0,NSRF-1 do begin
                     Phs_c[*, m, l, c, r] /= total(Phs_c[*, m, l, c, r] * weights) / 2.
                  endfor
               endfor
            endfor
         endif
      endfor

      print,''
      print,'All scattering calculations completed for each component'


;     **** Calculate the scattering parameters of the class for each required
;          effective radius

;     Define the arrays to hold the scattering parameters for the class as a
;     whole - these are the variables which are written out into the RT driver
;     files.

      Vavg    = fltarr(lutstr.NEfR)

;     At the reference wavelength
      Bext550 = fltarr(lutstr.NEfR)       ; Extinction coefficient
      w550    = fltarr(lutstr.NEfR)       ; Single scattering albedo
      g550    = fltarr(lutstr.NEfR)       ; Asymmetry parameter
      Phs550  = fltarr(NMom, lutstr.NEfR) ; Phase function
      AMom550 = fltarr(NMom, lutstr.NEfR) ; Legendre moments

;     For each channel
      BextRat = fltarr(NSRF, inststr.NChan, lutstr.NEfR) ; Ratio of Bext with that
;                                                          at the reference wavelength
      Bext    = fltarr(NSRF, inststr.NChan, lutstr.NEfR) ; Extinction coefficient
      w       = fltarr(NSRF, inststr.NChan, lutstr.NEfR) ; Single scattering albedo
      g       = fltarr(NSRF, inststr.NChan, lutstr.NEfR) ; Asymmetry parameter
      Phs     = fltarr(NMom, NSRF, inststr.NChan, lutstr.NEfR) ; Phase function
      AMom    = fltarr(NMom, NSRF, inststr.NChan, lutstr.NEfR) ; Legendre moments

      Vavg[*] = Vavg_c[0,*]

      for l=0,inststr.NChan-1 do begin
         for m=0,NSRF-1 do begin
            for r=0,lutstr.NEfR-1 do begin
;              Calculate 550 nm extinction coefficient
               if l eq 0 then begin
;                 Pre-calculate some factors that are used more than once
                  MratBext   = lut_MRat[*,r] * Bext550_c[*,r]
                  tMratBext  = total(MratBext)
                  MratBextw  = MratBext * w550_c[*,r]
                  tMratBextw = total(MratBextw)

                  Bext550[r] = tMratBext / total(lut_Mrat[*,r])
                  w550[r]    = tMratBextw / tMratBext
                  g550[r]    = total(MratBextw * g550_c[*,r]) / tMratBextw
                  for p=0,NMom-1 do $
                     Phs550[p,r] = total(MratBextw * Phs550_c[p,*,r]) / tMratBextw

;                 Calculate the Legendre moments for the phase function
                  legpexp, NMom, QV, weights, Phs550[*,r], Inlc, alc
                  AMom550[*,r] = alc / (2.0*findgen(NMom)+1.0)
               endif

;              Pre-calculate some factors that are used more than once
               MratBext   = lut_MRat[*,r] * Bext_c[m,l,*,r]
               tMratBext  = total(MratBext)
               MratBextw  = MratBext * w_c[m,l,*,r]
               tMratBextw = total(MratBextw)

               Bext[m,l,r]  = tMratBext / total(lut_Mrat[*,r])
               w[m,l,r]     = tMratBextw / tMratBext
               g[m,l,r]     = total(MratBextw * g_c[m,l,*,r]) / tMratBextw
               for p=0,NMom-1 do $
                  Phs[p,m,l,r] = total(MratBextw * Phs_c[p,m,l,*,r]) / tMratBextw

;              Calculate the Legendre moments for the phase function
               legpexp, NMom, QV, weights, Phs[*,m,l,r], Inlc, alc
               AMom[*,m,l,r] = alc / (2.0*findgen(NMom)+1.0)

;              Calculate the ratio of the extinction coefficient at the current
;              channel and 0.55 microns, allowing the spectral optical depth to be
;              calculated from the reference 0.55 micron value.
               BextRat[m,l,r] = Bext[m,l,r]/Bext550[r]
            endfor
         endfor
      endfor

;  **** Write the scattering parameters for the class as a whole for reuse
      save, NMom, Bext550, w550, g550, Phs550, AMom550, BextRat, Bext, w, g, $
            Phs, AMom, filename=out_path+'/scatfile.sav'

;  **** Create scattering properties output that can be feed in the RAL LUT
;       generation code
;     x = create_struct('Class', scatstr.lutname)
;     tmp1 = fltarr(inststr.NChan+1,lutstr.NEfR)
;     tmp1[0,*] = transpose(Bext550)
;     tmp1[1:*,*] = Bext
;     x = create_struct(x, 'KE', tmp1)
;     tmp1[0,*] = transpose(w550)*transpose(Bext550)
;     tmp1[1:*,*] = w*Bext
;     x = create_struct(x, 'KS', tmp1)
;     tmp1 = fltarr(inststr.NChan+1,lutstr.NEfR,NMom)
;     tmp1[0,*,*] = reform(transpose(AMom550),1,lutstr.NEfR,NMom)
;     tmp1[1:*,*,*] = transpose(AMom,[1,2,0])
;     x = create_struct(x, 'PLM', tmp1)
;     tmp1 = fltarr(inststr.NChan+1,lutstr.NEfR,NMom)
;     tmp1[0,*,*] = reform(transpose(Phs550),1,lutstr.NEfR,NMom)
;     tmp1[1:*,*,*] = transpose(Phs,[1,2,0])
;     x = create_struct(x, 'P', tmp1)
;     x = create_struct(x, 'WLS', [0.55, inststr.ChanWl], 'Re', lutstr.EfR)
;     x = create_struct(x, 'Theta', Ptheta*!radeg)

;     Output this structure using the save procedure
      save, x, filename=out_path+'/'+strupcase(inststr.name)+'_'+ $
            scatstr.outname+verstrng+'_scattering.str', /compress
   endelse


;  -----------------------------------------------------------------------------
;  Output scattering data into ORAC LUTs
;  -----------------------------------------------------------------------------

;  Write the Vavg LUT
   lutname = lutbase+verstrng+'_Vavg.sad'
   write_orac_lut_1d, lutname, lutstr.EfR, Vavg

;  Write the Bext LUT
   lutname = lutbase+verstrng+'_Bext.sad'
   write_orac_lut_2d, lutname, lutstr.AOD, lutstr.EfR, Bext

;  Write the Bext ratio LUT
   lutname = lutbase+verstrng+'_BextRat.sad'
   write_orac_lut_2d, lutname, lutstr.AOD, lutstr.EfR, BextRat

;  The rest of the LUTs are written out per channel

;  We want the center point of the SRF integration which, since the number of
;  points is odd, will be at the center of the channel's spectral interval.
   m = NSRF / 2

   for l=0,inststr.NChan-1 do begin
      ChStrng = 'Ch'+string(inststr.ChanNum[l], format=Chfmt)

;     Bext
      lutname = lutbase+'_Bext_'+ChStrng+verstrng+'.sad'
      BextLUT = fltarr(lutstr.NAOD, lutstr.NEfR)
      for a=0,lutstr.NAOD-1 do BextLUT[a,*] = Bext[m,l,*]
      write_orac_lut_2d, lutname, lutstr.AOD, lutstr.EfR, BextLUT, $
                         Wl=inststr.ChanWl[l], logAOD=lutstr.LAOD, $
                         logEfR=lutstr.LEfR

;     Bext ratio
      lutname = lutbase+'_BextRat_'+ChStrng+verstrng+'.sad'
      BextLUT = fltarr(lutstr.NAOD, lutstr.NEfR)
      for a=0,lutstr.NAOD-1 do BextLUT[a,*] = BextRat[m,l,*]
      write_orac_lut_2d, lutname, lutstr.AOD, lutstr.EfR, BextLUT, $
                         logAOD=lutstr.LAOD, logEfR=lutstr.LEfR

;     w
      lutname = lutbase+'_w_'+ChStrng+verstrng+'.sad'
      wLUT = fltarr(lutstr.NAOD, lutstr.NEfR)
      for a=0,lutstr.NAOD-1 do wLUT[a,*] = w[m,l,*]
      write_orac_lut_2d, lutname, lutstr.AOD, lutstr.EfR, wLUT, $
                         Wl=inststr.ChanWl[l], logAOD=lutstr.LAOD, $
                         logEfR=lutstr.LEfR

;     g
      lutname = lutbase+'_g_'+ChStrng+verstrng+'.sad'
      gLUT = fltarr(lutstr.NAOD, lutstr.NEfR)
      for a=0,lutstr.NAOD-1 do gLUT[a,*] = g[m,l,*]
      write_orac_lut_2d, lutname, lutstr.AOD, lutstr.EfR, gLUT, $
                         Wl=inststr.ChanWl[l], logAOD=lutstr.LAOD, $
                         logEfR=lutstr.LEfR
   endfor

;  **** Done with calculating scattering parameters, return if that is all that
;       is wanted

   print,''
   print,'Scattering parameters calculated for class '+lutbase+verstrng

   if keyword_set(scat_only) then begin
      return,0
   endif


;  -----------------------------------------------------------------------------
;  Run DISORT
;  -----------------------------------------------------------------------------

;  **** Setup the variables needed for the DISORT calls ****
   setup_disort, 60, NLayers, lutstr.NSat, lutstr.NAzi, NMom

;  **** Define the LUT table output variables themselves
   RFD  = fltarr(lutstr.NAOD,              lutstr.NEfR, inststr.NChan)
   TFD  = fltarr(lutstr.NAOD,              lutstr.NEfR, inststr.NChan)
   RD   = fltarr(lutstr.NAOD, lutstr.NSat, lutstr.NEfR, inststr.NChan)
   TD   = fltarr(lutstr.NAOD, lutstr.NSat, lutstr.NEfR, inststr.NChan)
   TB   = fltarr(lutstr.NAOD, lutstr.NSol, lutstr.NEfR, inststr.NChan)
   RFBD = fltarr(lutstr.NAOD, lutstr.NSol, lutstr.NEfR, inststr.NChan)
   TFBD = fltarr(lutstr.NAOD, lutstr.NSol, lutstr.NEfR, inststr.NChan)
   RBD  = fltarr(lutstr.NAOD, lutstr.NSat, lutstr.NSol, $
                 lutstr.NAzi, lutstr.NEfR, inststr.NChan)
   TBD  = fltarr(lutstr.NAOD, lutstr.NSat, lutstr.NSol, $
                 lutstr.NAzi, lutstr.NEfR, inststr.NChan)
   Em   = fltarr(lutstr.NAOD, lutstr.NSat, lutstr.NEfR, inststr.NChan)

   UUd  = fltarr(lutstr.NAOD, lutstr.NSat*2, lutstr.NSol, lutstr.NAzi, $
                 lutstr.NEfR, inststr.NChan)

;  **** Loop through the channels (and solar zenith angles) and run DISORT for
;       the beam and diffuse cases. Also produce the emissivity for the channels
;       that need it.

;  Define the Rayleigh scattering optical depth in each channel
   if keyword_set(no_rayleigh) then begin
      ColumnTauRay = fltarr(inststr.NChan)
      ColumnTauRay = replicate(1.e-6, inststr.NChan)
   endif else begin
      ColumnTauRay = (presstr.p[presstr.Nlevels-1] / 1013.0) / $
                     (117.03*inststr.ChanWl^4 - 1.316*inststr.ChanWl^2)
   endelse

   store_tau = fltarr(lutstr.NEfR, lutstr.NAOD, inststr.NChan)

   for l=0,inststr.NChan-1 do begin

      for m=0,NSRF-1 do begin
         print,'Running DISORT for channel '+strtrim(inststr.ChanNum[l],2)+ $
               ' (',strtrim(inststr.ChanWl[l],2),'um)'

;        Do we have a Gas optical depth profile for the current channel?
         if n_elements(gasdat) gt 0 then begin
            GasIndx = where(tag_names(gasstr) eq $
                            'C'+string(inststr.Channum[l],format='(I02)'))
            GasIndx = GasIndx[0]
         endif else GasIndx = -1
;        If we don't have gas OPD for this channel,then use zeros (i.e. no gas
;        absorption). Remember that the gas OPD is defined on the levels between
;        each atmospheric layer, not on the layers themselves.
         if GasIndx lt 0 then begin
            print, 'No gas optical depth profile. Assuming Rayleigh scattering only.'
            GasLvl = replicate(0.0,presstr.NLevels)
         endif else GasLvl = gasstr.(GasIndx).OPD

;        Calculate the cumulative Rayleigh optical depth at each level
         RayLvl = ColumnTauRay[l]*exp(-0.1188*presstr.H - 0.00116*presstr.H^2)

;        The optical depths of each layer from gas absorption and Rayleigh
;        scattering are defined as the difference between the optical depth at
;        the adjacent levels.
         TauGas = GasLvl[lindgen(NLayers)+1] - GasLvl[lindgen(NLayers)]
         TauRay = RayLvl[lindgen(NLayers)+1] - RayLvl[lindgen(NLayers)]

         for a=0,lutstr.NAOD-1 do begin

            for r=0,lutstr.NEfR-1 do begin
;              The optical depth from Aerosol is the desired total AODs for the
;              ORAC LUT * the relative AOD at each layer for this class * the
;              scaling factor relating AOD at this wavelength back to 550 nm.
               TauAer = lutstr.AOD[a] * AerRelTau * BextRat[m,l,r]
;              The optical depths are additive
               DTau = TauGas + TauRay + TauAer
               TotalTau = total(DTau)

               store_tau[r,a,l] = TotalTau

;              The single scattering albedo is weighted by optical depth.
;              NB. SSA for Rayleigh scattering = 1, and is effectively 0 for
;              gas absorption.
               SSALB = (TauRay + w[m,l,r]*TauAer) / DTau
;              Now check that we have no SSALB values over 1.0 (this can happen in
;              layers with no absorption due to rounding). DISORT has an internal
;              check for this and will exit with an error code if it fails.
               bd = where(SSALB gt 1.0)
               if bd[0] ge 0 then SSALB[bd] = 1.0;0.999999

;              Asymmetry parameter is only non-zero where we actually have aerosol
               ASYM = fltarr(NLayers)
               nonzero = where(TauAer gt 0.0)
               if nonzero[0] ge 0 then ASYM[nonzero] = g[m,l,r]

;              Now, use the GETMOM procedure (part of DISORT) to generate phase
;              function moments for the molecular scattering and then combine with
;              the aerosol moments generated earlier.
               PMom = fltarr(NMom, NLayers)
               for h=0,NLayers-1 do begin
                  if ASYM[h] eq 0.0 then GETMOM, 2, 0.0, NMom-1, PM $
                  else begin
                     GETMOM, 2, 0.0, NMom-1, mPM
                     PM = (mPM*TauRay[h] + AMom[*,m,l,r]*w[m,l,r]*TauAer[h]) / $
                          (TauRay[h] + w[m,l,r]*TauAer[h])
                  endelse
                  bd = where(PM gt 1.0)
                  if bd[0] ge 0 then PM[bd] = 1.0
                  PMom[*,h] = PM
               endfor

;              We are now ready to call DISORT. Call the fast diffuse calculation
;              first (errors and problems are more likely to turn up quickly that
;              way).

               print, 'Doing RT calculation for' + $
                      ' Channel: ' + string(inststr.ChanNum[l], $
                      inststr.ChanWl[l],format='(i3," (",f7.4,"um)")') + $
                      ', Point: ' + string(m,format='(i4)') + $
                      ', AOD: ' + string(lutstr.AOD[a],format='(e14.6)') + $
                      ', EfR: ' + string(lutstr.EfR[r],format='(f8.4)')
               print, ''

               print, '---------- DIFFUSE -----------'
               FBeam =   0.0           ; Direct beam intensity
               FIsot = 100.0           ; Isotropic illumination intensity
               UMu0  = cos(50.0*!dtor) ; A nominal value for beam zenith
               UTau  = [0.0, TotalTau] ; Define output layers (in terms of optical
                                       ; depth)
;              UMu is calculated to described downwelling as well as upwelling
;              radiance, as we also need the downwelling values
               UMu = fltarr(2*lutstr.NSat)
;              If 90 degrees is included in the list of satellite zeniths, take a
;              small angle off of it for the DISORT calculations in order to avoid
;              numerical problems.
               tmpSat = lutstr.Sat
               bd = where(tmpSat eq 90.0)
               if bd[0] ge 0 then tmpSat[bd] = 89.99
               UMu[0:lutstr.NSat-1] = -1.0* cos(tmpSat*!dtor)
               UMu[lutstr.NSat:2*lutstr.NSat-1] = $
                  -1.0*UMu[lutstr.NSat-lindgen(lutstr.NSat)-1]
;              DISORT will spit the dummy if there if abs(UMu)=1
               bd = where(abs(UMu) eq 1.0)
               if bd[0] ge 0 then UMu[bd] = 0.99999 * UMu[bd]/abs(UMu[bd])

               call_disort, DTau, SSAlb, PMom, UTau, UMu, lutstr.Azi, $
                            FBeam, UMu0, FISot, RFlDir, RFlDn, FlUp, $
                            dFdT, UAvg, UU, AlbMed, TrnMed

;              Generate the diffuse LUT variables
               RFD[a,r,l] += (100. * FlUp[0]  / (FIsot*!pi)) * ChanSRF[m,l]
               TFD[a,r,l] += (100. * RFlDn[1] / (FIsot*!pi)) * ChanSRF[m,l]
;              RD contains the Upwelling intensity
               RD[a,*,r,l] += $
                 (UU[2*lutstr.NSat-lindgen(lutstr.NSat)-1,0,0]) * ChanSRF[m,l]
;              TD contains the Downwelling intenisity
               TD[a,*,r,l] += $
                 (UU[              lindgen(lutstr.NSat),  1,0]) * ChanSRF[m,l]

               print, ''

;              If the channel has the emission flag set, calculate the
;              emissivity.
               if inststr.ChanEm[l] then begin
                  print, '---------- EMISSION ----------'
;                 Elisa's expression for emissivity....
;                 Em[a,*,r,l]  = 100.0*(1.0 - w(l,r)) * $
;                                (1.0 - exp(TotalTau*(-1.0/cos(lutstr.Sat*!dtor))))

;                 Use DISORT - I can't seem to get this to work correctly. For
;                 this we set both beam and diffuse input irradiances to zero
;                 and calculate emission across the window 1% of the nominal
;                 wavenumber (cm-1). Everything else is the same as diffuse
;                 calculations.
                  FBeam   = 0.0 ; Direct beam intensity
                  FIsot   = 0.0 ; Isotropic illumination intensity
                  wn      = 1e4 / inststr.ChanWl[l]
                  wnlo    = 0.995*wn
                  wnhi    = 1.005*wn
                  temp    = 250.0
                  incloud = where(TauAer gt 0.0,emnly)
                  emTau   = DTau[incloud]
                  emSSA   = SSAlb[incloud]
                  emPMo   = PMom[*,incloud]
                  emUTau  = [0.0, total(emTau)]

                  call_disort, emTau, emSSA, emPMo, emUTau, UMu, lutstr.Azi, $
                               FBeam, UMu0, FISot, RFlDir, RFlDn, FlUp, dFdT, $
                               UAvg, UU, AlbMed, TrnMed, /plank, wnlo=wnlo, $
                               wnhi=wnhi, temp=temp, nlayer=emnly

;                 Now we calculate the Plank emission across the wavelength
;                 interval.
                  BBE = PLKAVG(wnlo, wnhi, temp)

;                 Finally, combine to produce the emissivity
                  Em[a,*,r,l] += (100.0 * $
                     UU[2*lutstr.NSat-lindgen(lutstr.NSat)-1,0,0]/BBE) * ChanSRF[m,l]

                  print, ''
               endif

;              Now we loop over the solar zenith angles and do the direct beam
;              calculations. Note that this only needs to be done for channels with
;              a solar component to their signal.
               if inststr.ChanSol[l] then begin
                  print, '---------- DIRECT ----------'
                  for s=0,lutstr.NSol-1 do begin
;                    print, 'SZA: ', lutstr.Sol[s]

                     FBeam = 100.0 ; Direct beam intensity
                     FIsot =   0.0 ; Isotropic illumination intensity

;                    If a solar zenith angle of 90 degrees has been requested,
;                    alter it to something slightly smaller to prevent DISORT from
;                    crashing.
                     if lutstr.Sol[s] eq 90.0 then tmpSol = 89.99 $
                     else tmpSol = lutstr.Sol[s]
                     UMu0  = cos(tmpSol*!dtor)

                     call_disort, DTau, SSAlb, PMom, UTau, UMu, lutstr.Azi, $
                                  FBeam, UMu0, FISot, RFlDir, RFlDn, FlUp, $
                                  dFdT, UAvg, UU, AlbMed, TrnMed

;                    Generate the direct beam LUT variables
                     TB[a,s,r,l]   += (100. * RFlDir[1] / RFlDir[0]) * ChanSRF[m,l]
                     RFBD[a,s,r,l] += (100. * FlUp[0]   / RFlDir[0]) * ChanSRF[m,l]
                     TFBD[a,s,r,l] += (100. * RFlDn[1]  / RFlDir[0]) * ChanSRF[m,l]
                     for p=0,lutstr.NAzi-1 do begin
;                       Bug investigation.... is relative azimuth backwards?
;                       p2 = lutstr.NAzi - p - 1
                        p2 = p
;                       As with the diffuse case, RBD contains the upwelling
;                       intensity, while TBD contains the downwelling.
                        RBD[a,*,s,p2,r,l] += $
                           (UU[2*lutstr.NSat-lindgen(lutstr.NSat)-1,0,p] * !pi) * $
                           ChanSRF[m,l]
                        TBD[a,*,s,p2,r,l] += $
                           (UU[              lindgen(lutstr.NSat),  1,p] * !pi) * $
                           ChanSRF[m,l]
                        UUd[a,*,s,p2,r,l] = UU[*,0,p]
                     endfor
                  endfor

                  print, ''
               endif
            endfor ; End of EfR loop
         endfor ; End of AOD loop
         print,''
      endfor ; End of SRF loop
   endfor ; End of channel loop

   ; Normalize the RT operators wrt to the SRF.  Note 'sum' will be unity in
   ; monochomatic mode (when no SRFs were provided).
   for l=0,inststr.NChan-1 do begin
      sum = total(ChanSRF[*,l])

      RFD [*,*,l]       /= sum
      TFD [*,*,l]       /= sum
      RD  [*,*,*,l]     /= sum
      TD  [*,*,*,l]     /= sum
      TB  [*,*,*,l]     /= sum
      RFBD[*,*,*,l]     /= sum
      TFBD[*,*,*,l]     /= sum
      RBD [*,*,*,*,*,l] /= sum
      TBD [*,*,*,*,*,l] /= sum
      Em  [*,*,*,l]     /= sum
   endfor


;  save,/variables,filename='debug.sav',/compress


;  -----------------------------------------------------------------------------
;  Output reflectance, transmission and emission data into ORAC LUTs
;  -----------------------------------------------------------------------------

;  The rest of the LUTs are written out per channel
   for l=0,inststr.NChan-1 do begin
      ChStrng = 'Ch'+string(inststr.ChanNum[l], format=Chfmt)

;     The diffuse transmissions
      lutname = lutbase+'_TD_'+ChStrng+verstrng+'.sad'
      write_orac_lut_3d, lutname, lutstr.AOD, lutstr.Sat, lutstr.EfR, $
                         TD[*,*,*,l], Wl=inststr.ChanWl[l], $
                         Data2=TfD[*,*,l], logAOD=lutstr.LAOD, $
                         logEfR=lutstr.LEfR

;     The diffuse reflectances
      lutname = lutbase+'_RD_'+ChStrng+verstrng+'.sad'
      write_orac_lut_3d, lutname, lutstr.AOD, lutstr.Sat, lutstr.EfR, $
                         RD[*,*,*,l], Wl=inststr.ChanWl[l], $
                         Data2=RfD[*,*,l], logAOD=lutstr.LAOD, $
                         logEfR=lutstr.LEfR

;     Now check for a solar component to the channel, and output direct beam
;     transmissions and reflectances if needed.
      if inststr.ChanSol[l] then begin
;        The directional reflectance
         lutname = lutbase+'_RBD_'+ChStrng+verstrng+'.sad'
         write_orac_lut_5d, lutname, lutstr.AOD, lutstr.Sat, lutstr.Sol, $
                            lutstr.Azi, lutstr.EfR, RBD[*,*,*,*,*,l], $
                            Wl=inststr.ChanWl[l], Data2=RFBD[*,*,*,l], $
                            logAOD=lutstr.LAOD, logEfR=lutstr.LEfR

;        The bi-directional transmission and the diffuse-only transmission of
;        the beam
         lutname = lutbase+'_TBD_'+ChStrng+verstrng+'.sad'
         write_orac_lut_5d, lutname, lutstr.AOD, lutstr.Sat, lutstr.Sol, $
                            lutstr.Azi, lutstr.EfR, TBD[*,*,*,*,*,l], $
                            Wl=inststr.ChanWl[l], Data2=TFBD[*,*,*,l], $
                            logAOD=lutstr.LAOD, logEfR=lutstr.LEfR

;        The direct transmission of the direct beam
         lutname = lutbase+'_TB_'+ChStrng+verstrng+'.sad'
         write_orac_lut_3d, lutname, lutstr.AOD, lutstr.Sol, lutstr.EfR, $
                            TB[*,*,*,l], Wl=inststr.ChanWl[l], $
                            logAOD=lutstr.LAOD, logEfR=lutstr.LEfR
      endif

;     Finally, for those channels which require it, the Emission
      if inststr.ChanEm[l] then begin
         lutname = lutbase+'_EM_'+ChStrng+verstrng+'.sad'
         write_orac_lut_3d, lutname, lutstr.AOD, lutstr.Sat, $
                            lutstr.EfR, EM[*,*,*,l], Wl=inststr.ChanWl[l], $
                            logAOD=lutstr.LAOD, logEfR=lutstr.LEfR
      endif
   endfor


;  -----------------------------------------------------------------------------
;  Output termination timestamp.
;  -----------------------------------------------------------------------------

   openw, lun, lutbase + '_timestamp.txt', /get_lun
   printf, lun, strmid(timestamp(/utc), 0, 19) + 'Z'
   free_lun, lun


   print,'ORAC LUT generation completed for class '+lutbase+verstrng

   return,0
end
