;+
; function create_orac_aerosol_lut
;
; A routine FOR generating ORAC look-up tables (LUTs) from aerosol/cloud optical
; properties. The code does the following steps:
; * Reads in a series of driver text files
; * Performs the necessary scattering calculations to produce scattering
;   properties FOR the given aerosol/cloud class
; * Calls DISORT to perform the radiative transfer calculations
; * Writes the LUT "static application data" (SAD) files needed by ORAC
;
; INPUT ARGUMENTS:
; in_path     (string) The path in which the code can find the input file directories
; instfile     (string) The name of the instrument data driver file to use
; mmfile      (string) The name of the aerosol/cloud optical/microphysical
;                      properties  file to use
; lutfile      (string) The name of the LUT dimension and vertex driver file to
;                      use
; out_path    (string) The base bath containing the directory to save the generated LUTs
; atmospheres (string) Value '1' - '6' specifing the atmosphere to use
;
; INPUT KEYWORDS:
; channelID=strarr     Provide a list of channelIDs  to produce LUTs FOR.
;                      These must be a subset of the channelIDs listed in
;                      the instfile file. (By default the code will produce LUTs
;                      FOR all the channels defined in instfile).
; /no_rayleigh         Do not include Rayleigh scattering.
; /no_screen           By default the code uses string(13b) to prevent new-
;                      lines at the end of print-to-screen commands (FOR
;                      aesthetic reasons): these make a real mess IF output is
;                      piped to a file. Setting this keyword suppresses these
;                      control characters
; srf_quad=integer Set the quadrature for the integration across spectral response function.
;                      0 srf resolution
;                      1 single wavelength
;                      2 reduced srf resolution so that integral is accurate to ~ 0.1 %
; /reuse_scat          Reuse the scattering computations from a previous run
;                      stored in the file scatfile.sav created with the IDL
;                      save procedure. This file is saved whenever /reuse_scat
;                      is not set. And may only be used by a run with the same
;                      microphysical configuration and set of channels.
; /scat_only           Only compute the scattering properties skipping over
;                      the radiative transfer calculations. In this case the
;                      file scatfile.sav and the Bext, BextRat, w and g LUTS
;                      will be output WHEREas the reflection, transmission and
;                      emission LUTS will not.
; tmatrix_path=string  Path to the Dubovik T-Matrix LUT base directory.
;                      Required only when T-Matrix calculations are to to be
;                      made.
; version=string       Set a version string to place in the output LUT file
;                      names (defaults to nothing).
;
; RETURN VALUE:
; At the moment, IF the code finishes successfully, the value "0" will always be
; returned. Error/warning codes are a future option.
;
; OUTPUT ARGUMENTS:
; None
;
; HISTORY:
; 05/02/13, G Thomas: Original version.
; 08/02/13, G Thomas: Added channels keyword. Added code to make use of the
;                     log flag FOR AOD and EfR in lutfile files.
; 19/02/13, G Thomas: Debugging completed. Added check on solar component FOR
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
; 26/06/13, G Thomas: Altered the wavelength threshold FOR the force_k keyword
;                     to 2.0 microns, FOR producing LUTs FOR the SMASH project.
; XX/XX/15, G McGarragh: Add support FOR a gamma distribution as a single mode.
;    Useful FOR liquid water cloud.
; XX/XX/15, G McGarragh: Add support FOR a single mode of ice crystals according
;    Braun Baum et al.
; XX/XX/15, G McGarragh: Add support FOR a single mode of ice crystals according
;    Anthony Baran et al.
; XX/XX/15, G McGarragh: Change Bext output to be that actual Bext and not the
;    ratio with that of the reference wavelength.  Put the actual ratio output
;    into the BextRat LUT.  Output single scattering albedo into the 'w' LUT,
;    asymmetry parameter into the 'g' LUT and the average volume per particle
;    into the 'Vavg' LUT.
; XX/XX/15, G McGarragh: Add the /scat_only keyword to skip the time consuming
;    RT calculations IF the scattering output is all that is desired.
; XX/XX/15, G McGarragh: Add the /reuse_scat keyword to reuse the scattering
;    calculations from a previous run.  See the documentation FOR details.
; XX/XX/15, G McGarragh: The force_k keyword is now a vector WHERE the first
;    element is FOR the reference wavelength and the other elements are FOR the
;    channels.  Also added a corresponding force_n keyword.
; 08/06/16, G McGarragh: Add the /no_rayleigh keyword FOR the option to not
;    include Rayleigh scattering.  Useful FOR building LUTs FOR multilayer
;    algorithms.
; 31/08/16, G McGarragh: Add the tmatrix_path keyword as the path to the Dubovik
;    T-Matrix LUTs was fixed before.
; 12/10/16, G McGarragh: The phase functions interpolated from the Baran and
;    Baum ice crystal scattering properties must be normalized.
; 20/10/16, G McGarragh: Write the svn version and termination timestamp to
;    separate files.
; 21/10/16, G McGarragh: Copy the driver files to the output directory using the
;    LUT output base name.
; 14/04/18, G McGarragh: Add support FOR integration over channel spectral
;    response functions (SRFs).
; 21/04/18, G McGarragh: Change revision output from svn to git.
; 20/03/20, G Thomas: Added optical properties LUTs functionality.
; 20/07/20, RGG: Create V2 cdhf LUTS, 
; V2 Notes:  
; 1) *** Implicitly assumes there are less than 256 channels (need to add warning).
; 2) Creates LUT directory that mirrors driver file name.
; 26/08/23, RGG: Seperated aerosol & cloud LUT creation 
; 07/05/26  RGG: Removed keywords that were passed but never implemented 

; Begin the main LUT generation function
  function create_orac_aerosol_lut,in_path,     $
                           instfile,     $
                           mmfile,      $
                           lutfile,      $
                           out_path,    $
                           atmospheres, $                   
                           channelID     = channelID,     $
                           gas           = gas,           $
                           no_rayleigh   = no_rayleigh,   $
                           no_screen     = no_screen,     $
                           srf_quad      = srf_quad,      $
                           reuse_scat    = reuse_scat,    $
                           scat_only     = scat_only,     $
                           tmatrix_path  = tmatrix_path,  $
                           version       = version,       $
                           driver        = driver

;  -----------------------------------------------------------------------------
;  Test input and output files and directories
;  -----------------------------------------------------------------------------
   print, 'Reading input files...'

;  Check FOR the existence of the input file path and the required files
   ok = file_test(in_path, /directory, /read, /write)
   IF not ok then message, 'in_path not found, or is read/write protected'

;  Instrument definition file
   instdirfile = in_path+'/inst/'+instfile
   ok = file_test(instdirfile, /read)
   IF not ok then message, 'instdirfile not readable: ' + instdirfile

;  Microphysical definition file
   mmdirfile = in_path+'/microphysics/'+mmfile
   ok = file_test(mmdirfile, /read)
   IF not ok then message, 'mmdirfile not readable: ' + mmdirfile

;  LUT grid definition file
   lutdirfile = in_path +'/lut/'+lutfile
   ok = file_test(lutdirfile, /read)
   IF not ok then message, 'lutdirfile not readable: ' + lutdirfile

;  Pressure profile definition file
   Case Atmospheres of
     '0': begin
          atmfile = 'midsatm.dat'
          atmosphere_model = 'midlatitude summer'
          end
     '1': begin
          atmfile = 'tro.atm'
          atmosphere_model = 'tropical summer'
          end
     '2': begin
          atmfile = 'mls.atm'
          atmosphere_model = 'midlatitude summer'
          end
     '3': begin
          atmfile = 'mlw.atm'
          atmosphere_model = 'midlatitude winter'
          end
     '4': begin
          atmfile = 'sas.atm'
          atmosphere_model = 'subarctic  summer'
          end
     '5': begin
          atmfile = 'saw.atm'
          atmosphere_model = 'subarctic summer'
          end
     '6': begin
          atmfile = 'std.atm'
          atmosphere_model = 'US standard'
          end
   End
   atmdirfile = in_path+'/atm/'+atmfile
   ok = file_test(atmdirfile, /read)
   IF NOT ok then message, 'atmospheric file not readable: ' + atmdirfile

;  Finally, check that the output directory exists
   ok = file_test(out_path, /directory, /read, /write)
   IF NOT ok then print, 'out_path: '+out_path+ ' creating ..'
   FILE_MKDIR, out_path
   
;  -----------------------------------------------------------------------------
;  Convert Keywords to flags 
;  -----------------------------------------------------------------------------
   If KEYWORD_SET(no_Rayleigh) then Rayleigh_Flag = 0 else  Rayleigh_Flag = 1
   If KEYWORD_SET(gas) then gas_Flag = 1 else  Gas_Flag = 0

;  -----------------------------------------------------------------------------
;  Read in files
;  -----------------------------------------------------------------------------


;  **** Read the Instrument parameters file
   load_inststr, instdirfile, inststr, RequestedChannelID = ChannelID
   
;  **** Read the LUT parameters file
   load_lutstr, lutdirfile, inststr.max_sat_zenith, lutstr,/include_pressure
   
;  **** Read the spectral response functions for the relevant channels and generate integral quantities that depend upon the srf
   solar_spectrum_filename = in_path+'/sun/Gueymard2018.sssi'
   IF KEYWORD_SET(srf_quad) THEN QM = srf_quad else QM = 1
   load_srfstrarr,  inststr, solar_spectrum_filename, srfstrarr, QM, nwvl_max
   
;  **** Read the scattering parameters file
   load_mmdat, mmdirfile, mmstr

;  **** Read the atmospheric profile file
   load_atmstr, atmdirfile, Atmospheres, atmstr

;  **** If flagged Read the GAS OPD files
   IF Gas_Flag THEN begin
     gasdir = in_path+'/gas/'
     load_gasstr, atmospheres, inststr.platform, inststr.instrument, inststr.channelid,  gasdir, gasstr
   
;    Check that the Gas OPD values are on the same grid as the pressure profile
     FOR i=0,N_ELEMENTS(gasstr)-1 DO $
       IF ~ARRAY_EQUAL(atmstr.height, gasstr[i].height) THEN $
         MESSAGE, 'The pressure and gas OPD profiles must be on the same altitude grid'
   EndIf
   
;  -----------------------------------------------------------------------------
;  Report setup.
;  -----------------------------------------------------------------------------
   If (Inststr.View Gt 0) Then $
     print,'Dual View LUT calculation for substance ',mmstr.substance,' for the ', inststr.instrument,' instrument' $
   Else $
     print,'Single View LUT calculation for substance ',mmstr.substance,' for the ', inststr.instrument,' instrument'   
   print,'Output placed in: '+out_path
   print,atmosphere_model + ' (code = ',Atmospheres, ')'
   If (Gas_Flag) Then print,'Including gas absorption' else print,'No gas absorption'
   If (Rayleigh_Flag) Then print,'Including Rayleigh scattering' else print,'No Rayleigh scattering'

   print,'LUT Dimensions are:'
   print,'           Components ', strtrim(mmstr.NComp,2)
   print,'  Instrument channels ', strtrim(inststr.number_of_nadir_channels,2)
   print,'        Optical depth ', strtrim(lutstr.opd_n,2)
   print,'     Effective radius ', strtrim(lutstr.efr_n,2)
   print,'         Solar zenith ', strtrim(lutstr.soz_n,2)
   print,'    Instrument zenith ', strtrim(lutstr.saz_n,2)
   print,'     Relative azimuth ', strtrim(lutstr.raa_n,2)
   print,'            Pressures ', strtrim(lutstr.prs_n,2)
   print,'SRF quadrature method ', strtrim(QM,2)
   If (QM Eq 2) Then FOR i = 0,inststr.number_of_nadir_channels-1 DO print,'Channel ', i, ' Number of quadrature points: ',srfstrarr[i].nwvl
   

;  -----------------------------------------------------------------------------
; Create the output filename. ********** this has to be improved but will do for now
;  -----------------------------------------------------------------------------
  If (QM Eq 1) then MonoOrBand = 'm' else MonoOrBand ='b'
; substance (generic particle type)
;  liquid-water
;  water-ice
;  aerosol
;  volcanic-ash
; atmospheric model code 00 - 99
;  1X  scattering withing atmosphere including Rayleigh scattering and gasous absorption (generally used for aerosol)
;      X denotes the atmospheric gaseous model
;  00  scattering layer that excludes Rayleigh scattering (generally used for bottom level of multilevel cloud)
;  01  scattering layer that includes Rayleigh scattering for entire atmosphere (generally used for cloud)
  If KEYWORD_SET(gas) Then $
    Atmospheric_Model_Code = '1'+Atmospheres $
  else $
    If KEYWORD_SET(no_rayleigh) then $
      Atmospheric_Model_Code = '00' $
    else $
      Atmospheric_Model_Code = '01'

;  -----------------------------------------------------------------------------
;  Some miscellaneous setup.
;  -----------------------------------------------------------------------------

; Particle_model_code 3 digit string set in microphysical model definition file 
  Versions = '00'
  IF N_ELEMENTS(version) gt 0 then Versions    = string(Version,Format='(I2.2)')
  V2_LUT_Filename = out_path+'/'+ strlowcase(inststr.Platform)+'_'+strlowcase(inststr.Instrument)+'_'+MonoOrBand+'_'+strlowcase(mmstr.Substance)+'_a'+Atmospheric_Model_Code+'_p'+strlowcase(mmstr.shortname)+'_v'+versions+'.nc'
 
;  Check the size of the channel string needed...
   Chfmt = '(i0)'
;  IF max(inststr.ChannelID) ge 100 then Chfmt = '(i03)' $
;  ELSE Chfmt = '(i02)'

;  -----------------------------------------------------------------------------
;  Output git revision and copy the driver files to the output directory using
;  the LUT output base name.
;  -----------------------------------------------------------------------------
;   spawn, '; git --git-dir=' + file_dirname((routine_info('create_orac_lut', $     ; removed on switch to bash
;          /function, /source)).path) + '/.git rev-parse HEAD > ' + out_path + $    ; removed on switch to bash
;          '/git_revision.txt'                                                      ; removed on switch to bash
  
   FILE_COPY, driver,   out_path + '/', /OVERWRITE

;  -----------------------------------------------------------------------------
;  Interpolate the aerosol profile layers onto the atmos. pressure and gas OPD
;  layers, AND THE AEROSOL REFRACTIVE INDEX ONTO THE CHANNEL WAVELENGTHS.
;  -----------------------------------------------------------------------------

;  Firstly, we have to define the height of the layers, which lie between each
;  pressure level...
   nlayers = atmstr.nlevels -1
   hlayers = (atmstr.height[0:nlayers-1] + atmstr.height[1:nlayers]) / 2
   scatreltau = INTERPOL(mmstr.rext, mmstr.height, hlayers)
   scatreltau = scatreltau/total(scatreltau)

;  -----------------------------------------------------------------------------
;  Generate scattering properties either through calling a scattering code for
;  Mie or t-matrix or loading 'baum' or 'buran' or ice crystal properties, or by
;  reloading properties saved from a previous run.
;  -----------------------------------------------------------------------------

  IF ~KEYWORD_SET(reuse_scat) THEN BEGIN
    generate_scattering_properties,srfstrarr, scatoffset, nwvl_max, inststr, mmstr, lutstr, nmom, bext550, w550, g550, phs550, amom550, bextrat, bext, w, g, vavg, phs, amom,tmatrix_path=tmatrix_path,no_screen=no_screen
;   **** write the scattering parameters FOR the class as a whole FOR reuse
    SAVE, FILENAME = out_path + '/scatfile.sav', nmom, bext550, w550, g550, phs550, amom550, bextrat, bext, w, g, vavg, phs, amom
  ENDIF ELSE begin
;   **** read the scattering parameters for the class as a whole for reuse
    RESTORE, out_path + '/scatfile.sav'
  ENDELSE
  IF KEYWORD_SET(scat_only) THEN RETURN,0 ; end of scattering calulations. stop here if scat_only true, start here if reuse_scat true.
   
  

;  We want the center point of the SRF integration which, since the number of ************DODGY
;  points is odd, will be at the center of the channel's spectral interval.
   m = Nwvl_Max / 2
  
   BextRatOUT = reform(BextRat[m,*,*])
   BextOUT = reform(Bext[m,*,*])
    SSAOUT = reform(w[m,*,*])
    GOUT = reform(g[m,*,*])


   print,''
   print,'Scattering parameters calculated FOR class '+out_path+' Version: '+Versions

;  -----------------------------------------------------------------------------
;  Run DISORT
;  -----------------------------------------------------------------------------

;  **** Setup the variables needed FOR the DISORT calls ****
   setup_disort, 60, NLayers, lutstr.saz_n, lutstr.raa_n, NMom

;  **** Define the LUT table output variables themselves
   RFD  = FLTARR(inststr.Number_of_Channels, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n)
   TFD  = FLTARR(inststr.Number_of_Channels, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n)
   RD   = FLTARR(inststr.Number_of_Channels, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n, lutstr.saz_n)
   TD   = FLTARR(inststr.Number_of_Channels, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n, lutstr.saz_n)
   TB   = FLTARR(inststr.Number_of_Channels, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n)
   RFBD = FLTARR(inststr.Number_of_Channels, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n)
   TFBD = FLTARR(inststr.Number_of_Channels, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n)
   RBD  = FLTARR(inststr.Number_of_Channels, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n, lutstr.saz_n, lutstr.raa_n)
   Em   = FLTARR(inststr.Number_of_Channels, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n, lutstr.saz_n)

;  **** Loop through the channels (and solar zenith angles) and run DISORT FOR
;       the beam and diffuse cases. Also produce the emissivity FOR the channels
;       that need it.

   FOR l = 0, inststr.number_of_nadir_channels - 1 DO BEGIN
      print,'Running DISORT for channel '+string(inststr.channelid[l],Format='(I)')+ ' (',strtrim(srfstrarr[l].wvl_centre,2),' um)'
      FOR k = 0,lutstr.prs_n-1 DO BEGIN	 
;     Define the Rayleigh scattering optical depth in each channel
      IF KEYWORD_SET(no_rayleigh) THEN $
        columntauray = replicate(1.e-6, inststr.number_of_channels) $
      ELSE $
        columntauray = (atmstr.pressure[atmstr.nlevels-1] / lutstr.prs(k)) / (117.03*srfstrarr[*].wvl_centre ^4 - 1.316*srfstrarr[*].wvl_centre ^2)	  
      FOR m = 0, srfstrarr[l].nwvl - 1 DO BEGIN
 
;        Do we have a Gas optical depth profile for the current channel?
;        If we don't have gas OPD for this channel,then use zeros (i.e. no gas
;        absorption). Remember that the gas OPD is defined on the levels between
;        each atmospheric layer, not on the layers themselves.
 
         IF N_ELEMENTS(gasstr) gt 0 then begin
            GasIndx = (WHERE(Fix(Gasstr[*].ChannelID) eq inststr.ChannelID[l],Count)) [0] 
			If (Count Eq 0) Then Stop, 'Gas absorption channel mismatch'
            GasLvl = gasstr[GasIndx].Tau_Gas
         endif ELSE begin         
           print, 'No gas optical depth profile.'
           GasLvl = replicate(0.0,atmstr.NLevels)    
         endelse      


;        Calculate the cumulative Rayleigh optical depth at each level
         RayLvl = ColumnTauRay[l]*exp(-0.1188*atmstr.Height - 0.00116*atmstr.Height^2)

;        The optical depths of each layer from gas absorption and Rayleigh
;        scattering are defined as the difference between the optical depth at
;        the adjacent levels.
         TauGas = GasLvl[lindgen(NLayers)+1] - GasLvl[lindgen(NLayers)]
         TauRay = RayLvl[lindgen(NLayers)+1] - RayLvl[lindgen(NLayers)] 


         FOR a = 0,lutstr.opd_n-1 DO BEGIN
           FOR r = 0,lutstr.efr_n-1 DO BEGIN

;              THE OPTICAL DEPTH FROM AEROSOL IS THE DESIRED TOTAL AODS FOR THE
;              ORAC LUT * THE RELATIVE AOD AT EACH LAYER FOR THIS CLASS * THE
;              SCALING FACTOR RELATING AOD AT THIS WAVELENGTH BACK TO 550 NM.
               tauscat = lutstr.opd[a] * scatreltau * bextrat[m,l,r]
;              the optical depths are additive
               dtau = taugas + tauray + tauscat
               totaltau = total(dtau)

;              The single scattering albedo is weighted by optical depth.
;              NB. SSA FOR Rayleigh scattering = 1, and is effectively 0 FOR
;              gas absorption.
               SSALB = (TauRay + w[m,l,r]*tauscat) / DTau
;              Now check that we have no SSALB values over 1.0 (this can happen in
;              layers with no absorption due to rounding). DISORT has an internal
;              check FOR this and will exit with an error code IF it fails.
               bd = WHERE(SSALB gt 1.0)
               IF bd[0] ge 0 then SSALB[bd] = 1.0;0.999999

;              ASYMMETRY PARAMETER IS ONLY NON-ZERO WHERE WE ACTUALLY HAVE AEROSOL
               ASYM = FLTARR(NLayers)
               nonzero = WHERE(tauscat gt 0.0)
               IF nonzero[0] ge 0 then ASYM[nonzero] = g[m,l,r]

;              NOW, USE THE GETMOM PROCEDURE (PART OF DISORT) TO GENERATE PHASE
;              FUNCTION MOMENTS FOR THE MOLECULAR SCATTERING AND THEN COMBINE WITH
;              THE AEROSOL MOMENTS GENERATED EARLIER.
               PMom = FLTARR(NMom, NLayers)
               FOR h=0,NLayers-1 do begin
                  IF ASYM[h] eq 0.0 then GETMOM, 2, 0.0, NMom-1, PM $
                  ELSE begin
                     GETMOM, 2, 0.0, NMom-1, mPM
                     PM = (mPM*TauRay[h] + AMom[*,m,l,r]*w[m,l,r]*tauscat[h]) / $
                          (TauRay[h] + w[m,l,r]*tauscat[h])
                  endelse
                  bd = WHERE(PM gt 1.0)
                  IF bd[0] ge 0 then PM[bd] = 1.0
                  PMom[*,h] = PM
               ENDFOR
							 
;*****							 printf,dlut, inststr.channelid[l],lutstr.opd[a],lutstr.efr[r],DTau, SSAlb, PMom
;*****							 print, inststr.channelid[l],lutstr.opd[a],lutstr.efr[r],DTau, SSAlb, PMom

;              We are now ready to call DISORT. Call the fast diffuse calculation
;              first (errors and problems are more likely to turn up quickly that
;              way).

               print, 'RT calculation for' + $
                      ' Channel: ' + string(inststr.channelid[l],Format='(I2)'), string(srfstrarr[l].wvl_centre,format='(" (",f6.3,"um)")') + $
                   ', SRF Point: ' + string(m,format='(i3)') + $
                         ', Tau: ' + string(lutstr.opd[a],format='(g10.2)') + $
                         ', EfR: ' + string(lutstr.efr[r],format='(f5.1)')
 ;              print, ''

 ;              print, '---------- DIFFUSE -----------'
               FBeam =   0.0           ; Direct beam intensity
               FIsot = 100.0           ; Isotropic illumination intensity
               UMu0  = cos(50.0*!dtor) ; A nominal value for beam zenith
               UTau  = [0.0, TotalTau] ; Define output layers (in terms of optical
                                       ; depth)
;              UMu is calculated to described downwelling as well as upwelling
;              radiance, as we also need the downwelling values
               UMu = FLTARR(2*lutstr.saz_n)
;              IF 90 degrees is included in the list of satellite zeniths, take a
;              small angle off of it FOR the DISORT calculations in order to avoid
;              numerical problems.
               tmpSat = lutstr.saz
               bd = WHERE(tmpSat eq 90.0)
               IF bd[0] ge 0 then tmpSat[bd] = 89.99
               UMu[0:lutstr.saz_n-1] = -1.0* cos(tmpSat*!dtor)
               UMu[lutstr.saz_n:2*lutstr.saz_n-1] = $
                  -1.0*UMu[lutstr.saz_n-lindgen(lutstr.saz_n)-1]
;              DISORT will spit the dummy IF there IF abs(UMu)=1
               bd = WHERE(abs(UMu) eq 1.0)
               IF bd[0] ge 0 then UMu[bd] = 0.99999 * UMu[bd]/abs(UMu[bd])

               call_disort, DTau, SSAlb, PMom, UTau, UMu, lutstr.raa, $
                            FBeam, UMu0, FISot, RFlDir, RFlDn, FlUp, $
                            dFdT, UAvg, UU, AlbMed, TrnMed

;              Generate the diffuse LUT variables
               RFD[l, k,  r, a] += (100. * FlUp[0]  / (FIsot*!pi)) * srfstrarr[l].val[m]
               TFD[l, k,  r, a] += (100. * RFlDn[1] / (FIsot*!pi)) * srfstrarr[l].val[m]
;              RD contains the Upwelling intensity
               RD[l, k,  r, a, *] += (UU[2*lutstr.saz_n-lindgen(lutstr.saz_n)-1,0,0]) * srfstrarr[l].val[m]
;              TD contains the Downwelling intensity without the direct beam

; V7 vs V8
              IF inststr.Solar_Channel_Flag[l] then $
						    TD[l, k,  r, a, *] += (UU[lindgen(lutstr.saz_n),  1,0] -100*exp(totaltau/UMu[0:lutstr.saz_n-1])) * srfstrarr[l].val[m] $
							Else $
  							TD[l, k,  r, a, *] += (UU[lindgen(lutstr.saz_n),  1,0]) * srfstrarr[l].val[m]

;               print, ''

;              IF the channel has the emission flag set, calculate the
;              emissivity.
               IF inststr.Thermal_Channel_Flag[l] then begin
;                  print, '---------- EMISSION ----------'
;                 Elisa's expression FOR emissivity....
;                 Em[l, k,  r, a, *]  = 100.0*(1.0 - w(l,r)) * $
;                                (1.0 - exp(TotalTau*(-1.0/cos(lutstr.saz*!dtor))))

;                 Use DISORT - I can't seem to get this to work correctly. FOR
;                 this we set both beam and diffuse input irradiances to zero
;                 and calculate emission across the window 1% of the nominal
;                 wavenumber (cm-1). Everything ELSE is the same as diffuse
;                 calculations.
                  FBeam   = 0.0 ; Direct beam intensity
                  FIsot   = 0.0 ; Isotropic illumination intensity
                  wn      = 1e4 / srfstrarr[L].wvl_centre
                  wnlo    = 0.995*wn
                  wnhi    = 1.005*wn
                  temp    = 270.0 ; hack from 250
                  incloud = WHERE(tauscat gt 0.0,emnly)
                  emTau   = DTau[incloud]
                  emSSA   = SSAlb[incloud]
                  emPMo   = PMom[*,incloud]
                  emUTau  = [0.0, total(emTau)]

                  call_disort, emTau, emSSA, emPMo, emUTau, UMu, lutstr.raa, $
                               FBeam, UMu0, FISot, RFlDir, RFlDn, FlUp, dFdT, $
                               UAvg, UU, AlbMed, TrnMed, /plank, wnlo=wnlo, $
                               wnhi=wnhi,temp=replicate(temp,emnly+1), nlayer=emnly

;                 Now we calculate the Plank emission across the wavelength
;                 interval.
                  BBE = PLKAVG(wnlo, wnhi, temp)

;                 Finally, combine to produce the emissivity
                  Em[l, k,  r, a, *] += (100.0 *  UU[2*lutstr.saz_n-lindgen(lutstr.saz_n)-1,0,0]/BBE) * srfstrarr[l].val[m]

;                  print, ''
               ENDIF

;              Now we loop over the solar zenith angles and do the direct beam
;              calculations. Note that this only needs to be done FOR channels with
;              a solar component to their signal.
               IF inststr.Solar_Channel_Flag[l] then begin
 ;                print, '---------- DIRECT ----------'
                  FOR s=0,lutstr.soz_n-1 do begin
;                    print, 'SZA: ', lutstr.soz[s]

                     FBeam = 100.0 ; Direct beam intensity
                     FIsot =   0.0 ; Isotropic illumination intensity

;                    IF a solar zenith angle of 90 degrees has been requested,
;                    alter it to something slightly smaller to prevent DISORT from
;                    crashing.
                     IF lutstr.soz[s] eq 90.0 then tmpSol = 89.99 $
                     ELSE tmpsol = lutstr.soz[s]
                     umu0  = cos(tmpsol*!dtor)

                     call_disort, dtau, ssalb, pmom, utau, umu, lutstr.raa, $
                                  fbeam, umu0, fisot, rfldir, rfldn, flup, $
                                  dfdt, uavg, uu, albmed, trnmed

;                    Generate the direct beam LUT variables
                     tb[l, k, r,a,s]   += (100. * rfldir[1] / rfldir[0]) * srfstrarr[l].val[m]
                     rfbd[l, k, r,a,s] += (100. * flup[0]   / rfldir[0]) * srfstrarr[l].val[m]
                     tfbd[l, k, r,a,s] += (100. * rfldn[1]  / rfldir[0]) * srfstrarr[l].val[m]
                     FOR p=0,lutstr.raa_n-1 do begin
;                       Reverse azimuth to ORAC convention
                        p2 = lutstr.raa_n - p - 1 ; WARNING this means raa must be evenly spaced from 0 to 180 otherwise the reversal doesn't make sense                        
;                       As with the diffuse case, RBD contains the upwelling
;                       intensity, while TBD contains the downwelling.
                        rbd[l, k,  r, a, s, *, p2] += (uu[2*lutstr.saz_n-lindgen(lutstr.saz_n)-1,0,p] * !pi) * srfstrarr[l].val[m]
                     ENDFOR
										 
                  ENDFOR

 ;                 print, ''
               ENDIF
            ENDFOR ; End of EfR loop						
         ENDFOR ; End of AOD loop
 ;        print,''
       ENDFOR ; End of SRF loop  	
     ENDFOR ; End of pressure loop 	 
   ENDFOR ; End of channel loop
	
;*****close,dlut

; Normalize the RT operators wrt to the SRF.  Note 'sum' will be unity in
; monochomatic mode (when no SRFs were provided).
  FOR l=0,inststr.Number_of_nadir_Channels-1 do begin
      sum = total(srfstrarr[l].val[*])
      RFD [l,*,*,*]       /= sum
      TFD [l,*,*,*]       /= sum
      RD  [l,*,*,*,*]     /= sum
      TD  [l,*,*,*,*]     /= sum
      TB  [l,*,*,*,*]     /= sum
      RFBD[l,*,*,*,*]     /= sum
      TFBD[l,*,*,*,*]     /= sum
      RBD [l,*,*,*,*,*,*] /= sum
      Em  [l,*,*,*,*]     /= sum
  ENDFOR
 
; Replicate nadir view to forward view (ie as though instrument has twice the number of channels)
  If ( inststr.View Gt 0) then begin              
    FOR l=0, inststr.number_of_nadir_channels - 1 do begin
      RFD [inststr.number_of_nadir_channels +l,*,*,*]       = RFD [l,*,*,*]
      TFD [inststr.number_of_nadir_channels +l,*,*,*]       = TFD [l,*,*,*] 
      RD  [inststr.number_of_nadir_channels +l,*,*,*,*]     = RD  [l,*,*,*,*] 
      TD  [inststr.number_of_nadir_channels +l,*,*,*,*]     = TD  [l,*,*,*,*]  
      TB  [inststr.number_of_nadir_channels +l,*,*,*,*]     = TB  [l,*,*,*,*]
      RFBD[inststr.number_of_nadir_channels +l,*,*,*,*]     = RFBD[l,*,*,*,*]
      TFBD[inststr.number_of_nadir_channels +l,*,*,*,*]     = TFBD[l,*,*,*,*]
      RBD [inststr.number_of_nadir_channels +l,*,*,*,*,*,*] = RBD [l,*,*,*,*,*,*]
      Em  [inststr.number_of_nadir_channels +l,*,*,*,*]     = Em  [l,*,*,*,*]
    ENDFOR
;   Rebuild instrument structure to account for slant channels  Note that all dual instrument devices
;   have on-board callibration so use (rua, rub and ruc) not (rgu and rou) 
    inststr = {  instrument_filename: inststr.instrument_filename, $
                           platform : inststr.platform,$
                         instrument : inststr.instrument,$
			     instrument_version : inststr.instrument_version, $
                     max_sat_zenith : inststr.max_sat_zenith,$    
                 Number_of_Channels : inststr.Number_of_Channels,$
                          ChannelID : [inststr.ChannelID           , inststr.ChannelID + inststr.view],$                        
                 Solar_Channel_Flag : [inststr.Solar_Channel_Flag  , inststr.Solar_Channel_Flag],$
                 Mixed_Channel_Flag : [inststr.Mixed_Channel_Flag  , inststr.Mixed_Channel_Flag],$
               Thermal_Channel_Flag : [inststr.Thermal_Channel_Flag, inststr.Thermal_Channel_Flag],$
                           srf_file : [inststr.srf_file            , inststr.srf_file],$
                              oldf0 : [inststr.oldf0               , inststr.oldf0],$                         
                              oldf1 : [inststr.oldf1               , inststr.oldf1],$                         
                            oldnefr : [inststr.oldnefr             , inststr.oldnefr],$                         
                             oldwvn : [inststr.oldwvn              , inststr.oldwvn],$                         
                              oldb1 : [inststr.oldb1               , inststr.oldb1],$                         
                              oldb2 : [inststr.oldb2               , inststr.oldb2],$                         
                              oldt1 : [inststr.oldt1               , inststr.oldt1],$                         
                              oldt2 : [inststr.oldt2               , inststr.oldt2],$                         
                            oldnebt : [inststr.oldnebt             , inststr.oldnebt],$ 
                                rua : [inststr.rua                 , inststr.rua],$
                                rub : [inststr.rub                 , inststr.rub],$
                                ruc : [inststr.ruc                 , inststr.ruc],$
                              refbt : [inststr.refbt               , inststr.refbt],$
                               nedt : [inststr.nedt                , inststr.nedt] }
;   extend srfstrarr to cover slant channels   
    srfstrarr = [srfstrarr,srfstrarr] 
  EndIf
 
 
; -----------------------------------------------------------------------------
; Create LUT.
; -----------------------------------------------------------------------------

  IF (File_Test(V2_LUT_Filename)) then print,'Info: Over-writing ' + V2_LUT_Filename else  print,'Info: Creating ' + V2_LUT_Filename

  write_v2_lut, V2_LUT_Filename, lutstr, inststr, srfstrarr, Vavg, bextout, bextratout, SSAOUT, GOUT, TD, TfD, RD, RfD, RBD = RBD, RfBD = RfBD, TfBD = TfBd, TB = TB, EM = EM,/include_pressure

;  -----------------------------------------------------------------------------
;  Output termination timestamp.
;  -----------------------------------------------------------------------------
  openw, lun, out_path + '/timestamp.txt', /get_lun
  printf, lun, strmid(timestamp(/utc), 0, 19) + 'Z'
  free_lun, lun

  print,'ORAC LUT generation completed for substance '+mmstr.substance

  return,0
end
