  pro read_V2_LUT,  V2filename,  V2_LUT
  
  fid = ncdf_open(V2filename)

  Solar_Channels_Exist   = Boolean(NCDF_DIMID(fid,'solar_channels')  +1)
  Mixed_Channels_Exist   = Boolean(NCDF_DIMID(fid,'mixed_channels')  +1)
  Thermal_Channels_Exist = Boolean(NCDF_DIMID(fid,'thermal_channels')+1)
  
  ncdf_varget, fid, ncdf_varid(fid,'instrument_filename') , instrument_filename
  ncdf_varget, fid, ncdf_varid(fid,'platform')            , platform
  ncdf_varget, fid, ncdf_varid(fid,'instrument')          , instrument
  ncdf_varget, fid, ncdf_varid(fid,'max_sat_zenith')      , max_sat_zenith
  ncdf_varget, fid, ncdf_varid(fid,'number_of_channels')  , number_of_channels
  ncdf_varget, fid, ncdf_varid(fid,'SRF_file')            , srf_file
  
  ncdf_varget, fid, ncdf_varid(fid,'channel_id')          ,  channel_id
  ncdf_varget, fid, ncdf_varid(fid,'central_wavelength')  , central_wavelength
  ncdf_varget, fid, ncdf_varid(fid,'central_wavenumber')  , central_wavenumber
  ncdf_varget, fid, ncdf_varid(fid,'solar_channel_flag')  , solar_channel_flag
  ncdf_varget, fid, ncdf_varid(fid,'mixed_channel_flag')  , mixed_channel_flag
  ncdf_varget, fid, ncdf_varid(fid,'thermal_channel_flag'), thermal_channel_flag
 
  If (Solar_Channels_Exist) then begin
    ncdf_varget, fid, ncdf_varid(fid,'solar_channel_id')  , solar_channel_id
    ncdf_varget, fid, ncdf_varid(fid,'oldf0')             , oldf0
    ncdf_varget, fid, ncdf_varid(fid,'oldf1')             , oldf1
    ncdf_varget, fid, ncdf_varid(fid,'F0')                , F0

    snr_id = ncdf_varid(fid,'snr')
    If (snr_id ne -1) then ncdf_varget, fid, snr_id,  snr $
	else begin
      rgu_id = ncdf_varid(fid,'rgu')
      If (rgu_id ne -1) then begin
        ncdf_varget, fid, rgu_id                          , rgu
        ncdf_varget, fid, ncdf_varid(fid,'rou')           , rou
      endif else begin
        ncdf_varget, fid, ncdf_varid(fid,'rua')           , rua
        ncdf_varget, fid, ncdf_varid(fid,'rub')           , rub
        ncdf_varget, fid, ncdf_varid(fid,'rua')           , rua
 
      endelse
    endelse
  EndIf
  ncdf_varget, fid, ncdf_varid(fid,'oldnefr')             , oldnefr
  If (Mixed_Channels_Exist) then $
    ncdf_varget, fid, ncdf_varid(fid,'mixed_channel_id')  , mixed_channel_id 
  If (Thermal_Channels_Exist) then begin
    ncdf_varget, fid, ncdf_varid(fid,'thermal_channel_id'), thermal_channel_id  
    
    ncdf_varget, fid, ncdf_varid(fid,'refbt')             , refbt  
    ncdf_varget, fid, ncdf_varid(fid,'nedt')              , nedt  
    ncdf_varget, fid, ncdf_varid(fid,'B1')                , B1  
    ncdf_varget, fid, ncdf_varid(fid,'B2')                , B2  
    ncdf_varget, fid, ncdf_varid(fid,'T1')                , T1  
    ncdf_varget, fid, ncdf_varid(fid,'T2')                , T2  
    ncdf_varget, fid, ncdf_varid(fid,'oldwvn')            , oldwvn  
    ncdf_varget, fid, ncdf_varid(fid,'oldb1')             , oldb1  
    ncdf_varget, fid, ncdf_varid(fid,'oldb2')             , oldb2  
    ncdf_varget, fid, ncdf_varid(fid,'oldt1')             , oldt1  
    ncdf_varget, fid, ncdf_varid(fid,'oldt2')             , oldt2  
    ncdf_varget, fid, ncdf_varid(fid,'oldnebt')           , oldnebt  
    ncdf_varget, fid, ncdf_varid(fid,'thermal_channel_id'), thermal_channel_id  
  EndIf
  
  ncdf_varget, fid, ncdf_varid(fid,'optical_depth')          ,  optical_depth
  ncdf_varget, fid, ncdf_varid(fid,'effective_radius')          ,  effective_radius
  ncdf_varget, fid, ncdf_varid(fid,'satellite_zenith')          ,  satellite_zenith
  ncdf_varget, fid, ncdf_varid(fid,'solar_zenith')          ,  solar_zenith
  ncdf_varget, fid, ncdf_varid(fid,'relative_azimuth')          ,  relative_azimuth
  

  LUT_ID = ncdf_varid(fid,'average_volume_per_particle')
  ncdf_varget, fid, LUT_ID, Vavg     
  LUT_ID = ncdf_varid(fid,'extinction_coefficient')
  ncdf_varget, fid, LUT_ID, Bext     
  LUT_ID = ncdf_varid(fid,'extinction_coefficient_ratio')
  ncdf_varget, fid, LUT_ID, BextRat          
  LUT_ID = ncdf_varid(fid,'single_scatter_albedo')
  ncdf_varget, fid, LUT_ID, ssa   
  LUT_ID = ncdf_varid(fid,'asymmetry_parameter')
  ncdf_varget, fid, LUT_ID, asym    

   LUT_ID = ncdf_varid(fid,'T_dv')
   ncdf_varget, fid, LUT_ID, TD
   LUT_ID = ncdf_varid(fid,'T_dd')
   ncdf_varget, fid, LUT_ID, TfD
   LUT_ID = ncdf_varid(fid,'R_dv')
   ncdf_varget, fid, LUT_ID, RD
   LUT_ID = ncdf_varid(fid,'R_dd')
   ncdf_varget, fid, LUT_ID, RfD
   LUT_ID = ncdf_varid(fid,'T_00')
   ncdf_varget, fid, LUT_ID, TB
   LUT_ID = ncdf_varid(fid,'R_0v')
   ncdf_varget, fid, LUT_ID, RBD
   LUT_ID = ncdf_varid(fid,'T_0d')
   ncdf_varget, fid, LUT_ID, TfBD
   LUT_ID = ncdf_varid(fid,'R_0d')
   ncdf_varget, fid, LUT_ID, RfBD
   LUT_ID = ncdf_varid(fid,'E_md')
   ncdf_varget, fid, LUT_ID, EM

; All done, close the file.
  ncdf_close, fid
  
  V2_LUT ={instrument_filename: string(instrument_filename), $
                      platform: string(platform), $
                    instrument: string(instrument), $
                max_sat_zenith: max_sat_zenith, $
            number_of_channels: number_of_channels, $
                      srf_file: string(srf_file)       , $                        
                    channel_id: string(channel_id), $
            central_wavelength: central_wavelength, $
            central_wavenumber: central_wavenumber, $
           solar_channel_flag: solar_channel_flag, $
           mixed_channel_flag: mixed_channel_flag, $
         thermal_channel_flag: thermal_channel_flag, $
               optical_depths: N_elements(optical_depth), $
                optical_depth: optical_depth, $
              effective_radii: N_elements(effective_radius), $  
             effective_radius: effective_radius, $  
            satellite_zeniths: N_elements(satellite_zenith), $	   
             satellite_zenith: satellite_zenith, $	   
                solar_zeniths: N_elements(solar_zenith), $
                 solar_zenith: solar_zenith, $
            relative_azimuths: N_elements(relative_azimuth), $   
             relative_azimuth: relative_azimuth, $   
		                 Vavg: Vavg, $
		                 Bext: Bext, $
		              BextRat: BextRat, $
		                    w: ssa, $
		                    g: asym, $
                        oldf0: oldf0  , $
                        oldf1: oldf1, $
	              	       F0: F0, $	
              	      oldnefr:oldnefr, $
              	        refbt: refbt  , $
              	         nedt: nedt , $ 
              	           B1: B1 , $ 
              	           B2: B2, $  
                           T1: T1, $  
                           T2: T2, $  
	                   oldwvn: oldwvn  , $
	                    oldb1: oldb1 , $ 
	                    oldb2: oldb2 , $ 
	                    oldt1: oldt1 , $ 
	                    oldt2: oldt2, $  
	                  oldnebt: oldnebt , $
                           TD: TD, $
                          TfD: TfD, $
                           RD: RD, $ 
                          RfD: RfD, $
                           TB: TB, $
	                      RBD: RBD, $
	                     TfBD: TfBD, $
	                     RfBD: RfBD, $
	                       EM: EM  }
end

