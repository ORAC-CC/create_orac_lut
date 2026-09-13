;pro value_check, X
;  Q= where(X Lt 0, Count)
;  If Count Gt 0 Then begin
;    print,'Zeroing '+string(Count)+' array elements'
;  endif
;end

pro write_v2_lut, V2_LUT_Filename, lutstr, inststr, srfstrarr, Vavg, BextOut, BextRatOut, SSAOUT, GOUT, TD, TfD, RD, RfD, RBD = RBD, RfBD= RfBD, TfBD = TfBd, TB = TB, EM = EM, include_pressure=include_pressure

  NumberofSolarChannels = Total(inststr.Solar_Channel_Flag)
  Solar_Channels_Exist  = boolean(NumberofSolarChannels)
  IF (Solar_Channels_Exist) THEN Solar_Index = WHERE(inststr.Solar_Channel_Flag)
  
  NumberofThermalChannels = Total(inststr.Thermal_Channel_Flag)
  Thermal_Channels_Exist  = boolean(NumberofThermalChannels)
  IF (Thermal_Channels_Exist) THEN Thermal_Index = WHERE(inststr.Thermal_Channel_Flag)

  NumberofMixedChannels = Total(inststr.Mixed_Channel_Flag)
  Mixed_Channels_Exist  = boolean(NumberofMixedChannels)
  IF (Mixed_Channels_Exist) THEN Mixed_Index = WHERE(inststr.Mixed_Channel_Flag)
 
 ; snr sut = 1
 ; rgu sut = 2
 ; rua sut = 3
 
  dontcare=where('SNR' eq tag_names(inststr), sut)
  If (not sut) then begin
     dontcare=where('RUA' eq tag_names(inststr), match)
     sut = match + 2
  endif

  file_id = ncdf_create(V2_LUT_Filename, /CLOBBER, /NETCDF4_FORMAT)
   
  ncdf_control, 0, /verbose
  
; Define the dimensions to be used 
  opd_dim = ncdf_dimdef(file_id, 'optical_depth'   , lutstr.opd_n)
  efr_dim = ncdf_dimdef(file_id, 'effective_radius', lutstr.efr_n)
  saz_dim = ncdf_dimdef(file_id, 'satellite_zenith', lutstr.saz_n)
  soz_dim = ncdf_dimdef(file_id, 'solar_zenith'    , lutstr.soz_n)
  raa_dim = ncdf_dimdef(file_id, 'relative_azimuth', lutstr.raa_n)
  if (keyword_set(include_pressure)) then prs_dim = ncdf_dimdef(file_id, 'surface_pressure', lutstr.prs_n)
  chn_dim = ncdf_dimdef(file_id, 'channels'        , inststr.Number_of_Channels)
  str_dim = NCDF_DIMDEF(file_id, 'length'          , MAX(STRLEN(inststr.srf_file)))
  st1_dim = ncdf_dimdef(file_id, 'st1'             , STRLEN(inststr.instrument_filename))
  st2_dim = ncdf_dimdef(file_id, 'st2'             , STRLEN(inststr.platform)) 
  st3_dim = ncdf_dimdef(file_id, 'st3'             , STRLEN(inststr.instrument))
  st4_dim = ncdf_dimdef(file_id, 'st4'             , STRLEN(inststr.instrument_version))
  
  IF (Solar_Channels_Exist  ) THEN soc_dim = ncdf_dimdef(file_id,'solar_channels'  , NumberofSolarChannels)
  IF (Thermal_Channels_Exist) THEN thc_dim = ncdf_dimdef(file_id,'thermal_channels', NumberofThermalChannels)  
  IF (Mixed_Channels_Exist)   THEN mxc_dim = ncdf_dimdef(file_id,'mixed_channels'  , NumberofMixedChannels)  
  
; instrument properties            
  instrument_filename_id  = ncdf_vardef(file_id,'instrument_filename',[st1_dim], /char) 
  platform_id             = ncdf_vardef(file_id,'platform',[st2_dim], /char) 
  instrument_id           = ncdf_vardef(file_id,'instrument',[st3_dim], /char) 
  instrument_version_id   = ncdf_vardef(file_id,'instrument_version',[st4_dim], /char) 
  max_sat_zenith_id       = ncdf_vardef(file_id,'max_sat_zenith', /float) 
    ncdf_attput, file_id, max_sat_zenith_id,'units','degrees'
    ncdf_attput, file_id, max_sat_zenith_id,'valid_range',[0.0,90.]
  number_of_channels_id   = ncdf_vardef(file_id,'number_of_channels', /long) 
    ncdf_attput, file_id, number_of_channels_id,'units','dimensionless'
    ncdf_attput, file_id, number_of_channels_id,'valid_range',[0,inststr.number_of_channels]
    
; channel specfic instrument properties    
  srf_file_id  = ncdf_vardef(file_id,'SRF_file', [str_dim,chn_dim],/char) 
    ncdf_attput, file_id, srf_file_id,'long_name','file containing the spectral response for the channel'

  channel_id_id = ncdf_vardef(file_id,'channel_id', [chn_dim],/short)
    ncdf_attput, file_id, channel_id_id,'long_name','Instrument channel identifier'
    ncdf_attput, file_id, channel_id_id,'units','dimensionless'           
    ncdf_attput, file_id, channel_id_id,'valid_range',[0, inststr.number_of_channels]  

  central_wavelength_id = ncdf_vardef(file_id,'central_wavelength', [chn_dim],/float) 
    ncdf_attput, file_id, central_wavelength_id,'long_name','effective central wavelength for the channel'
    ncdf_attput, file_id, central_wavelength_id,'units','microns'
    ncdf_attput, file_id, central_wavelength_id,'valid_range',[0.0,(machar()).xmax]
   
  central_wavenumber_id = ncdf_vardef(file_id,'central_wavenumber', [chn_dim],/float) 
    ncdf_attput, file_id, central_wavenumber_id,'long_name','effective central wavenumber for the channel'
    ncdf_attput, file_id, central_wavenumber_id,'units','cm^{-1}'
    ncdf_attput, file_id, central_wavenumber_id,'valid_range',[0.0,(machar()).xmax]
   
  solar_channel_flag_id = ncdf_vardef(file_id,'solar_channel_flag', [chn_dim],/short)
    ncdf_attput, file_id, solar_channel_flag_id,'long_name','Flag set to 1 if channel measures reflected solar radiation otherwise 0'
    ncdf_attput, file_id, solar_channel_flag_id,'units','dimensionless'           
    ncdf_attput, file_id, solar_channel_flag_id,'valid_range',[0, 1]    

  mixed_channel_flag_id = ncdf_vardef(file_id,'mixed_channel_flag', [chn_dim],/short)
    ncdf_attput, file_id, mixed_channel_flag_id,'long_name','Flag set to 1 if channel measures reflected solar and emitted infrared radiation otherwise 0'
    ncdf_attput, file_id, mixed_channel_flag_id,'units','dimensionless'           
    ncdf_attput, file_id, mixed_channel_flag_id,'valid_range',[0, 1]  
  
  thermal_channel_flag_id = ncdf_vardef(file_id,'thermal_channel_flag', [chn_dim],/short)
    ncdf_attput, file_id, thermal_channel_flag_id,'long_name','Flag set to 1 if channel measures emitted infrared radiation otherwise 0'
    ncdf_attput, file_id, thermal_channel_flag_id,'units','dimensionless'           
    ncdf_attput, file_id, thermal_channel_flag_id,'valid_range',[0, 1]  

  IF (Solar_Channels_Exist  ) THEN Begin
    soc_id = ncdf_vardef(file_id,'solar_channel_id', [soc_dim],/short)
      ncdf_attput, file_id, soc_id,'long_name','Instrument channel identifier for solar channels'
      ncdf_attput, file_id, soc_id,'units','dimensionless'
      ncdf_attput, file_id, soc_id,'valid_range',[0, inststr.number_of_channels]
     
    oldf0_id = ncdf_vardef(file_id,'oldf0', [soc_dim],/float) 
      ncdf_attput, file_id, oldf0_id,'long_name','SAD file F0'
      ncdf_attput, file_id, oldf0_id,'units','W/(m^2 um)'
      ncdf_attput, file_id, oldf0_id,'valid_range',[0.0,(machar()).xmax]
    oldf1_id = ncdf_vardef(file_id,'oldf1', [soc_dim],/float) 
      ncdf_attput, file_id, oldf1_id,'long_name','SAD file F1'  
      ncdf_attput, file_id, oldf1_id,'units','W/(m^2 um)'
      ncdf_attput, file_id, oldf1_id,'valid_range',[0.0,(machar()).xmax] 
      
    F0_id  = ncdf_vardef(file_id,'F0', [soc_dim],/float)
      ncdf_attput, file_id, F0_id,'long_name','in-band solar radiance'
      ncdf_attput, file_id, F0_id,'units','W/(m^2 um)'
      ncdf_attput, file_id, F0_id,'valid_range',[0.0,(machar()).xmax]    
 
    Case sut of
    1: begin
       snr_id = ncdf_vardef(file_id,'snr', [soc_dim],/float)
        ncdf_attput, file_id, snr_id,'long_name','signal-to-noise ratio'
        ncdf_attput, file_id, snr_id,'units','dimensionless'
        ncdf_attput, file_id, snr_id,'valid_range',[0.0,(machar()).xmax]   
       end
    2: begin
      rgu_id = ncdf_vardef(file_id,'rgu', [soc_dim],/float)
        ncdf_attput, file_id, rgu_id,'long_name','radiance gain uncertainty'
        ncdf_attput, file_id, rgu_id,'units','W/(m^2 sr um count)'
        ncdf_attput, file_id, rgu_id,'valid_range',[0.0,(machar()).xmax]
      rou_id = ncdf_vardef(file_id,'rou', [soc_dim],/float)
        ncdf_attput, file_id, rou_id,'long_name','radiance offset uncertainty'
        ncdf_attput, file_id, rou_id,'units','W/(m^2 sr um)'
        ncdf_attput, file_id, rou_id,'valid_range',[0.0,(machar()).xmax]    
       end
    3: begin
       rua_id = ncdf_vardef(file_id,'rua', [soc_dim],/float)
        ncdf_attput, file_id, rua_id,'long_name','radiance uncertainty coefficient a'
        ncdf_attput, file_id, rua_id,'units','dimensionless'
      rub_id = ncdf_vardef(file_id,'rub', [soc_dim],/float)
        ncdf_attput, file_id, rub_id,'long_name','radiance uncertainty coefficient b'
        ncdf_attput, file_id, rub_id,'units','[W/(m^2 sr um)]^(1/2)'
      ruc_id = ncdf_vardef(file_id,'ruc', [soc_dim],/float)
        ncdf_attput, file_id, ruc_id,'long_name','radiance uncertainty coefficient c'
        ncdf_attput, file_id, ruc_id,'units','[W/(m^2 sr um)]'   
       end
    EndCase
; ***** TEMPORARY FOR BACK COMPATIBILITY ******  
    oldnefr_id = ncdf_vardef(file_id,'oldnefr', [chn_dim],/float) 
      ncdf_attput, file_id, oldnefr_id,'long_name','SAD file NeFr (measurement uncertainty)'
      ncdf_attput, file_id, oldnefr_id,'units','W/(m^2 sr um)'
      ncdf_attput, file_id, oldnefr_id,'valid_range',[0.0,(machar()).xmax] 
; *********************************************
  ENDIF
  
    IF (Mixed_Channels_Exist) THEN begin
    mxc_id = ncdf_vardef(file_id,'mixed_channel_id', [mxc_dim],/short)
      ncdf_attput, file_id, mxc_id,'long_name','Instrument channel identifier for mixed channels'
      ncdf_attput, file_id, mxc_id,'units','dimensionless'
      ncdf_attput, file_id, mxc_id,'valid_range',[0, inststr.number_of_channels]   
  ENDIF
 

  IF (Thermal_Channels_Exist) THEN begin
    thc_id = ncdf_vardef(file_id,'thermal_channel_id', [thc_dim],/short)
      ncdf_attput, file_id, thc_id,'long_name','Instrument channel identifier for thermal channels'
      ncdf_attput, file_id, thc_id,'units','dimensionless'
      ncdf_attput, file_id, thc_id,'valid_range',[0, inststr.number_of_channels]  
    refbt_id  = ncdf_vardef(file_id,'refbt', [thc_dim],/float)
      ncdf_attput, file_id, refbt_id,'long_name','reference brightness temperate at neDT has been calculated'
      ncdf_attput, file_id, refbt_id,'units','K'
      ncdf_attput, file_id, refbt_id,'valid_range',[100, 500]     
    nedt_id  = ncdf_vardef(file_id,'nedt', [thc_dim],/float)
      ncdf_attput, file_id, nedt_id,'long_name','noise equivalent delta temperature'
      ncdf_attput, file_id, nedt_id,'units','K'
      ncdf_attput, file_id, nedt_id,'valid_range',[0,(machar()).xmax]     
    B1_id  = ncdf_vardef(file_id,'B1', [thc_dim],/float)
    B2_id  = ncdf_vardef(file_id,'B2', [thc_dim],/float)
    T1_id  = ncdf_vardef(file_id,'T1', [thc_dim],/float)
    T2_id  = ncdf_vardef(file_id,'T2', [thc_dim],/float)
; ***** TEMPORARY FOR BACK COMPATIBILITY ******  
    oldwvn_id = ncdf_vardef(file_id,'oldwvn', [chn_dim],/float) 
      ncdf_attput, file_id, oldwvn_id,'long_name','SAD file Wvn'
      ncdf_attput, file_id, oldwvn_id,'units','1/cm'
      ncdf_attput, file_id, oldwvn_id,'valid_range',[0.0,(machar()).xmax]   
    oldb1_id = ncdf_vardef(file_id,'oldb1', [chn_dim],/float) 
      ncdf_attput, file_id, oldb1_id,'long_name','SAD file B1' 
    oldb2_id = ncdf_vardef(file_id,'oldb2', [chn_dim],/float) 
      ncdf_attput, file_id, oldb2_id,'long_name','SAD file B2'
    oldt1_id = ncdf_vardef(file_id,'oldt1', [chn_dim],/float) 
      ncdf_attput, file_id, oldt1_id,'long_name','SAD file T1'
    oldt2_id = ncdf_vardef(file_id,'oldt2', [chn_dim],/float) 
      ncdf_attput, file_id, oldt2_id,'long_name','SAD file T2'
    oldnebt_id = ncdf_vardef(file_id,'oldnebt', [chn_dim],/float) 
      ncdf_attput, file_id, oldnebt_id,'long_name','SAD file NeBT'
      ncdf_attput, file_id, oldnebt_id,'units','K'
      ncdf_attput, file_id, oldnebt_id,'valid_range',[0.0,(machar()).xmax]
; *********************************************
  ENDIF
	
; MICROPHYSICAL properties
		
  vavg_id = ncdf_vardef(file_id,'average_volume_per_particle', [efr_dim])  
    ncdf_attput, file_id, vavg_id,'long_name','average volume per particle'
    ncdf_attput, file_id, vavg_id,'units','to be investigated'
    ncdf_attput, file_id, vavg_id,'valid_range',[0.0,(machar()).xmax]
    
  bext_id = ncdf_vardef(file_id,'extinction_coefficient', [chn_dim, efr_dim]) 
    ncdf_attput, file_id, bext_id,'long_name','volume extinction coefficient'
    ncdf_attput, file_id, bext_id,'units','to be investigated'
    ncdf_attput, file_id, bext_id,'valid_range',[0.0,(machar()).xmax]
   
  bextrat_id =  ncdf_vardef(file_id,'extinction_coefficient_ratio', [chn_dim, efr_dim ]) 
    ncdf_attput, file_id, bextrat_id,'long_name','ratio of volume extinction coefficient to the volume extinction coefficient at 550 nm'
    ncdf_attput, file_id, bextrat_id,'units','dimensionless'
    ncdf_attput, file_id, bextrat_id,'valid_range',[0.0,(machar()).xmax]
   
  wid = ncdf_vardef(file_id,'single_scatter_albedo', [chn_dim, efr_dim]) 
    ncdf_attput, file_id, wid,'long_name','single scatter albedo'
    ncdf_attput, file_id, wid,'units','dimensionless'
    ncdf_attput, file_id, wid,'valid_range',[0.0,1.0]
   
  gid = ncdf_vardef(file_id,'asymmetry_parameter', [chn_dim, efr_dim]) 
    ncdf_attput, file_id, gid,'long_name','asymmetry parameter'
    ncdf_attput, file_id, gid,'units','dimensionless'
    ncdf_attput, file_id, gid,'valid_range',[-1.0,1.0]	
	

; LUT axes 
  opd_id = ncdf_vardef(file_id,'optical_depth', [opd_dim],/float)
    ncdf_attput, file_id, opd_id,'long_name','optical depth'
    ncdf_attput, file_id, opd_id,'spacing',lutstr.opd_spacing 
    ncdf_attput, file_id, opd_id,'units','dimensionless'
    ncdf_attput, file_id, opd_id,'valid_range',[0.0,(machar()).xmax]
   
  efr_id = ncdf_vardef(file_id,'effective_radius', [efr_dim],/float)
    ncdf_attput, file_id, efr_id,'long_name','particle effective radius'
    ncdf_attput, file_id, efr_id,'spacing',lutstr.efr_spacing 
    ncdf_attput, file_id, efr_id,'units','microns'
    ncdf_attput, file_id, efr_id,'valid_range',[0.0,(machar()).xmax]
   
  saz_id = ncdf_vardef(file_id,'satellite_zenith', [saz_dim],/float)
    ncdf_attput, file_id, saz_id,'long_name','satellite zenith angle'
    ncdf_attput, file_id, saz_id,'spacing',lutstr.saz_spacing 
    ncdf_attput, file_id, saz_id,'units','degrees'
    ncdf_attput, file_id, saz_id,'valid_range',[0.0, 180.0]

  soz_id = ncdf_vardef(file_id,'solar_zenith', [soz_dim],/float)
    ncdf_attput, file_id, soz_id,'long_name','solar zenith angle'
    ncdf_attput, file_id, soz_id,'spacing',lutstr.soz_spacing 
    ncdf_attput, file_id, soz_id,'units','degrees'
    ncdf_attput, file_id, soz_id,'valid_range',[0.0, 180.0]

  raa_id = ncdf_vardef(file_id,'relative_azimuth', [raa_dim],/float)
    ncdf_attput, file_id, raa_id,'long_name','satellite azimuth relative to the Sun'
    ncdf_attput, file_id, raa_id,'spacing',lutstr.raa_spacing 
    ncdf_attput, file_id, raa_id,'units','degrees'
    ncdf_attput, file_id, raa_id,'valid_range',[0.0, 180.0]
  if (keyword_set(include_pressure)) then begin
    prs_id = ncdf_vardef(file_id,'surface_pressure', [prs_dim],/float)
      ncdf_attput, file_id, prs_id,'long_name','surface pressure'
      ncdf_attput, file_id, prs_id,'spacing',lutstr.prs_spacing 
      ncdf_attput, file_id, prs_id,'units','hPa'
      ncdf_attput, file_id, prs_id,'valid_range',[900., 1100.]
  endif



; LUT Not 100% sure on the long names nor the valid range.  I think some of these can be bigger than 1
  if (keyword_set(include_pressure)) then $
    TDid = ncdf_vardef(file_id,'T_dv',  [chn_dim, prs_dim, efr_dim, opd_dim, saz_dim],/float) $; IT_dv = ITd 
  else $
    TDid = ncdf_vardef(file_id,'T_dv',  [chn_dim, efr_dim, opd_dim, saz_dim],/float) ; IT_dv = ITd  
      ncdf_attput, file_id, TDid,'long_name','diffuse transmission of direct light'
      ncdf_attput, file_id, TDid,'units','dimensionless'
      ncdf_attput, file_id, TDid,'valid_range',[0.0, 1.0]
  if (keyword_set(include_pressure)) then $
    TfDid = ncdf_vardef(file_id,'T_dd',  [chn_dim, prs_dim, efr_dim, opd_dim],/float)  $        ; IT_dd = ITfd 
  else $
    TfDid = ncdf_vardef(file_id,'T_dd',  [chn_dim, efr_dim, opd_dim],/float)          ; IT_dd = ITfd 
      ncdf_attput, file_id, TfDid,'long_name','diffuse transmission'
      ncdf_attput, file_id, TfDid,'units','dimensionless'
      ncdf_attput, file_id, TfDid,'valid_range',[0.0, 1.0]
  if (keyword_set(include_pressure)) then $
    RDid = ncdf_vardef(file_id,'R_dv',  [chn_dim, prs_dim, efr_dim, opd_dim, saz_dim],/float) $; IR_dv = IRd
  else $
    RDid = ncdf_vardef(file_id,'R_dv',  [chn_dim, efr_dim, opd_dim, saz_dim],/float) ; IR_dv = IRd
      ncdf_attput, file_id, RDid,'long_name','direct reflection of diffuse light'
      ncdf_attput, file_id, RDid,'units','dimensionless'
      ncdf_attput, file_id, RDid,'valid_range',[0.0, 1.0]
  if (keyword_set(include_pressure)) then $
    RFDid = ncdf_vardef(file_id,'R_dd',  [chn_dim, prs_dim, efr_dim, opd_dim],/float)  $        ; IR_dd = IRfd
  else $
    RFDid = ncdf_vardef(file_id,'R_dd',  [chn_dim, efr_dim, opd_dim],/float)          ; IR_dd = IRfd
      ncdf_attput, file_id, RFDid,'long_name','diffuse reflection of diffuse light'
      ncdf_attput, file_id, RFDid,'units','dimensionless'
      ncdf_attput, file_id, RFDid,'valid_range',[0.0, 1.0]
  IF (Solar_Channels_Exist) THEN begin
    if (keyword_set(include_pressure)) then $
      RBDid = ncdf_vardef(file_id,'R_0v', [soc_dim, prs_dim, efr_dim, opd_dim, soz_dim, saz_dim, raa_dim],/float) $; IR_0v = IRbd
    else $
      RBDid = ncdf_vardef(file_id,'R_0v', [soc_dim, efr_dim, opd_dim, soz_dim, saz_dim, raa_dim],/float) ; IR_0v = IRbd
       ncdf_attput, file_id, RBDid,'long_name','bi-directional reflectance'
       ncdf_attput, file_id, RBDid,'units','dimensionless'
       ncdf_attput, file_id, RBDid,'valid_range',[0.0, 1.0]
    if (keyword_set(include_pressure)) then $
      RFBDid = ncdf_vardef(file_id,'R_0d', [soc_dim, prs_dim, efr_dim, opd_dim, soz_dim],/float)    $               ; IR_0d = IRfbd
    else $
      RFBDid = ncdf_vardef(file_id,'R_0d', [soc_dim, efr_dim, opd_dim, soz_dim],/float)                   ; IR_0d = IRfbd
       ncdf_attput, file_id, RFBDid,'long_name','diffuse reflectance of direct beam'
       ncdf_attput, file_id, RFBDid,'units','dimensionless'
       ncdf_attput, file_id, RFBDid,'valid_range',[0.0, 1.0]
    if (keyword_set(include_pressure)) then $
      TFBDid = ncdf_vardef(file_id,'T_0d', [soc_dim, prs_dim, efr_dim, opd_dim, soz_dim],/float)  $                 ; IT_0d = ITfbd
    else $
      TFBDid = ncdf_vardef(file_id,'T_0d', [soc_dim, efr_dim, opd_dim, soz_dim],/float)                   ; IT_0d = ITfbd
       ncdf_attput, file_id, TFBDid,'long_name','diffuse transmission of diffuse light'
       ncdf_attput, file_id, TFBDid,'units','dimensionless'
       ncdf_attput, file_id, TFBDid,'valid_range',[0.0, 1.0] 
    if (keyword_set(include_pressure)) then $
      TBid = ncdf_vardef(file_id,'T_00', [soc_dim, prs_dim, efr_dim, opd_dim, soz_dim],/float)  $                   ; IT_00 = ITb
    else $
      TBid = ncdf_vardef(file_id,'T_00', [soc_dim, efr_dim, opd_dim, soz_dim],/float)                     ; IT_00 = ITb
       ncdf_attput, file_id, TBid,'long_name','direct transmission'
       ncdf_attput, file_id, TBid,'units','dimensionless'
       ncdf_attput, file_id, TBid,'valid_range',[0.0, 1.0]
  ENDIF  
     
  IF (Thermal_Channels_Exist) THEN BEGIN
    if (keyword_set(include_pressure)) then $
      EMid = ncdf_vardef(file_id,'E_md',   [thc_dim, prs_dim, efr_dim, opd_dim, saz_dim],/float) $
    else $
      EMid = ncdf_vardef(file_id,'E_md',   [thc_dim, efr_dim, opd_dim, saz_dim],/float) 
     ncdf_attput, file_id, EMid,'long_name','diffuse emissivity'
     ncdf_attput, file_id, EMid,'units','dimensionless'
     ncdf_attput, file_id, EMid,'valid_range',[0.0, 1.0]
  ENDIF
        
  ncdf_control, file_id, /endef
	
; output the axis values   
  ncdf_varput, file_id, opd_id, lutstr.opd
  ncdf_varput, file_id, efr_id, lutstr.efR
  ncdf_varput, file_id, saz_id, lutstr.saz
  ncdf_varput, file_id, soz_id, lutstr.soz
  ncdf_varput, file_id, raa_id, lutstr.raa
  if (keyword_set(include_pressure)) then ncdf_varput, file_id, prs_id, lutstr.prs
 
 ; instrument properties    
  ncdf_varput, file_id, instrument_filename_id, inststr.instrument_filename
  ncdf_varput, file_id, platform_id, inststr.platform
  ncdf_varput, file_id, instrument_id, inststr.instrument
  ncdf_varput, file_id, instrument_version_id, inststr.instrument_version
  ncdf_varput, file_id, max_sat_zenith_id, inststr.max_sat_zenith
  ncdf_varput, file_id, number_of_channels_id, inststr.number_of_channels
  ncdf_varput, file_id, srf_file_id, inststr.srf_file
  ncdf_varput, file_id, channel_id_id, inststr.ChannelID   
  ncdf_varput, file_id, central_wavelength_id, srfstrarr.wvl_centre
  ncdf_varput, file_id, central_wavenumber_id, srfstrarr.wvn_centre
  ncdf_varput, file_id, solar_channel_flag_id, inststr.solar_channel_flag
  ncdf_varput, file_id, mixed_channel_flag_id, inststr.mixed_channel_flag
  ncdf_varput, file_id, thermal_channel_flag_id, inststr.thermal_channel_flag  
  IF (Solar_Channels_Exist)   THEN begin
    ncdf_varput, file_id, soc_id, inststr.ChannelID(Solar_Index)
    ncdf_varput, file_id, oldf0_id, inststr.oldf0(Solar_Index)
    ncdf_varput, file_id, oldf1_id, inststr.oldf1(Solar_Index)
    ncdf_varput, file_id, F0_id, srfstrarr(Solar_Index).F0
    Case sut of
    1: ncdf_varput, file_id, snr_id, inststr.snr(Solar_Index) 
    2: begin
         ncdf_varput, file_id, rgu_id, inststr.rgu(Solar_Index)
         ncdf_varput, file_id, rou_id, inststr.rou(Solar_Index)    
       end
    3: begin
         ncdf_varput, file_id, rua_id, inststr.rua(Solar_Index)  
         ncdf_varput, file_id, rub_id, inststr.rub(Solar_Index)  
         ncdf_varput, file_id, ruc_id, inststr.ruc(Solar_Index)    
      end
    EndCase
; ***** TEMPORARY FOR BACK COMPATIBILITY ******  
    ncdf_varput, file_id, oldnefr_id, inststr.oldnefr(Solar_Index)
; *********************************************
  ENDIF
  IF (Mixed_Channels_Exist)   THEN begin
    ncdf_varput, file_id, mxc_id, inststr.ChannelID(Mixed_Index)
  ENDIF
  IF (Thermal_Channels_Exist) THEN BEGIN 
    ncdf_varput, file_id, thc_id, inststr.ChannelID(Thermal_Index)
    ncdf_varput, file_id, refbt_id, inststr.refbt(Thermal_Index)
    ncdf_varput, file_id, nedt_id, inststr.nedt(Thermal_Index)
    ncdf_varput, file_id, B1_id, srfstrarr(Thermal_Index).B1
    ncdf_varput, file_id, B2_id, srfstrarr(Thermal_Index).B2
    ncdf_varput, file_id, T1_id, srfstrarr(Thermal_Index).T1
    ncdf_varput, file_id, T2_id, srfstrarr(Thermal_Index).T2
; ***** TEMPORARY FOR BACK COMPATIBILITY ******  
    ncdf_varput, file_id,  oldwvn_id, inststr.oldwvn(Thermal_Index)
    ncdf_varput, file_id,   oldb1_id, inststr.oldb1(Thermal_Index)
    ncdf_varput, file_id,   oldb2_id, inststr.oldb2(Thermal_Index)
    ncdf_varput, file_id,   oldt1_id, inststr.oldt1(Thermal_Index)
    ncdf_varput, file_id,   oldt2_id, inststr.oldt2(Thermal_Index)
    ncdf_varput, file_id, oldnebt_id, inststr.oldnebt(Thermal_Index)
; *********************************************
  ENDIF
  
;  -----------------------------------------------------------------------------
;  Output reflectance, transmission and emission data into ORAC LUT
;  -----------------------------------------------------------------------------
 ; find a tidier way of  /100 
  ncdf_varput, file_id,   TDid, TD /100
  ncdf_varput, file_id,  TfDid, TfD/100
  ncdf_varput, file_id,   RDid, RD /100
  ncdf_varput, file_id,  RfDid, RfD/100
  IF (Solar_Channels_Exist) THEN begin
    if (keyword_set(include_pressure)) then begin
      ncdf_varput, file_id,  RBDid,  RBD[Solar_Index,*,*,*,*,*,*]/100
      ncdf_varput, file_id, RfBDid, RfBD[Solar_Index,*,*,*,*]    /100
      ncdf_varput, file_id, TfBDid, TfBD[Solar_Index,*,*,*,*]    /100
      ncdf_varput, file_id,   TBid,   TB[Solar_Index,*,*,*,*]    /100
	endif else begin
      ncdf_varput, file_id,  RBDid,  RBD[Solar_Index,*,*,*,*,*]/100
      ncdf_varput, file_id, RfBDid, RfBD[Solar_Index,*,*,*]    /100
      ncdf_varput, file_id, TfBDid, TfBD[Solar_Index,*,*,*]    /100
      ncdf_varput, file_id,   TBid,   TB[Solar_Index,*,*,*]    /100	
	endelse
  ENDIF
  IF (Thermal_Channels_Exist) THEN $
    if (keyword_set(include_pressure)) then $
      ncdf_varput, file_id,   EMid,   EM[Thermal_Index,*,*,*,*]/100 $
	else  $
      ncdf_varput, file_id,   EMid,   EM[Thermal_Index,*,*,*]/100 
	
  
;  -----------------------------------------------------------------------------
;  Output scattering data into ORAC LUT
;  -----------------------------------------------------------------------------

;  Write the Vavg LUT
   ncdf_varput, file_id, vavg_id, Vavg
   ncdf_varput, file_id, bext_id, BextOut
   ncdf_varput, file_id, bextrat_id, BextRatOut
   ncdf_varput, file_id, wid, SSAOUT   
   ncdf_varput, file_id, gid, GOUT
   
; All done, close the file.
  ncdf_close, file_id
end