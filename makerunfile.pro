; make run file to create script to form one or multiple luts

; THESE ARE THE MASTER LISTS DO NOT DELETE
  ; MM = ['aerosol_a70',$
        ; 'aerosol_a75',$
        ; 'aerosol_a76',$
        ; 'aerosol_a77',$
        ; 'aerosol_a78',$  
        ; 'aerosol_a79',$  
        ; 'biomass_a',$
        ; 'biomass_b',$
        ; 'biomass_c',$
        ; 'biomass_d',$
        ; 'biomass_australia',$
        ; 'biomass_indonesia',$
        ; 'liquid-water_old',$
        ; 'liquid-water_stg',$
        ; 'liquid-water_240',$
        ; 'liquid-water_253',$
        ; 'liquid-water_263',$
        ; 'liquid-water_273',$
        ; 'sulphuric-acid_300',$  	
        ; 'volcanic-ash_ctn',$
        ; 'volcanic-ash_ey1',$
        ; 'volcanic-ash_spr', $
        ; 'water-ice_sph',$
        ; 'water-ice_agg',$
        ; 'water-ice_ghm',$
        ; 'water-ice_src']
 
  ; platform = 'aqua'
  ; instrument = 'modis'  

  ; platform = 'earthcare'
  ; instrument = 'msi'
  
  ; platform = 'envisat'
  ; instrument = 'aatsr'   
   
  ; platform = 'ers1'
  ; instrument = 'atsr'
  
  ; platform = 'ers2'
  ; instrument = 'atsr2'
 
  ; platform = 'fengyun-4a'
  ; instrument = 'agri'
  
  ; platform = 'goes-16'
  ; instrument = 'abi'
  	
  ; platform = 'himawari-8'
  ; instrument = 'ahi'   

  ; platform = 'himawari-9'
  ; instrument = 'ahi' 
      
  ; platform = 'meteosat-8'
  ; instrument = 'seviri'

  ; platform = 'meteosat-9'
  ; instrument = 'seviri'
   
  ; platform = 'meteosat-10'
  ; instrument = 'seviri'
  
  ; platform = 'meteosat-11'
  ; instrument = 'seviri'

  ; platform = 'meteosat-12'
  ; instrument = 'fci'  

  ; platform = 'noaa-5'
  ; instrument = 'avhrr'
  
  ; platform = 'noaa-6'
  ; instrument = 'avhrr'
  
  ; platform = 'noaa-19'
  ; instrument = 'avhrr'

  ; platform = 'noaa-20'
  ; instrument = 'viirs'
  
  ; platform = 'sentinel-3a'
  ; instrument = 'slstr' 
   
  ; platform = 'sentinel-3b'
  ; instrument = 'slstr'  
    
  ; platform = 'terra'
  ; instrument = 'modis'    

   ; forward_model = 'aerosol' 
   ; forward_model = 'cloud' 
   
  ; control flags specify (unusual) options used within the lut creation code  
  ;  0 = Band Centre value
  ;  1 = Band average value at the resolution of the spectral response function
  ;  2 = Band average value at optimised resolution of the spectral response function
  ;  4 = scat_only  
  ;  8 = reuse_scat
  ;  16 = rayleigh
  ;  32 = gas
  ; to activate sum the control numbers to get control. It is usually zero

;##############################################################  
; Creates runfile for an instrument and selection of microphysical models   
; parameters that must be set
  
  platform   = 'sentinel-3a'
  instrument = 'slstr' 
  
  Job_Sheet = [{control: 0 ,forward_model:'aerosol', microphysical_model:'aerosol_a70',   LUT_type: 'aerosol'}, $
               {control: 0, forward_model:'cloud',   microphysical_model:'water-ice_sph', LUT_type: 'ice_cloud'} ]
  Job_Sheet = [{control: 0, forward_model:'aerosol', microphysical_model:'aerosol_a70',   LUT_type: 'aerosol_test'}]
  
  Versions = '21'     ; this is the output LUT version - change when new paramaters are used for LUT generation (or code output is 'changed')
;#########################################

  file = platform+'_'+instrument+'_run'
  openw,lun,file,/get_lun
  
  printf,lun,'#!/bin/bash'
  printf,lun,'source setup_oraclut_env.sh'
  printf,lun,'echo'
  For J = 0,N_elements(Job_Sheet) -1  do begin
    printf,lun,'export ORAC_LUT_DRIVER_FILE="${ORAC_LUT_INPUT_ROOT_DIR}/driver/' + platform + '_' + instrument + '_' + Job_Sheet(J).forward_model + '.driver"'
    case Job_Sheet(J).forward_model of         
      'aerosol': cmd = 'idl -e "create_orac_lut_aerosol, ' + strtrim(string(Job_Sheet(J).control), 2) + ', ''' + Job_Sheet(J).microphysical_model + ''', ''' + Job_Sheet(J).LUT_type + ''', ''' + Versions + '''"'
      'cloud':   cmd = 'idl -e "create_orac_lut_cloud, '   + strtrim(string(Job_Sheet(J).control), 2) + ', ''' + Job_Sheet(J).microphysical_model + ''', ''' + Job_Sheet(J).LUT_type + ''', ''' + Versions + '''"'
     EndCase      
    printf, lun, cmd
  EndFor
  printf,lun,'echo'
  close,lun

end