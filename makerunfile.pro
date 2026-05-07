; make run file to create multiple luts

; THIS IS THE MASTER LIST DO NOT DELETE
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
 
  ; platform = 'ers1'
  ; instrument = 'atsr'
  
  ; platform = 'ers2'
  ; instrument = 'atsr2'
 
  ; platform = 'fengyun-4a'
  ; instrument = 'agri'
  
  ; platform = 'meteosat-8'
  ; instrument = 'seviri'

  ; platform = 'meteosat-9'
  ; instrument = 'seviri'
   
  ; platform = 'meteosat-10'
  ; instrument = 'seviri'
  
  ; platform = 'aqua'
  ; instrument = 'modis' 
  
  ; platform = 'terra'
  ; instrument = 'modis'  
	
  ; platform = 'himawari-8'
  ; instrument = 'ahi' 
  
  ; platform = 'noaa-5'
  ; instrument = 'avhrr'
  
  ; platform = 'noaa-6'
  ; instrument = 'avhrr'
  
  ; platform = 'noaa-19'
  ; instrument = 'avhrr'

  ; platform = 'meteosat-11'
  ; instrument = 'seviri'

  ; platform = 'envisat'
  ; instrument = 'aatsr'   
   
  ; platform = 'earthcare'
  ; instrument = 'msi'

  ; platform = 'sentinel-3b'
  ; instrument = 'slstr'  
   
  ; platform = 'goes-16'
  ; instrument = 'abi'
  
  ; platform = 'himawari-9'
  ; instrument = 'ahi' 
    
  ; platform = 'meteosat-12'
  ; instrument = 'fci'  
  
  ; platform = 'meteosat-10'
  ; instrument = 'seviri'  

  ; platform = 'noaa-20'
  ; instrument = 'viirs'
  
MM = ['aerosol_a70', 'water-ice_sph']

  platform = 'sentinel-3a'
  instrument = 'slstr'    
  
  Test = 1
  Versions = '21'
  
  file = platform+'_'+instrument+'_run'
  openw,lun,file,/get_lun
  
  printf,lun,'#!/bin/bash'

  If (Test) Then Begin
   printf,lun,'#*** WARN USER THAT THIS IS A TEST RUN ***'
   printf,lun,'echo'
   printf,lun,'echo'
   printf,lun,'echo'
   printf,lun,'echo "****************************************"'
   printf,lun,'echo "*** TEST TEST TEST TEST TEST TEST ***"'
   printf,lun,'echo "*** THIS IS A TEST LUT GENERATION  ***"'
   printf,lun,'echo "****************************************"'
   printf,lun,'echo'
   printf,lun,'echo'
   printf,lun,'echo'
  EndIf
  printf,lun,'source setup_oraclut_env.sh'
  printf,lun,'echo'
  printf,lun,'echo'
  Old_Material=''
  
  For srf_quad= 1,1 do begin  
    Case srf_quad of
    1:printf,lun,'# monochromatic calculation'
    2:printf,lun,'# band  calculation'
    EndCase
    
    For M = 0,N_elements(MM) -1  do begin
      material = (strsplit(mm(m),'_', /extract)) [0]
      If (Material Ne OLd_Material) then begin
        printf,lun,'# '+material + ' section'
        old_material = material
      endif

      case material of
        'aerosol'       : FM = 'aerosol'
        'biomass'       : FM = 'cloud'
        'liquid-water'  : FM = 'cloud'
        'sulphuric-acid': FM = 'cloud'
        'volcanic-ash'  : FM = 'cloud'
        'water-ice'     : FM = 'cloud' 
      EndCase
  
	  Case FM of
	  'aerosol': printf,lun,'export CREATE_ORAC_AEROSOL_LUT_DRIVER="input_files/driver/'+platform+'_'+instrument+'_'+FM+'.driver"'
	  'cloud': printf,lun,'export CREATE_ORAC_CLOUD_LUT_DRIVER="input_files/driver/'+platform+'_'+instrument+'_'+FM+'.driver"'
	  EndCase 	  
	  Case FM of
	  'aerosol': printf,lun,'echo "using ..." $CREATE_ORAC_AEROSOL_LUT_DRIVER'
	  'cloud':   printf,lun,'echo "using ..." $CREATE_ORAC_CLOUD_LUT_DRIVER'
	  EndCase 
      
           
      mminstruction = ",mmfile='"+mm(m)+".mm'"
      
      case material of
        'aerosol'       : lutinstruction = ",lutfile='aerosol"
        'biomass'       : lutinstruction = ",lutfile='biomass-cloud"
        'liquid-water'  : lutinstruction = ",lutfile='liquid-water-cloud" 
        'sulphuric-acid': lutinstruction = ",lutfile='sulphuric-acid-cloud"
        'volcanic-ash'  : lutinstruction = ",lutfile='ash-plume"
        'water-ice'     : lutinstruction = ",lutfile='ice-cloud"          
      EndCase    
      
      If (Test) Then $
        lutinstruction = lutinstruction+"_test.lut'" $
       Else $
        lutinstruction = lutinstruction+".lut'"
		
	  Case FM of        
       'aerosol':  printf,lun,'idl -e "create_orac_aerosol_lut_wrapper,srf_quad='+string(srf_quad,format='(I1)')+mminstruction+lutinstruction+",tmatrix_path='/network/aopp/matin/eodg/shared/dubovik_tmatrix/',atmospheres=2,gas=1,version="+versions+'"'
      'cloud': begin
         printf,lun,'idl -e "create_orac_cloud_lut_wrapper,srf_quad='+string(srf_quad,format='(I1)')+mminstruction+lutinstruction+",tmatrix_path='/network/group/aopp/eodg/shared/dubovik_tmatrix/',version="+versions+'"'
        end
       endcase	   
       printf,lun
    EndFor
  EndFor 
  close,lun

end