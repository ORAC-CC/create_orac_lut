; make run file to create multiple luts

; THIS IS THE MASTER LIST DO NOT DELETE
  MM = ['aerosol_a70',$
        'aerosol_a75',$
        'aerosol_a76',$
        'aerosol_a77',$
        'aerosol_a78',$  
        'aerosol_a79',$  
        'biomass_a',$
        'biomass_b',$
        'biomass_c',$
        'biomass_d',$
        'biomass_australia',$
        'biomass_indonesia',$
        'liquid-water_old',$
        'liquid-water_stg',$
        'liquid-water_240',$
        'liquid-water_253',$
        'liquid-water_263',$
        'liquid-water_273',$
        'sulphuric-acid_300',$  	
        'volcanic-ash_ctn',$
        'volcanic-ash_ey1',$
        'volcanic-ash_spr', $
        'water-ice_sph',$
        'water-ice_agg',$
        'water-ice_ghm',$
        'water-ice_src']

 ; 

;MM = ['sulphuric-acid_htha','volcanic-ash_htha'] ; Himawari only

;  MM = ['biomass_australia_minus','biomass_australia_plus']

; aerosol only	
MM = ['aerosol_a70',$
        'aerosol_a75',$
        'aerosol_a76',$
        'aerosol_a77',$
        'aerosol_a79',$  
        'biomass_a',$
        'biomass_b',$
        'biomass_c',$
        'biomass_d',$
        'biomass_australia',$
        'biomass_indonesia']


MM = ['aerosol_a70',$
        'aerosol_a75',$
        'aerosol_a76',$
        'aerosol_a77',$
        'aerosol_a78',$
        'aerosol_a79'] 
	
; cloud only		
  MM = ['liquid-water_old',$
        'liquid-water_stg',$
        'liquid-water_240',$
        'liquid-water_253',$
        'liquid-water_263',$
        'liquid-water_273',$
        'water-ice_sph',$
        'water-ice_agg',$
        'water-ice_ghm',$
        'water-ice_src']
MM = ['aerosol_a70',$
        'aerosol_a75',$
        'aerosol_a76',$
        'aerosol_a77',$
        'aerosol_a79'] 
; ice only		
  MM = ['water-ice_sph',$
        'water-ice_agg',$
        'water-ice_ghm',$
        'water-ice_src']
MM = ['liquid-water_old', 'water-ice_ghm']  
  MM = ['liquid-water_stg',$
        'liquid-water_ap02',$
        'liquid-water_ap04',$
        'liquid-water_ap06',$
        'liquid-water_ap08',$
        'liquid-water_ap10']

MM = ['biomass_australia_minus', 'biomass_australia','biomass_australia_plus']


MM = ['liquid-water_old','liquid-water_stg']  

  MM = ['liquid-water_stg',$
        'liquid-water_240',$
        'liquid-water_253',$
        'liquid-water_263',$
        'liquid-water_273',$
        'water-ice_sph',$
        'water-ice_agg',$
        'water-ice_ghm',$
        'water-ice_src']

 ; MM = ['aerosol_a70','aerosol_a76','aerosol_a79'] 
;MM = [ 'aerosol_a70','aerosol_a75','aerosol_a76','aerosol_a77','aerosol_a78','aerosol_a79'] 
 

 ; MM = ['liquid-water_stg','aerosol_a70','aerosol_a75','water-ice_agg']
; aerosol only	
;MM = [  'biomass_australia',    'biomass_indonesia']

   MM = ['liquid-water_old']  

   
  platform = 'ers1'
  instrument = 'atsr'
  
  platform = 'ers2'
  ;instrument = 'atsr2'
 
  platform = 'fengyun-4a'
  instrument = 'agri'

  
  platform = 'meteosat-8'
  instrument = 'seviri'

  platform = 'meteosat-9'
  instrument = 'seviri'
   
  platform = 'meteosat-10'
  instrument = 'seviri'
  
    platform = 'aqua'
  instrument = 'modis' 
  
   platform = 'terra'
  instrument = 'modis'  
	
  platform = 'himawari-8'
  instrument = 'ahi' 
  
  platform = 'noaa-5'
  instrument = 'avhrr'
  
  platform = 'noaa-6'
  instrument = 'avhrr'
  
  platform = 'noaa-19'
  instrument = 'avhrr'

  platform = 'meteosat-11'
  instrument = 'seviri'

  platform = 'envisat'
  instrument = 'aatsr'   
   
  platform = 'earthcare'
  instrument = 'msi'

  
  platform = 'sentinel-3a'
  instrument = 'slstr'  
  
  platform = 'sentinel-3b'
  instrument = 'slstr'  
   
  platform = 'goes-16'
  instrument = 'abi'
  
  platform = 'himawari-9'
  instrument = 'ahi' 
  
    
   platform = 'meteosat-12'
  instrument = 'fci'  

    platform = 'noaa-20'
  instrument = 'viirs'
  
    platform = 'meteosat-10'
  instrument = 'seviri'  

  
  
  
  Test = 1
  Versions = '21'
  
  file = platform+'_'+instrument+'_run'
  openw,lun,file,/get_lun
  
  printf,lun,'#!/bin/bash'
  printf,lun,'source initpath.bash'
; normally use gfortran but here have used intel
  printf,lun,'module load intel-compilers/2022'


  If (Test) Then printf,lun,'#*** TEST RUN ***'
  
  Old_Material=''
  
  For srf_quad= 2,2 do begin  
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
      printf,lun,'echo "making ..." $CREATE_ORAC_LUT_DRIVER'
      
      
      mminstruction = ",mmfile='"+mm(m)+".mm'"
      
      case material of
        'aerosol'       : lutinstruction = ",lutfile='aerosol"
        'biomass'       : lutinstruction = ",lutfile='biomass-cloud"
;        'liquid-water'  : lutinstruction = ",lutfile='liquid-water-cloud-MODIS" 
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
       'aerosol':  printf,lun,'idl -e "create_orac_aerosol_lut_wrapper,srf_quad='+string(srf_quad,format='(I1)')+mminstruction+lutinstruction+",tmatrix_path='/network/group/aopp/eodg/shared/dubovik_tmatrix/',atmospheres=2,gas=1,version="+versions+'"'
      'cloud': begin
         printf,lun,'idl -e "create_orac_cloud_lut_wrapper,srf_quad='+string(srf_quad,format='(I1)')+mminstruction+lutinstruction+",tmatrix_path='/network/group/aopp/eodg/shared/dubovik_tmatrix/',version="+versions+'"'
 ;        if (material eq 'liquid-water' or mm(m) eq 'volcanic-ash_htha') then printf,lun,'idl -e "create_orac_cloud_lut_wrapper,srf_quad='+string(srf_quad,format='(I1)')+mminstruction+lutinstruction+ ",no_rayleigh=1,reuse_scat=1,tmatrix_path='/network/group/aopp/eodg/shared/dubovik_tmatrix/',version="+versions+'"'
       end
       endcase	   
       printf,lun
    EndFor
  EndFor 
  close,lun

end