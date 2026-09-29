pro interpret_rttov_filter, rttov_filter, platform, instrument, channelid
; convert the rttov spectral response function name into the instrument and platform names used in ORAC
; ie OSCAR convention for platform and instrument naming.  Channels in resolution then ascending order of wavelength.
; atsr, atsr-2 and aatsr RTTOV filter fiiles break the convention and go IR->Vis.  The channel number is reversed here.

  words = strsplit(rttov_filter,'_',/extract)
  RTTOV_platform = words[1]
  RTTOV_number = words[2]
  RTTOV_instrument=words[3]
  
; Generally we can use the rttov filename to create a new name  
  platform = RTTOV_platform
  instrument = RTTOV_instrument
  channelid = strmid(words[5],2,2)

; but there are exceptions covered below  
  Case RTTOV_platform of
  'ers': begin 
           channelvl=long(channelid) ; needed to reverse channel designation to OSCAR/ORAC standard
           case RTTOV_number of
             '1': begin 
                    platform = 'ers-1'
                    instrument= 'atsr'
                    If (channelid Lt 5) then $
                      channelid = string(5-channelvl,Format='(I2.2)') $
                    Else begin     
                      instrument= 'atsr-dt'    
                      channelid = '04'                      
                    endelse                    
                  end
             '2': begin
                    platform = 'ers-2'
                    instrument= 'atsr-2'
                    channelid = string(8-channelvl,Format='(I2.2)')
                  end
             endcase
         end
  'envisat': begin
               channelvl=long(channelid) ; needed to reverse channel designation to OSCAR/ORAC standard
               case RTTOV_instrument of
                 'atsr': instrument= 'aatsr'
                 'atsr-shifted':  instrument= 'aatsr-shifted'
               endcase
               channelid = string(8-channelvl,Format='(I2.2)')
             end

  'eos': begin
           case RTTOV_number of
             '1': platform ='terra'
             '2': platform ='aqua'
           endcase
         end
  'fy4': platform = 'fengyun-4a'
  'goes': platform=RTTOV_platform+'-'+RTTOV_number
  'himawari': platform=RTTOV_platform+'-'+RTTOV_number
  'metop': begin
             case RTTOV_number of
               '1': platform = 'metop-b'  ; note unexpected order 1 is B and 2 is A
               '2': platform = 'metop-a'
               '3': platform = 'metop-c'
             endcase
           end
  'msg': begin
             case RTTOV_number of
               '1': platform = 'meteosat-8' 
               '2': platform = 'meteosat-9'
               '3': platform = 'meteosat-10'
               '4': platform = 'meteosat-11'
             endcase
           end
  'mtg': begin
             case RTTOV_number of
               '1': platform = 'meteosat-12' 
             endcase
           end
  'noaa': platform=RTTOV_platform+'-'+RTTOV_number
  'sentinel2': begin
                case RTTOV_number of
                  '1': platform = 'sentinel-2a'
                  '2': platform = 'sentinel-2b'
                endcase
              end
  'sentinel3': begin
                case RTTOV_number of
                  '1': platform = 'sentinel-3a'
                  '2': platform = 'sentinel-3b'
                endcase
              end
  else:
  endcase
end
;+
; pro create_orac_lut_calc_gas_opd
;
; This procedure uses the gas extinction profiles, produced by
; create_orac_lut_run_modtran.pro, and instrument channel filter
; response functions, to calculate the 'gas opd' profiles for the
; channels of a range of instruments. These profiles are then output
; in a form suitable for use in the ORAC LUT generation code.
;
; INPUT ARGUMENTS
; inst      The instrument name to produce data for (e.g. MODIS)
; plat      The platform name for the above instrument (e.g. Terra)
; filterdir The directory in which filter response function data files
;           can be found. These must have a particular format and
;           conform to a particular naming convention - see
;           appropriate reading routine for details
; KEYWORDS
; outdir    Allows an output directory to be specified. By default,
;           output is written to the current working dir.
; tapedir   Allows the directory in the MODTRAN output files can be
;           found. By default, they are expected to be in the current
;           working dir.
; HISTORY
; 2018/02/09 Gareth Thomas: Rewritten from the old calc_gas_opd, which
;           has been kicking about the Oxford EODG group for about 15
;           years.
; 2018/02/22 Gareth Thomas: Updated to support using new version
;           create_orac_lut_read_chfilter, which now expects RTTOV
;           format response function files.
;           Also improved writing header lines, and changed to a .gas
;           file suffix.
; 2020/08/17 RGG: Rewrite to get information from files rttov files rather than define here. In some case, modify plat and inst to conform to OSCAR naming. https://www.wmo-sat.info/oscar/satellites
;- 
pro create_orac_lut_calc_gas_opd, tapedir, filterdir, rttov_filter, outdir

  For Atmosphere = 0,6 do begin
    A = string(Atmosphere,Format='(I1)')
    restore, tapedir+'/modtran_height_A'+A+'.sav' ; restores atmospheres and height
    nlevels = n_elements(heights)
    tau = fltarr(nlevels)
  
    interpret_rttov_filter, rttov_filter, platform, instrument, channelid
  
    outfile = 'ModtranGasOpd_A'+A+'_'+platform+'_'+instrument+'_'+'ch'+channelid+'.gas'   
  
    print, rttov_filter+' ---> '+outfile
  
    outname = outdir+'/'+outfile
             
;   Read in the channel spectral response function    
    filename = filterdir+'/'+rttov_filter
    read_srfstr, filename, srfstr
    wvl_centre = string(srfstr.wvl_centre ,Format='(F6.3)')
 
;   Transmission files are stored in 1 cm-1 resolution
;   Calaulate the wavenumber range to use
    wn1 = min(srfstr.wvn) 
    nwn = ceil(max(srfstr.wvn)) - floor(min(srfstr.wvn)) + 1
  
    openw,outlun,outname,/get_lun
      printf,outlun,'# MODTRAN 3.5 v1.1 produced atmospheric gas optical depth (without'
      printf,outlun,'# scattering for '+instrument+' channel '+channelid+' ('+wvl_centre+' um).'
      printf,outlun,'*atmosphere'
      printf,outlun, Atmospheres
      printf,outlun,'*instrument'
      printf,outlun,instrument
      printf,outlun,'*channelid'
      printf,outlun,channelid
      printf,outlun,'*nlevels'
      printf,outlun, nlevels

      tauold=999.
      for j = 0,nlevels-2 do begin ; Top level has 0 absorption by definition
        modfile=tapedir+'/tape7_'+A+'_'+strtrim(fix(Heights[j]),2)
        read_modtran, modfile, wn1, nwn, wvn, trans
;       interpolate response onto trans (response smoother)
        result =  interpol(srfstr.srf,srfstr.wvn,wvn)
        tau[j] = (-int_tabulated(wvn,result*alog(trans),/double)/int_tabulated(wvn,result,/double)) > 0
        if tau[j] ge tauold then tau[j]=0
        tauold=tau[j]
      endfor
    
;     Define the output format
;     1234567890123
;     XX0.12345e-78
      format = '(f5.1,e13.5)'
      for j = nlevels-1,0,-1 do begin
        printf,outlun,float(heights[j]),tau[j], format=format
      endfor
    close,outlun
    free_lun,outlun
  EndFor
END
