pro load_srfstrarr, inststr,  solar_spectrum_filename, srfstrarr, QM, nwvl_max
   
;  For each requested channel determine integrated quantities as well as set the quadrature for the LUT calculation   

;  If srf_quad is set it specifies the number of spectral integration points to use in LUT calculation.  
;  if value is set to 1 then use weighted centre wavelength and monochromatic rt. 
;  if the value is not set use the raw resolution of srf values 

;  nwvl_max:   maximum number of srf integration points
;  srf.wvl: 2-d array (nwvl_max, inststr.number_of_nadir_channels) of the wavelengths for each integration
;             point for each channel
;  srf.val:   2-d array (nwvl_max, inststr.number_of_nadir_channels) of the spectral response function values
;             for each integration point for each channel



;  work out the maximum number of spectral points that wil be used
  nwvl_max = 1
  FOR i = 0, inststr.number_of_nadir_channels - 1 DO BEGIN
    read_srfstr, 'input_files/srf/'+ inststr.srf_file[i], srfstr

    Case QM OF
      0: nwvl  = srfstr.nwvl
      1: nwvl = 1 
      2: BEGIN
          fe = 0.001  ; ie integrals must agree to better than .1 %
          segment,srfstr.wvl,srfstr.srf,x,y,fe,minn=12
          nwvl  = N_ELEMENTS(x)      
         END
    EndCase
    nwvl_max = nwvl_max > nwvl
  ENDFOR 
  
  load_solar_spectrum,solar_spectrum_filename,sssi
   
  srfstrarr = Replicate({nwvl: 0L , $
                   wvl_centre: 0.0, $
                   wvn_centre: 0.0, $
                           F0: 0.0, $     
                           B1: 0.0, $
                           B2: 0.0, $
                           T1: 0.0, $
                           T2: 0.0, $
                          wvl: FltARR(nwvl_max) , $
                          val: FltARR(nwvl_max)}, inststr.number_of_nadir_channels)
                           
  FOR i = 0,inststr.number_of_nadir_channels-1 DO BEGIN

    read_srfstr, 'input_files/srf/'+ inststr.srf_file[i], srfstr
   

    srfstrarr[i].wvl_centre  = srfstr.wvl_centre
    srfstrarr[i].wvn_centre  = srfstr.wvn_centre

;   Calcuate solar constant for channel (note different units for SW and thermal channels)  
    if (srfstr.wvl_centre Gt 3) then begin
      wvn=1d4/sssi.wvl
      srfstrarr[i].F0 =  (simpson_integral(srfstr.wvn,interpol(1d8*sssi.val/wvn^2,wvn, srfstr.wvn) *srfstr.srf)/  simpson_integral(srfstr.wvn,srfstr.srf))/!pi
    endif else $
      srfstrarr[i].F0 =  (simpson_integral(srfstr.wvl,interpol(sssi.val*10,sssi.wvl, srfstr.wvl) *srfstr.srf)/  simpson_integral(srfstr.wvl,srfstr.srf))/!pi
      
    bbconstants,srfstr.filename,b1,b2,t1,t2	
    srfstrarr[i].B1 = b1
    srfstrarr[i].B2 = b2
    srfstrarr[i].T1 = t1
    srfstrarr[i].T2 = t2

    Case QM OF
      0:BEGIN
         srfstrarr[i].nwvl  = srfstr.nwvl
         srfstrarr[i].wvl[0:srfstr.nwvl-1] = srfstr.wvl
         srfstrarr[i].val[0:srfstr.nwvl-1] = srfstr.srf
        END
      1:BEGIN
         srfstrarr[i].nwvl = 1 
         srfstrarr[i].wvl[0] = srfstrarr[i].wvl_centre
         srfstrarr[i].val[0] = 1 
        END
      2:BEGIN
         fe = 0.001  ; ie integrals must agree to better than .1 %
         segment,srfstr.wvl,srfstr.srf,x,y,fe,minn=12
         srfstrarr[i].nwvl  = N_ELEMENTS(x)
         srfstrarr[i].wvl[0:srfstrarr[i].nwvl-1] = x
         srfstrarr[i].val[0:srfstrarr[i].nwvl-1] = y   
      END
    EndCase
  ENDFOR
end