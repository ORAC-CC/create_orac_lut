 pro xx, x, a, f, pder

   b1=a(0)
   b2=a(1)  
   t1=a(2)   
   t2=a(3)
   
   teff=x*t2+t1
   c= exp(b2/teff)
   f = b1/(c-1)

 end
 
 function planck_wvn, Wavenumber, Temperature
; PURPOSE: Evaluate  planck function 
; Inputs:  WAVEnumber   wavenumbers(s) in cm-1 where Planck function is to be evaluated
;	         Temperature at specified wavelength points
; Outputs: Radiance in units mW m-2 str-1 (cm-1)-1

 c1=1.191043d-5
 c2=1.4387769d0

 return, c1*Wavenumber^3/(exp(c2*Wavenumber/Temperature)-1)
 
 end

pro load_old,inststr,channel,F0,F1,B1,B2,T1,T2
; herritage values from sad files


  case inststr.instrument of
    'aatsr': begin
      case channel of
        1: begin 
             F0 = 18.25 & F1 = 0.598 & B1 = 0 & B2 = 0 & T1 = 0 & T2 = 0 
           end
        2: begin 
             F0 = 18.25 & F1 = 0.598 & B1 = 0 & B2 = 0 & T1 = 0 & T2 = 0
           end
        3: begin 
             F0 = 18.25 & F1 = 0.598 & B1 = 0 & B2 = 0 & T1 = 0 & T2 = 0
           end
        4: begin 
             F0 = 18.25 & F1 = 0.598 & B1 = 0 & B2 = 0 & T1 = 0 & T2 = 0
           end
        5: begin 
             F0 =  5.39 & F1 = 0.177 & B1 = 228023.02   & B2 = 3849.0896 & T1 = 1.7269503   & T2 = 0.99634349
           end
        6: begin 
             F0 =  0    & F1 = 0     & B1 =   9323.4507 & B2 = 1326.0422 & T1 = 0.1033266   & T2 = 0.999253
           end
        7: begin 
             F0 =  0    & F1 = 0     & B1 =   6821.7295 & B2 = 1194.8956 & T1 = 0.055009452 & T2 = 0.99928154 
           end
        else: stop,'load_old in load_srfstrarr: unknown channel'
      endcase
      end
    else: stop,'load_old in load_srfstrarr: unknown instrument'
  endcase
end
 
pro load_srfstrarr, inststr,  solar_spectrum_filename, srfstrarr, nwvl_max, N_srf_points = n_srf_points
   
;  For each requested channel determine integrated quantities as well as set the quadrature for the LUT calculation   

;  If n_srf_points is set it specifies the number of spectral integration points to use in LUT calculation.  
;  if value is set to 1 then use weighted centre wavelength and monochromatic rt. 
;  if the value is not set use the raw resolution of srf values 

;  nwvl_max:   maximum number of srf integration points
;  srf.wvl: 2-d array (nwvl_max, inststr.number_of_nadir_channels) of the wavelengths for each integration
;             point for each channel
;  srf.val:   2-d array (nwvl_max, inststr.number_of_nadir_channels) of the spectral response function values
;             for each integration point for each channel


   IF KEYWORD_SET(n_srf_points) then $
     nwvl_max = n_srf_points $
   ELSE begin
;    work out the maximum number of spectral points  needed
     nwvl_max = 1
     FOR i = 0, inststr.number_of_nadir_channels - 1 DO BEGIN
        read_srfstr, 'input_files/srf/'+ inststr.srf_file[i], srfstr
        nwvl_max = nwvl_max > srfstr.nwvl
     ENDFOR
   ENDELSE
   
  load_solar_spectrum,solar_spectrum_filename,sssi
   
   srfstrarr = Replicate({nwvl: 0L , $
                    wvl_centre: 0.0, $
                    wvn_centre: 0.0, $
                            F0: 0.0, $
                            F1: 0.0, $
                            B1: 0.0, $
                            B2: 0.0, $
                            T1: 0.0, $
                            T2: 0.0, $
                           wvl: FltARR(nwvl_max) , $
                           val: FltARR(nwvl_max)}, inststr.number_of_nadir_channels)
                           

   hacky=[0.560000, 0.66000, 0.862000, 1.59400,3.74200,10.8570,12.0510,0.560000, 0.66000, 0.862000, 1.59400,3.74200,10.8570,12.0510]
   variation = .0328
   T0 = 200d0
   dt = .1d0
   tn = 1400
   Temperature=T0+dt*indgen(Tn)
   r=fltarr(tn)  
   weights=replicate(1.0,tn)
   y=[0,0,0,0]
   cc=0

   FOR i = 0,inststr.number_of_nadir_channels-1 DO BEGIN
     read_srfstr, 'input_files/srf/'+ inststr.srf_file[i], srfstr
     srfstrarr[i].wvl_centre  = srfstr.wvl_centre
     srfstrarr[i].wvn_centre  = srfstr.wvn_centre

;    Calcuate solar constant for channel (only matters for joint ir/solar channels)  
     wvn=1d4/sssi.wvl
     srfstrarr[i].F0 =  (int_tabulated(srfstr.wvn,interpol(1d8*sssi.val/wvn^2,wvn, srfstr.wvn) *srfstr.srf)/  int_tabulated(srfstr.wvn,srfstr.srf))/!pi
     srfstrarr[i].F1 =  srfstrarr[i].F0  * variation 

;    calculate planck function fit for channel    (only matters for IR channels ) 
;      a priori constants
      case 1 of
        (srfstr.wvl_centre gt 3) and (srfstr.wvl_centre lt 4) :a =[209052.19 , 3804.5361, 0.76911219, 1.0115517]
        (srfstr.wvl_centre gt 4) and (srfstr.wvl_centre lt 6) :a =[155731.42 , 3618.3402, 0.94714021, 1.0773587]
        (srfstr.wvl_centre gt 6) and (srfstr.wvl_centre lt 7) :a =[ 42890.197, 2962.5137, 2.1166978 , 1.3437637] 
        (srfstr.wvl_centre gt 7) and (srfstr.wvl_centre lt 8) :a = [ 30325.608, 2801.9281, 0.35307431, 1.4256249] 
        (srfstr.wvl_centre gt 8) and (srfstr.wvl_centre lt 9) :a = [ 18928.043, 2600.6576, 0.39913111, 1.5483174] 
        (srfstr.wvl_centre gt 9) and (srfstr.wvl_centre lt 10) :a =[ 13236.544, 2454.7584, 0.18254054, 1.6467873] 
        (srfstr.wvl_centre gt 10) and (srfstr.wvl_centre lt 11) :a =[ 9694.3846, 2335.1605, 0.86201090, 1.7366491]
        (srfstr.wvl_centre gt 11) and (srfstr.wvl_centre lt 12) :a = [ 7659.0934, 2246.4495, 0.81786035, 1.8095527] 
        (srfstr.wvl_centre gt 12) and (srfstr.wvl_centre lt 13) :a =[ 6661.0251, 2196.6090, 0.56743168, 1.8519080] 
        (srfstr.wvl_centre gt 13) and (srfstr.wvl_centre lt 15) :a =[ 4798.6471, 2082.6820, 0.28632467, 1.9606028]
        else : a = [0,0,0,0]
      endcase
      if (srfstr.wvl_centre gt 3) and (srfstr.wvl_centre lt 15) then begin
        for j = 0,tn-1 do begin
          r(j) = int_tabulated(srfstr.wvn, planck_wvn(srfstr.wvn, Temperature(j))*srfstr.srf)/  int_tabulated(srfstr.wvn,srfstr.srf) ; mappring from temperature to radiance
        endfor
;       find fit for channel
        yfit = curvefit(temperature,r,weights,a,function_name='xx',/noderivative,/double,itmaX=1000,tol=1d-3)
      endif
      srfstrarr[i].B1 = A(0)
      srfstrarr[i].B2 = A(1)
      srfstrarr[i].T1 = A(2)
      srfstrarr[i].T2 = A(3)
      
      print,' Calculated:', srfstrarr[i].F0,srfstrarr[i].F1,srfstrarr[i].B1,srfstrarr[i].B2,srfstrarr[i].T1,srfstrarr[i].T2
       If  inststr.instrument eq 'aatsr' then begin
         srfstrarr[i].wvl_centre =hacky[inststr.channelid[i]-1]
         srfstrarr[i].wvn_centre = 1e4/srfstrarr[i].wvl_centre
         load_old,inststr, i+1,F0,F1,B1,B2,T1,T2 
         srfstrarr[i].F0 = F0
         srfstrarr[i].F1 = F1  
         srfstrarr[i].B1 = B1
         srfstrarr[i].B2 = B2
         srfstrarr[i].T1 = T1
         srfstrarr[i].T2 = T2    
         print,'Using Hacky:', srfstrarr[i].F0,srfstrarr[i].F1,srfstrarr[i].B1,srfstrarr[i].B2,srfstrarr[i].T1,srfstrarr[i].T2
     ENDIF       

     IF KEYWORD_SET(n_srf_points) THEN BEGIN
       srfstrarr[i].nwvl = nwvl_max   
       IF (nwvl_max eq 1) THEN BEGIN
         srfstrarr[i].wvl[0] = srfstrarr[i].wvl_centre
         srfstrarr[i].val[0] = 1     
       ENDIF ELSE BEGIN
         chanwl1 = srfstr.wvl[0]
         chanwl2 = srfstr.wvl[-1]
         srfstrarr[i].wvl[*] = chanwl1 + FINDGEN(nwvl_max) * (chanwl2 - chanwl1) / (nwvl_max - 1)
         srfstrarr[i].val[*] = INTERPOL(srfstr.srf, srfstr.wvl, srfstrarr[i].wvl[*])
       ENDELSE
     ENDIF ELSE BEGIN
       srfstrarr[i].nwvl  = srfstr.nwvl
       srfstrarr[i].wvl[0:srfstr.nwvl-1] = srfstr.wvl
       srfstrarr[i].val[0:srfstr.nwvl-1] = srfstr.srf
     ENDELSE
   ENDFOR
  return
end