; programme to convert V1 Luts to a V2 single file

  platform='seviri'

  instrument='MSG4'
  table ='WAT'
  lutfile ='liquid-water-cloud.lut'
  Substance = 'liquid-water'
  MonoOrBand ='m'
  Atmospheric_Model_Code='01'
  shortname ='old'
  versions='00'
  

  indir = '/network/aopp/apres/ORAC_LUTS/'+platform+'/'+instrument+'/'+table+'/'
  
  in_path= 'input_files'
  case instrument of
  'MSG1':instrument_number = '8'
  'MSG2':instrument_number = '9'
  'MSG3':instrument_number = '10'
  'MSG4':instrument_number = '11'
  endcase
  
  sadbasefilename = indir+strupcase(platform)+'-'+instrument
  sadtablefilename = sadbasefilename+'_'+table

; read the Instrument definition file
; NOTE this is a dummy read - all relevant parameters are overwritten
  instfile ='meteosat-'+instrument_number+'_'+platform+'_v1.inst'
  instdirfile = in_path+'/inst/'+instfile
  ok = file_test(instdirfile, /read)
  IF not ok then message, 'instdirfile not readable: ' + instdirfile
; create structure (values will be overwriten)
; **** Read the Instrument parameters file
  load_inststr, instdirfile, inststr
  
; **** Read the Instrument spectral response file assuming QM = 1 & nwvl_max =10
; NOTE this is a dummy read - all relevant parameters are overwritten
  QM = 1
  Nwvl_max = 10
  solar_spectrum_filename = in_path+'/sun/Thekaekara1973.sssi'
  ok = file_test(solar_spectrum_filename , /read)
  IF not ok then message, 'solar spectrum file not readable: ' + solar_spectrum_filename
  load_srfstrarr,  inststr, solar_spectrum_filename, srfstrarr, QM, nwvl_max
  
;  **** Create LUT structure
 
;  Open the RBD file and get the LUT dimensions   
   rbd_filename = sadtablefilename+'_RBD_Ch1.sad'
   read_lut5d, rbd_filename, lut1, lut2
  
   LUTstr = { Opd_N      : lut1.opticaldepths,       Efr_N: lut1.effectiveradii , Soz_N       : lut1.sunzeniths, Saz_N      : lut1.satzeniths, Raa_N      : lut1.relativeazimuths, $
              OPD_Spacing: 'uneven_logarithmic'   , Efr_Spacing: 'uneven_linear'     , Soz_Spacing : 'uneven_linear', Saz_Spacing: 'uneven_linear', Raa_Spacing: 'uneven_linear',  $
              opd        : lut1.opticaldepth ,         EFR: lut1.effectiveradius, SOz         : lut1.sunzenith , Saz        : lut1.satzenith , Raa        : lut1.relativeazimuth    }
   
;  define & create Vavg LUT   
   Vavg = Fltarr(lutstr.efr_n)
   V1filename = sadtablefilename+'_Vavg.sad'
   read_lut1d, V1filename, LUTVavg  
   Vavg = LUTVavg.LTV 
;  define property LUTs    
   BextOUT = Fltarr(inststr.Number_of_nadir_Channels,lutstr.opd_n, lutstr.efr_n)        ; opd a redundant dimension remove before write
   BextRatOUT = Fltarr(inststr.Number_of_nadir_Channels,lutstr.opd_n, lutstr.efr_n)     ; opd a redundant dimension remove before write
   SSAOUT = Fltarr(inststr.Number_of_nadir_Channels,lutstr.opd_n, lutstr.efr_n)         ; opd a redundant dimension remove before write   
   GOUT =Fltarr(inststr.Number_of_nadir_Channels,lutstr.opd_n, lutstr.efr_n)            ; opd a redundant dimension remove before write   
   
;  **** Define the LUT table output variables themselves
   RFD  = FLTARR(inststr.Number_of_Channels, lutstr.efr_n, lutstr.opd_n)
   TFD  = FLTARR(inststr.Number_of_Channels, lutstr.efr_n, lutstr.opd_n)
   RD   = FLTARR(inststr.Number_of_Channels, lutstr.efr_n, lutstr.opd_n, lutstr.saz_n)
   TD   = FLTARR(inststr.Number_of_Channels, lutstr.efr_n, lutstr.opd_n, lutstr.saz_n)
   TB   = FLTARR(inststr.Number_of_Channels, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n)
   RFBD = FLTARR(inststr.Number_of_Channels, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n)
   TFBD = FLTARR(inststr.Number_of_Channels, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n)
   RBD  = FLTARR(inststr.Number_of_Channels, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n, lutstr.saz_n, lutstr.raa_n)
   Em   = FLTARR(inststr.Number_of_Channels, lutstr.efr_n, lutstr.opd_n, lutstr.saz_n)  
   

; Loop through channels populating the optical property and rt tables
  For I = 0,inststr.Number_of_Channels-1 do begin
    channel = I + 1
	
    suffix = '_Ch'+string(channel,format='(I0)')+'.sad' 
    sadfilename =sadbasefilename+suffix
    read_sad,sadfilename,sadvalues
	srfstrarr[i].WVN_CENTRE =sadvalues.wvn
    srfstrarr[i].WVL_CENTRE = sadvalues.descriptor
    srfstrarr[i].WVL[0]= sadvalues.descriptor
	

    V1filename = sadtablefilename+'_Bext'+suffix
    read_lut2d, V1filename, LUTbext, /Wavelength
	BextOUT(I,*,*) = LUTbext.LTV
   
    V1filename = sadtablefilename+'_BextRat'+suffix
    read_lut2d, V1filename, LUTbextrat 
	BextRatOUT(I,*,*) = LUTBextrat.LTV
      
    V1filename = sadtablefilename+'_w'+suffix
    read_lut2d, V1filename, LUTw, /Wavelength
    SSAOUT(I,*,*) = LUTw.LTV
	
	V1filename = sadtablefilename+'_g'+suffix
    read_lut2d, V1filename, LUTg, /Wavelength   
	GOUT(I,*,*) = LUTg.LTV 	
	
	V1filename = sadtablefilename+'_RD'+suffix
    read_lut3db, V1filename, LUTRD, LUTRfd 
    For ii = 0,lutstr.efr_n-1 do $
      For jj = 0, lutstr.opd_n-1 do $
        For kk = 0, lutstr.saz_n-1 do $
          RD(I,ii,jj,kk) = LUTRD.LTV(jj,kk,ii)
	RFD(I,*,*) = Transpose(LUTRfd.LTV)
	
	V1filename = sadtablefilename+'_TD'+suffix
	read_lut3db, V1filename, LUTTD, LUTTfd    
    For ii = 0,lutstr.efr_n-1 do $
      For jj = 0, lutstr.opd_n-1 do $
        For kk = 0, lutstr.saz_n-1 do $
          TD(I,ii,jj,kk) = LUTTD.LTV(jj,kk,ii)
	TFD(I,*,*) = Transpose(LUTTfd.LTV)

	If (Sadvalues.solarflag Eq 1 ) Then begin
      srfstrarr[i].F0 =sadvalues.solarconstant[0]	
    
	  V1filename = sadtablefilename+'_RBD'+suffix	
      read_lut5d, V1filename, LUTRBD, LUTRfBD  
      For ii = 0, lutstr.efr_n-1 do $
        For jj = 0, lutstr.opd_n-1 do $
          For kk = 0, lutstr.soz_n-1 do $
            For ll = 0, lutstr.saz_n-1 do $
              For mm = 0, lutstr.raa_n-1 do $
                RBD(I,ii,jj,kk,ll,mm) = LUTRbD.LTV(jj,ll,kk,lutstr.raa_n-1 -mm,ii)	 ; reverse azimuth
	   
      For ii = 0, lutstr.efr_n-1 do $
        For jj = 0, lutstr.opd_n-1 do $
          For kk = 0, lutstr.soz_n-1 do $
	        RFBD(I,ii,jj,kk) =LUTRfBD.LTV(jj,kk,ii)
	   
	  V1filename = sadtablefilename+'_TBD'+suffix
	  read_lut5d, V1filename, LUTTBD, LUTTfBD		 
	  For ii = 0, lutstr.efr_n-1 do $
        For jj = 0, lutstr.opd_n-1 do $
          For kk = 0, lutstr.soz_n-1 do $
	        TFBD(I,ii,jj,kk) =LUTTfBD.LTV(jj,kk,ii)
	  
	  V1filename = sadtablefilename+'_TB'+suffix
	  read_lut3da, V1filename, LUTTB
	  TB(I,*,*,*) = LUTTB.LTV
      For ii = 0,lutstr.efr_n-1 do $
        For jj = 0, lutstr.opd_n-1 do $
          For kk = 0, lutstr.saz_n-1 do $
            TB(I,ii,jj,kk) = LUTTB.LTV(jj,kk,ii)  
	endif	
	If (Sadvalues.thermalflag Eq 1 ) Then begin
	  srfstrarr[i].B1=sadvalues.B1
      srfstrarr[i].B2=sadvalues.B2
      srfstrarr[i].T1=sadvalues.T1
      srfstrarr[i].T2=sadvalues.T2
	  
	  inststr.nedt(i)=sadvalues.nebt
	  
	
	  V1filename = sadtablefilename+'_EM'+suffix
	  read_lut3da, V1filename, LUTEM	  
	  For ii = 0,lutstr.efr_n-1 do $
        For jj = 0, lutstr.opd_n-1 do $
          For kk = 0, lutstr.saz_n-1 do $
            EM(I,ii,jj,kk) = LUTEM.LTV(jj,kk,ii)  
	endif	
  endfor
 
; remove redundant dimensions
  BextOUT    = reform(BextOUT[*,0,*])
  BextRatOUT = reform(BextRatOUT[*,0,*])
  SSAOUT     = reform(SSAOUT[*,0,*])
  GOUT       = reform(GOUT[*,0,*])
  
; write out the table to current directory
  V2_LUT_Filename = strlowcase(inststr.Platform)+'_'+strlowcase(inststr.Instrument)+'_'+MonoOrBand+'_'+strlowcase(Substance)+'_a'+Atmospheric_Model_Code+'_p'+strlowcase(shortname)+'_v'+versions+'.nc'
  write_v2_lut, V2_LUT_Filename, lutstr, inststr, srfstrarr, Vavg, BextOUT, BextRatOUT, SSAOUT, GOUT, TD, TfD, RD, RfD, RBD = RBD, RfBD = RfBD, TfBD = TfBd, TB = TB, EM = EM


end