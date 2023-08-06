
  rttovfilterdir='~/project-oraclut/create_orac_lut/input_files/srf'
  spawn,'ls -1 '+rttovfilterdir, rttov_filters

  q=where(strmid(rttov_filters,0,2) eq 'rt',count) ; select the filterfiles
solar_spectrum_filename = 'sun/Thekaekara1973.sssi'
;solar_spectrum_filename = 'sun/2000ASTM.sssi'
  load_solar_spectrum,solar_spectrum_filename,sssi
  
  variation = .0328
 
  for j =0,count-1 do begin 

    srf= rttovfilterdir+'/'+rttov_filters[q[j]]
    read_srfstr,srf,srfstr
    
    If srfstr.wvl_centre gt 3 and srfstr.wvl_centre lt 5 then begin

;    Calcuate solar constant for channel (only matters for joint ir/solar channels)  
     wvn=1d4/sssi.wvl
     F0 =  (int_tabulated(srfstr.wvn,interpol(1d8*sssi.val/wvn^2,wvn, srfstr.wvn) *srfstr.srf)/  int_tabulated(srfstr.wvn,srfstr.srf))/!pi
     F1 =  F0  * variation 
     print, rttov_filters[q[j]],srfstr.wvl_centre,f0,f1
     endif
  Endfor
end