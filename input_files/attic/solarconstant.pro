
  rttovfilterdir='/network/home/g/grainger/project-oraclut/create_orac_lut/input_files/srf'
  rttovfilterdir='/home/g/grainger/project-oraclut/create_orac_lut/input_files/srf'
  spawn,'ls -1 '+rttovfilterdir, rttov_filters

  q=where(strmid(rttov_filters,0,2) eq 'rt',count) ; select the filterfiles
  solar_spectrum_filename = ['sun/Gueymard2018.sssi','sun/Thekaekara1973.sssi', 'sun/2000ASTM.sssi',   'sun/MCebKur.sssi',	  'sun/MChKur.sssi',	    'sun/MODWherli_WMO.sssi',		  'sun/MthKur.sssi', 'sun/MoldKur.sssi', 'sun/MNewKur.sssi']

  
  variation = .0328
 

 
  fst = [284,300]


  Form='(A20,7(1x,F8.2))'
  K=0
      
  load_solar_spectrum,solar_spectrum_filename[K],sssi
      
    For J = 0,count-1 do begin
      srf= rttovfilterdir+'/'+rttov_filters[q(j)]
      read_srfstr,srf,srfstr
      centre=srfstr.wvl_centre
      F0V=(simpson_integral(srfstr.wvl,interpol(sssi.val*10,sssi.wvl, srfstr.wvl) *srfstr.srf)/  simpson_integral(srfstr.wvl,srfstr.srf))/!pi
      wvn=1d4/sssi.wvl
      F0I =  (simpson_integral(srfstr.wvn,interpol(1d8*sssi.val/wvn^2,wvn, srfstr.wvn) *srfstr.srf)/  simpson_integral(srfstr.wvn,srfstr.srf))/!pi    
     print, rttov_filters[q(j)],'Channel:',centre,f0v,f0i
 ;     if ((3 lt srfstr.wvl_centre) and (srfstr.wvl_centre lt 4.6)) then print, rttov_filters[q(j)],'Channel:',centre,f0i

    Endfor
end