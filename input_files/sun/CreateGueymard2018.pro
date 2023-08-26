; convert raw Gueymard 2018 spectrum into sssi file.

loadxy,'1-s2.0-S0038092X1830433X-mmc1.txt',0,1,x,y

openw,lun,'Gueymard2018.sssi',/get_lun
  printf,lun,"* superterrestrial solar spectral irradiance from C. Gueymard, 'Revised composite extraterrestrial spectrum based on recent solar irradiance observations', Solar Energy, 2018"
  printf,lun,"* wavelength in um irradiance in mW cm^{-2} um^{-1}"
  printf,lun,"* solar constant 136.11 mW cm^{-2}"
  for i = 0,n_elements(x)-1 do printf,lun, x(i)*1E-3,y(i)*100,format='(F9.4,1x,E10.4)'
 close,lun
 end