pro load_solar_spectrum,solar_spectrum_filename,sssi
  loadxy, solar_spectrum_filename,0,1,wvl,val
  
  sssi ={Name: File_Basename(solar_spectrum_filename), $ 
          wvl: wvl, $
          val: val}
  return
end