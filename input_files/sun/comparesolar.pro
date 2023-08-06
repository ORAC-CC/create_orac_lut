colour_ps,9

solar_spectrum_filename = ['Thekaekara1973.sssi','2000ASTM.sssi','Gueymard2018.sssi']

  for j =0,N_elements(solar_spectrum_filename)-1 do begin 
    load_solar_spectrum,solar_spectrum_filename(j),sssi
    if j eq 0 then plot_oo, sssi.wvl,sssi.val,xrange=[0.5,4.1] else $
    oplot, sssi.wvl,sssi.val,color=j+2
 
  Endfor
  
  solar_spectrum_filename = '2000ASTM.sssi'
      load_solar_spectrum,solar_spectrum_filename,sssi1  
  solar_spectrum_filename = 'Gueymard2018.sssi'
      load_solar_spectrum,solar_spectrum_filename,sssi2
      val = interpol(sssi2.val,sssi2.wvl,sssi1.wvl)
      
end