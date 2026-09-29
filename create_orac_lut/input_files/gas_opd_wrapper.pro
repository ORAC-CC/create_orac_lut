; programme to build gas transmission files 
; first run MODTRAN and create gas transmission at all wavelengths
; you need to run create_orac_lut_run_modtran (idl) in the modtran directory.
; It will create the spectral transmission files required in the output directory

  modtrantapedir='modtran3.5-v1.1/output'
  rttovfilterdir='srf'
  outdir='gas'

  spawn,'ls -1 '+rttovfilterdir, rttov_filters

  q=where(strmid(rttov_filters,0,2) eq 'rt',count) ; select the filterfiles
 q=where(strmid(rttov_filters,0,2) eq 'rt'  and strmid(rttov_filters,7,3) eq 'mtg'  ,count) ; select the filterfiles

  for i =0,count-1 do create_orac_lut_calc_gas_opd, modtrantapedir, rttovfilterdir, rttov_filters[q[i]],outdir 

end