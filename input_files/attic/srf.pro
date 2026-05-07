


  rttovfilterdir='~/project-oraclut/create_orac_lut/input_files/srf'
  spawn,'ls -1 '+rttovfilterdir, rttov_filters
  qq=where(strmid(rttov_filters,0,2) eq 'rt',count) ; select the filterfiles
 
  for j =0,count-1 do begin 

    srf= rttovfilterdir+'/'+rttov_filters[qq[j]]
    read_srfstr,srf,srfstr
    print,srfstr.filename,srfstr.wvl_centre
    xm = max(srfstr.srf)
    q= where(srfstr.srf gt xm/2)
    print, srfstr.wvn(q(0)),srfstr.srf(q(0))
    print,xm
    print, srfstr.wvn(q(-1)),srfstr.srf(q(-1))
    
  endfor

end