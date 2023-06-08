  colour_ps,8
  rttovfilterdir ='/home/grainger/project-oraclut/create_orac_lut/input_files/srf'

  spawn,'ls -1 '+rttovfilterdir, rttov_filters

  q=where(strmid(rttov_filters,0,2) eq 'rt',count) ; select the filterfiles
 
  for j =0,count-1 do begin 
; for j =285,285 do begin 
    print,rttov_filters[q[j]]
    srf= rttovfilterdir+'/'+rttov_filters[q[j]]
    read_srfstr,srf,srfstr
    plot,srfstr.wvn,srfstr.srf,title=rttov_filters[q[j]]
;   plot,srfstr.wvl,srfstr.srf,title=rttov_filters[q[j]]
    
    fe = 0.001  ; ie integrals must agree to better than .1 %
    segment,srfstr.wvn,srfstr.srf,x,y,fe,minn=12   
;   segment,srfstr.wvl,srfstr.srf,x,y,fe,minn=12   
  
      
   
    oplot,x,y,color=4
    oplot,x,y,color=4,psym=2,symsize =1

    print,j,n_elements(srfstr.wvn),n_elements(x),integrate_trapeziodal(srfstr.wvn,srfstr.srf),integrate_trapeziodal(x,y),abs(integrate_trapeziodal(srfstr.wvn,srfstr.srf)-integrate_trapeziodal(x,y))/integrate_trapeziodal(srfstr.wvn,srfstr.srf)
    empty
 ;  wait,0.5
 ;  if (n_elements(x) gt 30) then stop
  endfor
end