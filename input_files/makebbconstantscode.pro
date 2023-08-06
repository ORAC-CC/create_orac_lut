pro xx, x, a, f, pder
 
  c1=1.191043d-5
  c2=1.4387769d0

  nu=a(0)
  t1=a(1)  
  t2=a(2)  

  b1 = c1 * nu^3   
  b2 = c2 * nu
   
  teff=x*t2+t1
  c= exp(b2/teff)
  cm1 = c-1
  f = b1/(cm1)
   
  dn1 = 3*c1*nu^2/cm1 - b1*c2*c/(cm1*cm1*teff)
  dt1 = B1*B2 * C /(cm1*teff)^2
  dt2 = B1*B2*x*c/  (cm1*teff)^2   
   
   pder =[[dn1],[dt1],[dt2]]
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


  makepro = 0 ; TURN THIS ON TO MAKE FILE USED IN LUT CREATION
 

  c1=1.191043d-5  ; mW\,m$^{-2}$\,(cm$^{-1}$)$^{-4}$\,sr$^{-1}$
  c2=1.4387769d0  ; cm K
  rttovfilterdir='~/project-oraclut/create_orac_lut/input_files/srf'
  spawn,'ls -1 '+rttovfilterdir, rttov_filters

  q=where(strmid(rttov_filters,0,2) eq 'rt',count) ; select the filterfiles

; code to make black body fits used in orac
  
  T0 = 150d0
  dt = .1d0
  tn = 2001
  Temperature=T0+dt*indgen(Tn)
  r=fltarr(tn)  
  weights=replicate(1.0,tn)

  X="(A40,1x,'Lambda:',F5.2,1x,'B1:',F7.0,1x,'B2:',F5.0,1x,'nu0:',F5.0,1x,'T1:',F6.4,1x,'T2:',F6.4,1x,'Fit:',F6.4)"
  
  max_error = 0
  if (Makepro) then openw,LUN,'bbconstants.pro',/get_lun
  if (Makepro) then printf,lun,'pro bbconstants,filename,b1,b2,t1,t2'
  if (Makepro) then printf,lun,'  case filename of'
 
  for j =0,count-1 do begin 

    srf= rttovfilterdir+'/'+rttov_filters[q[j]]
    read_srfstr,srf,srfstr
    If srfstr.wvl_centre gt 3 then begin
      a=[1d4/srfstr.wvl_centre,0.001,1]
         
      for i = 0,tn-1 do begin
        r(i) = int_tabulated(srfstr.wvn, planck_wvn(srfstr.wvn, Temperature(i))*srfstr.srf)/  int_tabulated(srfstr.wvn,srfstr.srf)
      endfor
  ;    yfit = curvefit(temperature,r,weights,a,function_name='xx',/noderivative,/double,itmaX=1000,tol=1d-3)
       yfit = curvefit(temperature,r,weights,a,function_name='xx',/double)
      b1 = c1 * a[0]^3   
      b2 = c2 * a[0]
    endif else begin
      a=[0,0,0]
      b1 = 0
      b2 = 0
    Endelse
    If srfstr.wvl_centre gt 3 then begin
      rms =sqrt(total((yfit -r)^2)/tn)
      max_error = max_error > rms
      print,rttov_filters[q[j]],srfstr.wvl_centre, b1,b2, a[0],a[1],a[2] , rms,Format=X
    endif
    if (Makepro) then printf,lun,string(39B),srfstr.filename,string(39B), Format="( a1,a,a1, ': begin' )"
    if (Makepro) then printf,lun, b1,b2, a[1],a[2],Format="( '    b1=',G12.6,' & b2=',G12.6,' & t1=',G12.6,' & t2=',G12.6)"
;    print,      b1,b2, a[1],a[2],Format="( '    b1=',G12.6,' & b2=',G12.6,' & t1=',G12.6,' & t2=',G12.6)"
    if (Makepro) then printf,lun,'end'
  endfor
  if (Makepro) then printf,lun,' endcase'
  if (Makepro) then printf,lun,'end'
  if (Makepro) then close,lun
  if (Makepro) then free_lun,lun
print,'maximum error : ', max_error
end