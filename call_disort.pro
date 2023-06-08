; Check the input variables, define the remaining DISORT arguments and call
; DISORT.
;
; HISTORY:
; 21/06/13, G Thomas: Original version.
; 19/10/16, G McGarragh: Compress the layer stack by combining layers with
;    identical values for SSALB and PMOM.  This can lead to a significant speed
;    increase in some cases like when Rayleigh scattering is the only radiative
;    component in layers without cloud.

pro call_disort, DTAUC, SSALB, PMOM, UTAU, UMU, PHI, FBEAM, UMU0, FISOT, $
                 RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED, $
                 plank=plank, wnlo=wnlo, wnhi=wnhi, temp=temp, nlayer=nlayer
								 
   common DISORT_VARS, NLYR, NMOM, NTAU, NUMU, NPHI, NSTR, MAXCLY, MAXULV, $
                       MAXUMU, MAXPHI, MAXMOM, USRTAU, USRANG, PHI0, IBCND, $
                       LAMBER, ALBEDO, ONLYFL, ACCUR, PRNT, HEADER
;   print, 'Calling DISORT'

;  If the nlayer keyword has been set, we are only running on a sub-set of the
;  atmospheric layers...
   if keyword_set(nlayer) then begin
      NLYR1 = nlayer
      MAXCLY1 = nlayer
   endif else begin
      NLYR1 = NLYR
      MAXCLY1 = MAXCLY
   endelse

;  Do some array size checks, just to make sure all is as expected
   if n_elements(DTAUC)     ne NLYR1  then stop, 'DTAUC dimension missmatch'
   if n_elements(SSALB)     ne NLYR1  then stop, 'SSALB dimension missmatch'
   if n_elements(PMOM[*,0]) ne NMOM+1 then stop, 'PMOM dimension 1 missmatch'
   if n_elements(PMOM[0,*]) ne NLYR1  then stop, 'PMOM dimension 2 missmatch'
;  if n_elements(TEMPER)    ne NLYR   then stop, 'TEMPER dimension missmatch'
   if n_elements(UTAU)      ne NTAU   then stop, 'UTAU dimension missmatch'
   if n_elements(UMU)       ne NUMU   then stop, 'UMU dimension missmatch'
   if n_elements(PHI)       ne NPHI   then stop, 'PHI dimension missmatch'

;  Compress the layer stack by combining layers with identical values fr SSALB
;  and PMOM
   MAXCLY1 = NLYR1

   DTAUC1 = fltarr(NLYR1)
   SSALB1 = fltarr(NLYR1)
   PMOM1  = fltarr(NMOM+1,NLYR1)

   ii = 0
   DTAUC1[0]  = DTAUC[0]
   SSALB1[0]  = SSALB[0]
   PMOM1[*,0] = PMOM[*,0]
   for i = 1, NLYR1 - 1 do begin
      if (SSALB1[ii] eq SSALB[i]) and array_equal(PMOM1[*,ii], PMOM[*,i]) then begin
           DTAUC1[ii] += DTAUC[i]
      endif else begin
           ii++
           DTAUC1[ii] = DTAUC[i]
           SSALB1[ii] = SSALB[i]
           PMOM1[*,ii] = PMOM[*,i]
      endelse
   endfor

   NLYR1 = ii + 1

   TEMPER = fltarr(nlyr)

;  Set thermal related inputs based on whether the plank keyword was set or not
   if keyword_set(plank) then begin
      if (n_elements(wnlo) eq 0) or (n_elements(wnhi) eq 0)  then $
         message, 'Wavenumber must be specified for emission calculations!'
      if (n_elements(temp) eq 0) then $
         message, 'Temperature must be specified for emission calculations!'

      WVNMLO    = wnlo
      WVNMHI    = wnhi
      TEMPER[*] = temp
      BTEMP     = 0.0
      TTEMP     = 0.0
      TEMIS     = 0.0
   endif else begin
      PLANK     = 0
      WVNMLO    = 0.0
      WVNMHI    = 0.0
      TEMPER[*] = 0.0
      BTEMP     = 0.0
      TTEMP     = 0.0
      TEMIS     = 0.0
   endelse

					; stop
;  Call DISORT
   DISORT, NLYR1, DTAUC1, SSALB1, NMOM, PMOM1, TEMPER, WVNMLO, WVNMHI, USRTAU, $
           NTAU, UTAU, NSTR, USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0, $
           PHI0, FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, PLANK, ONLYFL, $
           ACCUR, PRNT, MAXCLY1, MAXULV, MAXUMU, MAXPHI, MAXMOM, RFLDIR, $
           RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED
end
