; Diagnostic shadow of create_orac_lut/call_disort.pro.
;
; The numerical statements through the DISORT call are deliberately identical
; to the V21 production wrapper.  For the one direct-beam call in the compact
; validation case, this routine saves the complete state on both sides of the
; DLM boundary and exits before create_orac_cloud_lut can write a LUT.

pro call_disort, DTAUC, SSALB, PMOM, UTAU, UMU, PHI, FBEAM, UMU0, FISOT, $
                 RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED, $
                 plank=plank, wnlo=wnlo, wnhi=wnhi, temp=temp, nlayer=nlayer

   common DISORT_VARS, NLYR, NMOM, NTAU, NUMU, NPHI, NSTR, MAXCLY, MAXULV, $
                       MAXUMU, MAXPHI, MAXMOM, USRTAU, USRANG, PHI0, IBCND, $
                       LAMBER, ALBEDO, ONLYFL, ACCUR, PRNT, HEADER

   if keyword_set(nlayer) then begin
      NLYR1 = nlayer
      MAXCLY1 = nlayer
   endif else begin
      NLYR1 = NLYR
      MAXCLY1 = MAXCLY
   endelse

   if n_elements(DTAUC)     ne NLYR1  then stop, 'DTAUC dimension missmatch'
   if n_elements(SSALB)     ne NLYR1  then stop, 'SSALB dimension missmatch'
   if n_elements(PMOM[*,0]) ne NMOM+1 then stop, 'PMOM dimension 1 missmatch'
   if n_elements(PMOM[0,*]) ne NLYR1  then stop, 'PMOM dimension 2 missmatch'
   if n_elements(UTAU)      ne NTAU   then stop, 'UTAU dimension missmatch'
   if n_elements(UMU)       ne NUMU   then stop, 'UMU dimension missmatch'
   if n_elements(PHI)       ne NPHI   then stop, 'PHI dimension missmatch'

   if keyword_set(plank) then begin
      if (n_elements(wnlo) eq 0) or (n_elements(wnhi) eq 0) then $
         message, 'Wavenumber must be specified for emission calculations!'
      if n_elements(temp) eq 0 then $
         message, 'Temperature must be specified for emission calculations!'

      DTAUC1 = DTAUC
      SSALB1 = SSALB
      PMOM1  = PMOM
      WVNMLO = wnlo
      WVNMHI = wnhi
      if n_elements(TEMP) ne NLYR1+1 then stop, 'TEMP dimension missmatch'
      TEMPER = temp
      BTEMP  = 0.0
      TTEMP  = 0.0
      TEMIS  = 0.0
   endif else begin
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
      PLANK  = 0
      TEMPER = fltarr(nlyr)
      WVNMLO = 0.0
      WVNMHI = 0.0
      TEMPER[*] = 0.0
      BTEMP  = 0.0
      TTEMP  = 0.0
      TEMIS  = 0.0
   endelse

   capture_this_call = (FBEAM eq 100.0) and $
      (abs(UMU0 - cos(80.0*!dtor)) lt 1.0e-7)
   if capture_this_call then begin
      input_file = getenv('ORACLUT_V21_DISORT_INPUT')
      output_file = getenv('ORACLUT_V21_DISORT_OUTPUT')
      if (input_file eq '') or (output_file eq '') then $
         message, 'ORACLUT V21 capture paths are not set'

      DTAUC_ORIGINAL = DTAUC
      SSALB_ORIGINAL = SSALB
      PMOM_ORIGINAL = PMOM
      TYPE_DTAUC = size(DTAUC1, /type)
      TYPE_SSALB = size(SSALB1, /type)
      TYPE_PMOM = size(PMOM1, /type)
      TYPE_TEMPER = size(TEMPER, /type)
      TYPE_UTAU = size(UTAU, /type)
      TYPE_UMU = size(UMU, /type)
      TYPE_PHI = size(PHI, /type)
      TYPE_PRNT = size(PRNT, /type)

      save, filename=input_file, DTAUC_ORIGINAL, SSALB_ORIGINAL, PMOM_ORIGINAL, $
         NLYR1, DTAUC1, SSALB1, NMOM, PMOM1, TEMPER, WVNMLO, WVNMHI, $
         USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, UMU, NPHI, PHI, IBCND, $
         FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, $
         PLANK, ONLYFL, ACCUR, PRNT, MAXCLY1, MAXULV, MAXUMU, MAXPHI, $
         MAXMOM, TYPE_DTAUC, TYPE_SSALB, TYPE_PMOM, TYPE_TEMPER, TYPE_UTAU, $
         TYPE_UMU, TYPE_PHI, TYPE_PRNT
   endif

   DISORT, NLYR1, DTAUC1, SSALB1, NMOM, PMOM1, TEMPER, WVNMLO, WVNMHI, USRTAU, $
           NTAU, UTAU, NSTR, USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0, $
           PHI0, FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, PLANK, ONLYFL, $
           ACCUR, PRNT, MAXCLY1, MAXULV, MAXUMU, MAXPHI, MAXMOM, RFLDIR, $
           RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED

   if capture_this_call then begin
      save, filename=output_file, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED
      print, 'ORACLUT V21 DISORT boundary capture complete'
      exit, status=0
   endif
end
