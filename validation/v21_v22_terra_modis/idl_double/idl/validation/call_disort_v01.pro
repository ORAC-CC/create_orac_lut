; Validation shadow of create_orac_lut/call_disort.pro for the double-precision
; DISORT study (validation/v21_v22_terra_modis/idl_double).
;
; Lines from the signature to the end of the thermal/compression block are a
; verbatim copy of the V21 production call_disort.pro (lines 11-92); the
; historical single-precision DISORT call follows unchanged and its results
; are returned, so the generator computes exactly the V21 LUT.
;
; In addition, when ORACLUT_DP_CAPTURE_DIR is set, every call is solved again
; by DISORT_DP (DISORTtoIDL_DP, 8-byte REAL kernel) at each stream count in
; ORACLUT_DP_NSTR_LIST, from copies of the inputs taken before the historical
; call (that call dithers SSALB = 1 in place), and once more at the V21 stream
; count to test repeatability.  Inputs, single-precision outputs, all
; double-precision outputs and timings are saved, one file per call.

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
   if n_elements(UTAU)      ne NTAU   then stop, 'UTAU dimension missmatch'
   if n_elements(UMU)       ne NUMU   then stop, 'UMU dimension missmatch'
   if n_elements(PHI)       ne NPHI   then stop, 'PHI dimension missmatch'

;  Set thermal related inputs based on whether the plank keyword was set or not
   if keyword_set(plank) then begin
      if (n_elements(wnlo) eq 0) or (n_elements(wnhi) eq 0)  then $
         message, 'Wavenumber must be specified for emission calculations!'
      if (n_elements(temp) eq 0) then $
         message, 'Temperature must be specified for emission calculations!'

     DTAUC1 = DTAUC
     SSALB1 = SSALB
     PMOM1  = PMOM
      WVNMLO    = wnlo
      WVNMHI    = wnhi
      if n_elements(TEMP)    ne NLYR1+1   then stop, 'TEMP dimension missmatch'
      TEMPER    = temp
      BTEMP     = 0.0
      TTEMP     = 0.0
      TEMIS     = 0.0
   endif else begin
;    Compress the layer stack by combining layers with identical values fr SSALB
;    and PMOM
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
   
   
   
      PLANK     = 0
      TEMPER = fltarr(nlyr)
      WVNMLO    = 0.0
      WVNMHI    = 0.0
      TEMPER[*] = 0.0
      BTEMP     = 0.0
      TTEMP     = 0.0
      TEMIS     = 0.0
   endelse

;  Pre-call copies: the historical kernel overwrites SSALB = 1 in place.
   DTAUC1_IN = DTAUC1
   SSALB1_IN = SSALB1
   PMOM1_IN  = PMOM1
   TEMPER_IN = TEMPER
   UTAU_IN   = UTAU
   UMU_IN    = UMU
   PHI_IN    = PHI

;  Call DISORT (historical single precision, V21)
   T0 = systime(/seconds)
   DISORT, NLYR1, DTAUC1, SSALB1, NMOM, PMOM1, TEMPER, WVNMLO, WVNMHI, USRTAU, $
           NTAU, UTAU, NSTR, USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0, $
           PHI0, FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, PLANK, ONLYFL, $
           ACCUR, PRNT, MAXCLY1, MAXULV, MAXUMU, MAXPHI, MAXMOM, RFLDIR, $
           RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED
   TIME_SP = systime(/seconds) - T0

   capture_dir = getenv('ORACLUT_DP_CAPTURE_DIR')
   if capture_dir eq '' then return
   EXISTING = file_search(capture_dir + '/call_*.sav', count=NEXISTING)
   call_index = long(NEXISTING) + 1L
   NSTR_LIST = long(strsplit(getenv('ORACLUT_DP_NSTR_LIST'), ',', /extract))
   NK = n_elements(NSTR_LIST)
   DP_RFLDIR = dblarr(NTAU, NK)
   DP_RFLDN  = dblarr(NTAU, NK)
   DP_FLUP   = dblarr(NTAU, NK)
   DP_DFDT   = dblarr(NTAU, NK)
   DP_UAVG   = dblarr(NTAU, NK)
   DP_UU     = dblarr(NUMU, NTAU, NPHI, NK)
   DP_ALBMED = dblarr(NUMU, NK)
   DP_TRNMED = dblarr(NUMU, NK)
   DP_TIME   = dblarr(NK)
   for K = 0, NK - 1 do begin
      T0 = systime(/seconds)
      DISORT_DP, NLYR1, DTAUC1_IN, SSALB1_IN, NMOM, PMOM1_IN, TEMPER_IN, WVNMLO, WVNMHI, $
                 USRTAU, NTAU, UTAU_IN, NSTR_LIST[K], USRANG, NUMU, UMU_IN, NPHI, PHI_IN, $
                 IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, $
                 PLANK, ONLYFL, ACCUR, PRNT, MAXCLY1, MAXULV, MAXUMU, MAXPHI, MAXMOM, $
                 D_RFLDIR, D_RFLDN, D_FLUP, D_DFDT, D_UAVG, D_UU, D_ALBMED, D_TRNMED
      DP_TIME[K] = systime(/seconds) - T0
      DP_RFLDIR[*, K] = D_RFLDIR
      DP_RFLDN[*, K]  = D_RFLDN
      DP_FLUP[*, K]   = D_FLUP
      DP_DFDT[*, K]   = D_DFDT
      DP_UAVG[*, K]   = D_UAVG
      DP_UU[*, *, *, K] = D_UU
      DP_ALBMED[*, K] = D_ALBMED
      DP_TRNMED[*, K] = D_TRNMED
   endfor
;  Repeatability: the V21 stream count again, in the same process.
   DISORT_DP, NLYR1, DTAUC1_IN, SSALB1_IN, NMOM, PMOM1_IN, TEMPER_IN, WVNMLO, WVNMHI, $
              USRTAU, NTAU, UTAU_IN, NSTR, USRANG, NUMU, UMU_IN, NPHI, PHI_IN, $
              IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, $
              PLANK, ONLYFL, ACCUR, PRNT, MAXCLY1, MAXULV, MAXUMU, MAXPHI, MAXMOM, $
              R_RFLDIR, R_RFLDN, R_FLUP, R_DFDT, R_UAVG, R_UU, R_ALBMED, R_TRNMED
   TYPE_DP_UU = size(D_UU, /type)
   CALL_INDEX_SAVED = call_index
   NSTR_V21 = NSTR
   save, filename=capture_dir + '/call_' + string(call_index, format='(I05)') + '.sav', $
      CALL_INDEX_SAVED, NLYR1, DTAUC1_IN, SSALB1_IN, NMOM, PMOM1_IN, TEMPER_IN, $
      WVNMLO, WVNMHI, USRTAU, NTAU, UTAU_IN, NSTR_V21, USRANG, NUMU, UMU_IN, NPHI, $
      PHI_IN, IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, $
      PLANK, ONLYFL, ACCUR, PRNT, MAXCLY1, RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, $
      ALBMED, TRNMED, TIME_SP, NSTR_LIST, DP_RFLDIR, DP_RFLDN, DP_FLUP, DP_DFDT, $
      DP_UAVG, DP_UU, DP_ALBMED, DP_TRNMED, DP_TIME, R_RFLDIR, R_RFLDN, R_FLUP, $
      R_DFDT, R_UAVG, R_UU, R_ALBMED, R_TRNMED, TYPE_DP_UU
end
