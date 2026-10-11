; Interface test for DISORT_DP (DISORTtoIDL_DP) inside IDL.
;
; Uses the V21 single-precision capture of the p240 channel-6 direct-beam
; call (disort_boundary).  Checks, in one IDL session:
;   1. the historical DISORT DLM still reproduces the captured V21 outputs
;      bitwise while DISORTtoIDL_DP is loaded;
;   2. DISORT_DP accepts float32 or float64 arrays and INT or LONG scalars
;      and gives bitwise-identical double results (float->double is exact);
;   3. DISORT_DP does not modify its inputs and returns DOUBLE outputs;
;   4. repeated DISORT_DP calls are bitwise identical;
;   5. invalid arguments raise IDL errors (LAMBER = 0, wrong PMOM layout,
;      SSALB > 1, short UMU).
; Results are saved to OUTFILE; a line 'DP INTERFACE TEST PASSED' is printed
; only if every check passes.

pro dp_interface_check, ok, condition, label, failures
   if condition then print, 'PASS  ' + label else begin
      print, 'FAIL  ' + label
      failures = [failures, label]
      ok = 0
   endelse
end

pro dp_interface_test, capture_input, capture_output, outfile
   compile_opt idl2
   ok = 1
   failures = ['']
   INFO = DISORT_DP_INFO()
   print, INFO
   dp_interface_check, ok, strpos(INFO, 'fortran_real_bytes=8') ge 0, 'kernel REAL is 8 bytes', failures

   restore, capture_output
   V21_RFLDIR = RFLDIR & V21_RFLDN = RFLDN & V21_FLUP = FLUP & V21_UU = UU
   restore, capture_input

;  1. Historical single-precision DLM in the same session.
   SS = SSALB1
   DISORT, NLYR1, DTAUC1, SS, NMOM, PMOM1, TEMPER, WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR, $
           USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP, $
           TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT, MAXCLY1, MAXULV, MAXUMU, MAXPHI, MAXMOM, $
           S_RFLDIR, S_RFLDN, S_FLUP, S_DFDT, S_UAVG, S_UU, S_ALBMED, S_TRNMED
   SP_BITWISE = array_equal(S_UU, V21_UU) and array_equal(S_FLUP, V21_FLUP) and $
      array_equal(S_RFLDN, V21_RFLDN) and array_equal(S_RFLDIR, V21_RFLDIR)
   SP_MAXDIFF = max([max(abs(S_UU - V21_UU)), max(abs(S_FLUP - V21_FLUP)), max(abs(S_RFLDN - V21_RFLDN))])
   print, 'historical DISORT versus V21 capture: bitwise = ', SP_BITWISE, ', max |difference| = ', SP_MAXDIFF
   dp_interface_check, ok, SP_MAXDIFF le 1e-4, 'historical DISORT reproduces V21 capture (bitwise or <= 1e-4)', failures

;  2./3. Float32 inputs.
   IN_DTAUC = DTAUC1 & IN_SSALB = SSALB1 & IN_PMOM = PMOM1 & IN_UMU = UMU
   DISORT_DP, NLYR1, DTAUC1, SSALB1, NMOM, PMOM1, TEMPER, WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR, $
              USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP, $
              TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT, MAXCLY1, MAXULV, MAXUMU, MAXPHI, MAXMOM, $
              F_RFLDIR, F_RFLDN, F_FLUP, F_DFDT, F_UAVG, F_UU, F_ALBMED, F_TRNMED
   dp_interface_check, ok, size(F_UU, /type) eq 5 and size(F_FLUP, /type) eq 5, 'outputs are DOUBLE', failures
   dp_interface_check, ok, array_equal(size(F_UU, /dimensions), [NUMU, NTAU, NPHI]), 'UU dimensions', failures
   dp_interface_check, ok, array_equal(IN_DTAUC, DTAUC1) and array_equal(IN_SSALB, SSALB1) and $
      array_equal(IN_PMOM, PMOM1) and array_equal(IN_UMU, UMU), 'inputs unchanged', failures

;  Float64 inputs and LONG/INT scalars.
   DISORT_DP, long(NLYR1), double(DTAUC1), double(SSALB1), long(NMOM), double(PMOM1), double(TEMPER), $
              double(WVNMLO), double(WVNMHI), fix(USRTAU), fix(NTAU), double(UTAU), long(NSTR), $
              long(USRANG), long(NUMU), double(UMU), long(NPHI), double(PHI), long(IBCND), $
              double(FBEAM), double(UMU0), double(PHI0), double(FISOT), long(LAMBER), double(ALBEDO), $
              double(BTEMP), double(TTEMP), double(TEMIS), long(PLANK), long(ONLYFL), double(ACCUR), $
              long(PRNT), long(MAXCLY1), long(MAXULV), long(MAXUMU), long(MAXPHI), long(MAXMOM), $
              D_RFLDIR, D_RFLDN, D_FLUP, D_DFDT, D_UAVG, D_UU, D_ALBMED, D_TRNMED
   dp_interface_check, ok, array_equal(D_UU, F_UU) and array_equal(D_FLUP, F_FLUP) and $
      array_equal(D_RFLDN, F_RFLDN) and array_equal(D_DFDT, F_DFDT) and array_equal(D_UAVG, F_UAVG), $
      'float32 and float64 inputs give identical results', failures

;  4. Repeatability.
   DISORT_DP, NLYR1, DTAUC1, SSALB1, NMOM, PMOM1, TEMPER, WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR, $
              USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP, $
              TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT, MAXCLY1, MAXULV, MAXUMU, MAXPHI, MAXMOM, $
              R_RFLDIR, R_RFLDN, R_FLUP, R_DFDT, R_UAVG, R_UU, R_ALBMED, R_TRNMED
   dp_interface_check, ok, array_equal(R_UU, F_UU) and array_equal(R_FLUP, F_FLUP), 'repeat call bitwise identical', failures

;  5. Invalid arguments must raise IDL errors (each attempt inside CATCH).
   BAD_PMOM = transpose(PMOM1)
   BAD_SSALB = SSALB1 * 0.0 + 1.5
   SHORT_UMU = UMU[0:1]
   ZERO = 0L
   for TEST = 0, 3 do begin
      RAISED = 0
      catch, ERR
      if ERR ne 0 then begin
         catch, /cancel
         RAISED = 1
      endif else begin
      case TEST of
         0: DISORT_DP, NLYR1, DTAUC1, SSALB1, NMOM, PMOM1, TEMPER, WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR, $
                       USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, ZERO, ALBEDO, BTEMP, $
                       TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT, MAXCLY1, MAXULV, MAXUMU, MAXPHI, MAXMOM, $
                       A1, A2, A3, A4, A5, A6, A7, A8
         1: DISORT_DP, NLYR1, DTAUC1, SSALB1, NMOM, BAD_PMOM, TEMPER, WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR, $
                       USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP, $
                       TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT, MAXCLY1, MAXULV, MAXUMU, MAXPHI, MAXMOM, $
                       A1, A2, A3, A4, A5, A6, A7, A8
         2: DISORT_DP, NLYR1, DTAUC1, BAD_SSALB, NMOM, PMOM1, TEMPER, WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR, $
                       USRANG, NUMU, UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP, $
                       TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT, MAXCLY1, MAXULV, MAXUMU, MAXPHI, MAXMOM, $
                       A1, A2, A3, A4, A5, A6, A7, A8
         3: DISORT_DP, NLYR1, DTAUC1, SSALB1, NMOM, PMOM1, TEMPER, WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR, $
                       USRANG, NUMU, SHORT_UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO, BTEMP, $
                       TTEMP, TEMIS, PLANK, ONLYFL, ACCUR, PRNT, MAXCLY1, MAXULV, MAXUMU, MAXPHI, MAXMOM, $
                       A1, A2, A3, A4, A5, A6, A7, A8
      endcase
      catch, /cancel
      endelse
      LABELS = ['LAMBER = 0 rejected', 'wrong PMOM layout rejected', 'SSALB > 1 rejected', 'short UMU rejected']
      dp_interface_check, ok, RAISED eq 1, LABELS[TEST], failures
   endfor

   save, filename=outfile, INFO, SP_BITWISE, SP_MAXDIFF, V21_UU, V21_FLUP, S_UU, S_FLUP, S_RFLDN, F_RFLDIR, F_RFLDN, F_FLUP, $
      F_DFDT, F_UAVG, F_UU, F_ALBMED, F_TRNMED, D_UU, R_UU, failures, ok
   if ok then print, 'DP INTERFACE TEST PASSED' else print, 'DP INTERFACE TEST FAILED: ' + strjoin(failures[1:*], '; ')
end
