; This library contains routines which generate the perturbed mixing ratios
; and/or component mode radii to generate the required aerosol/cloud class
; effective radii (create_range) ... and call the scattering code (Mie or T-
; matrix) to generate the scattering properties for each component, which are
; then combined and passed into DISORT.
;
; This routine are essentially copies of their equivalents from the
; old LUT code.
;
; HISTORY
; 21/06/13, G Thomas: Original version.


; procedure create_bwgp
;
; Wrapper procedure for the various scattering codes used in generating optical
; properties for ORAC LUTs. This procedure is very similar, but not identical
; to its name-sake in the old AOPP LUT generating code.
;
; INPUT ARGUMENTS:
;
; INPUT KEYWORDS:
;
; OUTPUT ARGUMENTS:
;
; HISTORY:
; 21/06/13, G Thomas: Original version.

; call parameters
;mmstr.distname[c], lut_Rm[c,r], mmstr.S[c], AerM[*,*,c], srfstrarr[*].wvl[*], QV, Bext1, w1, g1, Phs1, scode=scode, tmatrix_path=tmatrix_path, eps=epsvals, neps=nepsvals
;mmstr.distname[c], lut_Rm[c,r], mmstr.S[c], AerM550[c] , 0.55               , QV, Bext1, w1, g1, Phs1, scode=scode, tmatrix_path=tmatrix_path, eps=epsvals, neps=nepsvals, Vavg=Vavg1
  
pro create_bwgp, distname, Rm, S, RI, wl, Dqv, Bext, w, g, Phi, scode=scode, tmatrix_path=tmatrix_path, eps=eps, neps=neps, Vavg=Vavg

  IF n_elements(RI) ne n_elements(wl) THEN message, 'create_bwgp: Array size mismatch!'

  Inp = n_elements(Dqv)
  Inw = n_elements(wl)
  wn = 1.0 / wl
  Bext = fltarr(Inw)
  w = Bext
  g = Bext
  Phi = fltarr(Inp,Inw)

; check if we're using T-Matrix calculations
  dotmatrix = 0
  IF KEYWORD_SET(scode) THEN $
    IF (strlowcase(scode) eq 'tmatrix') THEN begin
      dotmatrix = 1 
      q = WHERE (wl lT 6, count)  ; only care about tmatrix calculations for wavelengths less thsan 6 um
;     Check that the user has passed the asymmetry information required for the Dubovik code.
      IF (count GT 0) AND ((n_elements(eps) eq 0) or (n_elements(neps) eq 0)) THEN message, 'Asymmetry parameters, "eps",  and their relative numbers, "neps", are required to use T-Matrix calculations' 
    ENDIF

  FOR i = 0, Inw-1 DO BEGIN ; LOOP OVER WAVELENGTH OR CHANNEL
    IF (wl(i) NE 0) THEN BEGIN ; If using native resolution spectral response functions the number of wavelengths 
;                                may be different between channels so ignore wl which are zero

;     check if we're using T-Matrix or Mie scattering...
;     note that we don't bother with non-sphericity for thermal-only channels because:
;     (a) Scattering should be of limited importance
;     (b) Dubovik's tabulation doesn't cover the range of likely refractive indices for such situations

      IF dotmatrix and (wl(i) lt 6) THEN BEGIN
;       Check for the absorption coefficient is within the Dubovik range (-0.5 to -0.0005). if it's outside, this range then:
;       * If it is greater than -0.0005, then print a warning message, but assign it to the extreme value.
;         This can often happen with essentially non-absorbing particles, so the change from -1e-6 to -5e-4 
;         will basically have no effect on the result.
;       * If it is less than -0.5, issue an error message and stop.

        IF imaginary(RI(i)) gt -5e-4 THEN BEGIN
          print,''
          message, /info, 'Warning: Imaginary RI altered from '+ strtrim(imaginary(RI[i]),2)+ ' to -5e-04, so as to lie within the Dubovik range.'
          RItmp = complex(float(RI[i]),-5e-4)
        ENDIF ELSE $
          IF imaginary(RI[i]) lt -0.5 THEN BEGIN
            print,''
            message, 'Imaginary RI is too negative for the Dubovik LUT range: k = '+strtrim(imaginary(RI[i]),2)+'; min value = -0.5'
          ENDIF ELSE $
            RItmp = RI[i]

        dubovik_lognormal_multiple_eps, tmatrix_path, 1.0, Rm, S, wn[i], RItmp, eps, neps, Dqv=Dqv, Bexttmp, Bscatmp, wtmp, gtmp, ph, /no_mie, /renorm_ph, /silent
        Phi[*,i] = ph
      ENDIF ELSE BEGIN
;       We're not using T-Matrix, so call the normal Mie scattering code
        mie_size_dist_new, distname, 1.0, [Rm, S, 0.001, 100.0], wn[i], RI[i], Dqv=Dqv, /dlm, xres=0.4, Bexttmp, Bscatmp, wtmp, gtmp, SPM, Vavg=Vavg
        Phi[*,i] = SPM[0,*]
	;	Bexttmp2=Bexttmp
	;	Bscatmp2=Bscatmp
	;	wtmp2=wtmp
	;	gtmp2=gtmp
	;	SPM2=SPM
	;	Vavg2=Vavg
   ; mie_size_dist_new, distname, 1.0, [Rm, S, 0.001, 100.0], wn[i], RI[i], Dqv=Dqv, /dlm, xres=0.4, Bexttmp2, Bscatmp2, wtmp2, gtmp2, SPM2, Vavg=Vavg2
	;	print,Bexttmp,Bexttmp2
	;	print,Bscatmp,Bscatmp2
	;	print,wtmp,wtmp2
	;	print,gtmp,gtmp2
	;	print,Vavg,Vavg2
      ENDELSE
    Bext[i] = Bexttmp
    w[i] = wtmp
    g[i] = gtmp
    ENDIF
  ENDFOR
END
