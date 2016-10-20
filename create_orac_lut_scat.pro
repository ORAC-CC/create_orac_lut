; This library contains routines which generate the perturbed mixing ratios
; and/or component mode radii to generate the required aerosol/cloud class
; effective radii (create_range) ... and call the scattering code (Mie or T-
; matrix) to generate the scattering properties for each component, which are
; then combined and passed into DISORT.
;
; Both of these routines are essentially copies of their equivalents from the
; old LUT code.
;
; HISTORY
; 21/06/13, G Thomas: Original version.


; procedure create_range
;
; Determine the mixing ratios/mode radii to achieve a specific effective radius
; from a 1, 2, 3 or 4 mode distribution.
;
; INPUT ARGUMENTS:
; Rat   (array) Relative ratios of each of the modes
; MR    (array) Mode radius of each of the modes
; Spd   (array) Spread of each of the modes
; Radii (array) Target effective radii
;
; INPUT KEYWORDS:
; None
;
; OUTPUT ARGUMENTS:
; Nrat  (array) Mixing ratios for modes
; NMR   (array) Mode radii of each of the three modes
;
; ------------------------------------------------------------------------------
; History of create_range.pro:
; (?? June 2005) Don Grainger
; (?? July 2005) Elisa Carboni: Changed the computation of mixing ratio between
;     Er(0) and Er(2) -> request sum of mixing = 1 and more...
; (15 July 2005) Elisa Carboni: Add 2 and 4 (mode) size distributions
; (07 June 2010) Gareth Thomas: Added verbose keyword (if not specified, text
;     output will be minimal)
; ------------------------------------------------------------------------------
;
; History:
; 04/02/13, G Thomas: Copied directly from "create_range.pro" used with the old
;                     AOPP LUT generation code.
; 25/07/14, G Thomas: Changed so that code deals with size modes (rather than
;                     assuming all components have different sizes).
; 06/02/15, G Thomas: Added old_run_style keyword, which makes the code behave
;                     as it did before the 25/07/14 fix.
; 17/02/15, G Thomas: Bug fix to the indexing of the Mode_Ratq array (which
;                     stores the mixing prior mixing ratio of all components in
;                     each size mode).

pro create_range, Rat, MR, Spd, Radii, Nrat_o, NMR_o, verbose=verbose, $
                  old_run_style=old_run_style

   if keyword_set(verbose) then message,/info, ' '

;  Calculate the effective radius of each of the modes
   Er = MR^3 * exp(4.5*alog(Spd)^2) / (Mr^2 * exp(2.0*alog(Spd)^2))

;  Calculate the effective radius of all of the modes
   Ea = total(Rat * Mr^3 * exp(4.5*alog(Spd)^2),/double) / $
        total(Rat * Mr^2 * exp(2.0*alog(Spd)^2),/double)

;  Check for different components with the same size (as with the CCI LUTs). The
;  ratio between components with identical size distributions are kept equal,
;  unless the old_run_style keyword has been set, and then we blindly assume all
;  components have different sizes and have been passed in size order.
   if keyword_set(old_run_style) then begin
      if keyword_set(verbose) then print,'old_run_style option selected - ' + $
         'all components assumed to have different sizes'
      nunq  = n_elements(Er)
      srtEr = lindgen(nunq)
      unqEr = srtEr
   endif else begin
      srtEr = sort(Er)
      unqEr = uniq(Er[srtEr]) ; Index to unique efr values in terms of srtEr...
      nunq  = n_elements(unqEr)
   endelse

   if keyword_set(verbose) then begin
      print,'Effective radius of input components    ',Er
      print,'Effective radius of input aerosol class ',Ea
      print,strtrim(nunq,2),' unique size modes detected'
   endif

   Erq = Er[srtEr[unqEr]]  ; Unique subsets of Effective radius
   Spq = Spd[srtEr[unqEr]] ;                   Log-normal spread
   MRq = MR[srtEr[unqEr]]  ;                   Mode radius

;  Calculate the mixing ratio of each size mode by summing the MR values for
;  each unique size
;  Also calculate the mixing ratio within each size mode, so that we can keep
;  this ratio constant in the output.
   Ratq = fltarr(nunq)
   if nunq ne N_elements(Rat) then begin
;     Figure out which size mode has the largest number of components in it and
;     use this to define the size of the size mode mixing ratio array.
      if keyword_set(verbose) then begin
         print,'Have multiple components in each size mode... '
      endif
      max_mode_count = unqEr[0]+1
      for i=1,nunq-1 do if (unqEr[i] - unqEr[i-1]) gt max_mode_count then $
         max_mode_count = unqEr[i] - unqEr[i-1]
      Mode_Ratq = fltarr(max_mode_count, nunq)

      Ratq[0] = total(Rat[srtEr[0:unqEr[0]]])
      Mode_Ratq[0:unqEr[0],0] = Rat[srtEr[0:unqEr[0]]] / Ratq[0]
      for i=1,nunq-1 do begin
         Ratq[i] = total(Rat[srtEr[unqEr[i-1]+1:unqEr[i]]])
         Mode_Ratq[0:unqEr[i]-unqEr[i-1]-1,i] = $
            Rat[srtEr[unqEr[i-1]+1:unqEr[i]]] / Ratq[i]
      endfor
   endif else begin
      if keyword_set(verbose) then begin
         print,'Each size mode corresponds to a single component... '
      endif
      Ratq = Rat
      Mode_Ratq = replicate(1.0,1,nunq)
   endelse

;  Define place holder output and intermediate arrays
   Nrat = Fltarr(nunq,N_elements(Radii)) ; Mode mixing ratios
   NMR  = Fltarr(nunq,N_elements(radii)) ; Mode radii

   g = Fltarr(nunq)
   f = Fltarr(nunq)

   for I = 0, N_elements(Radii)-1 do begin
      Nmr[*,I] = Mrq
      g = exp(2.0*alog(Spq)^2)
      f = exp(4.5*alog(Spq)^2)

      if (nunq eq 1) then begin
         Nrat[*,I] = 1.
         Nmr[*,I]  = Mrq[0]*Radii[I]/Erq[0]
      endif
      if (nunq eq 2) then begin
         if (Radii[I] Le Erq[0]) then begin
            Nrat[*,I] = [1., 0]
            Nmr[*,I]  = [Mrq[0]*Radii[I]/Erq[0] , Mrq[1]]
         endif
         if ((Radii[I] Gt Erq[0]) And (Radii[I] Lt Erq[1])) then begin
            den1 =  (Nmr[0,I]^2) * (Radii[I]*g[0]- Nmr[0,I] *f[0])
            den2 = -(Nmr[1,I]^2) * (Radii[I]*g[1]- Nmr[1,I] *f[1])

            NRat[0,I] = den2/(den1+den2)
            NRat[1,I] = 1-NRat[0,I]

            type_re = total(NRat[*,I] * Nmr[*,I]^3 * f) / $
                      total(NRat[*,I] * Nmr[*,I]^2 * g)
         endif
         if (Radii[I] Ge Erq[1]) then Begin
            Nrat[*,I] = [ 0 , 1]
            Nmr[*,I]  = [Mrq[0], Mrq[1] * Radii[I] / Erq[1] ]
         endif
      endif
      if (nunq eq 3) then begin
         if (Radii[I] Le Erq[0]) then begin
            Nrat[*,I] = [1.,                      0,      0     ]
            Nmr[*,I]  = [Mrq[0]*Radii[I]/Erq[0] , Mrq[1], Mrq[2]]

            type_re = total(Nrat[*,I] * Nmr[*,I]^3 * f) / $
                      total(Nrat[*,I] * Nmr[*,I]^2 * g)
            if keyword_set(verbose) then $
               print,'(a) Re: ',strtrim(Radii[I],2), $
                     ', New effective radius: ',strtrim(type_re,2)
         endif
         if ((Radii[I] Gt Erq[0]) And (Radii[I] Lt Ea)) then begin ; Case B
            den1 = (Ratq[2]/Ratq[1]) * Nmr[2,I]^2 * (Radii[I]*g[2]- Nmr[2,I]*f[2])
            den2 =   Nmr[1,I]^2 * (Radii[I]*g[1]- Nmr[1,I]*f[1])
            den3 = -(Nmr[0,I]^2 * (Radii[I]*g[0]- Nmr[0,I]*f[0])) * $
                    (Ratq[2]/Ratq[1]+1)

            NRat[1,I] = (-Nmr[0,I]^2 *(Radii[I]*g[0]- Nmr[0,I]*f[0])) / $
                        (den1+den2+den3)

            NRat[0,I] = 1- NRat[1,I] * ((Ratq[2]/Ratq[1])+1)
            NRat[2,I] = NRat[1,I] * (Ratq[2]/Ratq[1])

            type_re = total(Nrat[*,I] * Nmr[*,I]^3 * f) / $
                        total(Nrat[*,I] * Nmr[*,I]^2 * g)
            if keyword_set(verbose) then $
               print,'(b) Re: ',strtrim(Radii[I],2), $
                     ', New effective radius: ',strtrim(type_re,2)
         endif
         if ((Radii[I] Gt Ea) And (Radii[I] Lt Erq[2])) then begin ; Case A
            den1 = (Ratq[0]/Ratq[1]) * Nmr[0,I]^2 * (Radii[I]*g[0]- Nmr[0,I]*f[0])
            den2 =   Nmr[1,I]^2 * (Radii[I]*g[1]- Nmr[1,I]*f[1])
            den3 = -(Nmr[2,I]^2 * (Radii[I]*g[2]- Nmr[2,I]*f[2])) * $
                    (Ratq[0]/Ratq[1]+1)

            NRat[1,I] = (-Nmr[2,I]^2*(Radii[I]*g[2]- Nmr[2,I]*f[2])) / $
                        (den1+den2+den3)

            NRat[0,I] = NRat[1,I] * (Ratq[0]/Ratq[1])
            NRat[2,I] = 1- NRat[1,I] * ((Ratq[0]/Ratq[1])+1)

            type_re = total(Nrat[*,I] * Nmr[*,I]^3 * f) / $
                      total(Nrat[*,I] * Nmr[*,I]^2 * g)
            if keyword_set(verbose) then $
               print,'(c) Re: ',strtrim(Radii[I],2), $
                     ', New effective radius: ',strtrim(type_re,2)
         endif
         if (Radii[I] Ge Erq[2]) then Begin
            Nrat[*,I] = [0.    , 0    , 1]
            Nmr[*,I]  = [Mrq[0], Mrq[1], Mrq[2] * Radii[I] / Erq[2]]

            type_re = total(Nrat[*,I] * Nmr[*,I]^3 * f) / $
                      total(Nrat[*,I] * Nmr[*,I]^2 * g)
            if keyword_set(verbose) then $
               print,'(d) Re: ',strtrim(Radii[I],2), $
                     ', New effective radius: ',strtrim(type_re,2)
         endif
      endif
      if (nunq eq 4) then begin
         if (Radii[I] Le Erq[0]) then begin
             Nrat[*,I] = [1., 0.,0.,0.]
             Nmr[*,I]  = [Mrq[0]*Radii[I]/Erq[0] , Mrq[1],Mrq[2], Mrq[3]]
         endif
         if ((Radii[I] Gt Erq[0]) And (Radii[I] Lt Ea)) then begin
            den1 = Nmr[3,I]^2 * (Radii[I]*g[3]- Nmr[3,I]*f[3])
            den2 = (Ratq[1]/Ratq[3]) *  Nmr[1,I]^2 * (Radii[I]*g[1]- Nmr[1,I]*f[1])
            den3 = (Ratq[2]/Ratq[3]) *  Nmr[2,I]^2 * (Radii[I]*g[2]- Nmr[2,I]*f[2])
            den4 = -(1+(Ratq[1]/Ratq[3])+(Ratq[2]/Ratq[3])) * $
                    Nmr[0,I]^2 * (Radii[I]*g[0]- Nmr[0,I]*f[0])

            NRat[3,I] = -(Nmr[0,I]^2 * (Radii[I]*g[0]- Nmr[0,I]*f[0])) / $
                         (den1+den2+den3+den4)

            NRat[1,I] = (Ratq[1]/Ratq[3]) * NRat[3,I]
            NRat[2,I] = (Ratq[2]/Ratq[3]) * NRat[3,I]
            NRat[0,I] = 1- NRat[3,I]- NRat[2,I]-  NRat[1,I]
         endif
         if ((Radii[I] Gt Ea) And (Radii[I] Lt Erq[3])) then begin
            den1 = Nmr[0,I]^2 * (Radii[I]*g[0]- Nmr[0,I]*f[0])
            den2 = (Ratq[1]/Ratq[0]) *  Nmr[1,I]^2 * (Radii[I]*g[1]- Nmr[1,I]*f[1])
            den3 = (Ratq[2]/Ratq[0]) *  Nmr[2,I]^2 * (Radii[I]*g[2]- Nmr[2,I]*f[2])
            den4 = -(1+(Ratq[1]/Ratq[0])+(Ratq[2]/Ratq[0])) * $
                    Nmr[3,I]^2 * (Radii[I]*g[3]- Nmr[3,I]*f[3])

            NRat[0,I] = -(Nmr[3,I]^2 * (Radii[I]*g[3]- Nmr[3,I]*f[3])) / $
                         (den1+den2+den3+den4)

            NRat[1,I] = (Ratq[1]/Ratq[0]) * NRat[0,I]
            NRat[2,I] = (Ratq[2]/Ratq[0]) * NRat[0,I]
            NRat[3,I] = 1- NRat[0,I]- NRat[2,I]- NRat[1,I]
         endif
         if (Radii[I] Ge Erq[3]) then Begin
            Nrat[*,I] = [0., 0., 0., 1]
            Nmr[*,I]  = [Mrq[0], Mrq[1],Mrq[2], Mrq[3] * Radii[I] / Erq[3] ]
         endif
      endif
   endfor

;  Now, if we have multiple components in each size mode, then we need to expand
;  the output arrays to all components.
   if (nunq ne n_elements(Rat)) then begin
      NRat_o = Fltarr(n_elements(Rat),N_elements(Radii))
      NMR_o  = Fltarr(n_elements(Rat),N_elements(Radii))
      for j=0,unqEr[0] do begin
;        Reapply the mixing ratio for the components within this size mode
         NRat_o[srtEr[j],*] = NRat[0,*] * Mode_Ratq[j,0]
         NMR_o [srtEr[j],*] = NMR[0,*]
      endfor
      for i=1,nunq-1 do begin
         for j=unqEr[i-1]+1,unqEr[i] do begin
            NRat_o[srtEr[j],*] = NRat[i,*] * Mode_Ratq[j-unqEr[i-1]-1,i]
            NMR_o [srtEr[j],*] = NMR[i,*]
         endfor
      endfor
   endif else begin
      NRat_o = NRat
      NMR_o  = NMR
   endelse

   if keyword_set(verbose) then message,/info, 'Returning'
end


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

pro create_bwgp, distname, Rm, S, RI, wl, Dqv, Bext, w, g, Phi, $
                 scode=scode, tmatrix_path=tmatrix_path, eps=eps, $
                 neps=neps, Vavg=Vavg

  if n_elements(RI) ne n_elements(wl) then $
     message, 'create_bwgp: Array size mismatch!'

   Inp = n_elements(Dqv)
   Inw = n_elements(wl)
   wn = 1.0 / wl
   Bext = fltarr(Inw)
   w = Bext
   g = Bext
   Phi = fltarr(Inp,Inw)

;  Check if we're using T-Matrix calculations
   if n_elements(scode) gt 0 then if strlowcase(scode) eq 'tmatrix' then $
      dotmatrix = 1 else dotmatrix = 0

   for i = 0, Inw-1 do begin ; LOOP OVER WAVELENGTH OR CHANNEL

;     Check if we're using T-Matrix or Mie scattering...
;     Note that we don't bother with non-sphericity for thermal-only channels
;     because:
;     (a) Scattering should be of limited importance
;     (b) Dubovik's tabulation doesn't cover the range of likely refractive
;         indices for such situations
      if dotmatrix and (wl(i) lt 6) then begin
;        Check that the path to the T-Matrix LUTs has been passed.
         if (not keyword_set(tmatrix_path)) then $
            message, 'Must specifiy the path to the Dubovik T-Matrix LUT ' + $
                     'base directory (keyword tmatrix_path) when T-Matrix ' + $'
                     'calculations are required'

;        Check that the user has passed the asymmetry information required for
;        the Dubovik code.
         if (n_elements(eps) eq 0) or (n_elements(neps) eq 0) then $
            message, 'Asymmetry parameters, "eps",  and their relative ' + $
                     'numbers, "neps", are required to use T-Matrix calculations'

;        TEMPORARY TESTING
;         eps[*]  = 1.0
;         neps[*] = 1.0 / float(n_elements(neps))
;        ******

;        Check for the absorption coefficient to make sure it's within the
;        Dubovik range (-0.5 to -0.0005). if it's outside, this range then:
;        * If it is greater than -0.0005, then print a warning message, but
;          assign it to the extreme value.
;          This can often happen with essentially non-absorbing particles, so
;          the change from -1e-6 to -5e-4 will basically have no effect on the
;          result.
;        * If it is less than -0.5, issue an error message and stop.
         if imaginary(RI(i)) gt -5e-4 then begin
            print,''
            message, /info, 'Warning: Imaginary RI altered from '+ $
                     strtrim(imaginary(RI[i]),2)+ $
                     ' to -5e-04, so as to lie within the Dubovik LUT range.'
            RItmp = complex(float(RI[i]),-5e-4)
         endif else if imaginary(RI[i]) lt -0.5 then begin
            print,''
            message, 'Imaginary RI is too negative for the Dubovik LUT range:'+ $
                     'k = '+strtrim(imaginary(RI[i]),2)+'; min value = -0.5'
         endif else RItmp = RI[i]

         dubovik_lognormal_multiple_eps, tmatrix_path, 1.0, Rm, S, wn[i], $
                                         RItmp, eps, neps, Dqv=Dqv, Bexttmp, $
                                         Bscatmp, wtmp, gtmp, ph, /no_mie, $
                                         /renorm_ph, /silent
         Phi[*,i] = ph
      endif else begin
;        We're not using T-Matrix, so call the normal Mie scattering code
         mie_size_dist, distname, 1.0, [Rm, S, 0.001, 100.0], wn[i], RI[i], $
                        Dqv=Dqv, /dlm, xres=0.4, Bexttmp, Bscatmp, wtmp, $
                        gtmp, SPM, Vavg=Vavg
         Phi[*,i] = SPM[0,*]
      endelse

      Bext[i] = Bexttmp
      w[i] = wtmp
      g[i] = gtmp

   endfor
end
