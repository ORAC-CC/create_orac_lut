pro generate_scattering_properties,srfstrarr,scatoffset, nwvl_max,  inststr, mmstr, lutstr, nmom, bext550, w550, g550, phs550, amom550, bextrat, bext, w, g, vavg, phs, amom,tmatrix_path=tmatrix_path,no_screen=no_screen

;  IF the no_screen keyword has been set, we suppress the control character
;  used to prevent a new-line FOR some print statements (see below). This is
;  useful IF the output is being piped into a file, FOR example (WHERE the
;  control character just messes up the formatting).
   IF keyword_set(no_screen) then begin
      newlinechar = ' '
      newlinestr  = ')'
   ENDIF ELSE begin
      newlinechar = string(13b)
      newlinestr  = ',$)'
   endelse


;  Interpolate the components refractive index values onto the 0.55 micron
;  reference wavelength and the instrument channels.
   AerM550 = complexarr(mmstr.NComp)
   AerM = complexarr(Nwvl_Max, inststr.Number_of_nadir_Channels, mmstr.NComp)

;  Which member of the mmstr structure is the first component with refractive
;  index/asymmetry information substructure?
   scatoffset = WHERE(tag_names(mmstr) eq 'COMP1')
   scatoffset = scatoffset[0]
   IF scatoffset eq -1 then message, "Can't find COMP1 substructure in mmstr"

   IF mmstr.(scatoffset).code ne '' then begin
      FOR c=0,mmstr.NComp-1 do begin
         cc = scatoffset + c
         AerM550[c]  = interpol(mmstr.(cc).Cm, mmstr.(cc).wl, 0.55)
         AerM[*,*,c] = interpol(mmstr.(cc).Cm, mmstr.(cc).wl, srfstrarr[*].wvl[*])
;         print, 0.55, srfstrarr[*].wvl[*]
;         print, AerM550[c],AerM[*,*,c]
      ENDFOR

;     IF the force_n keyword has been specified, replace real RI values with
;     those given in force_n.
      IF N_ELEMENTS(force_n) gt 0 then begin
         IF force_n[0] ne ' ' then begin
            FOR c=0,mmstr.NComp-1 do begin
               AerM550[c] = AerM550[c] - complex(float(AerM550[c]), 0.0) + $
                            complex(float(force_n[0]), 0.0)
            ENDFOR
         ENDIF

         FOR i=0,inststr.Number_of_nadir_Channels-1 do begin
            IF force_n[i+1] ne ' ' then begin
               FOR c=0,mmstr.NComp-1 do begin
                  AerM[*,i,c] = AerM[*,i,c] - complex(float(AerM[*,i,c]), 0.0) + $
                                complex(float(force_n[i+1]), 0.0)
               ENDFOR
            ENDIF
         ENDFOR
      ENDIF

;     IF the force_k keyword has been specified, replace imaginary RI values
;     with those given in force_k.
      IF N_ELEMENTS(force_k) gt 0 then begin
         IF force_k[0] ne ' ' then begin
            FOR c=0,mmstr.NComp-1 do begin
               AerM550[c] = AerM550[c] - complex(0.0,imaginary(AerM550[c])) + $
                            complex(0.0,float(force_k[0]))
            ENDFOR
         ENDIF

         FOR i=0,inststr.Number_of_nadir_Channels-1 do begin
            IF force_k[i+1] ne ' ' then begin
               FOR c=0,mmstr.NComp-1 do begin
                  AerM[*,i,c] = AerM[*,i,c] - complex(0.0,imaginary(AerM[*,i,c])) + $
                                complex(0.0,float(force_k[i+1]))
               ENDFOR
            ENDIF
         ENDFOR
      ENDIF
   ENDIF


;   Generate the range of component mixing-ratios/mode radii required to provide the required effective radii
    PRINT, mmstr.comptype[0]
    IF mmstr.comptype[0] EQ 'opac' OR mmstr.comptype[0] EQ 'user' THEN BEGIN
      PRINT, mmstr.distname[0]
      CASE mmstr.distname[0] OF
         'log_normal': create_range, mmstr.mrat, mmstr.rm, mmstr.s, lutstr.efr, lut_mrat, lut_rm ; create_range returns (ncomp, nefr) arrays of mixing ratio and rm
     'modified_gamma': BEGIN
                         lut_MRat = DBLARR(N_ELEMENTS(mmstr.mrat),N_ELEMENTS(lutstr.efr))
                         lut_Rm   = DBLARR(N_ELEMENTS(mmstr.mrat),N_ELEMENTS(lutstr.efr))
                         lut_MRat[0,*] = mmstr.mrat
                         lut_Rm  [0,*] = lutstr.efr   
                       END
      ENDCASE
    ENDIF ELSE BEGIN
         lut_MRat = DBLARR(N_ELEMENTS(mmstr.MRat),N_ELEMENTS(lutstr.EfR))
         lut_MRat[0,*] = mmstr.MRat
    ENDELSE



;     **** Generate the quadrature points FOR the scattering phase function

;     Check IF the NMom keyword has been set, IF it hasn't we use the default
;     value of 1000.
      IF N_ELEMENTS(n_theta) eq 0 then begin
         NMom = 1000
;        x = 2. * !pi * 240. / .47;
;        NMom = fix(2 * (x + 4.05 * x^(1./3.) + 8))
      ENDIF ELSE begin
         NMom = n_theta
      endelse

;     The quadrature procedure gives us our phase function angles
      quadrature, 'g', NMom, Abscissas, Weights
;     Note that QV = cos(scattering_angle)
      QV0 =  1.0 ; theta = 0
      QV1 = -1.0 ; theta = 180
      QV = ((QV1-QV0)*Abscissas + (QV0+QV1)) / 2d0
      PTheta = acos(QV)

;     **** call the scattering code for the required range of components and mode radii

;     define the arrays which hold the scattering parameters
      vavg_c    = FLTARR(mmstr.ncomp, lutstr.efr_n) ; average volume per particle

      ; at the reference wavelength
      bext550_c = FLTARR(mmstr.ncomp, lutstr.efr_n) ; extinction coefficient
      w550_c    = FLTARR(mmstr.ncomp, lutstr.efr_n) ; single scattering albedo
      g550_c    = FLTARR(mmstr.ncomp, lutstr.efr_n) ; asymmetry parameter
      phs550_c  = FLTARR(nmom, mmstr.ncomp, lutstr.efr_n) ; phase function

      ; at the individual channels
      bext_c    = FLTARR(nwvl_max, inststr.Number_of_nadir_Channels, mmstr.ncomp, lutstr.efr_n)
      w_c       = FLTARR(nwvl_max, inststr.Number_of_nadir_Channels, mmstr.ncomp, lutstr.efr_n)
      g_c       = FLTARR(nwvl_max, inststr.Number_of_nadir_Channels, mmstr.ncomp, lutstr.efr_n)
      phs_c     = FLTARR(nmom, nwvl_max, inststr.Number_of_nadir_Channels, mmstr.ncomp, lutstr.efr_n)

      FOR c=0,mmstr.NComp-1 do begin
         IF mmstr.comptype[0] eq 'opac' or mmstr.comptype[0] eq 'user' then begin
            FOR r=0,lutstr.efr_n-1 do begin
               print,' Performing scattering calculations for component ' + $
                     mmstr.compname[c] + ', EfR: ', lutstr.EfR[r];, newlinechar, $
                    ; format='(A,f8.4,A'+newlinestr

;              Note that the mode radii of each component don't usually change
;              from one effective radius to the next (it's the mixing ratio that
;              changes). Here we check IF the component has changed size and
;              only call the scattering code IF it has.
               calculated = 0
               IF r gt 0 then begin
                  IF (lut_Rm[c,r] eq lut_Rm[c,r-1]) and (Bext_c[0,0,c,r-1] ne 0) then calculated = 1
               ENDIF
               IF calculated then begin
                  Vavg_c[c,r]     = Vavg_c[c,r-1]
                  Bext550_c[c,r]  = Bext550_c[c,r-1]
                  w550_c[c,r]     = w550_c[c,r-1]
                  g550_c[c,r]     = g550_c[c,r-1]
                  Phs550_c[*,c,r] = Phs550_c[*,c,r-1]

                  Bext_c[*,*,c,r]  = Bext_c[*,*,c,r-1]
                  w_c[*,*,c,r]     = w_c[*,*,c,r-1]
                  g_c[*,*,c,r]     = g_c[*,*,c,r-1]
                  Phs_c[*,*,*,c,r] = Phs_c[*,*,*,c,r-1]
               ENDIF ELSE begin ; This is a new mode radius, do the calculation
                  IF lut_MRat[c,r] gt 0 then begin
                     cc = scatoffset + c

;                    IF the Mie keyword has been set, we use Mie scattering FOR
;                    all components, regardless of the driver settings.
                     IF keyword_set(mie) then scode = 'mie' $
                     ELSE scode = mmstr.(cc).code

;                    Set up the values of eps and neps, which only exist in the
;                    mmstr structure IF tmatrix scattering is to be used.
                     IF strlowcase(scode) eq 'tmatrix' then begin
                        epsvals = mmstr.(cc).eps
                        nepsvals = mmstr.(cc).neps
                     ENDIF ELSE begin
                        epsvals = 0
                        nepsvals = 0
                     endelse

;                    Calculate Bext at 550 nm (the reference wavelength) and  Vavg (average volume per particle).
                     Vavg1 = 0.
                     create_bwgp, mmstr.distname[c], lut_Rm[c,r], mmstr.S[c], AerM550[c], 0.55, QV, Bext1, w1, g1, Phs1, scode=scode, tmatrix_path=tmatrix_path, eps=epsvals, neps=nepsvals, Vavg=Vavg1
                     Vavg_c[c,r]     = Vavg1
                     Bext550_c[c,r]  = Bext1
                     w550_c[c,r]     = w1
                     g550_c[c,r]     = g1
                     Phs550_c[*,c,r] = Phs1

;                    calculate bext, w (single scatter albedo), g (asymmetry parameter) and phs (phase function) for each instrument channel.
                     create_bwgp, mmstr.distname[c], lut_Rm[c,r], mmstr.S[c], AerM[*,*,c], srfstrarr[*].wvl[*], QV, Bext1, w1, g1, Phs1, scode=scode, tmatrix_path=tmatrix_path, eps=epsvals, neps=nepsvals
                     Bext_c[*,*,c,r]  = reform(Bext1, Nwvl_Max, inststr.Number_of_nadir_Channels)
                     w_c[*,*,c,r]     = reform(w1,    Nwvl_Max, inststr.Number_of_nadir_Channels)
                     g_c[*,*,c,r]     = reform(g1,    Nwvl_Max, inststr.Number_of_nadir_Channels)
                     Phs_c[*,*,*,c,r] = reform(Phs1,  NMom, Nwvl_Max, inststr.Number_of_nadir_Channels)
                  ENDIF
               endelse
            ENDFOR

         ENDIF ELSE IF mmstr.comptype[c] eq 'baran' then begin
            init_baran, mmstr.compname[c], baran

            read_baran, baran, 0.55, lutstr.EfR, Bext1, w1, g1, $
                        PTheta * 180. / !pi, Phs1
            Bext550_c[c,*]  = Bext1
            w550_c[c,*]     = w1
            g550_c[c,*]     = g1
            Phs550_c[*,c,*] = Phs1

            ; The Baran optical property data set FOR thermal channels does not
            ; have g in which case it sets g to zero but g is used below *only*
            ; to distinguish particle layers from layers without particles so we
            ; set g to a nonzero value IF it is zero.
            g550_c[WHERE(g_c eq 0)] = -999.

            ; Normalize
            FOR r=0,lutstr.efr_n-1 do begin
               Phs550_c[*, c, r] /= total(Phs550_c[*, c, r] * weights) / 2.
            ENDFOR

            read_baran, baran, srfstrarr[*].wvl[*], lutstr.EfR, Bext1, w1, g1, $
                        PTheta * 180. / !pi, Phs1
            Bext_c[*,*,c,*]  = reform(Bext1, Nwvl_Max, inststr.Number_of_nadir_Channels, lutstr.efr_n)
            w_c[*,*,c,*]     = reform(w1,    Nwvl_Max, inststr.Number_of_nadir_Channels, lutstr.efr_n)
            g_c[*,*,c,*]     = reform(g1,    Nwvl_Max, inststr.Number_of_nadir_Channels, lutstr.efr_n)
            Phs_c[*,*,*,c,*] = reform(Phs1,  NMom, Nwvl_Max, inststr.Number_of_nadir_Channels, lutstr.efr_n)

            g_c[WHERE(g_c eq 0)] = -999.

            ; Normalize
            FOR r=0,lutstr.efr_n-1 do begin
               FOR l=0,inststr.Number_of_nadir_Channels-1 do begin
                  FOR m=0,Nwvl_Max-1 do begin
                     Phs_c[*, m, l, c, r] /= total(Phs_c[*, m, l, c, r] * weights) / 2.
                  ENDFOR
               ENDFOR
            ENDFOR

         ENDIF ELSE IF mmstr.comptype[c] eq 'baum' then begin
            read_baum_lambda, mmstr.compname[c], 0.55, lutstr.EfR, Bext1, $
                              w1, g1, PTheta * 180. / !pi, Phs1
            Bext550_c[c,*]  = Bext1
            w550_c[c,*]     = w1
            g550_c[c,*]     = g1
            Phs550_c[*,c,*] = Phs1

            ; Normalize
            FOR r=0,lutstr.efr_n-1 do begin
               Phs550_c[*, c, r] /= total(Phs550_c[*, c, r] * weights) / 2.
            ENDFOR

            ; Choose between spectral or instrument/channel specific properties
            IF mmstr.compname2[c] eq '' then begin
               read_baum_lambda, mmstr.compname[c], srfstrarr[*].wvl[*], lutstr.EfR, $
                                 Bext1, w1, g1, PTheta * 180. / !pi, Phs1
            ENDIF ELSE begin
               read_baum_channel, mmstr.compname2[c], fix(channels), lutstr.EfR, $
                                  Bext1, w1, g1, PTheta * 180. / !pi, Phs1
            endelse
            Bext_c[*,*,c,*]  = reform(Bext1, Nwvl_Max, inststr.Number_of_nadir_Channels, lutstr.efr_n)
            w_c[*,*,c,*]     = reform(w1,    Nwvl_Max, inststr.Number_of_nadir_Channels, lutstr.efr_n)
            g_c[*,*,c,*]     = reform(g1,    Nwvl_Max, inststr.Number_of_nadir_Channels, lutstr.efr_n)
            Phs_c[*,*,*,c,*] = reform(Phs1,  NMom, Nwvl_Max, inststr.Number_of_nadir_Channels, lutstr.efr_n)

            ; Normalize
            FOR r=0,lutstr.efr_n-1 do begin
               FOR l=0,inststr.Number_of_nadir_Channels-1 do begin
                  FOR m=0,Nwvl_Max-1 do begin                     ; ****** check should be nwvls for this channel?
                     Phs_c[*, m, l, c, r] /= total(Phs_c[*, m, l, c, r] * weights) / 2.
                  ENDFOR
               ENDFOR
            ENDFOR
         ENDIF
      ENDFOR

      print,''
      print,'All scattering calculations completed for each component'


;     **** Calculate the scattering parameters of the class for each required
;          effective radius

;     Define the arrays to hold the scattering parameters for the class as a
;     whole - these are the variables which are written out into the RT driver
;     files.

      Vavg    = FLTARR(lutstr.efr_n)

;     At the reference wavelength
      Bext550 = FLTARR(lutstr.efr_n)       ; Extinction coefficient
      w550    = FLTARR(lutstr.efr_n)       ; Single scattering albedo
      g550    = FLTARR(lutstr.efr_n)       ; Asymmetry parameter
      Phs550  = FLTARR(NMom, lutstr.efr_n) ; Phase function
      AMom550 = FLTARR(NMom, lutstr.efr_n) ; Legendre moments

;     FOR each channel
      BextRat = FLTARR(Nwvl_Max, inststr.Number_of_nadir_Channels, lutstr.efr_n) ; Ratio of Bext with that
;                                                          at the reference wavelength
      Bext    = FLTARR(Nwvl_Max, inststr.Number_of_nadir_Channels, lutstr.efr_n) ; Extinction coefficient
      w       = FLTARR(Nwvl_Max, inststr.Number_of_nadir_Channels, lutstr.efr_n) ; Single scattering albedo
      g       = FLTARR(Nwvl_Max, inststr.Number_of_nadir_Channels, lutstr.efr_n) ; Asymmetry parameter
      Phs     = FLTARR(NMom, Nwvl_Max, inststr.Number_of_nadir_Channels, lutstr.efr_n) ; Phase function
      AMom    = FLTARR(NMom, Nwvl_Max, inststr.Number_of_nadir_Channels, lutstr.efr_n) ; Legendre moments

      Vavg[*] = Vavg_c[0,*]

      FOR l=0,inststr.Number_of_nadir_Channels-1 do begin
         FOR m=0,Nwvl_Max-1 do begin                      ;*** CHECK SHOULD BE n wav for this chaennel
            FOR r=0,lutstr.efr_n-1 do begin
;              Calculate 550 nm extinction coefficient
               IF l eq 0 then begin
;                 Pre-calculate some factors that are used more than once
                  MratBext   = lut_MRat[*,r] * Bext550_c[*,r]
                  tMratBext  = total(MratBext)
                  MratBextw  = MratBext * w550_c[*,r]
                  tMratBextw = total(MratBextw)

                  Bext550[r] = tMratBext / total(lut_Mrat[*,r])
                  w550[r]    = tMratBextw / tMratBext
                  g550[r]    = total(MratBextw * g550_c[*,r]) / tMratBextw
                  FOR p=0,NMom-1 do $
                     Phs550[p,r] = total(MratBextw * Phs550_c[p,*,r]) / tMratBextw

;                 Calculate the Legendre moments FOR the phase function
                  legpexp, NMom, QV, weights, Phs550[*,r], Inlc, alc
                  AMom550[*,r] = alc / (2.0*findgen(NMom)+1.0)
               ENDIF

;              Pre-calculate some factors that are used more than once
               MratBext   = lut_MRat[*,r] * Bext_c[m,l,*,r]
               tMratBext  = total(MratBext)
               MratBextw  = MratBext * w_c[m,l,*,r]
               tMratBextw = total(MratBextw)

               Bext[m,l,r]  = tMratBext / total(lut_Mrat[*,r])
               w[m,l,r]     = tMratBextw / tMratBext
               g[m,l,r]     = total(MratBextw * g_c[m,l,*,r]) / tMratBextw
               FOR p=0,NMom-1 do $
                  Phs[p,m,l,r] = total(MratBextw * Phs_c[p,m,l,*,r]) / tMratBextw

;              Calculate the Legendre moments FOR the phase function
               legpexp, NMom, QV, weights, Phs[*,m,l,r], Inlc, alc
               AMom[*,m,l,r] = alc / (2.0*findgen(NMom)+1.0)

;              Calculate the ratio of the extinction coefficient at the current
;              channel and 0.55 microns, allowing the spectral optical depth to be
;              calculated from the reference 0.55 micron value.
               BextRat[m,l,r] = Bext[m,l,r]/Bext550[r]
            ENDFOR
         ENDFOR
      ENDFOR
end