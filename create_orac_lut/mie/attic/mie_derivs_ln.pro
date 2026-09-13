; Copyright (C) 1998-2017 University of Oxford
;
; This source code is licensed under the GNU General Public License (GPL),
; Version 3.  See the file COPYING for more details.


;+
; NAME:
;     mie_derivs_ln
;
; PURPOSE:
;     Calculates the scattering parameters and their analytical derivatives
;     (wrt the parameters of the distribution) of a log normal distribution of
;     spherical particles.
;
;     A the derivation of expressions for the analytical derivatives of Mie
;     scattering terms is covered by:
;     Grainger, R.G., J. Lucas, G.E. Thomas, G. Ewan, "The Calculation of Mie
;     Derivatives", Submitted to Appl. Opt., 2004.
;
; CATEGORY:
;     EODG Mie routines
;
; CALLING SEQUENCE:
;     mie_derivs_ln, N, Rm, S, wavenumber, Cm [, Dqv=Dqv], Bext, Bsca, $
;     dBextdN, dBextdRm, dBextdS, dBscadN, dBscadRm, dBscadS, [, i1] [, i2] $
;     [, di1dN] [, di1dRm] [, di1dS] [, di2dN] [, di2dRm] [, di2dS]
;
; INPUTS:
;     N:          Number density of the particle distribution
;     Rm:         Median radius of the particle distribution (microns)
;     S:          Spread of the distribution, such that the standard deviation
;                 of ln(r) is ln(S)
;     wavenumber: Wavenumber of light, defined as 1/wavelength. A positive
;                 scalar whose units match those of Rm
;     Cm:         Complex refractive index of the particle(s)
;
; KEYWORD PARAMETERS:
;     Dqv:        Cosines of scattering angles at which to compute the intensity
;                 functions etc
;
; OUTPUTS:
;     Bext:       Volume extinction coefficient
;     Bsca:       Volume scattering coefficient
;     dBextdN:    Derivative of the extinction coefficient wrt the number density
;     dBextdRm:   Derivative of the extinction coefficient wrt the mean radius
;     dBextdS:    Derivative of the extinction coefficient wrt the spread
;     dBscadN:    Derivative of the scattering coefficient wrt the number density
;     dBscadRm:   Derivative of the scattering coefficient wrt the mean radius
;     dBscadS:    Derivative of the scattering coefficient wrt the
;                 spread
;     i1:         First intensity function - intensity of light polarised in
;                 the plane perpendicular to the directions
;                 of incident light propagation and observation. Only
;                 calculated if Dqv is specified.
;     i2:         Second intensity function - intensity of light polarised in
;                 the plane parallel to the directions of
;                 incident light propagation and observation. Only
;                 calculated if Dqv is specified.
;     di1dN:      Derivatives of the first intensity function wrt the number density
;     di1dRm:     Derivatives of the first intensity function wrt the mean radius
;     di1dS:      Derivatives of the first intensity function wrt the spread
;     di2dN:      Derivatives of the second intensity function wrt the number density
;     di2dRm:     Derivatives of the second intensity function wrt the mean radius
;     di2dS:      Derivatives of the second intensity function wrt the spread
;
;     NB. The values of involving the intensity functions are arrays of the same
;     dimension as Dqv and are only calculated if Dqv is specified.
;
; KEYWORD OUTPUTS:
;
; RESTRICTIONS:
;     Note, this procedure calls the mie_single and quadrature procedures.
;
; MODIFICATION HISTORY:
;     G. Thomas, Sep 2003: mie_derivs_ln.pro
;     G. Thomas, Nov 2003: minor bug fixes
;     G. Thomas, Feb 2004: Explicit double precision added throughout and header
;         updated quadrature. Also added npts and info keywords.
;-

pro mie_derivs_ln, N, Rm, S, wavenumber, Cm, Dqv=Dqv, Bext, Bsca, $
                   dBextdN, dBextdRm, dBextdS, dBscadN, dBscadRm, dBscadS, $
                   i1, i2, di1dN, di1dRm, di1dS, di2dN, di2dRm, di2dS

    Common miedervln, absc, wght

;   Create vectors for size integration
    Tq = gauss_cvf(0.999D0)
    Rl = exp(alog(Rm)+Tq*alog(S))  & Ru = exp(alog(Rm)-Tq*alog(S)+alog(4))
    if 2D0 * !dpi * Ru * wavenumber ge 12000 then begin
        Ru = 11999D0 / ( 2D0 * !dpi * wavenumber )
        message,/continue,'Warning:  Radius upper bound truncated to avoid '+ $
              'size parameter overflow.'
    endif

    if not keyword_set(npts) then begin
;      Accurate calculation requires 0.1 step size in x
       Npts = (long(2D0 * !dpi * (ru-rl) * wavenumber/0.1)) > 200
    endif

;   quadrature on the radii
    if n_elements(wght) ne Npts then quadrature,'T',Npts,absc,wght

    shift_quadrature,absc,wght,Rl,Ru,R,Wghtr

    Tmp = exp(-0.5D0*(alog(R/Rm) / alog(S))^2) / (sqrt(2D0*!dpi) * alog(S) * R)

    W1 = N * Tmp

    Dx = 2D0 * !dpi * R * wavenumber
    if keyword_set(info) then info = { Npts    : Npts, $
                                       MinSize : Dx[0], $
                                       MaxSize : Dx[Npts-1] }

;   Create Mie variables
    Dqxt=dblarr(npts)
    Dqsc=dblarr(npts)

    Dx = 2d0 * !dpi * R * wavenumber
;   If an array of cos(theta) is provided, calculate phase function
    if n_elements(Dqv) gt 0 then begin
        Inp = n_elements(Dqv)
        i1 = dblarr(Inp) & i2 = i1
        di1dN = i1 &  di2dN = i1
        di1dRm = i1 & di2dRm = i1
        di1dS = i1 & di2dS = i1
        Mie_single, Dx, Cm, Dqv=Dqv, Dqxt, Dqsc, Dqbk, Dg, Xs1, Xs2, Dph
    endif else begin
        inp = 1
        Dqv = 1d0
        Mie_single,  Dx, Cm, Dqxt, Dqsc, Dqbk, Dg
    endelse

    lnRRm = alog(R/Rm) ;Precalculate for speed

    Bext = Total(wghtr * W1 * DQxt * !dpi * R^2) ; Extinction
    dBextdN = Bext / N
    dBextdRm = Total(wghtr * W1 * DQxt * lnRRm * !dpi * R^2) / (alog(S)^2 * Rm)
    dBextdS = (Total(wghtr * W1 * DQxt * lnRRm^2 * !dpi * R^2) / alog(S)^2 $
             - Total(wghtr * W1 * DQxt * !dpi * R^2)) / (S * alog(S))

    Bsca = Total(wghtr * W1 * DQsc * !dpi * R^2) ;Scattering
    dBscadN = Total(wghtr * W1 * DQsc * !dpi * R^2) / N
    dBscadRm = Total(wghtr * W1 * DQsc * lnRRm * !dpi * R^2) / (alog(S)^2 * Rm)
    dBscadS = Total(wghtr * W1 * DQsc * lnRRm^2 * !dpi * R^2) / (S * alog(S)^3) $
            - Total(wghtr * W1 * DQsc * !dpi * R^2) / (S * alog(S))

    if n_elements(Dqv) gt 0 then $ ; Intensity functions
        for i =0,Inp-1 do begin
            i1(i) = Total(wghtr * W1 * real_part(Xs1(i,*)*conj(Xs1(i,*))))
            di1dN(i) = i1(i) / N
            di1dRm(i) = Total(wghtr * W1 * lnRRm * real_part(Xs1(i,*)*conj(Xs1(i,*)))) $
                        / (alog(S)^2 * Rm)
            di1dS(i) = Total(wghtr * W1 * lnRRm * real_part(Xs1(i,*)*conj(Xs1(i,*)))) $
                       / (S * alog(S)^3) $
                       - Total(wghtr * W1 * real_part(Xs1(i,*)*conj(Xs1(i,*)))) $
                       / (S * alog(S))

            i2(i) = Total(wghtr * W1 * real_part(Xs2(i,*)*conj(Xs2(i,*))))
            di2dN(i) = i2(i) / N
            di2dRm(i) = Total(wghtr * W1 * lnRRm * real_part(Xs2(i,*)*conj(Xs2(i,*)))) $
                        / (alog(S)^2 * Rm)
            di2dS(i) = Total(wghtr * W1 * lnRRm * real_part(Xs2(i,*)*conj(Xs2(i,*)))) $
                       / (S * alog(S)^3) $
                       - Total(wghtr * W1 * real_part(Xs2(i,*)*conj(Xs2(i,*)))) $
                       / (S * alog(S))
    endfor
end
