; Copyright (C) 1998-2017 University of Oxford
;
; This source code is licensed under the GNU General Public License (GPL),
; Version 3.  See the file COPYING for more details.


;+
; NAME:
;     mie_uoc_d
;
; PURPOSE:
;     Calculates the scattering parameters of a series of particles using the
;     Mie scattering theory.
;
; CATEGORY:
;     EODG Mie routines
;
; CALLING SEQUENCE:
;     mie_uoc_d, Dx, Cm, Dqv [, Inp=Inp] Xs1, Xs2, Dqxt, Dqsc, Dqbk, Dg, Dph
;
; INPUTS:
;     Dx:   Particle size parameter(s)
;     Cm:   Complex refractive index of the particle(s)
;     Dqv:  Cosines of scattering angles at which to compute the intensity
;           functions etc.
;
; KEYWORD PARAMETERS:
;     Inp:  Number of scattering angles at which to calculate intensity
;           functions etc (default = 1)
;
; OUTPUTS:
;     Xs1:  First amplitude function - amplitude of light polarised in the plane
;           perpendicular to the directions of incident light propagation and
;           observation.  Xs1 is a complex array of the same dimension as Dqv.
;     Xs2:  Second amplitude function - amplitude of light polarised in the
;           plane parallel to the directions of incident light propagation and
;           observation. Xs2 is a complex array of the same dimension as Dqv.
;     Dqxt: Extinction efficiency
;     Dqsc: Scattering efficiency
;     Dqbk: Backscattering efficiency
;     Dg:   Asymmetry parameter
;     Dph:  Phase function - an array of the same dimension as Dqv.
;
; OPTIONAL OUTPUTS:
;
; KEYWORD OUTPUTS:
;
; RESTRICTIONS:
;
; MODIFICATION HISTORY:
;     G. Thomas, 1998: mie_uoc.pro (translation of mieint.f to IDL)
;     R. Grainger, 2001: mie_uoc_d.pro (Added support for arrays of particle
;        sizes and included calculation of phase function)
;     G. Thomas, Sept 2003: (Put into EODG routines format)
;-

pro mie_uoc_d, Dx, Cm, Dqv, Inp=Inp, Xs1, Xs2, Dqxt, Dqsc, Dqbk, Dg, Dph

  Imaxx = 12000
  RIMax = 2.5
  Itermax = Imaxx * RIMax
  Imaxnp = 1100 ; Change this as required
  Sizes = n_elements(Dx)
  if (keyword_set(Inp) then begin
    Inp = 1
    Dqv = 0
  endif

  Xs1 = complexarr(Inp,Sizes)
  Xs2 = complexarr(Inp,Sizes)
  Dqxt = dblarr(Sizes)
  Dqsc = dblarr(Sizes)
  Dqbk = dblarr(Sizes)
  Dg = dblarr(Sizes)
  Dph = dblarr(Inp,Sizes)

  for Size = 0, Sizes - 1 do begin

    if (Dx(Size) gt Imaxx) then message, 'Error: Size Parameter Overflow in Mie'
    Ir = 1.D0 / Cm
    Y =  Dx(Size) * Cm

    if (Dx(Size) lt 0.02) then NStop = 2 else begin
      case 1 OF
        (Dx(Size) le 8.0)    : NStop = Dx(Size) + 4.00*Dx(Size)^(1./3.) + 2.0
        (Dx(Size) lt 4200.0) : NStop = Dx(Size) + 4.05*Dx(Size)^(1./3.) + 2.0
        else                 : NStop = Dx(Size) + 4.00*Dx(Size)^(1./3.) + 2.0
      endcase
    end
    NmX = fix(max([NStop,abs(Y)]) + 15.)
    D = dcomplexarr(Nmx+1)

    for N = Nmx-1,1,-1 do begin
      A1 = (N+1) / Y
      D(N) = A1 - 1/(A1+D(N+1))
    end

    Sm = dcomplexarr(Inp)
    Sp = dcomplexarr(Inp)
    Pi0 = dcomplexarr(Inp)
    Pi1 = dcomplex(replicate(1.D0,Inp),replicate(0.D0,Inp))

    Psi0 = cos(Dx(Size))
    Psi1 = sin(Dx(Size))
    Chi0 =-sin(Dx(Size))
    Chi1 = cos(Dx(Size))
    Xi0 = dcomplex(Psi0,Chi0)
    Xi1 = dcomplex(Psi1,Chi1)

    Dg(Size) = 0.D0
    Dqsc(Size) = 0.D0
    Dqxt(Size) = 0.D0
    Tnp1 = 1

    for N = 1,Nstop do begin
      DN = double(N)
      Tnp1 = Tnp1 + 2
      Tnm1 = Tnp1 - 2
      A2 = Tnp1 / (DN*(DN+1.D0))
      Turbo = (DN+1.D0) / DN
      Rnx = DN/Dx(Size)
      Psi = double(Tnm1)*Psi1/Dx(Size) - Psi0
      Chi = Tnm1*Chi1/Dx(Size)       - Chi0
      Xi = dcomplex(Psi,Chi)
      A = ((D[N]*Ir+Rnx)*Psi-Psi1) / ((D[N]*Ir+Rnx)*  Xi-  Xi1)
      B = ((D[N]*Cm+Rnx)*Psi-Psi1) / ((D[N]*Cm+Rnx)*  Xi-  Xi1)
      Dqxt(Size) = Tnp1 *      double(A + B)          + Dqxt(Size)
      Dqsc(Size) = Tnp1 * double(A*conj(A) + B*conj(B)) + Dqsc(Size)
      if (N gt 1) then Dg(Size) = Dg(Size) $
                           + (dN*dN - 1) * double(ANM1*conj(A) + BNM1 * conj(B)) / dN $
                           + TNM1 * double(ANM1*conj(BNM1)) / (dN*dN - dN)
      Anm1 = A
      Bnm1 = B

      S = Dqv * Pi1
      T = S - Pi0
      Taun = N*T - Pi0
      Sp = (A2 * (A + B)) * (Pi1 + Taun) + Sp
      Sm = (A2 * (A - B)) * (Pi1 - Taun) + Sm
      Pi0 = Pi1
      Pi1 = S + T*Turbo

      Psi0 = Psi1
      Psi1 = Psi
      Chi0 = Chi1
      Chi1 = Chi
      Xi1 = dcomplex(Psi1,Chi1)

    end ; for Nstop

    if (Dg(Size) gt 0) then Dg(Size) = 2 * Dg(Size) / Dqsc(Size)
    Dqsc(Size) =  2 * Dqsc(Size) / Dx(Size)^2
    Dqxt(Size) =  2 * Dqxt(Size) / Dx(Size)^2

    Xs1(*,Size) = (Sp + Sm) / 2
    Xs2(*,Size) = (Sp - Sm) / 2
    Dph(*,Size) = 2 * double(Xs1(*,Size)*conj(Xs1(*,Size)) + $
                             Xs2(*,Size)*conj(Xs2(*,Size))) / (Dx(Size)^2 * Dqsc(Size))
    Dqbk(Size) =  4 * abs(Xs1(Inp-1)^2) / Dx(Size)^2
  end

  return
end
