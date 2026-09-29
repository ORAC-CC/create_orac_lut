; Copyright (C) 1998-2017 University of Oxford
;
; This source code is licensed under the GNU General Public License (GPL),
; Version 3.  See the file COPYING for more details.


pro scat_phase_matrix_exp, n_l, qx, qw, F, n_l2, coefs
;+
; NAME:
;     scat_phase_matrix_exp
;
; PURPOSE:
;     Expand the given single scattering phase matrix in terms of generalized
;     spherical functions in the form required by radiative transfer solvers
;     such as VLIDORT or XRTM.
;
; CATEGORY:
;     EODG Mie routines
;
; CALLING SEQUENCE:
;
; INPUTS:
;     n_l: Maximum order to compute.
;     qx:  Array of quadrature points (scattering angle cosines) for the
;          integration of the input scattering phase matrix F.
;     qw:  Array of quadrature weights (sum(qw) = unity) for the integration of
;          the input scattering phase matrix F. n_elements(qx) = n_elements(qx).
;     F:   2D array (6, n_elements(qx)) of the elements of the single scattering
;          phase matrix for each quadrature point where
;          F[1,*] = F11[*]
;          F[2,*] = F22[*]
;          F[3,*] = F33[*]
;          F[4,*] = F44[*]
;          F[5,*] = F12[*]
;          F[6,*] = F34[*]
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;     coefs: 2D array (6, n_l) of expansions coefficients.
;
; KEYWORD OUTPUTS:
;
; RESTRICTIONS:

; MODIFICATION HISTORY:
;     G. McGarragh, 28 Jul 2015: Taken from LMie and translated from C.
;-

    n_q = n_elements(qx)

    b1 = dblarr(6)
    b2 = dblarr(6)

    P = dblarr(4, n_l)

    coefs = dblarr(6, n_l)
    coefs[*,*] = 0.d

    for i = 0, n_q - 1 do begin
        gen_sph_funcs, n_l, qx[i], P0, P1, P2, P3

        b1[*] = F[*,i] * qw[i]
        a = b1[1]
        b1[1] = a + b1[2]
        b1[2] = a - b1[2]

        coefs[0,*] += P0 * b1[0]
        coefs[1,*] += P2 * b1[1]
        coefs[2,*] += P3 * b1[2]
        coefs[3,*] += P0 * b1[3]
        coefs[4,*] += P1 * b1[4]
        coefs[5,*] += P1 * b1[5]
    endfor

    for j = 0, n_l - 1 do begin
        a = coefs[1,j]
        coefs[1,j] = (a + coefs[2,j]) / 2.d
        coefs[2,j] = (a - coefs[2,j]) / 2.d
    endfor

    for i = 0, n_l - 1 do begin
        a = (2 * i + 1) / 2.d
        coefs[*,i] *= a
    endfor
end
