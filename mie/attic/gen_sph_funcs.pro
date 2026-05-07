; Copyright (C) 1998-2017 University of Oxford
;
; This source code is licensed under the GNU General Public License (GPL),
; Version 3.  See the file COPYING for more details.


;+
; NAME:
;     gen_sph_funcs
;
; PURPOSE:
;     Computes generalized spherical functions required for the expantion of a
;     single scattering phase matrix in the form required by radiative transfer
;     solvers such as VLIDORT or XRTM.
;
; CATEGORY:
;     EODG Mie routines
;
; CALLING SEQUENCE:
;     gen_sph_funcs, n_l, mu, p00, p0p2, p2p2, p2m2
;
; INPUTS:
;     n_l: Maximum order to compute.
;     mu:  Array of scattering angle (0 - 180) cosines for which to compute
;          functions for.
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;     p00:  2D array (n_elements(mu), n_l) of P^l_{0,0} values where l is order.
;     p0p2: 2D array (n_elements(mu), n_l) of P^l_{0,2} values where l is order.
;     p2p2: 2D array (n_elements(mu), n_l) of P^l_{2,2} values where l is order.
;     p2m2: 2D array (n_elements(mu), n_l) of P^l_{2,-2} values where l is order.
;
; OPTIONAL OUTPUTS:
;
; KEYWORD OUTPUTS:
;
; RESTRICTIONS:
;
; MODIFICATION HISTORY:
;     G. McGarragh, 28 Jul 2015: Taken from LMie and translated from C.
;         Originally based on Mischenko 1991.
;-

pro gen_sph_funcs, n_l, mu, p00, p0p2, p2p2, p2m2

    n_mu = n_elements(mu)

    da_0p2 = .25d * sqrt(6.d)

    p00 = dblarr(n_mu, n_l)
    p0p2 = dblarr(n_mu, n_l)
    p2p2 = dblarr(n_mu, n_l)
    p2m2 = dblarr(n_mu, n_l)

    for i = 0, n_mu - 1 do begin
        mu2 = mu[i]*mu[i]

        if n_l gt 0 then begin
            p00[i,0] = 1.d
            p0p2[i,0] = 0.d
            p2p2[i,0] = 0.d
            p2m2[i,0] = 0.d
        endif

        if n_l gt 1 then begin
            p00[i,1] = mu[i]
            p0p2[i,1] = 0.d
            p2p2[i,1] = 0.d
            p2m2[i,1] = 0.d
        endif

        if n_l gt 2 then begin
            p00[i,2] = (3.d * mu2 - 1.d) / 2.d
            p0p2[i,2] = da_0p2 * (mu2 - 1.d)

            da = (1.d + mu[i])
            p2p2[i,2] = 0.25d * da*da

            da = (1.d - mu[i])
            p2m2[i,2] = 0.25d * da*da
        endif
    endfor

    for j = 2, n_l - 1 - 1 do begin
        j_2    = j * j
        jp1    = j + 1.d
        jjp1   = j * jp1
        jp1_2  = jp1 * jp1
        twojp1 = 2 * j + 1.d

        da_0p2 = sqrt(j_2 - 4.d)
        db_0p2 = sqrt(jp1_2 - 4.d)

        da_2x2 = twojp1 * 4.d
        db_2x2 = jp1 * (j_2 - 4.d)
        dc_2x2 = j * (jp1_2 - 4.d)

        for i = 0, n_mu - 1 do begin
            da = twojp1 * mu[i]
            p00 [i,j+1] = (da * p00[i,j] - j * p00[i,j-1]) / jp1

            p0p2[i,j+1] = (da * p0p2[i,j] - da_0p2 * p0p2[i,j-1]) / db_0p2

            da = twojp1 * jjp1 * mu[i]
            p2p2[i,j+1] = ((da - da_2x2) * p2p2[i,j] - db_2x2 * p2p2[i,j-1]) / dc_2x2

            p2m2[i,j+1] = ((da + da_2x2) * p2m2[i,j] - db_2x2 * p2m2[i,j-1]) / dc_2x2
         endfor
    endfor
end
