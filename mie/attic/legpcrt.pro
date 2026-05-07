; Copyright (C) 1998-2017 University of Oxford
;
; This source code is licensed under the GNU General Public License (GPL),
; Version 3.  See the file COPYING for more details.


pro legpcrt, Inlc, lc, Inp, qv, phase

; Inlc  = Number of Legendre coefficient
; lc    = Legendre coefficients
; Inp   = Number of quadrature value=number of angles
; qv    = cos(angles)
; phase = Recomputed phase function

;   Imaxnp = 1100
    Imaxnp = 20000

    phase = dblarr(Inp)
    lpn   = dblarr(Imaxnp)
    lpmn1 = dblarr(Imaxnp)
    lpmn2 = dblarr(Imaxnp)

    for i = 0, Inp - 1 do begin
        N = 0
        lpmn2(i) = 1d0
        phase(i) = lc(N)

        N = 1
        lpmn1(i) = qv(i)
        phase(i) = phase(i)+lc(N)*lpmn1(i)
        for N = 2, Inlc - 1 do begin
            lpn(i) = ((2*N-1)*qv(i)*lpmn1(i)-(N-1)*lpmn2(i))/N
            phase(i) = phase(i)+lc(N)*lpn(i)
            lpmn2(i) = lpmn1(i)
            lpmn1(i) = lpn(i)
        endfor
    endfor

    return
end
