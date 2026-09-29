; Copyright (C) 1998-2017 University of Oxford
;
; This source code is licensed under the GNU General Public License (GPL),
; Version 3.  See the file COPYING for more details.


pro legpexp, Inp, qv, qw, phase, Inlc, lc

; Expands function as a Legendre series.  The function is assumed to have been
; evaluated at Legendre quadrature points - Aqv.  The evaluation stops when the
; number of terms exceeds the number of quadrature points or when the absolute
; value of the Legendre coefficient is less than E-5.
;
; Converted from the "Alegpexp" fortran subroutine written by Don Grainger.
;
; Inp             = Number of points
; qv              = Legendre point
; qw              = Legendre weight
; phase           = Input (phase) function
; lc              = Legendre coefficients
; lpnm, lpnm, lpn = Legendre polynomials

;   Imaxnp = 1100
    Imaxnp = 20000

    lc = dblarr(Inp)

    if Inp gt Imaxnp then stop, 'Error in legpexp: Too many quadrature points'

    lc(0) = total(phase*qw)/2d0
    lc(1) = 3d0 * total(phase*qv*qw)/2d0
    lpnm2 = replicate(1d0,Inp)
    lpnm1 = qv
    n=2

    while n lt Inp do begin
;       calculate the nth Legendre polynomial
        lpn = (double(2*n-1)/n) * qv * lpnm1 - (double(n-1)/n) * lpnm2
        lpnm2 = lpnm1
        lpnm1 = lpn

;       integrate up Legendre coefficient
        lc(n) = (2*n+1) * total(phase * lpn * qw)/2
        if abs(lc(n)) lt 1d-9 then goto, j20
        n = n+1
    endwhile

j20:Inlc = n-1

end
