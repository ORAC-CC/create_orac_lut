; Copyright (C) 1998-2017 University of Oxford
;
; This source code is licensed under the GNU General Public License (GPL),
; Version 3.  See the file COPYING for more details.


; See main subroutine for comments
Function BesselZero, s

    A1 =    -0.60792710185402662866327677925836d0
    A2 =     6.159589352810601113491669960271217d-0002
    A3 =    -0.981080301612647885671934857251079d0
    A4 =    11.7507817547326698965409187559573d0
    A5 = -2731.24978203593727776707267899459d0

    b = 4d0*s + 1d0
    BslZ = 0.25d0 * !DPi * b * (1d0 + A1/(b^2) + A2/(b^4) + A3/(b^6) + A4/(b^8) + A5/(b^10))
    Return, BslZ
End ; Function BesselZero

Function FirstGuess, QuadType, N, Term

    Case strupcase(Quadtype) of
    'G': FG = cos(!dpi *(term-0.25d0)/(n+0.5d0))
    'R': FG = (cos(!dpi *(term-0.25d0)/(n+0.5d0)) + cos(!dpi *(term-0.25d0)/(n-0.5d0))) / 2
    'L': FG = cos(BesselZero(term)/sqrt((n-0.5d0)^2 + (0.25d0 - 1d0/!pi^2)))
    EndCase
    Return, FG
End ; Function FirstGuess

Function NewtonG, QuadType, N, X

    Case strupcase(Quadtype) of
    'G': begin
        Legendre, n, x, Pl, Pm, Pn
        Nwt = (1d0 - x^2)*Pn / (Double(n) * (Pm - x*Pn))
    End
    'R': begin
        Legendre, n, x, Pl, Pm, Pn
        Nwt = (1d0 + x)*(Pm+Pn) / ( ((n-1)*Pl - (n-1)*x*Pm +n*Pm - n*x*Pn)/(1-x) - (Pm+Pn) )
    End
    'L': begin
        Legendre, n-1, x, Pl, Pm, Pn
        Nwt = (1d0 - x^2)*(Pm-x*Pn) / (n*Pl + 2d0*(1d0-n)* x * Pm + (x^2*(n-1d0) - 1d0)*Pn)
    End
    EndCase
    Return, Nwt
End ; Function Newton

Pro Legendre, nn, x, Pl, Pm, Pn
    Case NN of
    0: Pn = 1D0 ; Pl and Pm are undefined
    1: begin    ; Pl is undefined
        Pm = 1D0
        Pn = X
    End
    Else: begin ; NN  >= 2
        Pl = 1d0  ; value doesn't matter
        Pm = 1d0
        Pn = x
        For n = 2,nn do begin
            in = 1d0/Double(n)
            Pl = Pm
            Pm = Pn
            Pn = (2d0-in) * x * Pm - (1d0-in) * Pl
        EndFor
    EndElse
    EndCase
    Return
End ; Procedure Legrendre


; Begin main subroutine

Pro Quadrature, Quadtype, NPts, Abscissa, Weight
; Assign N abscissa and weights for integration on the interval [-1,1] for a
; variety of quadrature types:
; S Simpson's
; T Trapezium
; G Gaussian
; R Radau
; L Lobarto
; HISTORY
; R. Grainger, 7 Jun 2005 : Enlarged from GET's translation of RGG FORTRAN
; G. Thomas, 10 Jun 2005 : Fixed Trapezium quadrature (it was providing results
;     which were a factor of 2 too small).

    N = Double(NPts)        ; Double version of Npts
    Abscissa = dblarr(NPts) ; Quadrature points
    Weight   = dblarr(NPts) ; and weights
    Sigfig   = 14           ; Minimum precision of Double data type in IDL

    Case strupcase(Quadtype) of
    'T': begin
        Abscissa       = -1D0+2*Dindgen(NPts)/(N-1)
        Weight(*)      =  2D0/(N-1)
        Weight(0)      =  1D0/(N-1)
        Weight(Npts-1) =  1D0/(N-1)
        Zeros = 1     ; Prevents unnecessary calculation of Legendre moments etc.
    End
    'S': begin
        Stop,'Simpson Method not yet implemented, sorry'
;       If (Npts Mod 2  Eq 0) Then stop,'Error in quadrature:N must be even',N
    End
    Else: begin
        If (NPts gt 2000) then  Stop,'Error in quadrature: Too many quadrature points',n
        Case strupcase(Quadtype) of
        'G': begin
          Term = 0
          Zeros = NPts
        End
        'R': begin
          Abscissa(Npts-1) = -1D0
          Weight(NPts-1)   =  2D0/(N*N)
          Term  = 0
          Zeros = Npts-1
        End
        'L': begin
          Abscissa(0) = 1
          Weight(0) = 2d0/(n*(n-1))
          Abscissa(NPts-1) = -1
          Weight(NPts-1) = 2d0/(n*(n-1))
          Term = 1
          Zeros = NPts-2
        End
        Else: Stop,'Error in quadrature: Invalid quadrature type'
        EndCase

        For Zero=1,Zeros do begin
            lx = 19262d0   ; Any number, Don's birthday apparently.
            x = Firstguess(QuadType,NPts,Zero)
            While Abs(x - lx) gt 1d1^(-SigFig) do begin
              lx = x
              x = x - NewtonG(QuadType,NPts,x)
            EndWhile
            Abscissa(Term)=x
            Legendre,n,x,Pl,Pm,Pn
            Case strupcase(Quadtype) of
            'G': Weight(Term) = 2d0 * (1d0-X*X) / (N*Pm)^2
            'R': Weight(Term) = (1d0-x) / (n*Pm)^2
            'L': Weight(Term) = 2d0 / (Double(n*(n-1))*Pm^2)
            EndCase
            Term = Term + 1
        EndFor
        Abscissa = Reverse(Abscissa)
        Weight = Reverse(Weight)
    EndElse
    EndCase

    Return
End ; Of main procedure Quadrature
