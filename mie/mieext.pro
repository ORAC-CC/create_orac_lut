; Copyright (C) 1998-2017 University of Oxford
;
; This source code is licensed under the GNU General Public License (GPL),
; Version 3.  See the file COPYING for more details.


pro mieext, Npts, Dx, Cm, Dqxt, Dqsc, Dg

    Imaxx = 2600
    Itermax = 3500
    Imaxnp = 1100 ; Change this as required

    Dqxt(*) = 0D0
    Dqsc(*) = 0D0
    Dg(*) = 0D0

    for I = 0, Npts - 1 do begin
        if (Dx(I) gt Imaxx) then message, 'Error: Size Parameter Overflow in Mie'
        Ir = 1.D0 / Cm
        Y =  Dx(I) * Cm

        if (Dx(I) lt 0.02) then NStop = 2 else $
          begin
            case 1 OF
            (Dx(I) le 8)    : NStop = Dx(I) + 4.00*Dx(I)^(1./3.) + 2.0
            (Dx(I) lt 4200) : NStop = Dx(I) + 4.05*Dx(I)^(1./3.) + 2.0
            else            : NStop = Dx(I) + 4.00*Dx(I)^(1./3.) + 2.0
           	endcase
          end
        NmX = fix(max([NStop,abs(Y)]) + 15.)
        D = dcomplexarr(Nmx+1)

        for N = Nmx-1,1,-1 do $
          begin
            A1 = (N+1) / Y
            D[N] = A1 - 1/(A1+D[N+1])
          end

        Psi0 = cos(Dx(I))
        Psi1 = sin(Dx(I))
        Chi0 =-sin(Dx(I))
        Chi1 = cos(Dx(I))
        Xi0 = dcomplex(Psi0,Chi0)
        Xi1 = dcomplex(Psi1,Chi1)

        Tnp1 = 1

        for N = 1,Nstop do $
          begin
            DN = double(N)
            Tnp1 = Tnp1 + 2
            Tnm1 = Tnp1 - 2
            A2 = Tnp1 / (DN*(DN+1D0))
            Turbo = (DN+1D0) / DN
            Rnx = DN/Dx(I)
            Psi = double(Tnm1)*Psi1/Dx(I) - Psi0
            Chi = Tnm1*Chi1/Dx(I)       - Chi0
            Xi = dcomplex(Psi,Chi)
            A = ((D[N]*Ir+Rnx)*Psi-Psi1) / ((D[N]*Ir+Rnx)*  Xi-  Xi1)
            B = ((D[N]*Cm+Rnx)*Psi-Psi1) / ((D[N]*Cm+Rnx)*  Xi-  Xi1)
            Dqxt(I) = Tnp1 * double(A + B) + Dqxt(I)
            Dqsc(I) = Tnp1 * double(A*conj(A) + B*conj(B)) + Dqsc(I)
            if (N gt 1) then Dg(I) = Dg(I) $
                                 + (dN*dN - 1) * double(ANM1*conj(A) + BNM1 * conj(B)) / dN $
                                 + TNM1 * double(ANM1*conj(BNM1)) / (dN*dN - dN)
            Anm1 = A
            Bnm1 = B
            APB = A2 * (A + B)
            AMB = A2 * (A - B)


            Psi0 = Psi1
            Psi1 = Psi
            Chi0 = Chi1
            Chi1 = Chi
            Xi1 = dcomplex(Psi1,Chi1)
          end; for Nstop

        if (Dg(I) gt 0) then Dg(I) = 2 * Dg(I) / Dqsc(I)
        Dqsc(I) =  2 * Dqsc(I) / Dx(I)^2
        Dqxt(I) =  2 * Dqxt(I) / Dx(I)^2
    endfor
end
