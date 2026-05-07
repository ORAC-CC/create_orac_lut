! Copyright (C) 1998-2017 University of Oxford
!
! This source code is licensed under the GNU General Public License (GPL),
! Version 3.  See the file COPYING for more details.


      Real * 4 Function MIEEXTIDL(Argc, Argv)
!     IDL interface to use MieExt.  Function returns non-zero
!     value to indicate an error has occured.

!     INPUT/OUTPUT
      Integer * 8 Argc
      Integer * 8 Argv(*)
!     LOCAL
      Integer * 4 Error

      Error = 0
      Call MieExt(%Val(Argv(1)), %Val(Argv(2)), %Val(Argv(3)),
     1           %Val(Argv(4)), %Val(Argv(5)), %Val(Argv(6)), Error)
      MieExtIDL = Error

      Return
      End

      Subroutine MieExt(Npts, DxA, SCm, DqxtA, DqscA, DgA, Error)
      Integer * 2  Imaxx
      Parameter (Imaxx = 12000)
      Real * 4     RIMax          ! largest real part of refractive index
      Parameter (RIMax = 2.5)
      Integer * 2  Itermax
      Parameter (Itermax = Imaxx*RiMax) ! must be large enough to cope with the largest possible nmx = x *real(scm) + 15
                                        ! or nmx =  Dx + 4.05*Dx**(1./3.) + 2.0
!     INPUT
      Integer*4    Npts
      Real * 8	   DxA(Npts)
      Complex      SCm
!     OUTPUT
      Real * 8	   DqxtA(Npts)
      Real * 8	   DqscA(Npts)
      Real * 8	   DgA(Npts)
      Integer * 4  Error
!     LOCAL
      Real * 8	   Dx
      Real * 8	   Dqxt
      Real * 8	   Dqsc
      Integer * 2  I
      Integer * 2  NStop
      Integer * 2  NmX
      Integer * 4  N	! N*N > 32767 ie N > 181
      Real * 8	   Chi,Chi0,Chi1
      Real * 8	   APsi,APsi0,APsi1
      Real * 16	   Psi,Psi0,Psi1
      Complex * 8  Cm
      Complex * 8  Ir
      Complex * 16 A,ANM1,APB
      Complex * 16 B,BNM1,AMB
      Complex * 16 D(Itermax)
      Complex * 16 Xi,Xi0,Xi1
      Complex * 16 Y
!     ACCELERATOR VARIABLES
      Integer * 2  Tnp1
      Integer * 2  Tnm1
      Real * 16    Dn
      Real * 8	   Rnx
      Real * 8	   Turbo
      Real * 8	   A2
      Complex * 16 A1


      Do Ipt = 1, Npts
        Dx = DxA(Ipt)

        If (Dx.Gt.Imaxx) Then
          Error = 1
          Return
        EndIf
        Cm = SCm
        Ir = 1 / Cm
        Y =  Dx * Cm
        If (Dx.Lt.0.02) Then
	  NStop = 2
        Else
	  If (Dx.Le.8.0) Then
	    NStop = Dx + 4.00*Dx**(1./3.) + 2.0
	  Else
	    If (Dx.Lt. 4200.0) Then
	      NStop = Dx + 4.05*Dx**(1./3.) + 2.0
	    Else
	      NStop = Dx + 4.00*Dx**(1./3.) + 2.0
	    End If
	  End If
        End If
        NmX = Max(Real(NStop),Real(Abs(Y))) + 15.
        If (Nmx .gt. Itermax) then
          Stop 'Fatal Error NMX Too Small'
        End if
        D(NmX) = Dcmplx(0,0)
        Do N = Nmx-1,1,-1
          A1 = (N+1) / Y
	  D(N) = A1 - 1/(A1+D(N+1))
        End Do
        Psi0 = Cos(Dx)
        Psi1 = Sin(Dx)
        Chi0 =-Sin(Dx)
        Chi1 = Cos(Dx)
        APsi0 = Psi0
        APsi1 = Psi1
        Xi0 = Dcmplx(APsi0,Chi0)
        Xi1 = Dcmplx(APsi1,Chi1)
        Dqsc = 0
        Dqxt = 0
        Dg = 0
        Tnp1 = 1
        Do N = 1,Nstop
	  DN = N
	  Tnp1 = Tnp1 + 2
	  Tnm1 = Tnp1 - 2
	  A2 = Tnp1 / (DN*(DN+1D0))
	  Turbo = (DN+1D0) / DN
	  Rnx = DN/Dx
	  Psi = Dble(Tnm1)*Psi1/Dx - Psi0
	  APsi = Psi
	  Chi = Tnm1*Chi1/Dx       - Chi0
	  Xi = Dcmplx(APsi,Chi)
	  A = ((D(N)*Ir+Rnx)*APsi-APsi1) / ((D(N)*Ir+Rnx)*  Xi-  Xi1)
	  B = ((D(N)*Cm+Rnx)*APsi-APsi1) / ((D(N)*Cm+Rnx)*  Xi-  Xi1)
	  Dqxt = Tnp1 *      Dble(A + B)          + Dqxt
	  Dqsc = Tnp1 * (A*Conjg(A) + B*Conjg(B)) + Dqsc
	  If (N.Gt.1)
     &      Dg = Dg + (dN*dN - 1) * Dble(ANM1*Conjg(A) + BNM1 * Conjg(B)) / dN + TNM1 * Dble(ANM1*Conjg(BNM1)) / (dN*dN - dN)
	  Anm1 = A
	  Bnm1 = B
	  Psi0 = Psi1
	  Psi1 = Psi
	  Apsi1 = Psi1
	  Chi0 = Chi1
	  Chi1 = Chi
	  Xi1 = Dcmplx(APsi1,Chi1)
        End Do
        If (Dg .GT.0) DgA(Ipt) = 2 * Dg / Dqsc
        DqscA(Ipt) =  2 * Dqsc / Dx**2
        DqxtA(Ipt) =  2 * Dqxt / Dx**2
      End Do
      Return
      End
