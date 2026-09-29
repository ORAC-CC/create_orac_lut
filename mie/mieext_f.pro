; Copyright (C) 1998-2017 University of Oxford
;
; This source code is licensed under the GNU General Public License (GPL),
; Version 3.  See the file COPYING for more details.


pro mieext_f, Npts, Dx, Cm, Dqxt, Dqsc, Dg

; Multi architecture Mie routine.
;
; INPUT
; Npts, I4,       Number of size parameters in <tt>Dx</tt>
; Dx,   R8(Npts), Particle size parameter(s)
; Cm,   C8,       Complex refractive index of the particles
;
; OUTPUT
; Dqxt, R8(Npts), Extinction efficiency
; Dqsc, R8(Npts), Scattering efficiency
; Dg,   R8(Npts), Asymmetry parameter

    case !Version.arch of
        'x86'  : status = CALL_EXTERNAL('/home/crun/eodg/idl/mie/mieext_x86.so'  , $
                                        'mieextidl_', Npts, Dx, Cm, Dqxt, Dqsc, Dg)
        'alpha': status = CALL_EXTERNAL('/home/crun/eodg/idl/mie/mieext_alpha.so', $
                                        'mieextidl_', Npts, Dx, Cm, Dqxt, Dqsc, Dg)
        else: Mieext, Npts, Dx, Cm, Dqxt, Dqsc, Dg
    endcase
end
