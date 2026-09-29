"""IDL-structured Python routines for ORAC LUT generation.

Each module in this package is the Python counterpart of one legacy IDL
routine of the coherent working tree (``current_local_create_orac_lut``):

   load_inststr.py                    <- load_inststr.pro
   load_lutstr.py                     <- load_lutstr.pro
   load_srfstrarr.py                  <- load_srfstrarr.pro, read_srfstr.pro,
                                         load_solar_spectrum.pro, bbconstants.pro,
                                         integrate_trapeziodal.pro, simpson_integral.pro
   load_mmdat.py                      <- load_mmdat.pro, read_ri.pro
   load_atmstr.py                     <- load_atmstr.pro
   load_gasstr.py                     <- load_gasstr.pro
   interpol.py                        <- IDL built-in INTERPOL (profile interpolation)
   generate_scattering_properties.py  <- generate_scattering_properties.pro, create_range.pro
   create_bwgp.py                     <- create_bwgp.pro, mie_size_dist_new.pro,
                                         quadrature.pro, shift_quadrature.pro, legpexp.pro
   dubovik.py                         <- dubovik_* kernel readers and log-normal routines
   setup_disort.py                    <- setup_disort.pro
   call_disort.py                     <- call_disort.pro
   write_v2_lut.py                    <- write_v2_lut.pro

The main programs create_orac_cloud_lut() and create_orac_aerosol_lut() live in
the top-level file create_orac_luts.py, next to create_orac_cloud_lut.pro and
create_orac_aerosol_lut.pro.  IDL structures are represented by
types.SimpleNamespace so that ``inststr.channelid`` reads as in IDL.

The numerical kernels (Mie: mie/dlm-code/mieint.f; DISORT: create_orac_lut/
disort2/src) are reached through the already validated wrappers
oraclut.optics.legacy_mie and oraclut.radiative_transfer.legacy_disort.
"""

from .load_inststr import load_inststr
from .load_lutstr import load_lutstr, lut_quadrature
from .load_srfstrarr import load_srfstrarr, read_srfstr, load_solar_spectrum, bbconstants
from .segment import segment, int_tabulated
from .load_mmdat import load_mmdat, read_ri
from .load_atmstr import load_atmstr, read_presdat
from .load_gasstr import load_gasstr, read_gasstr
from .interpol import interpol
from .generate_scattering_properties import generate_scattering_properties, create_range
from .baum import BaumTable, read_baum_lambda
from .create_bwgp import create_bwgp, mie_size_dist_new, quadrature, shift_quadrature, legpexp
from .dubovik import (DubovikFilenames, DubovikGrid, DubovikKext, DubovikScatteringMatrix,
                      VALID_EPS, dubovik_lognormal, dubovik_lognormal_multiple_eps,
                      read_dubovik_filenames, read_dubovik_grid, read_dubovik_kernel_kext,
                      read_dubovik_kernel_scatt_matrix)
from .setup_disort import setup_disort
from .call_disort import call_disort
from .write_v2_lut import write_v2_lut

__all__ = [
   "load_inststr", "load_lutstr", "lut_quadrature", "load_srfstrarr", "read_srfstr",
   "load_solar_spectrum", "bbconstants", "segment", "int_tabulated", "load_mmdat", "read_ri", "load_atmstr",
   "read_presdat", "load_gasstr", "read_gasstr", "interpol",
   "generate_scattering_properties", "create_range", "BaumTable", "read_baum_lambda", "create_bwgp", "mie_size_dist_new",
   "quadrature", "shift_quadrature", "legpexp", "setup_disort", "call_disort", "write_v2_lut",
   "DubovikFilenames", "DubovikGrid", "DubovikKext", "DubovikScatteringMatrix", "VALID_EPS",
   "dubovik_lognormal", "dubovik_lognormal_multiple_eps", "read_dubovik_filenames",
   "read_dubovik_grid", "read_dubovik_kernel_kext", "read_dubovik_kernel_scatt_matrix",
]
