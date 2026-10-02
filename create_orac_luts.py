#!/usr/bin/env python
"""create_orac_luts.py - generate an ORAC V2 look-up table the way the IDL did.

   python create_orac_luts.py runs/<run-file>.run

This is the Python counterpart of the legacy IDL programs

   makerunfile_v2.pro                       -> the run file (runs/*.run)
   create_orac_cloud_lut_wrapper.pro        -> read_runfile() + main()
   create_orac_cloud_lut.pro                -> create_orac_cloud_lut()
   create_orac_aerosol_lut_wrapper.pro      -> read_runfile() + main()
   create_orac_aerosol_lut.pro              -> create_orac_aerosol_lut()

and it is meant to be read side by side with them.  The two main functions
keep the IDL sequence of operations, loop hierarchy, branch structure and
variable names (inststr, lutstr, srfstrarr, mmstr, atmstr, gasstr, scatreltau,
bextrat, tauscat, dtau, ssalb, pmom, umu, utau, rfd, tfd, rd, td, tb, rfbd, tfbd,
rbd, em ...).  The routines they call live in src/oraclut/idl_mirror/, one
Python module per IDL routine (see docs/idl_python_mapping.md).

The scientific data flow is

   instrument (.inst), LUT grid (.lut), SRFs + solar spectrum, microphysics (.mm)
   + refractive indices (.ri), atmosphere (.atm), gas optical depths (.gas)
      -> generate_scattering_properties: size distribution -> Mie -> bulk
         Bext, w, g, phase function -> Legendre moments (per channel, effective radius)
      -> for every channel / SRF point / (surface pressure) / optical depth /
         effective radius: layer dtau, ssalb, pmom  ->  call_disort
         (diffuse, thermal emission, direct beam for each solar zenith)
      -> RT operators RFD, TFD, RD, TD, TB, RFBD, TFBD, RBD, EM
      -> write_v2_lut (NetCDF).

Numerical results are those of the validated Python implementation
(src/oraclut/pipeline.py); the Mie and DISORT kernels are the preserved
Fortran sources.  Wherever the IDL expression has been simplified or reordered
to keep the products bitwise identical to the validated ones, the IDL form is
quoted in a comment beginning "IDL:".

Versions and numerical changes.  Three things are recorded separately:

1. LUT product / grid version (the run-file ``version``, e.g. V23).  It
   identifies a product set and its grids (e.g. Grid B); it is not a
   description of the source code.
2. Microphysical numerical integration.  Up to source revision 9d663e9 (the
   revision that produced the V23 reference LUTs) every Mie size distribution
   was integrated over the IDL's fixed 0.001-100 um.  Since the 2026-10
   integration-limit change, liquid-water modified-gamma components stop at
   the first node of the same legacy radius lattice at or beyond
   3.5 x effective radius (beyond 100 um when necessary); other components are
   unchanged (generate_scattering_properties.radius_upper_factor).
3. Legendre / moment calculation.  Up to source revision c6ad545 NMom = 1000
   (the IDL's fixed value) was the Gauss-Legendre order, the number of
   Legendre coefficients and the number of DISORT moments for every phase
   function.  Since the 2026-10 Legendre change, each size-distribution-
   averaged Mie phase function is sampled on a Gauss-Legendre order above its
   polynomial-degree bound and keeps the expansion length L given by King's
   criterion (Grainger 1990, section 4.5), accepted only if the series
   reproduces the directly calculated phase function to six significant
   figures; DISORT receives those L moments (padded to NSTR + 1).  Baum and
   T-matrix (tabulated) classes keep the fixed nmom expansion: their Legendre
   convergence is not solved (src/oraclut/idl_mirror/legendre_expansion.py).

The exact source revision used for each validation comparison is recorded in
validation/REPORT_lut_numerics_development.md.
"""

import ast
import contextlib
import datetime
import errno
import functools
import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from oraclut.idl_mirror import (   # noqa: E402  (import after sys.path is set)
   call_disort, generate_scattering_properties, interpol, load_atmstr, load_gasstr,
   load_inststr, load_lutstr, load_mmdat, load_srfstrarr, setup_disort, write_v2_lut,
)
from oraclut.idl_mirror.generate_scattering_properties import uses_adaptive_legendre   # noqa: E402
from oraclut.radiative_transfer.legacy_disort import getmom, plkavg   # noqa: E402  (DISORT GETMOM / PLKAVG)


@contextlib.contextmanager
def private_workdir():
   """Yield a unique temporary directory for one LUT calculation."""

   root = os.environ.get("TMPDIR")
   root_path = Path(root) if root and Path(root).is_dir() else None
   job_id = os.environ.get("SLURM_JOB_ID", "local")
   with tempfile.TemporaryDirectory(prefix=f"oraclut-{job_id}-", dir=root_path) as directory:
      yield Path(directory)


def _with_private_workdir(function):
   """Give direct and command-line LUT calls an isolated working directory."""

   @functools.wraps(function)
   def wrapped(*args, **kwargs):
      requested = kwargs.pop("work_path", None)
      if requested is not None:
         work_path = Path(requested)
         work_path.mkdir(parents=True, exist_ok=True)
         return function(*args, work_path=work_path, **kwargs)
      with private_workdir() as work_path:
         return function(*args, work_path=work_path, **kwargs)

   return wrapped


def _replicate_dual_view(inststr, srfstrarr, rt_arrays, optical_arrays):
   """Apply the legacy IDL forward-view replication after the nadir RT loops."""
   if inststr.view <= 0:
      return inststr, srfstrarr, rt_arrays, optical_arrays

   nadir = inststr.number_of_nadir_channels
   vector_names = (
      "solar_channel_flag", "mixed_channel_flag", "thermal_channel_flag",
      "oldf0", "oldf1", "oldnefr", "oldwvn", "oldb1", "oldb2",
      "oldt1", "oldt2", "oldnebt", "snr", "rgu", "rou", "rua",
      "rub", "ruc", "refbt", "nedt",
   )
   for name in vector_names:
      values = getattr(inststr, name)
      setattr(inststr, name, np.concatenate((values, values)))
   inststr.channelid = np.concatenate((inststr.channelid, inststr.channelid + inststr.view))
   inststr.srf_file = list(inststr.srf_file) + list(inststr.srf_file)
   inststr.number_of_channels = 2 * nadir
   srfstrarr = list(srfstrarr) + list(srfstrarr)

   # IDL copies every RT operator along its leading channel dimension.
   rt_arrays = tuple(np.concatenate((values, values), axis=0) for values in rt_arrays)
   # The legacy IDL calculation does not copy the optical-property arrays when
   # it creates the forward-view channels.  write_v2_lut writes the nadir
   # values and leaves the remaining channel slots at the NetCDF fill value.
   # Keep the arrays at their calculated (nadir-channel) size; the writer
   # performs that legacy padding at serialization time.
   return inststr, srfstrarr, rt_arrays, optical_arrays


def _print_execution_configuration(
   *, forward_model, inststr, mmstr, mmfile, lutfile, lutstr, qm, atmospheres,
   gas_flag, rayleigh_flag, scat_only, nstreams, nmom, version,
   srfstrarr, output,
):
   """Print the resolved configuration and loaded SRF counts for one task."""

   srf_description = {
      1: "centre-wavelength / monochromatic treatment",
      2: "segmented SRF integration",
   }.get(int(qm), "unsupported/unknown")
   print("ORAC LUT calculation")
   print(f"Array task:           {os.environ.get('SLURM_ARRAY_TASK_ID', 'single')}")
   print(f"Platform:             {inststr.platform}")
   print(f"Instrument:           {inststr.instrument}")
   print(f"Forward model:        {forward_model}")
   print(f"Microphysics:         {mmfile}")
   print(f"LUT definition:       {lutfile}")
   print(f"Version:              {version if version is not None else '00'}")
   print(f"Channels:             {np.asarray(inststr.channelid).tolist()}")
   print(f"SRF quadrature:       {qm} ({srf_description})")
   print(f"Atmosphere:           {atmospheres}")
   print(f"Gas absorption:       {'on' if gas_flag else 'off'}")
   print(f"Rayleigh scattering:  {'on' if rayleigh_flag else 'off'}")
   print(f"Scattering only:      {'yes' if scat_only else 'no'}")
   print(f"DISORT streams:       {nstreams}")
   if nmom is None:
      print("Legendre moments:     adaptive (King's criterion on each averaged Mie phase function)")
   else:
      print(f"Legendre moments:     {nmom} (fixed; Baum / T-matrix tabulated phase functions)")
   print(f"Output:               {output}")
   print("LUT dimensions:")
   print(f"  Optical depth:       {lutstr.opd_n}")
   print(f"  Effective radius:    {lutstr.efr_n}")
   print(f"  Solar zenith:        {lutstr.soz_n}")
   print(f"  Satellite zenith:    {lutstr.saz_n}")
   print(f"  Relative azimuth:    {lutstr.raa_n}")
   if hasattr(lutstr, "prs_n"):
      print(f"  Pressure:            {lutstr.prs_n}")
   print("SRF integration setup:")
   print(f"  SRF quadrature method: {qm}")
   for channel, srf in zip(np.asarray(inststr.channelid).tolist(), srfstrarr):
      print(f"  Channel {channel}: SRF integration points: {srf.nwvl}")


def _publish_lut(private_file, final_file):
   """Validate and atomically publish a completed private NetCDF product."""

   from netCDF4 import Dataset

   private_file = Path(private_file)
   final_file = Path(final_file)
   lock_file = final_file.with_name(f".{final_file.name}.publish.lock")
   try:
      lock_fd = os.open(lock_file, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
      os.close(lock_fd)
   except FileExistsError as exc:
      raise FileExistsError(f"Another job is publishing {final_file}") from exc
   try:
      if final_file.exists():
         raise FileExistsError(
            f"Refusing to overwrite existing LUT {final_file}; choose a new version/output or remove it deliberately"
         )
      # Opening the closed file catches truncated/invalid NetCDF output before
      # it becomes visible in the shared production directory.
      with Dataset(private_file, "r"):
         pass
      try:
         # This is atomic when scratch and output share a filesystem.
         os.replace(private_file, final_file)
      except OSError as exc:
         if exc.errno != errno.EXDEV:
            raise
         # SLURM scratch is commonly a different filesystem from the repository.
         # Copy only after validation to a unique staging file in the destination
         # directory, then rename within that filesystem.  The staging file is
         # removed on every exit path and is never a production product.
         staging_fd, staging_name = tempfile.mkstemp(prefix=f".{final_file.name}.", suffix=".tmp", dir=final_file.parent)
         os.close(staging_fd)
         staging = Path(staging_name)
         try:
            with private_file.open("rb") as source, staging.open("wb") as destination:
               shutil.copyfileobj(source, destination)
               destination.flush()
               os.fsync(destination.fileno())
            os.replace(staging, final_file)
         finally:
            staging.unlink(missing_ok=True)
   finally:
      lock_file.unlink(missing_ok=True)


def _scattering_cache_path(out_path, work_path, reuse_scat):
   """Return the private or explicitly persistent scattering-cache path."""

   return Path(out_path) / "scatfile.npz" if reuse_scat else Path(work_path) / "scatfile.npz"


# ==============================================================================
# Run file: the Python counterpart of makerunfile_v2.pro plus the driver file
# read by create_orac_*_lut_wrapper.pro.  Every scientifically important
# setting of a LUT calculation is set in one place.
# ==============================================================================

# Every key a run file must set; nothing is defaulted silently.
REQUIRED_RUN_KEYS = (
   "platform", "instrument", "forward_model", "in_path", "instfile", "mmfile", "lutfile",
   "atmospheres", "channelid", "srf_quad", "nstreams", "version", "out_path",
)
# Keys with the IDL keyword defaults (not set = not present, as in IDL).
# nmom is deprecated: it was the fixed Legendre expansion of the IDL and of
# the Python generator up to c6ad545, and is now accepted only for Baum /
# T-matrix (tabulated) classes.
OPTIONAL_RUN_KEYS = {"gas": 0, "no_rayleigh": 0, "reuse_scat": 0, "scat_only": 0, "tmatrix_path": None, "nmom": None}


def _check_legendre_configuration(mmstr, nmom):
   """Mie classes use the adaptive expansion (nmom obsolete); tabulated classes need nmom."""

   if not uses_adaptive_legendre(mmstr):
      if nmom is None:
         raise ValueError("nmom is required for Baum / T-matrix (tabulated) phase functions; their Legendre "
                          "convergence is not covered by the adaptive Mie expansion")
      return
   if nmom is not None:
      raise ValueError("nmom is obsolete for Mie size distributions: the Legendre expansion length is determined "
                       "from each averaged phase function; remove nmom from the run file")

# IDL create_orac_*_lut.pro: Case Atmospheres of ... (MODTRAN model codes)
ATMOSPHERE_FILES = {
   "0": ("midsatm.dat", "midlatitude summer"),
   "1": ("tro.atm", "tropical summer"),
   "2": ("mls.atm", "midlatitude summer"),
   "3": ("mlw.atm", "midlatitude winter"),
   "4": ("sas.atm", "subarctic  summer"),
   "5": ("saw.atm", "subarctic summer"),
   "6": ("std.atm", "US standard"),
}


def read_runfile(runfile):
   """Read a run file of ``name = value`` lines (comments start with ';' or '#').

   Values use Python/IDL literal syntax: strings in quotes, integers, and
   channel lists such as [1, 9].  Missing required keys and unknown keys are
   errors.
   """

   run = {}
   for number, raw in enumerate(Path(runfile).read_text().splitlines(), start=1):
      line = raw.split(";", 1)[0].split("#", 1)[0].strip()
      if line == "":
         continue
      if "=" not in line:
         raise ValueError(f"{runfile}:{number}: expected 'name = value', got {raw!r}")
      name, value = (part.strip() for part in line.split("=", 1))
      name = name.lower()
      if name not in REQUIRED_RUN_KEYS and name not in OPTIONAL_RUN_KEYS:
         raise ValueError(f"{runfile}:{number}: unknown setting {name!r}. Recognised settings: "
                          + ", ".join(REQUIRED_RUN_KEYS + tuple(OPTIONAL_RUN_KEYS)))
      if name in run:
         raise ValueError(f"{runfile}:{number}: {name!r} is set twice")
      try:
         run[name] = ast.literal_eval(value)
      except (ValueError, SyntaxError) as exc:
         raise ValueError(f"{runfile}:{number}: cannot read the value of {name!r}: {value!r}") from exc
   missing = [name for name in REQUIRED_RUN_KEYS if name not in run]
   if missing:
      if "channelid" in missing:
         raise ValueError(
            f"{runfile}: channelid must be specified explicitly in the run file, for example "
            "channelid = [1] or channelid = [1, 4, 9]"
         )
      raise ValueError(f"{runfile}: required settings missing: {', '.join(missing)}")
   for name, default in OPTIONAL_RUN_KEYS.items():
      run.setdefault(name, default)
   if run["forward_model"] not in ("cloud", "aerosol"):
      raise ValueError(f"{runfile}: forward_model must be 'cloud' or 'aerosol', got {run['forward_model']!r}")
   if str(run["atmospheres"]) not in ATMOSPHERE_FILES:
      raise ValueError(f"{runfile}: atmospheres must be a MODTRAN code 0-6, got {run['atmospheres']!r}")
   channels = run["channelid"]
   if not isinstance(channels, (list, tuple)) or not channels:
      raise ValueError(
         f"{runfile}: channelid must be specified explicitly in the run file and contain at least one "
         "channel, for example channelid = [1] or channelid = [1, 4, 9]"
      )
   try:
      run["channelid"] = [int(channel) for channel in channels]
   except (TypeError, ValueError) as exc:
      raise ValueError(f"{runfile}: channelid must be a list of integer channel IDs, for example [1, 9]") from exc
   if len(set(run["channelid"])) != len(run["channelid"]):
      raise ValueError(f"{runfile}: channelid contains duplicate channel IDs: {run['channelid']}")
   return run


# ==============================================================================
# create_orac_cloud_lut: Python counterpart of create_orac_cloud_lut.pro
# ==============================================================================

@_with_private_workdir
def create_orac_cloud_lut(in_path, instfile, mmfile, lutfile, out_path, atmospheres,
                          channelid=None, gas=0, no_rayleigh=0, srf_quad=None, reuse_scat=0, scat_only=0,
                          tmatrix_path=None, version=None, driver=None, nstreams=60, nmom=None, work_path=None):
   """Generate an ORAC LUT with the cloud formulation (discrete particle layer).

   Arguments follow the IDL function; ``nstreams`` (60) is the value hard-wired
   in the IDL setup_disort.  ``nmom`` was the IDL's fixed number of Legendre
   moments (1000, generate_scattering_properties.pro).  It is obsolete for Mie
   classes, whose expansion length comes from each averaged phase function,
   and is required only for Baum / T-matrix (tabulated) classes
   (_check_legendre_configuration).
   """

   # -----------------------------------------------------------------------------
   # Test input and output files and directories
   # -----------------------------------------------------------------------------
   print("Reading input files...")
   in_path = Path(in_path)
   if not in_path.is_dir():
      raise FileNotFoundError("in_path not found: " + str(in_path))

   # Instrument definition file
   instdirfile = in_path / "inst" / instfile
   if not instdirfile.is_file():
      raise FileNotFoundError("instdirfile not readable: " + str(instdirfile))

   # Microphysical definition file
   mmdirfile = in_path / "microphysics" / mmfile
   if not mmdirfile.is_file():
      raise FileNotFoundError("mmdirfile not readable: " + str(mmdirfile))

   # LUT grid definition file
   lutdirfile = in_path / "lut" / lutfile
   if not lutdirfile.is_file():
      raise FileNotFoundError("lutdirfile not readable: " + str(lutdirfile))

   # Pressure profile definition file
   atmospheres = str(atmospheres)
   atmfile, atmosphere_model = ATMOSPHERE_FILES[atmospheres]
   atmdirfile = in_path / "atm" / atmfile
   if not atmdirfile.is_file():
      raise FileNotFoundError("atmospheric file not readable: " + str(atmdirfile))

   # Finally, check that the output directory exists
   out_path = Path(out_path)
   if not out_path.is_dir():
      print("out_path: " + str(out_path) + " creating ..")
   out_path.mkdir(parents=True, exist_ok=True)
   work_path = Path(work_path)

   # -----------------------------------------------------------------------------
   # Convert keywords to flags
   # -----------------------------------------------------------------------------
   rayleigh_flag = 0 if no_rayleigh else 1
   gas_flag = 1 if gas else 0

   # -----------------------------------------------------------------------------
   # Read in files
   # -----------------------------------------------------------------------------

   # **** Read the instrument parameters file
   inststr = load_inststr(instdirfile, requestedchannelid=channelid)

   # **** Read the LUT parameters file
   lutstr = load_lutstr(lutdirfile, inststr.max_sat_zenith)

   # **** Read the spectral response functions for the relevant channels and
   #      generate the integral quantities that depend upon the SRF
   solar_spectrum_filename = in_path / "sun" / "Gueymard2018.sssi"
   qm = srf_quad if srf_quad else 1
   srfstrarr, nwvl_max = load_srfstrarr(inststr, solar_spectrum_filename, qm, in_path)

   # **** Read the scattering parameters file
   mmstr = load_mmdat(mmdirfile, in_path)
   _check_legendre_configuration(mmstr, nmom)

   # **** Read the atmospheric profile file
   atmstr = load_atmstr(atmdirfile, atmospheres)

   # **** If flagged read the gas OPD files
   gasstr = None
   if gas_flag:
      gasdir = in_path / "gas"
      gasstr = load_gasstr(atmospheres, inststr.platform, inststr.instrument, inststr.channelid, gasdir)
      # Check that the gas OPD values are on the same grid as the pressure profile
      for onegasstr in gasstr:
         if not np.array_equal(atmstr.height, onegasstr.height):
            raise ValueError("The pressure and gas OPD profiles must be on the same altitude grid")

   # -----------------------------------------------------------------------------
   # Report setup.
   # -----------------------------------------------------------------------------
   if inststr.view > 0:
      print("Dual View LUT calculation for substance ", mmstr.substance, " for the ", inststr.instrument, " instrument")
   else:
      print("Single View LUT calculation for substance ", mmstr.substance, " for the ", inststr.instrument, " instrument")
   print("Output placed in: " + str(out_path))
   print(atmosphere_model + " (code = ", atmospheres, ")")
   print("Including gas absorption" if gas_flag else "No gas absorption")
   print("Including Rayleigh scattering" if rayleigh_flag else "No Rayleigh scattering")

   print("LUT Dimensions are:")
   print("           Components ", mmstr.ncomp)
   print("  Instrument channels ", inststr.number_of_nadir_channels)
   print("        Optical depth ", lutstr.opd_n)
   print("     Effective radius ", lutstr.efr_n)
   print("         Solar zenith ", lutstr.soz_n)
   print("    Instrument zenith ", lutstr.saz_n)
   print("     Relative azimuth ", lutstr.raa_n)
   print("SRF quadrature method ", qm)

   # -----------------------------------------------------------------------------
   # Create the output filename
   # -----------------------------------------------------------------------------
   monoorband = "m" if qm == 1 else "b"
   # atmospheric model code: 1X Rayleigh + gaseous absorption for MODTRAN model X
   # (generally aerosol); 00 no Rayleigh (bottom layer of a multilayer cloud);
   # 01 Rayleigh for the whole atmosphere, no gas (generally cloud)
   if gas:
      atmospheric_model_code = "1" + atmospheres
   elif no_rayleigh:
      atmospheric_model_code = "00"
   else:
      atmospheric_model_code = "01"

   versions = "00"
   if version is not None:
      versions = f"{int(version):02d}"
   v2_lut_filename = out_path / (inststr.platform.lower() + "_" + inststr.instrument.lower() + "_" + monoorband
                                 + "_" + mmstr.substance.lower() + "_a" + atmospheric_model_code
                                 + "_p" + mmstr.shortname.lower() + "_v" + versions + ".nc")

   _print_execution_configuration(
      forward_model="cloud", inststr=inststr, mmstr=mmstr, mmfile=mmfile,
      lutfile=lutfile, lutstr=lutstr, qm=qm, atmospheres=atmospheres,
      gas_flag=gas_flag, rayleigh_flag=rayleigh_flag, scat_only=scat_only,
      nstreams=nstreams, nmom=nmom, version=version, srfstrarr=srfstrarr,
      output=v2_lut_filename,
   )

   # Preserve legacy driver provenance only in private scratch.  A Python
   # ``runs/*.run`` file is the input configuration and must remain under
   # ``runs/``; copying any driver into the LUT output directory pollutes
   # ``luts/`` and creates another concurrent-job collision point.
   if driver is not None and Path(driver).suffix.lower() != ".run":
      shutil.copy(driver, work_path)

   # -----------------------------------------------------------------------------
   # Interpolate the particle profile layers onto the atmospheric pressure (and
   # gas OPD) layers.
   # -----------------------------------------------------------------------------

   # Firstly, define the height of the layers, which lie between each pressure level
   nlayers = atmstr.nlevels - 1
   hlayers = (atmstr.height[0:nlayers] + atmstr.height[1:nlayers + 1]) / np.float32(2)
   scatreltau = interpol(mmstr.rext, mmstr.height, hlayers)
   scatreltau = scatreltau / np.sum(scatreltau, dtype=np.float32)

   # -----------------------------------------------------------------------------
   # Generate scattering properties by calling the scattering code (Mie), or by
   # reloading properties saved from a previous run.
   # -----------------------------------------------------------------------------
   # A normal calculation writes only to its private work area.  Explicit
   # reuse retains the legacy persistent-cache contract and reads the cache
   # supplied in out_path without modifying it.
   scatfile = _scattering_cache_path(out_path, work_path, reuse_scat)
   if not reuse_scat:
      (lmom, bext550, w550, g550, phs550, amom550, bextrat, bext, w, g, vavg, phs, amom) = \
         generate_scattering_properties(srfstrarr, nwvl_max, inststr, mmstr, lutstr, nmom, tmatrix_path=tmatrix_path)
      # **** write the scattering parameters for the class as a whole for reuse
      np.savez(scatfile, lmom=lmom, bext550=bext550, w550=w550, g550=g550, phs550=phs550, amom550=amom550,
               bextrat=bextrat, bext=bext, w=w, g=g, vavg=vavg, phs=phs, amom=amom)
   else:
      # **** read the scattering parameters for the class as a whole for reuse
      with np.load(scatfile) as saved:
         if "lmom" not in saved:
            raise ValueError(f"{scatfile} holds fixed-nmom scattering properties from before the adaptive Legendre "
                             "change; recalculate them")
         lmom = saved["lmom"]
         bext550, w550, g550, phs550, amom550 = (saved[k] for k in ("bext550", "w550", "g550", "phs550", "amom550"))
         bextrat, bext, w, g, vavg, phs, amom = (saved[k] for k in ("bextrat", "bext", "w", "g", "vavg", "phs", "amom"))
   if scat_only:
      return 0

   # We want the centre point of the SRF integration which, since the number of
   # points is odd, is at the centre of the channel's spectral interval.
   m = nwvl_max // 2
   bextratout = bextrat[m, :, :]
   bextout = bext[m, :, :]
   ssaout = w[m, :, :]
   gout = g[m, :, :]

   print("")
   print("Scattering parameters calculated for class " + str(out_path) + " Version: " + versions)

   # -----------------------------------------------------------------------------
   # Run DISORT
   # -----------------------------------------------------------------------------

   # **** The variables needed for the DISORT calls are set up for each phase
   #      function inside the loops below (setup_disort), because the number
   #      of Legendre moments differs from one phase function to the next.

   # **** Define the LUT table output variables themselves (IDL FLTARR(channels, efr, opd, ...))
   f32 = np.float32
   # The DISORT loops calculate nadir channels only.  Dual-view instruments
   # expose 2*N channels in the instrument structure, but the legacy IDL
   # allocates these working operators for N nadir channels and replicates them
   # exactly once after the loops.
   nchan = inststr.number_of_nadir_channels
   rfd = np.zeros((nchan, lutstr.efr_n, lutstr.opd_n), f32)
   tfd = np.zeros((nchan, lutstr.efr_n, lutstr.opd_n), f32)
   rd = np.zeros((nchan, lutstr.efr_n, lutstr.opd_n, lutstr.saz_n), f32)
   td = np.zeros((nchan, lutstr.efr_n, lutstr.opd_n, lutstr.saz_n), f32)
   tb = np.zeros((nchan, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n), f32)
   rfbd = np.zeros((nchan, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n), f32)
   tfbd = np.zeros((nchan, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n), f32)
   rbd = np.zeros((nchan, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n, lutstr.saz_n, lutstr.raa_n), f32)
   em = np.zeros((nchan, lutstr.efr_n, lutstr.opd_n, lutstr.saz_n), f32)

   # **** Loop through the channels (and solar zenith angles) and run DISORT for
   #      the beam and diffuse cases.  Also produce the emissivity for the
   #      channels that need it.

   # Define the Rayleigh scattering optical depth in each channel
   wvl_centre = np.asarray([s.wvl_centre for s in srfstrarr], f32)
   if no_rayleigh:
      columntauray = np.full(nchan, f32(1.0e-6), f32)
   else:
      columntauray = (atmstr.pressure[atmstr.nlevels - 1] / f32(1013.0)) / \
                     (f32(117.03) * wvl_centre ** f32(4.0) - f32(1.316) * wvl_centre ** f32(2.0))

   # Molecular (Rayleigh) phase moments come from the DISORT GETMOM procedure
   # (IDL: GETMOM, 2, 0.0, NMom-1 inside the layer loop), called below with the
   # number of moments of each particle phase function.

   # Indices of the upwelling view directions in UU (IDL: UU[2*saz_n-lindgen(saz_n)-1, ...])
   up = slice(2 * lutstr.saz_n - 1, lutstr.saz_n - 1, -1)

   for l in range(inststr.number_of_nadir_channels):
      print("Running DISORT for channel " + str(inststr.channelid[l]) + " (" + str(srfstrarr[l].wvl_centre) + " um)")
      for m in range(srfstrarr[l].nwvl):

         # Do we have a gas optical depth profile for the current channel?  If
         # not, use zeros (no gas absorption).  Remember that the gas OPD is
         # defined on the levels between each atmospheric layer, not on the layers.
         if gasstr is not None:
            gasindx = [i for i in range(len(gasstr)) if int(gasstr[i].channelid) == int(inststr.channelid[l])]
            if len(gasindx) == 0:
               raise ValueError("Gas absorption channel mismatch")
            gaslvl = gasstr[gasindx[0]].tau_gas
         else:
            print("No gas optical depth profile.")
            gaslvl = np.zeros(atmstr.nlevels, f32)

         # Calculate the cumulative Rayleigh optical depth at each level
         raylvl = columntauray[l] * np.exp(f32(-0.1188) * atmstr.height - f32(0.00116) * atmstr.height ** f32(2.0))

         # The optical depths of each layer from gas absorption and Rayleigh
         # scattering are the differences between adjacent levels.
         taugas = gaslvl[1:nlayers + 1] - gaslvl[0:nlayers]
         tauray = raylvl[1:nlayers + 1] - raylvl[0:nlayers]

         for a in range(lutstr.opd_n):
            for r in range(lutstr.efr_n):

               # THE OPTICAL DEPTH FROM THE PARTICLES IS THE DESIRED TOTAL OPTICAL DEPTH
               # FOR THE ORAC LUT * THE RELATIVE OPTICAL DEPTH AT EACH LAYER FOR THIS CLASS *
               # THE SCALING FACTOR RELATING THE OPTICAL DEPTH AT THIS WAVELENGTH BACK TO 550 NM.
               tauscat = lutstr.opd[a] * scatreltau * bextrat[m, l, r]
               # the optical depths are additive
               dtau = taugas + tauray + tauscat
               totaltau = np.sum(dtau, dtype=f32)

               # The single scattering albedo is weighted by optical depth.
               # NB. SSA for Rayleigh scattering = 1, and is effectively 0 for gas
               # absorption.  (Layers with dtau = 0 are given SSALB = 0.)
               ssalb = np.divide(tauray + w[m, l, r] * tauscat, dtau, out=np.zeros_like(dtau), where=dtau != 0)
               # Check that we have no SSALB values over 1.0 (this can happen in
               # layers with no absorption due to rounding): DISORT would exit.
               ssalb[ssalb > 1.0] = 1.0

               # THE ASYMMETRY PARAMETER IS ONLY NON-ZERO WHERE WE ACTUALLY HAVE PARTICLES
               asym = np.zeros(nlayers, f32)
               asym[tauscat > 0.0] = g[m, l, r]

               # NOW USE THE GETMOM PROCEDURE (PART OF DISORT) TO GENERATE PHASE FUNCTION
               # MOMENTS FOR THE MOLECULAR SCATTERING AND COMBINE THEM WITH THE PARTICLE
               # MOMENTS GENERATED EARLIER, WEIGHTED BY SCATTERING OPTICAL DEPTH.
               # The number of Legendre moments of this phase function (determined from
               # the averaged phase function, generate_scattering_properties).  DISORT
               # needs NMOM >= NSTR (delta-M uses PMOM(NSTR)), so a shorter expansion is
               # padded with zero moments, below the termination threshold.
               nlmom = max(int(lmom[m, l, r]), nstreams + 1)
               disort_vars = setup_disort(nstreams, nlayers, lutstr.saz_n, lutstr.raa_n, nlmom)
               rayleigh_pm = getmom(2, 0.0, nlmom - 1)
               particle_pm = np.zeros(nlmom, f32)
               particle_pm[:min(nlmom, amom.shape[0])] = amom[:min(nlmom, amom.shape[0]), m, l, r]
               pmom = np.zeros((nlmom, nlayers), f32, order="F")
               for h in range(nlayers):
                  if asym[h] == 0.0:
                     pm = rayleigh_pm                                          # Rayleigh scattering only
                  else:
                     pm = (rayleigh_pm * tauray[h] + particle_pm * w[m, l, r] * tauscat[h]) / \
                          (tauray[h] + w[m, l, r] * tauscat[h])
                     if np.any(pm > 1.0):
                        pm = pm.copy()
                        pm[pm > 1.0] = 1.0
                        print("Unusual warning: phase moment  gt 1")
                  pmom[:, h] = pm

               # We are now ready to call DISORT.  Call the fast diffuse calculation
               # first (errors and problems are more likely to turn up quickly that way).
               print("RT calculation for Channel: " + str(inststr.channelid[l]) + f" ({srfstrarr[l].wvl_centre:6.3f}um)"
                     + ", SRF Point: " + str(m) + ", Tau: " + str(lutstr.opd[a]) + ", EfR: " + str(lutstr.efr[r]))

               # ---------- DIFFUSE -----------
               fbeam = 0.0                                                    # Direct beam intensity
               fisot = 100.0                                                  # Isotropic illumination intensity
               umu0 = f32(np.cos(np.deg2rad(f32(50.0))))                      # A nominal value for beam zenith
               utau = np.asarray([0.0, totaltau], f32)                        # Output layers (in optical depth)
               # UMu describes downwelling as well as upwelling radiance.  If 90
               # degrees is among the satellite zeniths take a small angle off it,
               # and keep |UMu| < 1, to avoid numerical problems in DISORT.
               tmpsat = lutstr.saz.copy()
               tmpsat[tmpsat == 90.0] = f32(89.99)
               umu = np.zeros(2 * lutstr.saz_n, f32)
               umu[0:lutstr.saz_n] = -f32(1.0) * np.cos(np.deg2rad(tmpsat)).astype(f32)
               umu[lutstr.saz_n:2 * lutstr.saz_n] = -f32(1.0) * umu[lutstr.saz_n - 1::-1]
               bd = np.abs(umu) == 1.0
               umu[bd] = f32(0.99999) * umu[bd] / np.abs(umu[bd])

               rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed = call_disort(
                  disort_vars, dtau, ssalb, pmom, utau, umu, lutstr.raa, fbeam, umu0, fisot)

               # Generate the diffuse LUT variables.
               # IDL: RFD += (100.*FlUp[0]/(FIsot*!pi)) * val[m]; with FIsot = 100 the factor
               # cancels and the validated Python evaluates FlUp[0]/!pi directly.
               rfd[l, r, a] += (flup[0] / f32(np.pi)) * srfstrarr[l].val[m]
               tfd[l, r, a] += (rfldn[1] / f32(np.pi)) * srfstrarr[l].val[m]
               # RD contains the upwelling intensity
               rd[l, r, a, :] += uu[up, 0, 0] * srfstrarr[l].val[m]
               # TD contains the downwelling intensity without the direct beam
               if inststr.solar_channel_flag[l]:
                  td[l, r, a, :] += (uu[0:lutstr.saz_n, 1, 0] - f32(100) * np.exp(totaltau / umu[0:lutstr.saz_n])) * srfstrarr[l].val[m]
               else:
                  td[l, r, a, :] += uu[0:lutstr.saz_n, 1, 0] * srfstrarr[l].val[m]

               # If the channel has the emission flag set, calculate the emissivity.
               if inststr.thermal_channel_flag[l]:
                  # ---------- EMISSION ----------
                  # Use DISORT: both beam and diffuse irradiances zero, emission across
                  # a window of 1% of the nominal wavenumber, otherwise as the diffuse case.
                  fbeam = 0.0                                                 # Direct beam intensity
                  fisot = 0.0                                                 # Isotropic illumination intensity
                  wn = f32(1e4) / srfstrarr[l].wvl_centre
                  wnlo = f32(0.995) * wn
                  wnhi = f32(1.005) * wn
                  temp = 250.0
                  incloud = np.flatnonzero(tauscat > 0.0)
                  emnly = incloud.size
                  if emnly > 0:
                     emtau = dtau[incloud]
                     emssa = ssalb[incloud]
                     empmo = np.asfortranarray(pmom[:, incloud])
                     emutau = np.asarray([0.0, np.sum(emtau, dtype=f32)], f32)

                     rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed = call_disort(
                        disort_vars, emtau, emssa, empmo, emutau, umu, lutstr.raa, fbeam, umu0, fisot,
                        plank=True, wnlo=wnlo, wnhi=wnhi, temp=np.full(emnly + 1, temp, f32), nlayer=emnly)
                     # Now calculate the Planck emission across the wavelength interval.
                     bbe = plkavg(float(wnlo), float(wnhi), temp)
                     # Finally, combine to produce the emissivity.
                     # IDL: Em += (100.0*UU[up,0,0]/BBE) * val[m], divided by 100 again in
                     # write_v2_lut; the validated Python keeps the fraction UU/BBE throughout.
                     em[l, r, a, :] += (uu[up, 0, 0] / f32(bbe)) * srfstrarr[l].val[m]

               # Now loop over the solar zenith angles and do the direct beam
               # calculations.  Only needed for channels with a solar component.
               if inststr.solar_channel_flag[l]:
                  # ---------- DIRECT ----------
                  for s in range(lutstr.soz_n):
                     fbeam = 100.0                                            # Direct beam intensity
                     fisot = 0.0                                              # Isotropic illumination intensity
                     # A solar zenith of 90 degrees is altered slightly to prevent DISORT crashing
                     tmpsol = f32(89.99) if lutstr.soz[s] == 90.0 else lutstr.soz[s]
                     umu0 = f32(np.cos(np.deg2rad(tmpsol)))

                     rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed = call_disort(
                        disort_vars, dtau, ssalb, pmom, utau, umu, lutstr.raa, fbeam, umu0, fisot)

                     # Generate the direct beam LUT variables
                     tb[l, r, a, s] += (f32(100.0) * rfldir[1] / rfldir[0]) * srfstrarr[l].val[m]
                     rfbd[l, r, a, s] += (f32(100.0) * flup[0] / rfldir[0]) * srfstrarr[l].val[m]
                     tfbd[l, r, a, s] += (f32(100.0) * rfldn[1] / rfldir[0]) * srfstrarr[l].val[m]
                     for p in range(lutstr.raa_n):
                        # Reverse azimuth to the ORAC convention (raa must be evenly spaced 0-180)
                        p2 = lutstr.raa_n - p - 1
                        # As with the diffuse case, RBD contains the upwelling intensity
                        rbd[l, r, a, s, :, p2] += (uu[up, 0, p] * f32(np.pi)) * srfstrarr[l].val[m]
            # End of EfR loop
         # End of AOD loop
      # End of SRF loop
   # End of channel loop

   # Normalise the RT operators with respect to the SRF.  Note 'sum' is unity in
   # monochromatic mode (srf_quad = 1).
   for l in range(inststr.number_of_nadir_channels):
      srfsum = np.sum(srfstrarr[l].val, dtype=f32)
      rfd[l] /= srfsum
      tfd[l] /= srfsum
      rd[l] /= srfsum
      td[l] /= srfsum
      tb[l] /= srfsum
      rfbd[l] /= srfsum
      tfbd[l] /= srfsum
      rbd[l] /= srfsum
      em[l] /= srfsum

   # Replicate nadir view to forward view exactly as the legacy IDL does after
   # completing the nadir radiative-transfer loops.
   inststr, srfstrarr, rt_arrays, optical_arrays = _replicate_dual_view(
      inststr, srfstrarr,
      (rfd, tfd, rd, td, tb, rfbd, tfbd, rbd, em),
      (bextout, bextratout, ssaout, gout),
   )
   rfd, tfd, rd, td, tb, rfbd, tfbd, rbd, em = rt_arrays
   bextout, bextratout, ssaout, gout = optical_arrays

   # -----------------------------------------------------------------------------
   # Create LUT.
   # -----------------------------------------------------------------------------
   private_lut_filename = work_path / v2_lut_filename.name
   print("Info: Creating " + str(v2_lut_filename))

   write_v2_lut(private_lut_filename, lutstr, inststr, srfstrarr, vavg, bextout, bextratout, ssaout, gout,
                td, tfd, rd, rfd, rbd=rbd, rfbd=rfbd, tfbd=tfbd, tb=tb, em=em)
   _publish_lut(private_lut_filename, v2_lut_filename)

   # -----------------------------------------------------------------------------
   # Output termination timestamp.
   # -----------------------------------------------------------------------------
   (work_path / "timestamp.txt").write_text(datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ") + "\n")

   print("ORAC LUT generation completed for substance " + mmstr.substance)
   return 0


# ==============================================================================
# create_orac_aerosol_lut: Python counterpart of create_orac_aerosol_lut.pro
#
# The differences from the cloud formulation are exactly those of the IDL:
#   * the LUT grid has a surface-pressure dimension (load_lutstr /include_pressure);
#   * the operators carry that extra dimension (index k);
#   * the Rayleigh column optical depth is scaled to each LUT surface pressure
#     inside the pressure loop instead of to 1013 hPa;
#   * the emission temperature is 270 K instead of 250 K;
#   * write_v2_lut is called with include_pressure.
# ==============================================================================

@_with_private_workdir
def create_orac_aerosol_lut(in_path, instfile, mmfile, lutfile, out_path, atmospheres,
                            channelid=None, gas=0, no_rayleigh=0, srf_quad=None, reuse_scat=0, scat_only=0,
                            tmatrix_path=None, version=None, driver=None, nstreams=60, nmom=None, work_path=None):
   """Generate an ORAC LUT with the aerosol formulation (particles through the column).

   Arguments as create_orac_cloud_lut (including the meaning of ``nmom``).
   """

   # -----------------------------------------------------------------------------
   # Test input and output files and directories
   # -----------------------------------------------------------------------------
   print("Reading input files...")
   in_path = Path(in_path)
   if not in_path.is_dir():
      raise FileNotFoundError("in_path not found: " + str(in_path))

   # Instrument definition file
   instdirfile = in_path / "inst" / instfile
   if not instdirfile.is_file():
      raise FileNotFoundError("instdirfile not readable: " + str(instdirfile))

   # Microphysical definition file
   mmdirfile = in_path / "microphysics" / mmfile
   if not mmdirfile.is_file():
      raise FileNotFoundError("mmdirfile not readable: " + str(mmdirfile))

   # LUT grid definition file
   lutdirfile = in_path / "lut" / lutfile
   if not lutdirfile.is_file():
      raise FileNotFoundError("lutdirfile not readable: " + str(lutdirfile))

   # Pressure profile definition file
   atmospheres = str(atmospheres)
   atmfile, atmosphere_model = ATMOSPHERE_FILES[atmospheres]
   atmdirfile = in_path / "atm" / atmfile
   if not atmdirfile.is_file():
      raise FileNotFoundError("atmospheric file not readable: " + str(atmdirfile))

   # Finally, check that the output directory exists
   out_path = Path(out_path)
   if not out_path.is_dir():
      print("out_path: " + str(out_path) + " creating ..")
   out_path.mkdir(parents=True, exist_ok=True)
   work_path = Path(work_path)

   # -----------------------------------------------------------------------------
   # Convert keywords to flags
   # -----------------------------------------------------------------------------
   rayleigh_flag = 0 if no_rayleigh else 1
   gas_flag = 1 if gas else 0

   # -----------------------------------------------------------------------------
   # Read in files
   # -----------------------------------------------------------------------------

   # **** Read the instrument parameters file
   inststr = load_inststr(instdirfile, requestedchannelid=channelid)

   # **** Read the LUT parameters file (with the surface-pressure grid)
   lutstr = load_lutstr(lutdirfile, inststr.max_sat_zenith, include_pressure=True)

   # **** Read the spectral response functions for the relevant channels and
   #      generate the integral quantities that depend upon the SRF
   solar_spectrum_filename = in_path / "sun" / "Gueymard2018.sssi"
   qm = srf_quad if srf_quad else 1
   srfstrarr, nwvl_max = load_srfstrarr(inststr, solar_spectrum_filename, qm, in_path)

   # **** Read the scattering parameters file
   mmstr = load_mmdat(mmdirfile, in_path)
   _check_legendre_configuration(mmstr, nmom)

   # **** Read the atmospheric profile file
   atmstr = load_atmstr(atmdirfile, atmospheres)

   # **** If flagged read the gas OPD files
   gasstr = None
   if gas_flag:
      gasdir = in_path / "gas"
      gasstr = load_gasstr(atmospheres, inststr.platform, inststr.instrument, inststr.channelid, gasdir)
      # Check that the gas OPD values are on the same grid as the pressure profile
      for onegasstr in gasstr:
         if not np.array_equal(atmstr.height, onegasstr.height):
            raise ValueError("The pressure and gas OPD profiles must be on the same altitude grid")

   # -----------------------------------------------------------------------------
   # Report setup.
   # -----------------------------------------------------------------------------
   if inststr.view > 0:
      print("Dual View LUT calculation for substance ", mmstr.substance, " for the ", inststr.instrument, " instrument")
   else:
      print("Single View LUT calculation for substance ", mmstr.substance, " for the ", inststr.instrument, " instrument")
   print("Output placed in: " + str(out_path))
   print(atmosphere_model + " (code = ", atmospheres, ")")
   print("Including gas absorption" if gas_flag else "No gas absorption")
   print("Including Rayleigh scattering" if rayleigh_flag else "No Rayleigh scattering")

   print("LUT Dimensions are:")
   print("           Components ", mmstr.ncomp)
   print("  Instrument channels ", inststr.number_of_nadir_channels)
   print("        Optical depth ", lutstr.opd_n)
   print("     Effective radius ", lutstr.efr_n)
   print("         Solar zenith ", lutstr.soz_n)
   print("    Instrument zenith ", lutstr.saz_n)
   print("     Relative azimuth ", lutstr.raa_n)
   print("            Pressures ", lutstr.prs_n)
   print("SRF quadrature method ", qm)

   # -----------------------------------------------------------------------------
   # Create the output filename
   # -----------------------------------------------------------------------------
   monoorband = "m" if qm == 1 else "b"
   if gas:
      atmospheric_model_code = "1" + atmospheres
   elif no_rayleigh:
      atmospheric_model_code = "00"
   else:
      atmospheric_model_code = "01"

   versions = "00"
   if version is not None:
      versions = f"{int(version):02d}"
   v2_lut_filename = out_path / (inststr.platform.lower() + "_" + inststr.instrument.lower() + "_" + monoorband
                                 + "_" + mmstr.substance.lower() + "_a" + atmospheric_model_code
                                 + "_p" + mmstr.shortname.lower() + "_v" + versions + ".nc")

   _print_execution_configuration(
      forward_model="aerosol", inststr=inststr, mmstr=mmstr, mmfile=mmfile,
      lutfile=lutfile, lutstr=lutstr, qm=qm, atmospheres=atmospheres,
      gas_flag=gas_flag, rayleigh_flag=rayleigh_flag, scat_only=scat_only,
      nstreams=nstreams, nmom=nmom, version=version, srfstrarr=srfstrarr,
      output=v2_lut_filename,
   )

   # Preserve legacy driver provenance only in private scratch; see the cloud
   # formulation above.
   if driver is not None and Path(driver).suffix.lower() != ".run":
      shutil.copy(driver, work_path)

   # -----------------------------------------------------------------------------
   # Interpolate the aerosol profile layers onto the atmospheric pressure and
   # gas OPD layers.
   # -----------------------------------------------------------------------------

   # Firstly, define the height of the layers, which lie between each pressure level
   nlayers = atmstr.nlevels - 1
   hlayers = (atmstr.height[0:nlayers] + atmstr.height[1:nlayers + 1]) / np.float32(2)
   scatreltau = interpol(mmstr.rext, mmstr.height, hlayers)
   scatreltau = scatreltau / np.sum(scatreltau, dtype=np.float32)

   # -----------------------------------------------------------------------------
   # Generate scattering properties by calling the scattering code (Mie), or by
   # reloading properties saved from a previous run.
   # -----------------------------------------------------------------------------
   scatfile = _scattering_cache_path(out_path, work_path, reuse_scat)
   if not reuse_scat:
      (lmom, bext550, w550, g550, phs550, amom550, bextrat, bext, w, g, vavg, phs, amom) = \
         generate_scattering_properties(srfstrarr, nwvl_max, inststr, mmstr, lutstr, nmom, tmatrix_path=tmatrix_path)
      # **** write the scattering parameters for the class as a whole for reuse
      np.savez(scatfile, lmom=lmom, bext550=bext550, w550=w550, g550=g550, phs550=phs550, amom550=amom550,
               bextrat=bextrat, bext=bext, w=w, g=g, vavg=vavg, phs=phs, amom=amom)
   else:
      # **** read the scattering parameters for the class as a whole for reuse
      with np.load(scatfile) as saved:
         if "lmom" not in saved:
            raise ValueError(f"{scatfile} holds fixed-nmom scattering properties from before the adaptive Legendre "
                             "change; recalculate them")
         lmom = saved["lmom"]
         bext550, w550, g550, phs550, amom550 = (saved[k] for k in ("bext550", "w550", "g550", "phs550", "amom550"))
         bextrat, bext, w, g, vavg, phs, amom = (saved[k] for k in ("bextrat", "bext", "w", "g", "vavg", "phs", "amom"))
   if scat_only:
      return 0

   # We want the centre point of the SRF integration
   m = nwvl_max // 2
   bextratout = bextrat[m, :, :]
   bextout = bext[m, :, :]
   ssaout = w[m, :, :]
   gout = g[m, :, :]

   print("")
   print("Scattering parameters calculated for class " + str(out_path) + " Version: " + versions)

   # -----------------------------------------------------------------------------
   # Run DISORT
   # -----------------------------------------------------------------------------

   # **** The variables needed for the DISORT calls are set up for each phase
   #      function inside the loops below (setup_disort), because the number
   #      of Legendre moments differs from one phase function to the next.

   # **** Define the LUT table output variables themselves (IDL FLTARR(channels, prs, efr, opd, ...))
   f32 = np.float32
   # The DISORT loops calculate nadir channels only.  Dual-view instruments
   # expose 2*N channels in the instrument structure, but the legacy IDL
   # allocates these working operators for N nadir channels and replicates them
   # exactly once after the loops.
   nchan = inststr.number_of_nadir_channels
   rfd = np.zeros((nchan, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n), f32)
   tfd = np.zeros((nchan, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n), f32)
   rd = np.zeros((nchan, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n, lutstr.saz_n), f32)
   td = np.zeros((nchan, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n, lutstr.saz_n), f32)
   tb = np.zeros((nchan, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n), f32)
   rfbd = np.zeros((nchan, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n), f32)
   tfbd = np.zeros((nchan, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n), f32)
   rbd = np.zeros((nchan, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n, lutstr.soz_n, lutstr.saz_n, lutstr.raa_n), f32)
   em = np.zeros((nchan, lutstr.prs_n, lutstr.efr_n, lutstr.opd_n, lutstr.saz_n), f32)

   # **** Loop through the channels (and solar zenith angles) and run DISORT for
   #      the beam and diffuse cases.  Also produce the emissivity for the
   #      channels that need it.

   wvl_centre = np.asarray([s.wvl_centre for s in srfstrarr], f32)
   up = slice(2 * lutstr.saz_n - 1, lutstr.saz_n - 1, -1)

   for l in range(inststr.number_of_nadir_channels):
      print("Running DISORT for channel " + str(inststr.channelid[l]) + " (" + str(srfstrarr[l].wvl_centre) + " um)")
      for k in range(lutstr.prs_n):
         # Define the Rayleigh scattering optical depth in each channel, scaled
         # to the LUT surface pressure lutstr.prs[k]
         if no_rayleigh:
            columntauray = np.full(nchan, f32(1.0e-6), f32)
         else:
            columntauray = (atmstr.pressure[atmstr.nlevels - 1] / lutstr.prs[k]) / \
                           (f32(117.03) * wvl_centre ** f32(4.0) - f32(1.316) * wvl_centre ** f32(2.0))
         for m in range(srfstrarr[l].nwvl):

            # Do we have a gas optical depth profile for the current channel?
            if gasstr is not None:
               gasindx = [i for i in range(len(gasstr)) if int(gasstr[i].channelid) == int(inststr.channelid[l])]
               if len(gasindx) == 0:
                  raise ValueError("Gas absorption channel mismatch")
               gaslvl = gasstr[gasindx[0]].tau_gas
            else:
               print("No gas optical depth profile.")
               gaslvl = np.zeros(atmstr.nlevels, f32)

            # Calculate the cumulative Rayleigh optical depth at each level
            raylvl = columntauray[l] * np.exp(f32(-0.1188) * atmstr.height - f32(0.00116) * atmstr.height ** f32(2.0))

            # Layer optical depths from gas absorption and Rayleigh scattering
            taugas = gaslvl[1:nlayers + 1] - gaslvl[0:nlayers]
            tauray = raylvl[1:nlayers + 1] - raylvl[0:nlayers]

            for a in range(lutstr.opd_n):
               for r in range(lutstr.efr_n):

                  # THE OPTICAL DEPTH FROM AEROSOL IS THE DESIRED TOTAL AOD FOR THE ORAC LUT *
                  # THE RELATIVE AOD AT EACH LAYER FOR THIS CLASS * THE SCALING FACTOR RELATING
                  # AOD AT THIS WAVELENGTH BACK TO 550 NM.
                  tauscat = lutstr.opd[a] * scatreltau * bextrat[m, l, r]
                  # the optical depths are additive
                  dtau = taugas + tauray + tauscat
                  totaltau = np.sum(dtau, dtype=f32)

                  # The single scattering albedo is weighted by optical depth.
                  ssalb = np.divide(tauray + w[m, l, r] * tauscat, dtau, out=np.zeros_like(dtau), where=dtau != 0)
                  ssalb[ssalb > 1.0] = 1.0

                  # ASYMMETRY PARAMETER IS ONLY NON-ZERO WHERE WE ACTUALLY HAVE AEROSOL
                  asym = np.zeros(nlayers, f32)
                  asym[tauscat > 0.0] = g[m, l, r]

                  # The number of Legendre moments of this phase function (determined from
                  # the averaged phase function, generate_scattering_properties).  DISORT
                  # needs NMOM >= NSTR (delta-M uses PMOM(NSTR)), so a shorter expansion is
                  # padded with zero moments, below the termination threshold.
                  nlmom = max(int(lmom[m, l, r]), nstreams + 1)
                  disort_vars = setup_disort(nstreams, nlayers, lutstr.saz_n, lutstr.raa_n, nlmom)
                  rayleigh_pm = getmom(2, 0.0, nlmom - 1)
                  particle_pm = np.zeros(nlmom, f32)
                  particle_pm[:min(nlmom, amom.shape[0])] = amom[:min(nlmom, amom.shape[0]), m, l, r]
                  # Molecular and aerosol phase moments combined per layer
                  pmom = np.zeros((nlmom, nlayers), f32, order="F")
                  for h in range(nlayers):
                     scattau = tauray[h] + w[m, l, r] * tauscat[h]
                     if asym[h] == 0.0 or scattau == 0.0:
                        pm = rayleigh_pm
                     else:
                        pm = (rayleigh_pm * tauray[h] + particle_pm * w[m, l, r] * tauscat[h]) / scattau
                        if np.any(pm > 1.0):
                           pm = pm.copy()
                           pm[pm > 1.0] = 1.0
                     pmom[:, h] = pm

                  print("RT calculation for Channel: " + str(inststr.channelid[l]) + f" ({srfstrarr[l].wvl_centre:6.3f}um)"
                        + ", SRF Point: " + str(m) + ", Tau: " + str(lutstr.opd[a]) + ", EfR: " + str(lutstr.efr[r]))

                  # ---------- DIFFUSE -----------
                  fbeam = 0.0
                  fisot = 100.0
                  umu0 = f32(np.cos(np.deg2rad(f32(50.0))))
                  utau = np.asarray([0.0, totaltau], f32)
                  tmpsat = lutstr.saz.copy()
                  tmpsat[tmpsat == 90.0] = f32(89.99)
                  umu = np.zeros(2 * lutstr.saz_n, f32)
                  umu[0:lutstr.saz_n] = -f32(1.0) * np.cos(np.deg2rad(tmpsat)).astype(f32)
                  umu[lutstr.saz_n:2 * lutstr.saz_n] = -f32(1.0) * umu[lutstr.saz_n - 1::-1]
                  bd = np.abs(umu) == 1.0
                  umu[bd] = f32(0.99999) * umu[bd] / np.abs(umu[bd])

                  rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed = call_disort(
                     disort_vars, dtau, ssalb, pmom, utau, umu, lutstr.raa, fbeam, umu0, fisot)

                  # Generate the diffuse LUT variables (IDL: 100.*FlUp[0]/(FIsot*!pi), FIsot = 100)
                  rfd[l, k, r, a] += (flup[0] / f32(np.pi)) * srfstrarr[l].val[m]
                  tfd[l, k, r, a] += (rfldn[1] / f32(np.pi)) * srfstrarr[l].val[m]
                  # RD contains the upwelling intensity
                  rd[l, k, r, a, :] += uu[up, 0, 0] * srfstrarr[l].val[m]
                  # TD contains the downwelling intensity without the direct beam
                  if inststr.solar_channel_flag[l]:
                     td[l, k, r, a, :] += (uu[0:lutstr.saz_n, 1, 0] - f32(100) * np.exp(totaltau / umu[0:lutstr.saz_n])) * srfstrarr[l].val[m]
                  else:
                     td[l, k, r, a, :] += uu[0:lutstr.saz_n, 1, 0] * srfstrarr[l].val[m]

                  # If the channel has the emission flag set, calculate the emissivity.
                  if inststr.thermal_channel_flag[l]:
                     # ---------- EMISSION ----------
                     fbeam = 0.0
                     fisot = 0.0
                     wn = f32(1e4) / srfstrarr[l].wvl_centre
                     wnlo = f32(0.995) * wn
                     wnhi = f32(1.005) * wn
                     temp = 270.0                                             # IDL: "hack from 250"
                     incloud = np.flatnonzero(tauscat > 0.0)
                     emnly = incloud.size
                     if emnly > 0:
                        emtau = dtau[incloud]
                        emssa = ssalb[incloud]
                        empmo = np.asfortranarray(pmom[:, incloud])
                        emutau = np.asarray([0.0, np.sum(emtau, dtype=f32)], f32)

                        rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed = call_disort(
                           disort_vars, emtau, emssa, empmo, emutau, umu, lutstr.raa, fbeam, umu0, fisot,
                           plank=True, wnlo=wnlo, wnhi=wnhi, temp=np.full(emnly + 1, temp, f32), nlayer=emnly)
                        # Planck emission across the wavelength interval
                        bbe = plkavg(float(wnlo), float(wnhi), temp)
                        # Emissivity as a fraction (IDL: 100*UU/BBE here, /100 in write_v2_lut)
                        em[l, k, r, a, :] += (uu[up, 0, 0] / f32(bbe)) * srfstrarr[l].val[m]

                  # Direct beam calculations for each solar zenith (solar channels only)
                  if inststr.solar_channel_flag[l]:
                     # ---------- DIRECT ----------
                     for s in range(lutstr.soz_n):
                        fbeam = 100.0
                        fisot = 0.0
                        tmpsol = f32(89.99) if lutstr.soz[s] == 90.0 else lutstr.soz[s]
                        umu0 = f32(np.cos(np.deg2rad(tmpsol)))

                        rfldir, rfldn, flup, dfdt, uavg, uu, albmed, trnmed = call_disort(
                           disort_vars, dtau, ssalb, pmom, utau, umu, lutstr.raa, fbeam, umu0, fisot)

                        # Generate the direct beam LUT variables
                        tb[l, k, r, a, s] += (f32(100.0) * rfldir[1] / rfldir[0]) * srfstrarr[l].val[m]
                        rfbd[l, k, r, a, s] += (f32(100.0) * flup[0] / rfldir[0]) * srfstrarr[l].val[m]
                        tfbd[l, k, r, a, s] += (f32(100.0) * rfldn[1] / rfldir[0]) * srfstrarr[l].val[m]
                        for p in range(lutstr.raa_n):
                           # Reverse azimuth to the ORAC convention
                           p2 = lutstr.raa_n - p - 1
                           rbd[l, k, r, a, s, :, p2] += (uu[up, 0, p] * f32(np.pi)) * srfstrarr[l].val[m]
               # End of EfR loop
            # End of AOD loop
         # End of SRF loop
      # End of pressure loop
   # End of channel loop

   # Normalise the RT operators with respect to the SRF ('sum' is unity for srf_quad = 1)
   for l in range(inststr.number_of_nadir_channels):
      srfsum = np.sum(srfstrarr[l].val, dtype=f32)
      rfd[l] /= srfsum
      tfd[l] /= srfsum
      rd[l] /= srfsum
      td[l] /= srfsum
      tb[l] /= srfsum
      rfbd[l] /= srfsum
      tfbd[l] /= srfsum
      rbd[l] /= srfsum
      em[l] /= srfsum

   # Replicate nadir view to forward view exactly as the legacy IDL does after
   # completing the nadir radiative-transfer loops.
   inststr, srfstrarr, rt_arrays, optical_arrays = _replicate_dual_view(
      inststr, srfstrarr,
      (rfd, tfd, rd, td, tb, rfbd, tfbd, rbd, em),
      (bextout, bextratout, ssaout, gout),
   )
   rfd, tfd, rd, td, tb, rfbd, tfbd, rbd, em = rt_arrays
   bextout, bextratout, ssaout, gout = optical_arrays

   # -----------------------------------------------------------------------------
   # Create LUT.
   # -----------------------------------------------------------------------------
   private_lut_filename = work_path / v2_lut_filename.name
   print("Info: Creating " + str(v2_lut_filename))

   write_v2_lut(private_lut_filename, lutstr, inststr, srfstrarr, vavg, bextout, bextratout, ssaout, gout,
                td, tfd, rd, rfd, rbd=rbd, rfbd=rfbd, tfbd=tfbd, tb=tb, em=em, include_pressure=True)
   _publish_lut(private_lut_filename, v2_lut_filename)

   # -----------------------------------------------------------------------------
   # Output termination timestamp.
   # -----------------------------------------------------------------------------
   (work_path / "timestamp.txt").write_text(datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ") + "\n")

   print("ORAC LUT generation completed for substance " + mmstr.substance)
   return 0


# ==============================================================================
# main: Python counterpart of the generated <platform>_<instrument>_run script
# and of create_orac_*_lut_wrapper.pro
# ==============================================================================

def run(runfile):
   """Read one run file and generate its LUT.  Returns the IDL-style status (0)."""

   runfile = Path(runfile)
   if not runfile.is_file():
      raise FileNotFoundError(f"run file not found: {runfile}")
   settings = read_runfile(runfile)

   # Relative paths in the run file are taken from the repository root (where
   # this program lives), so the run file reads the same from any directory.
   def from_root(path):
      path = Path(path)
      return path if path.is_absolute() else REPO_ROOT / path

   in_path = from_root(settings["in_path"])
   out_path = from_root(settings["out_path"])

   # The instrument file must belong to the platform / instrument named in the
   # run file (makerunfile_v2 selected the driver from platform + instrument).
   inststr = load_inststr(in_path / "inst" / settings["instfile"])
   if inststr.platform != settings["platform"].lower() or inststr.instrument != settings["instrument"].lower():
      raise ValueError(f"{runfile}: platform/instrument {settings['platform']}/{settings['instrument']} do not "
                       f"match {settings['instfile']} ({inststr.platform}/{inststr.instrument})")
   requested = settings["channelid"]
   valid = [int(channel) for channel in inststr.channelid]
   invalid = [channel for channel in requested if channel not in valid]
   if invalid:
      raise ValueError(
         f"{runfile}: requested channels {requested} include invalid channel(s) {invalid} for "
         f"{inststr.platform}/{inststr.instrument}; valid channels are {valid}"
      )

   print("making ... " + str(runfile))
   if settings["forward_model"] == "cloud":
      generate = create_orac_cloud_lut
   else:
      generate = create_orac_aerosol_lut
   status = generate(in_path, settings["instfile"], settings["mmfile"], settings["lutfile"], out_path,
                     settings["atmospheres"],
                     channelid=settings["channelid"], gas=settings["gas"], no_rayleigh=settings["no_rayleigh"],
                     srf_quad=settings["srf_quad"], reuse_scat=settings["reuse_scat"], scat_only=settings["scat_only"],
                     tmatrix_path=settings["tmatrix_path"],
                     version=settings["version"], driver=runfile,
                     nstreams=settings["nstreams"], nmom=settings["nmom"])
   if status != 0:
      print("create_orac_lut failed with code: " + str(status))
   return status


def main(argv=None):
   argv = sys.argv[1:] if argv is None else list(argv)
   if len(argv) != 1 or argv[0] in ("-h", "--help"):
      print("usage: python create_orac_luts.py runs/<run-file>.run")
      print("       (see runs/template.run for the settings; docs/idl_python_mapping.md for the IDL correspondence)")
      return 2
   try:
      return run(argv[0])
   except (FileNotFoundError, ValueError, NotImplementedError) as exc:
      # A clear one-line failure for a bad run file or an unported option,
      # instead of a traceback (the IDL wrapper printed the failure code).
      print("create_orac_luts.py: " + str(exc), file=sys.stderr)
      return 2


if __name__ == "__main__":
   raise SystemExit(main())
