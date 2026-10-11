"""Generate an ORAC LUT from the exact source files of its LUT release (release launcher, file version v03).

   python create_orac_luts_release_v01.py runs/<product>.run
   python create_orac_luts_release_v01.py --check <lut-release>      check a release's files; compute nothing

Two version numbers are involved and they are independent:

* the LUT release (V21 = the IDL reference tables, V22 = the Python
  reproduction of V21, ...), requested by the run file's ``version``;
* the version of each individual source file (create_bwgp_v01.py,
  create_bwgp_v02.py, ...), which advances only when that file changes.

A LUT release is the set of file versions listed, with SHA-256, in
``releases/lut_v<NN>.release`` or ``releases/lut_v<NN>_<label>.release``.  This launcher

1. reads ``version`` from the run file and opens the manifest of that release:
   the one named by the environment variable ORAC_LUT_RELEASE_MANIFEST if it is
   set (a file name in releases/ or a path), otherwise the only manifest that
   defines that release (it refuses if several do, naming them);
2. refuses to run unless every listed file exists with exactly the listed
   SHA-256 (the release is sealed);
3. imports the driver and every ``oraclut`` module from the files the manifest
   names (create_bwgp_v01.py is loaded as oraclut.idl_mirror.create_bwgp), so
   the modules keep importing each other by their ordinary names and each file
   can move to a new version on its own; any other ``oraclut`` module, and
   any versioned file name, is refused;
4. checks that the compiled Mie and DISORT libraries are the ones the
   manifest pins (building the single-precision ones privately if they are
   missing; the double-precision DISORT library,
   create_orac_lut/disort2_dp/libdisort_dp.so, is never built here), and that
   a manifest declaring ``disort_precision`` is run at that precision
   (ORAC_LUT_DISORT_PRECISION);
5. records every input file the calculation opens, and the Git state, at the
   start and again at the end;
6. runs the release's driver on the run file and, just before the finished
   product is published, writes the provenance into the product as NetCDF
   global attributes (``lut_release``, ``provenance_*``).  The scientific
   variables, dimensions and their attributes are not touched.

Version v02 differs from v01 only in the manifest selection (step 1), the
double-precision library entry and the precision check (step 4).  Version v03
adds one rule: a pinned ``config`` file (CODE_VERSION) that has since been
edited for a newer release is served to the run from its byte-identical
SourceVault copy in the same directory (CODE_VERSION_vNN), and the
substitution is recorded in the provenance, so sealed historical releases
remain runnable.  Neither changes anything numerical.  The importlib finder and
loader classes below are the minimum Python needs for step 3; they exist so
that per-file versions can be combined freely without editing the imports
inside the sealed files.

See SOURCE_VERSIONING.md.
"""

import ast
import datetime
import hashlib
import importlib
import importlib.abc
import importlib.util
import json
import os
import platform
import re
import shutil
import socket
import subprocess
import sys
import tempfile
from pathlib import Path

LAUNCHER_FILE = Path(__file__).resolve()
REPO_ROOT = LAUNCHER_FILE.parent
RELEASE_DIR = REPO_ROOT / "releases"
MANIFEST_SCHEMA = "oraclut-release-manifest-1"
PROVENANCE_SCHEMA = "oraclut-lut-provenance-1"
PUBLIC_REPOSITORY = "https://github.com/ORAC-CC/create_orac_lut"
ROLES = ("driver", "package", "module", "config", "adapter", "kernel", "library")
PYTHON_ROLES = ("driver", "package", "module")
VERSIONED_STEM = re.compile(r"_v\d\d$")
GUARDED_TOP_NAMES = ("oraclut", "create_orac_luts")
PROGRAM = LAUNCHER_FILE.name


class ReleaseRefusal(RuntimeError):
   """The release cannot run, or its provenance cannot be recorded truthfully."""


def refuse(message):
   raise ReleaseRefusal(message)


def sha256_bytes(data):
   return hashlib.sha256(data).hexdigest()


def sha256_file(path):
   digest = hashlib.sha256()
   with open(path, "rb") as handle:
      for block in iter(lambda: handle.read(1 << 20), b""):
         digest.update(block)
   return digest.hexdigest()


def utc_now():
   return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def relative(path):
   path = Path(path)
   try:
      return str(path.resolve().relative_to(REPO_ROOT))
   except ValueError:
      return str(path)


# ==============================================================================
# Release manifests
# ==============================================================================

def read_manifest(path):
   """Read a release manifest: ``key = value`` header lines and five-column entries."""

   path = Path(path)
   header = {}
   entries = []
   for number, raw in enumerate(path.read_text().splitlines(), start=1):
      line = raw.split("#", 1)[0].strip()
      if line == "":
         continue
      if "=" in line and len(line.split()) <= 3 and line.split()[1:2] == ["="]:
         name, value = (part.strip() for part in line.split("=", 1))
         if name in header:
            refuse(f"{path}:{number}: {name!r} is set twice")
         header[name] = value
         continue
      fields = line.split()
      if len(fields) != 5:
         refuse(f"{path}:{number}: expected 'role logical-name file-version sha256 path', got {raw!r}")
      role, logical, file_version, sha256, file = fields
      if role not in ROLES:
         refuse(f"{path}:{number}: unknown role {role!r}")
      if not re.fullmatch(r"[0-9a-f]{64}", sha256):
         refuse(f"{path}:{number}: malformed SHA-256 {sha256!r}")
      stem = Path(file).stem
      suffix = VERSIONED_STEM.search(stem)
      if file_version == "pinned":
         if suffix:
            refuse(f"{path}:{number}: {file} has a version suffix but is listed as pinned")
      elif not re.fullmatch(r"v\d\d", file_version) or not suffix or suffix.group(0)[1:] != file_version:
         refuse(f"{path}:{number}: file version {file_version!r} does not match the file name {file}")
      if role in PYTHON_ROLES and not (suffix and file.endswith(".py")):
         refuse(f"{path}:{number}: Python source {file} must be a versioned *_vNN.py file")
      entries.append({"role": role, "logical": logical, "file_version": file_version, "sha256": sha256, "path": file})
   if header.get("manifest_schema") != MANIFEST_SCHEMA:
      refuse(f"{path}: manifest_schema must be {MANIFEST_SCHEMA!r}, got {header.get('manifest_schema')!r}")
   try:
      header["lut_release"] = int(header["lut_release"])
   except (KeyError, ValueError):
      refuse(f"{path}: lut_release must be an integer")
   for key in ("logical", "path"):
      seen = [entry[key] for entry in entries]
      duplicates = sorted({item for item in seen if seen.count(item) > 1})
      if duplicates:
         refuse(f"{path}: duplicate {key} entries: {', '.join(duplicates)}")
   drivers = [entry for entry in entries if entry["role"] == "driver"]
   if len(drivers) != 1:
      refuse(f"{path}: exactly one driver entry is required, found {len(drivers)}")
   return {"path": path, "sha256": sha256_file(path), "header": header, "entries": entries}


def find_release_manifest(lut_release):
   """The manifest that defines LUT release ``lut_release`` (ORAC_LUT_RELEASE_MANIFEST, or the only one)."""

   chosen = os.environ.get("ORAC_LUT_RELEASE_MANIFEST", "").strip()
   if chosen:
      path = Path(chosen)
      if not path.is_absolute():
         path = RELEASE_DIR / chosen if (RELEASE_DIR / chosen).is_file() else REPO_ROOT / chosen
      if not path.is_file():
         refuse(f"ORAC_LUT_RELEASE_MANIFEST names {chosen}, which is not a file")
      manifest = read_manifest(path)
      if manifest["header"]["lut_release"] != lut_release:
         refuse(f"{path} defines V{manifest['header']['lut_release']}, but the run file asks for V{lut_release}")
      if not re.fullmatch(rf"lut_v{lut_release:02d}(_[A-Za-z0-9.-]+)?\.release", path.name):
         refuse(f"{path} defines V{lut_release} but should be named lut_v{lut_release:02d}[_<label>].release")
      return manifest
   matches = []
   for path in sorted(RELEASE_DIR.glob("lut_v*.release")):
      manifest = read_manifest(path)
      if manifest["header"]["lut_release"] == lut_release:
         matches.append(manifest)
   if not matches:
      refuse(f"no release manifest defines LUT release V{lut_release} in {RELEASE_DIR} "
             "(releases made before source-file versioning are recovered from their Git tags; see SOURCE_VERSIONING.md)")
   if len(matches) > 1:
      refuse(f"several manifests define LUT release V{lut_release}: " + ", ".join(str(m["path"]) for m in matches)
             + "; choose one with ORAC_LUT_RELEASE_MANIFEST=<file name in releases/>")
   manifest = matches[0]
   if not re.fullmatch(rf"lut_v{lut_release:02d}(_[A-Za-z0-9.-]+)?\.release", manifest["path"].name):
      refuse(f"{manifest['path']} defines V{lut_release} but should be named lut_v{lut_release:02d}[_<label>].release")
   return manifest


CONFIG_SUBSTITUTIONS = {}


def sourcevault_copy_with_hash(path, digest):
   """A same-directory SourceVault copy (name_vNN or stem_vNN.suffix) of ``path`` with SHA-256 ``digest``."""

   pattern = f"{path.stem}_v[0-9][0-9]{path.suffix}"
   for candidate in sorted(path.parent.glob(pattern)):
      if candidate.is_file() and sha256_file(candidate) == digest:
         return candidate
   return None


def verify_manifest_files(manifest, roles):
   """Refuse unless every listed file of the given roles exists with the listed SHA-256."""

   problems = []
   for entry in manifest["entries"]:
      if entry["role"] not in roles:
         continue
      path = REPO_ROOT / entry["path"]
      if not path.is_file():
         problems.append(f"MISSING  {entry['path']}")
      elif sha256_file(path) != entry["sha256"]:
         copy = sourcevault_copy_with_hash(path, entry["sha256"]) if entry["role"] == "config" else None
         if copy is not None:
            CONFIG_SUBSTITUTIONS[str(path)] = str(copy)
            continue
         problems.append(f"MODIFIED {entry['path']} (SHA-256 {sha256_file(path)[:16]}..., manifest {entry['sha256'][:16]}...)")
   if problems:
      refuse(f"LUT release V{manifest['header']['lut_release']} is sealed by {relative(manifest['path'])} but its files differ:\n   "
             + "\n   ".join(problems)
             + "\nA changed file must become a new file version (file_vNN.py) in a new LUT release; sealed files are never edited.")


# ==============================================================================
# Run file
# ==============================================================================

def read_run_file_version(run_file):
   """The run file's ``version`` (same comment and assignment rules as the driver's read_runfile)."""

   version = None
   for raw in Path(run_file).read_text().splitlines():
      line = raw.split(";", 1)[0].split("#", 1)[0].strip()
      if "=" not in line:
         continue
      name, value = (part.strip() for part in line.split("=", 1))
      if name.lower() == "version":
         try:
            version = ast.literal_eval(value)
         except (ValueError, SyntaxError):
            refuse(f"{run_file}: cannot read version = {value!r}")
   if not isinstance(version, int) or isinstance(version, bool):
      refuse(f"{run_file}: the run file must set an integer version (the LUT release), got {version!r}")
   return version


# ==============================================================================
# Loading modules from the release's files
# ==============================================================================

class ReleaseLoader(importlib.abc.Loader):
   """Execute exactly the verified bytes of one manifest file (no .pyc cache)."""

   def __init__(self, entry, data):
      self.entry = entry
      self.data = data
      self.path = str(REPO_ROOT / entry["path"])

   def create_module(self, spec):
      return None

   def exec_module(self, module):
      module.__cached__ = None
      code = compile(self.data, self.path, "exec", dont_inherit=True)
      exec(code, module.__dict__)

   def get_source(self, name):
      return self.data.decode("utf-8")

   def is_package(self, name):
      return self.entry["role"] == "package"


class ReleaseFinder(importlib.abc.MetaPathFinder):
   """Map each logical module name of the release to the file its manifest names."""

   def __init__(self, manifest):
      self.manifest = manifest
      self.by_name = {entry["logical"]: entry for entry in manifest["entries"] if entry["role"] in PYTHON_ROLES}
      self.loaded = {}

   def find_spec(self, name, path=None, target=None):
      top = name.split(".", 1)[0]
      if name in self.by_name:
         entry = self.by_name[name]
         data = (REPO_ROOT / entry["path"]).read_bytes()
         if sha256_bytes(data) != entry["sha256"]:
            raise ImportError(f"{entry['path']} changed after the release was verified (SHA-256 mismatch)")
         package = entry["role"] == "package"
         locations = [str((REPO_ROOT / entry["path"]).parent)] if package else None
         spec = importlib.util.spec_from_file_location(name, REPO_ROOT / entry["path"], loader=ReleaseLoader(entry, data),
                                                       submodule_search_locations=locations)
         spec.cached = None
         self.loaded[name] = entry
         return spec
      if top in GUARDED_TOP_NAMES or top.startswith("create_orac_luts"):
         release = self.manifest["header"]["lut_release"]
         if VERSIONED_STEM.search(name.rsplit(".", 1)[-1]):
            raise ImportError(f"{name}: versioned files are loaded only under their logical names "
                              f"(LUT release V{release}, {relative(self.manifest['path'])})")
         raise ImportError(f"{name} is not part of LUT release V{release} ({relative(self.manifest['path'])})")
      return None


def check_not_preimported():
   for name in sys.modules:
      if name.split(".", 1)[0] in GUARDED_TOP_NAMES:
         refuse(f"module {name} was imported before the release finder was installed; "
                f"run {PROGRAM} as the program, not from a Python session that already imported oraclut")


def verify_bindings(finder):
   """Every oraclut / driver module in this process must come from the release's own files."""

   problems = []
   for name, module in list(sys.modules.items()):
      if name.split(".", 1)[0] not in GUARDED_TOP_NAMES:
         continue
      entry = finder.by_name.get(name)
      loader = getattr(module, "__loader__", None)
      if entry is None or not isinstance(loader, ReleaseLoader) or loader.entry is not entry:
         problems.append(f"{name} ({getattr(module, '__file__', None)})")
      elif Path(module.__file__).resolve() != (REPO_ROOT / entry["path"]).resolve():
         problems.append(f"{name} ({module.__file__} instead of {entry['path']})")
   if problems:
      refuse("modules not loaded from the release manifest: " + ", ".join(problems))


# ==============================================================================
# Compiled kernel libraries
# ==============================================================================

def verify_libraries(manifest):
   """Check (and if missing, build privately) the Mie and DISORT libraries the manifest pins."""

   checkout = manifest["header"].get("library_checkout")
   if checkout != str(REPO_ROOT):
      refuse(f"{relative(manifest['path'])} pins kernel libraries built in {checkout}, but this checkout is {REPO_ROOT}; "
             "the DISORT library embeds its absolute source paths, so the pins must be re-established for a new location "
             "(see SOURCE_VERSIONING.md)")
   wrappers = {"mie-library": sys.modules["oraclut.optics.legacy_mie"],
               "disort-library": sys.modules["oraclut.radiative_transfer.legacy_disort"]}
   disort_module = sys.modules["oraclut.radiative_transfer.legacy_disort"]
   declared = manifest["header"].get("disort_precision")
   if declared is not None:
      if not hasattr(disort_module, "disort_precision"):
         refuse(f"{relative(manifest['path'])} declares disort_precision = {declared}, but its legacy_disort has no precision selection")
      selected = disort_module.disort_precision()
      if selected != declared:
         refuse(f"{relative(manifest['path'])} is a {declared}-precision release, but ORAC_LUT_DISORT_PRECISION selects {selected}")
   records = []
   for entry in manifest["entries"]:
      if entry["role"] != "library":
         continue
      if entry["logical"] == "disort-double-library":
         path = REPO_ROOT / entry["path"]
         if getattr(disort_module, "DOUBLE_LIBRARY", None) is None or disort_module.DOUBLE_LIBRARY.resolve() != path.resolve():
            refuse(f"the double-precision DISORT library would not be loaded from the pinned {entry['path']}")
         if not path.is_file():
            refuse(f"{entry['path']} is missing; build it with create_orac_lut/disort2_dp/build_dlm.sh")
         digest = sha256_file(path)
         if digest != entry["sha256"]:
            refuse(f"{entry['path']} has SHA-256 {digest[:16]}..., not the pinned {entry['sha256'][:16]}...")
         status = path.stat()
         records.append({"logical": entry["logical"], "path": entry["path"], "sha256": digest, "built_by_launcher": False,
                         "device": status.st_dev, "inode": status.st_ino})
         continue
      wrapper = wrappers.get(entry["logical"])
      if wrapper is None:
         refuse(f"unknown library entry {entry['logical']!r}")
      path = REPO_ROOT / entry["path"]
      if wrapper._library_path().resolve() != path.resolve():
         refuse(f"the {entry['logical']} would be loaded from {wrapper._library_path()} (ORACLUT_FORTRAN / ORACLUT_FORTRAN_TAG), "
                f"not the pinned {entry['path']}")
      built_here = False
      if not path.is_file():
         # Build in a private directory (private object files), check, then install atomically.
         private = Path(tempfile.mkdtemp(prefix=".release-build-", dir=path.parent))
         try:
            candidate = private / path.name
            wrapper._build_library(candidate)
            if sha256_file(candidate) != entry["sha256"]:
               refuse(f"a fresh build of {entry['path']} has SHA-256 {sha256_file(candidate)[:16]}..., not the pinned "
                      f"{entry['sha256'][:16]}... (different compiler or flags)")
            os.replace(candidate, path)
            built_here = True
         finally:
            shutil.rmtree(private, ignore_errors=True)
      digest = sha256_file(path)
      if digest != entry["sha256"]:
         refuse(f"{entry['path']} has SHA-256 {digest[:16]}..., not the pinned {entry['sha256'][:16]}...; it was built from other "
                "sources or with another toolchain.  Delete it and rerun to rebuild from the pinned sources.")
      status = path.stat()
      records.append({"logical": entry["logical"], "path": entry["path"], "sha256": digest, "built_by_launcher": built_here,
                      "device": status.st_dev, "inode": status.st_ino})
   return records


def mapped_shared_objects():
   """Paths of the shared objects mapped into this process, with device/inode as in /proc/self/maps."""

   mapped = {}
   try:
      with open("/proc/self/maps") as handle:
         for line in handle:
            fields = line.split()
            if len(fields) >= 6 and fields[5].startswith("/") and ".so" in fields[5]:
               mapped[fields[5]] = fields[4]
   except OSError:
      pass
   return mapped


# ==============================================================================
# Inputs actually read
# ==============================================================================

INPUTS = {"files": {}, "other": set(), "violations": [], "busy": False, "excluded": ()}


def excluded_roots():
   roots = {sys.prefix, sys.base_prefix, sys.exec_prefix, tempfile.gettempdir(), "/proc", "/sys", "/dev", "/etc", "/usr", "/lib", "/lib64"}
   if os.environ.get("TMPDIR"):
      roots.add(os.environ["TMPDIR"])
   return tuple(sorted({os.path.realpath(root) for root in roots if root}))


def record_input(path, via):
   """Record one file opened for reading: path, real path, SHA-256 at the time of opening, identity."""

   if INPUTS["busy"]:
      return
   INPUTS["busy"] = True
   try:
      path = os.fspath(path)
      if isinstance(path, bytes):
         path = os.fsdecode(path)
      absolute = os.path.abspath(path)
      real = os.path.realpath(absolute)
      if any(real == root or real.startswith(root + os.sep) for root in INPUTS["excluded"]):
         return
      if not os.path.isfile(real):
         return
      inside = real.startswith(str(REPO_ROOT) + os.sep)
      if inside and (("/__pycache__/" in real) or (real.startswith(str(REPO_ROOT / "src") + os.sep) and real.endswith(".py"))
                     or os.path.basename(real).startswith("create_orac_luts")):
         INPUTS["violations"].append(real)
         return
      if real in INPUTS["files"]:
         return
      status = os.stat(real)
      INPUTS["files"][real] = {"path": relative(absolute), "realpath": real, "sha256": sha256_file(real), "size": status.st_size,
                               "device": status.st_dev, "inode": status.st_ino, "mtime_ns": status.st_mtime_ns, "via": via}
   finally:
      INPUTS["busy"] = False


def audit_hook(event, args):
   if event != "open":
      return
   path, mode, flags = args
   if not isinstance(path, (str, bytes, os.PathLike)):
      return
   if mode is not None:
      if any(character in mode for character in "wax+"):
         return
   elif flags is not None and flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT):
      return
   record_input(path, "open")


def install_netcdf_recording():
   """Record NetCDF files opened for reading (the netCDF C library bypasses Python's open)."""

   import netCDF4

   class RecordingDataset(netCDF4.Dataset):
      def __init__(self, filename, mode="r", *args, **kwargs):
         if mode == "r":
            record_input(filename, "netCDF4.Dataset")
         super().__init__(filename, mode, *args, **kwargs)

   netCDF4.Dataset = RecordingDataset


def recheck_inputs():
   """Refuse if any recorded input changed while the calculation ran."""

   changed = []
   for real, item in INPUTS["files"].items():
      try:
         status = os.stat(real)
      except OSError:
         changed.append(f"{item['path']} (disappeared)")
         continue
      if (status.st_size, status.st_ino, status.st_mtime_ns) != (item["size"], item["inode"], item["mtime_ns"]):
         changed.append(item["path"])
   if changed:
      refuse("input files changed during the calculation: " + ", ".join(changed))
   if INPUTS["violations"]:
      refuse("the calculation read source files outside the release binding: " + ", ".join(sorted(set(INPUTS["violations"]))))


# ==============================================================================
# Git and runtime environment
# ==============================================================================

def git_executable():
   return os.environ.get("ORACLUT_GIT") or shutil.which("git") or "/usr/bin/git"


def git(*args):
   environment = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
   environment["LC_ALL"] = "C"
   completed = subprocess.run([git_executable(), "-C", str(REPO_ROOT), *args], capture_output=True, text=True,
                              timeout=120, env=environment)
   if completed.returncode != 0:
      raise OSError(completed.stderr.strip() or completed.stdout.strip() or f"git {' '.join(args)} failed")
   return completed.stdout.rstrip("\n")


def git_state(manifest, launcher_paths):
   """Commit, tree status and, for every release file, whether it is in that commit unchanged."""

   tracked_paths = [entry["path"] for entry in manifest["entries"] if entry["role"] != "library"] + launcher_paths
   try:
      commit = git("rev-parse", "HEAD")
      state = {"available": True, "commit": commit, "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
               "describe": git("describe", "--tags", "--long", "--always"),
               "status_porcelain": git("status", "--porcelain", "--untracked-files=all").splitlines()}
      try:
         state["origin"] = git("remote", "get-url", "origin")
      except OSError:
         state["origin"] = None
      committed = {}
      for line in git("ls-tree", "-r", "HEAD", "--", *tracked_paths).splitlines():
         meta, path = line.split("\t", 1)
         committed[path] = meta.split()[2]
      current = git("hash-object", "--", *tracked_paths).splitlines()
      state["files_in_commit"] = {path: committed.get(path) == blob for path, blob in zip(tracked_paths, current)}
      state["release_files_not_in_commit"] = sorted(path for path, ok in state["files_in_commit"].items() if not ok)
      state["working_tree"] = "MODIFIED" if state["status_porcelain"] else "CLEAN"
      return state
   except (OSError, subprocess.TimeoutExpired, ValueError, IndexError) as exc:
      return {"available": False, "detail": str(exc).splitlines()[0] if str(exc) else type(exc).__name__}


def first_line(command):
   try:
      completed = subprocess.run(command, capture_output=True, text=True, timeout=30)
      return (completed.stdout or completed.stderr).splitlines()[0]
   except (OSError, subprocess.TimeoutExpired, IndexError):
      return "unavailable"


def runtime_environment():
   """Interpreter, libraries, CPU and environment that can affect the numbers."""

   record = {"python": sys.version, "executable": sys.executable, "flags": repr(sys.flags), "sys_path": list(sys.path),
             "platform": platform.platform(), "hostname": socket.gethostname(), "cwd": os.getcwd()}
   try:
      with open("/proc/cpuinfo") as handle:
         record["cpu_model"] = next((line.split(":", 1)[1].strip() for line in handle if line.startswith("model name")), None)
   except OSError:
      record["cpu_model"] = None
   import numpy
   import netCDF4
   record["numpy"] = {"version": numpy.__version__, "file": numpy.__file__}
   try:
      from numpy._core._multiarray_umath import __cpu_features__
      record["numpy"]["cpu_features_enabled"] = sorted(name for name, enabled in __cpu_features__.items() if enabled)
   except ImportError:
      record["numpy"]["cpu_features_enabled"] = None
   record["netCDF4"] = {"version": netCDF4.__version__, "file": netCDF4.__file__,
                        "libnetcdf": netCDF4.__netcdf4libversion__, "hdf5": netCDF4.__hdf5libversion__}
   conda = {}
   for package in ("python", "numpy", "netcdf4", "libnetcdf", "hdf5", "libgfortran5", "libopenblas", "libgcc"):
      for meta in sorted(Path(sys.prefix, "conda-meta").glob(f"{package}-[0-9]*.json")):
         conda[package] = meta.stem
   record["conda_packages"] = conda
   compiler = os.environ.get("ORACLUT_FORTRAN", "gfortran")
   record["compilers"] = {"fortran": compiler, "fortran_path": shutil.which(compiler), "fortran_version": first_line([compiler, "--version"]),
                          "gcc_path": shutil.which("gcc"), "gcc_version": first_line(["gcc", "--version"])}
   prefixes = ("NPY_", "OPENBLAS_", "OMP_", "MKL_", "LD_", "ORACLUT_", "PYTHON", "SLURM_JOB_ID", "SLURM_JOB_NAME", "SLURM_ARRAY")
   record["environment"] = {key: value for key, value in sorted(os.environ.items()) if key.startswith(prefixes)}
   return record


# ==============================================================================
# NetCDF provenance attributes
# ==============================================================================

def source_table(manifest, finder, git_start):
   """The release's files: compact text and machine-readable entries."""

   in_commit = git_start.get("files_in_commit", {}) if git_start.get("available") else {}
   lines = ["file-version  sha256                                                            path   (role, logical name)"]
   entries = []
   for entry in manifest["entries"]:
      lines.append(f"{entry['file_version']:<12}  {entry['sha256']}  {entry['path']}   ({entry['role']}, {entry['logical']})")
      item = dict(entry)
      if entry["role"] in PYTHON_ROLES:
         item["loaded"] = entry["logical"] in finder.loaded
      if entry["role"] != "library":
         item["in_git_commit"] = in_commit.get(entry["path"]) if in_commit else None
      entries.append(item)
   return "\n".join(lines), entries


def write_attributes(netcdf_file, attributes):
   """Add global attributes (text as NC_STRING, via the release's own writer helper) and read them back."""

   from netCDF4 import Dataset
   put_attribute = sys.modules["oraclut.io.v2"].put_attribute
   with Dataset(netcdf_file, "r+") as dataset:
      for name in attributes:
         if name in dataset.ncattrs():
            refuse(f"{netcdf_file}: global attribute {name!r} already exists")
      for name, value in attributes.items():
         put_attribute(dataset, name, value)
   sys.modules["oraclut.io.string_attributes"].check_text_attributes_are_strings(netcdf_file)
   netcdf_c = sys.modules["oraclut.io.netcdf_c"]
   with netcdf_c.open_file(netcdf_file) as ncid:
      for name, value in attributes.items():
         if not isinstance(value, str):
            continue
         xtype, length = netcdf_c.inq_att(ncid, netcdf_c.NC_GLOBAL, name)
         stored = netcdf_c.get_att_string(ncid, netcdf_c.NC_GLOBAL, name, length)
         if xtype != netcdf_c.NC_STRING or stored != [value]:
            refuse(f"{netcdf_file}: global attribute {name!r} did not read back as the NC_STRING that was written")


def build_attributes(release, finished):
   """Every provenance attribute; compact human-readable ones first, the complete JSON record last."""

   import numpy

   manifest = release["manifest"]
   lut_release = manifest["header"]["lut_release"]
   driver = next(entry for entry in manifest["entries"] if entry["role"] == "driver")
   git_start = release["git_start"]
   git_end = finished["git_end"]
   table, entries = source_table(manifest, release["finder"], git_start)
   commit = git_start.get("commit", "unavailable") if git_start.get("available") else "unavailable"
   if git_start.get("available"):
      outside = git_start["release_files_not_in_commit"]
      git_summary = (f"working tree {git_start['working_tree']} at {commit}; "
                     + ("every release file is in this commit unchanged" if not outside else
                        f"{len(outside)} release file(s) are not in this commit unchanged: {', '.join(outside)}"))
   else:
      git_summary = f"Git metadata unavailable: {git_start.get('detail')}"
   if git_end.get("available") and git_end.get("commit") != git_start.get("commit"):
      git_summary += f"; HEAD moved to {git_end.get('commit')} during the calculation (the start state applies)"
   record = {
      "provenance_schema": PROVENANCE_SCHEMA,
      "lut_release": lut_release,
      "versioning_note": "LUT release numbers (V21, V22, ...) and source-file versions (_v01, _v02, ...) are independent; "
                         "each file's version advances only when that file changes.",
      "generated_utc": finished["generated_utc"],
      "started_utc": release["started_utc"],
      "generation_command": release["command"],
      "launcher": release["launcher"],
      "release_manifest": {"path": relative(manifest["path"]), "sha256": manifest["sha256"], "header": manifest["header"]},
      "driver": {"logical": driver["logical"], "path": driver["path"], "file_version": driver["file_version"], "sha256": driver["sha256"]},
      "source_files": entries,
      "compiled_libraries": release["libraries"],
      "compiled_libraries_mapped": finished["libraries_mapped"],
      "run_file": release["run_file"],
      "inputs_read": sorted(INPUTS["files"].values(), key=lambda item: item["path"]),
      "repository": {"public": PUBLIC_REPOSITORY, "origin": git_start.get("origin"), "checkout": str(REPO_ROOT)},
      "git_at_start": git_start,
      "git_at_end": git_end,
      "runtime": release["runtime"],
      "config_substitutions": release.get("config_substitutions", {}),
      "shared_objects_mapped": finished["shared_objects"],
   }
   return {
      "lut_release": f"V{lut_release}",
      "lut_release_version": numpy.int32(lut_release),
      "provenance_schema": PROVENANCE_SCHEMA,
      "provenance_note": "LUT release and source-file versions are independent: the files listed in provenance_source_files, "
                         "with their own versions and SHA-256, generated this LUT release.  Complete record: provenance_json.",
      "provenance_generated_utc": finished["generated_utc"],
      "provenance_generation_driver": f"{driver['path']} (logical name {driver['logical']}, file version {driver['file_version']}, "
                                      f"sha256 {driver['sha256']})",
      "provenance_launcher": f"{release['launcher']['path']} sha256 {release['launcher']['sha256']}",
      "provenance_release_manifest": f"{relative(manifest['path'])} sha256 {manifest['sha256']}",
      "provenance_generation_command": release["command"],
      "provenance_run_file": f"{release['run_file']['path']} sha256 {release['run_file']['sha256']}",
      "provenance_repository": f"{PUBLIC_REPOSITORY} (origin {git_start.get('origin')}; checkout {REPO_ROOT})",
      "provenance_git_commit": commit,
      "provenance_git_state": git_summary,
      "provenance_source_files": table,
      "provenance_json": json.dumps(record, indent=1, sort_keys=True, default=str),
   }


def publish_with_provenance(release, private_file, final_file):
   """Replacement for the driver's _publish_lut: embed the provenance, then publish as before."""

   try:
      # First judge what the calculation read, then stop recording: the
      # launcher's own checks below read files that are not inputs.
      recheck_inputs()
      INPUTS["busy"] = True
      verify_bindings(release["finder"])
      if sha256_file(release["run_file"]["realpath"]) != release["run_file"]["sha256"]:
         refuse(f"the run file {release['run_file']['path']} changed during the calculation")
      verify_manifest_files(release["manifest"], ("driver", "package", "module", "config", "adapter", "kernel", "library"))
      mapped = mapped_shared_objects()
      libraries_mapped = []
      for library in release["libraries"]:
         path = str((REPO_ROOT / library["path"]).resolve())
         libraries_mapped.append({"path": library["path"], "mapped": path in mapped,
                                  "maps_inode_matches": (path in mapped and str(library["inode"]) == mapped[path])})
         if path in mapped and str(library["inode"]) != mapped[path]:
            refuse(f"{library['path']} mapped in this process is not the file whose SHA-256 was checked")
      shared_objects = {}
      for path in sorted(mapped):
         try:
            shared_objects[path] = sha256_file(path)
         except OSError:
            shared_objects[path] = None
      finished = {"generated_utc": utc_now(), "git_end": git_state(release["manifest"], release["launcher_paths"]),
                  "libraries_mapped": libraries_mapped, "shared_objects": shared_objects}
      write_attributes(private_file, build_attributes(release, finished))
   except Exception as exc:
      # Keep the computed product so that the provenance can be completed without recomputing.
      keep = Path(final_file).with_name(Path(final_file).name + ".unpublished")
      try:
         shutil.copyfile(private_file, keep)
         where = f"; the computed product was kept as {keep}"
      except OSError:
         where = ""
      refuse(f"provenance could not be embedded ({exc}){where}")
   finally:
      INPUTS["busy"] = False
   release["original_publish"](private_file, final_file)
   print(f"Release provenance embedded: LUT release V{release['manifest']['header']['lut_release']}, "
         f"{relative(release['manifest']['path'])}, driver {release['manifest_driver']}", flush=True)


# ==============================================================================
# Main
# ==============================================================================

def substitute_config_files():
   """Serve sealed CODE_VERSION content from its SourceVault copy when the working file was edited later."""

   code_version = str(REPO_ROOT / "CODE_VERSION")
   if code_version not in CONFIG_SUBSTITUTIONS:
      return
   version_module = sys.modules.get("oraclut.version")
   if version_module is None:
      refuse("CODE_VERSION must be substituted but oraclut.version is not loaded")
   copy = Path(CONFIG_SUBSTITUTIONS[code_version])
   original = version_module.read_code_version

   def read_sealed_code_version(path=version_module.CODE_VERSION_FILE):
      if Path(path).resolve() == Path(code_version).resolve():
         return original(copy)
      return original(path)

   version_module.read_code_version = read_sealed_code_version
   print(f"Sealed CODE_VERSION served from {relative(copy)} (the working CODE_VERSION belongs to a later release)")


def check_release(lut_release):
   manifest = find_release_manifest(lut_release)
   verify_manifest_files(manifest, ROLES)
   print(f"LUT release V{lut_release}: {relative(manifest['path'])} (sha256 {manifest['sha256']})")
   for entry in manifest["entries"]:
      served = CONFIG_SUBSTITUTIONS.get(str(REPO_ROOT / entry["path"]))
      note = f"  (served from {relative(Path(served))})" if served else ""
      print(f"   OK  {entry['file_version']:<7} {entry['sha256'][:16]}  {entry['path']}{note}")
   print(f"all {len(manifest['entries'])} files match the manifest")
   return 0


def run_release(run_file_argument, argv):
   invoked_from = os.getcwd()
   run_file = Path(run_file_argument)
   if not run_file.is_absolute():
      run_file = Path(invoked_from) / run_file
   if not run_file.is_file():
      refuse(f"run file not found: {run_file_argument}")
   check_not_preimported()
   os.chdir(REPO_ROOT)
   started = utc_now()
   run_bytes = run_file.read_bytes()
   lut_release = read_run_file_version(run_file)
   manifest = find_release_manifest(lut_release)
   verify_manifest_files(manifest, ("driver", "package", "module", "config", "adapter", "kernel"))
   driver_entry = next(entry for entry in manifest["entries"] if entry["role"] == "driver")
   launcher_paths = [relative(LAUNCHER_FILE), relative(manifest["path"])]
   git_start = git_state(manifest, launcher_paths)

   install_netcdf_recording()
   finder = ReleaseFinder(manifest)
   sys.meta_path.insert(0, finder)
   driver = importlib.import_module(driver_entry["logical"])
   verify_bindings(finder)
   substitute_config_files()

   # The release's own run-file reader, version and run-option checks.
   settings = driver.read_runfile(run_file)
   if int(settings["version"]) != lut_release:
      refuse(f"{run_file}: version {settings['version']} does not select LUT release V{lut_release}")
   for option in ("reuse_scat", "scat_only"):
      if settings[option]:
         refuse(f"{run_file}: {option} = {settings[option]} is not allowed in a release run "
                "(cached optics carry no provenance; scat_only publishes no product)")
   out_path = Path(settings["out_path"])
   out_path = out_path if out_path.is_absolute() else REPO_ROOT / out_path
   validation_root = REPO_ROOT / "validation"
   is_validation_output = str(out_path.resolve()).startswith(str(validation_root.resolve()) + os.sep)
   if not is_validation_output:
      if not git_start.get("available"):
         refuse("products outside validation/ require Git metadata to show that the release's files are committed")
      if git_start["release_files_not_in_commit"]:
         refuse("products outside validation/ require every release file (and the launcher) to be committed unchanged; "
                "not committed: " + ", ".join(git_start["release_files_not_in_commit"]))

   libraries = verify_libraries(manifest)
   release = {
      "manifest": manifest, "finder": finder, "started_utc": started, "git_start": git_start, "libraries": libraries,
      "launcher_paths": launcher_paths, "manifest_driver": driver_entry["path"],
      "launcher": {"path": relative(LAUNCHER_FILE), "sha256": sha256_file(LAUNCHER_FILE), "manifest_schema": MANIFEST_SCHEMA},
      "command": " ".join([Path(sys.executable).name, relative(LAUNCHER_FILE), *argv]),
      "run_file": {"path": relative(run_file), "realpath": str(run_file.resolve()), "sha256": sha256_bytes(run_bytes),
                   "text": run_bytes.decode("utf-8", "replace"), "invoked_from": invoked_from},
      "runtime": runtime_environment(),
      "config_substitutions": {relative(Path(k)): relative(Path(v)) for k, v in CONFIG_SUBSTITUTIONS.items()},
      "original_publish": driver._publish_lut,
   }

   # A trial write of the provenance attributes before hours of calculation.
   trial_dir = Path(tempfile.mkdtemp(prefix="release-trial-"))
   try:
      from netCDF4 import Dataset
      trial = trial_dir / "trial.nc"
      with Dataset(trial, "w", format="NETCDF4"):
         pass
      trial_finished = {"generated_utc": utc_now(), "git_end": git_start, "libraries_mapped": [], "shared_objects": {}}
      write_attributes(trial, build_attributes(release, trial_finished))
   finally:
      shutil.rmtree(trial_dir, ignore_errors=True)

   print("ORAC LUT release")
   print("-" * 60)
   print(f"LUT release:       V{lut_release} ({relative(manifest['path'])}, sha256 {manifest['sha256'][:16]}...)")
   print(f"Driver:            {driver_entry['path']} (file version {driver_entry['file_version']})")
   print(f"Release files:     {len(manifest['entries'])} verified; {len(finder.loaded)} Python modules loaded from them")
   for library in libraries:
      print(f"Library:           {library['path']} sha256 {library['sha256'][:16]}..."
            + (" (built by this launcher)" if library["built_by_launcher"] else ""))
   if git_start.get("available"):
      print(f"Git:               {git_start['commit'][:12]} {git_start['working_tree']}; "
            f"{len(git_start['release_files_not_in_commit'])} release file(s) not in this commit")
   else:
      print(f"Git:               unavailable ({git_start.get('detail')})")
   print("-" * 60, flush=True)

   driver._publish_lut = lambda private_file, final_file: publish_with_provenance(release, private_file, final_file)
   INPUTS["excluded"] = excluded_roots()
   record_input(run_file, "run file")
   sys.addaudithook(audit_hook)
   return driver.main([str(run_file)])


def main(argv=None):
   argv = sys.argv[1:] if argv is None else list(argv)
   try:
      if len(argv) == 2 and argv[0] == "--check":
         return check_release(int(argv[1]))
      if len(argv) != 1 or argv[0].startswith("-"):
         print(f"usage: python {PROGRAM} runs/<product>.run")
         print(f"       python {PROGRAM} --check <lut-release>")
         print("       (ORAC_LUT_RELEASE_MANIFEST=<file in releases/> chooses among manifests of the same release)")
         return 2
      return run_release(argv[0], argv)
   except ReleaseRefusal as exc:
      print(f"{PROGRAM}: refused: {exc}", file=sys.stderr)
      return 2


if __name__ == "__main__":
   sys.exit(main())
