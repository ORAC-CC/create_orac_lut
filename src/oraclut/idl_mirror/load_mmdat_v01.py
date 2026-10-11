"""Python counterparts of load_mmdat.pro and read_ri.pro.

IDL:  pro load_mmdat, file, mmstr
      function read_ri, Filename

The microphysical (.mm) file gives the substance, short name, the relative
optical-depth profile with height, and one or more components, each with a
size distribution, scattering code, mixing ratio and refractive-index file.
The IDL stores the per-component refractive index in sub-structures Comp1,
Comp2, ...; here they are the list ``mmstr.comp`` (mmstr.comp[c].code,
.wl, .cm), which is what the IDL tag arithmetic scatoffset + c indexes.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np


def read_ri(filename):
   """Read a refractive-index (.ri) file: '#TAG = value' header, FORMAT gives the columns.

   Returns wavl (microns), wavn (cm^-1), n, dn, k, dk as double precision
   arrays, as the IDL does (k positive, as in the file).
   """

   filename = Path(filename)
   if filename.suffix == ".gz":
      import gzip
      text = gzip.open(filename, "rt").read()
   else:
      text = filename.read_text()
   lines = text.splitlines()

   ri = SimpleNamespace(description="", distributedby="", substance="", sampleform="",
                        temperature="", concentration="", reference="", doi="", source="",
                        contact="", comment="", format="")
   header = [line for line in lines if line.startswith("#")]
   values = [line for line in lines if not line.startswith("#") and line.strip() != ""]
   last_tag = None
   for line in header:
      if line.startswith("##"):
         # continuation of the previous tag
         if last_tag is not None:
            setattr(ri, last_tag, getattr(ri, last_tag) + " " + line[2:].strip())
      else:
         i = line.index("=")
         tagname = line[1:i].strip().lower()
         tagvalue = line[i + 1:].strip()
         if hasattr(ri, tagname):
            setattr(ri, tagname, tagvalue)
            last_tag = tagname
         else:
            print("Warning in read_ri: " + tagname.upper() + " tag unknown so ignored")

   fmt = ri.format.upper().split()
   count = len(fmt)
   data = np.asarray([[float(x) for x in line.split()[:count]] for line in values], dtype=np.float64)
   nvals = data.shape[0]
   ri.vals = nvals
   ri.wavn = np.zeros(nvals)
   ri.wavl = np.zeros(nvals)
   ri.n = np.zeros(nvals)
   ri.dn = np.zeros(nvals)
   ri.k = np.zeros(nvals)
   ri.dk = np.zeros(nvals)
   for i, name in enumerate(fmt):
      if name == "WAVN":
         ri.wavn = data[:, i]
      elif name == "WAVL":
         ri.wavl = data[:, i]
      elif name == "N":
         ri.n = data[:, i]
      elif name == "DN":
         ri.dn = data[:, i]
      elif name == "K":
         ri.k = data[:, i]
      elif name == "DK":
         ri.dk = data[:, i]
      else:
         raise ValueError(f"read_ri: unknown FORMAT column {name!r} in {filename}")
   # Convert wavl to wavn or vice versa
   if ri.wavl[0] == 0:
      ri.wavl = 1e4 / ri.wavn
   if ri.wavn[0] == 0:
      ri.wavn = 1e4 / ri.wavl
   return ri


def load_mmdat(file, in_path):
   """Read the microphysical definition ``file`` and return the structure ``mmstr``.

   ``in_path`` locates ``ri/``; the IDL hard-codes 'input_files/ri/'.
   """

   in_path = Path(in_path)
   lines = Path(file).read_text().splitlines()

   ncomp = 0
   substance = ""
   shortname = ""
   description = ""
   nlayer = 0
   height = np.zeros(0, dtype=np.float32)
   rext = np.zeros(0, dtype=np.float32)
   mrat = []
   distname = []
   comptype = []
   compname = []
   compname2 = []
   rm = []
   s = []
   comp = []

   # Comment lines start with "#", data descriptors start with "*"
   i = 0
   while i < len(lines) and lines[i].startswith("#"):
      i += 1

   # Extract the first word from the non-comment lines
   while i < len(lines) and lines[i].strip().lower() != "end":
      line = lines[i]
      chunks = line.split()
      if not chunks:
         i += 1
         continue
      label = chunks[0].lower()
      if label == "component":
         ncomp += 1
         code = ""
         has_ri = False
         has_as = False
         if len(chunks) == 1:
            raise ValueError("Scattering file components must specify the component type as the second field")
         comptype1 = chunks[1]
         compname21 = ""
         if comptype1.lower() == "user":
            if len(chunks) < 3:
               raise ValueError("Scattering file component type 'user' requires a component name as its first+ field.")
            compname1 = " ".join(chunks[2:])
         elif comptype1.lower() == "opac":
            raise NotImplementedError("OPAC components (rd_optdat_p) are not ported to Python")
         elif comptype1.lower() in ("baran", "baum"):
            compname1 = chunks[2]
            if len(chunks) == 4:
               compname21 = chunks[3]
         else:
            raise ValueError("Invalid component type: " + comptype1)
         i += 1
         # Now inspect each data descriptor and read the data accordingly.
         distname1 = ""
         rm1 = 0.0
         s1 = 0.0
         mrat1 = 0.0
         wl = n = k = None
         eps = neps = None
         while i < len(lines) and lines[i].startswith("*"):
            descriptor = lines[i][1:].strip().lower()
            if descriptor == "size":
               distname1 = lines[i + 1].strip()
               rm1, s1 = (float(x) for x in lines[i + 2].split()[:2])
               i += 3
            elif descriptor == "scattering code":
               code = lines[i + 1].strip()
               i += 2
            elif descriptor == "mixing ratio":
               mrat1 = float(lines[i + 1])
               i += 2
            elif descriptor == "refractive index":
               has_ri = True
               rifilename = lines[i + 1].strip()
               ri = read_ri(in_path / "ri" / rifilename)
               wl = ri.wavl
               n = ri.n
               k = -ri.k                                # IDL: k = -ri.k (absorbing => negative)
               i += 2
            elif descriptor == "aspect ratio":
               has_as = True
               nar = int(lines[i + 1])
               rows = np.asarray([[float(x) for x in lines[i + 2 + j].split()[:2]] for j in range(nar)], dtype=np.float32)
               eps = rows[:, 0]
               neps = rows[:, 1]
               i += 2 + nar
            else:
               raise ValueError(f"Unknown component data descriptor in {file}: {lines[i]!r}")
         mrat.append(mrat1)
         distname.append(distname1)
         comptype.append(comptype1)
         compname.append(compname1)
         compname2.append(compname21)
         rm.append(rm1)
         s.append(s1)
         # IDL: tmp = {code: code} (+ WL, CM = complex(n, k)) (+ eps, neps); mmstr.Comp<n> = tmp
         tmp = SimpleNamespace(code=code)
         if has_ri:
            tmp.wl = wl
            tmp.cm = n + 1j * k                         # IDL: complex(n, k); kept double here
         if has_as:
            tmp.eps = eps
            tmp.neps = neps
         comp.append(tmp)
      elif label == "profile":
         nalt = int(lines[i + 1])
         rows = np.asarray([[float(x) for x in lines[i + 2 + j].split()[:2]] for j in range(nalt)], dtype=np.float32)
         nlayer = nalt
         height = rows[:, 0]
         rext = rows[:, 1]
         i += 2 + nalt
      elif label == "substance":
         substance = "_".join(chunks[1:])
         i += 1
      elif label == "description":
         description = " ".join(chunks[1:])
         i += 1
      elif label == "shortname":
         shortname = " ".join(chunks[1:])
         i += 1
      else:
         print("Warning: unknown label line found and skipped: " + line)
         i += 1

   if ncomp == 0:
      raise ValueError(f"Microphysical file {file} defines no component")
   if substance == "" or shortname == "":
      raise ValueError(f"Microphysical file {file} must define 'substance' and 'shortname'")

   mmstr = SimpleNamespace(
      substance=substance, shortname=shortname, description=description, ncomp=ncomp,
      nlayer=nlayer, height=height, rext=rext,
      comp=comp,
      mrat=np.asarray(mrat, dtype=np.float64), distname=distname, comptype=comptype,
      compname=compname, compname2=compname2,
      rm=np.asarray(rm, dtype=np.float64), s=np.asarray(s, dtype=np.float64),
   )
   # Deliberate difference: mrat, rm and s are double precision (the validated
   # Python reads them so); the IDL keeps them as single precision floats.
   return mmstr
