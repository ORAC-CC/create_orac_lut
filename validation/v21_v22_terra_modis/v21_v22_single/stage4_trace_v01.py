"""Stage 4: stage-by-stage trace of the scattering-property arithmetic, uncorrected Python vs IDL.

For one effective radius of one compact case (default: case A, channel 6,
r_e = 29 um, the failing p240 node) the two arithmetic chains of
scattering_replica.py are evaluated side by side:

   "python": the uncorrected V22 Python arithmetic (reproduces the 17de4c0
             Python scatfile bitwise);
   "idl":    the IDL arithmetic (reproduces the IDL scatfile.sav bitwise in all
             four compact cases).

Only the final scatfile arrays were captured inside IDL; the intermediate
"idl" arrays are those of the replica that reproduces them, i.e. inferred IDL
intermediates.  For every boundary: dtype, shape, elements differing, max
|difference|, max ulp (in the array's own precision) and the first differing
element.  The Mie kernels are also compared with identical (IDL) inputs.

Usage: python stage4_trace.py [case] [radius_index]   -> results/stage4_trace_<case>_r<i>.json
"""

import json
import sys

import numpy as np

import scattering_replica as sr


def ulp_distance(a, b):
   a = np.asarray(a)
   b = np.asarray(b)
   if a.dtype == np.float32:
      ia, ib = a.ravel().view(np.int32).astype(np.int64), b.astype(np.float32).ravel().view(np.int32).astype(np.int64)
      ia = np.where(ia < 0, -(ia & 0x7FFFFFFF), ia)
      ib = np.where(ib < 0, -(ib & 0x7FFFFFFF), ib)
      return np.abs(ia - ib)
   ia = a.astype(np.float64).ravel().view(np.int64)
   ib = b.astype(np.float64).ravel().view(np.int64)
   big = np.int64(0x7FFFFFFFFFFFFFFF)
   ia = np.where(ia < 0, -(ia & big), ia).astype(object)
   ib = np.where(ib < 0, -(ib & big), ib).astype(object)
   return np.array([abs(x - y) for x, y in zip(ia, ib)], dtype=object)


def stage(name, python, idl, operation):
   python, idl = np.asarray(python), np.asarray(idl)
   row = {"stage": name, "dtype": str(idl.dtype), "shape": list(idl.shape), "operation": operation}
   if python.shape != idl.shape:
      row.update(status="shape differs", python_shape=list(python.shape))
      return row
   if np.iscomplexobj(idl):
      d = (python.real != idl.real) | (python.imag != idl.imag)
      u = max(int(np.max(ulp_distance(idl.real, python.real))), int(np.max(ulp_distance(idl.imag, python.imag))))
   else:
      d = python != idl
      u = int(np.max(ulp_distance(idl, python))) if idl.size else 0
   row.update(bitwise=not bool(np.any(d)), differing=int(np.count_nonzero(d)), max_abs=float(np.max(np.abs(python - idl))) if idl.size else 0.0,
              max_ulp=u)
   if np.any(d):
      e = int(np.argmax(d.ravel()))
      row["first"] = {"index": e, "idl": repr(idl.ravel()[e]), "python": repr(python.ravel()[e])}
   return row


def main():
   case = sys.argv[1] if len(sys.argv) > 1 else "A"
   r = int(sys.argv[2]) if len(sys.argv) > 2 else 1
   out_p, det_p, pre_p = sr.scattering(case, sr.PYTHON)
   out_i, det_i, pre_i = sr.scattering(case, sr.IDL)
   rows = [stage("refractive index at 0.55 um (AerM550)", pre_p["aerm550"], pre_i["aerm550"], "INTERPOL of COMPLEX(n,-k) (IDL) vs double interp (Python)"),
           stage("refractive index at the channel (AerM)", pre_p["aerm"], pre_i["aerm"], "as above"),
           stage("phase-function cosines QV", pre_p["qv"], pre_i["qv"], "quadrature.pro 'G' Newton iteration (IDL) vs numpy leggauss"),
           stage("Gauss-Legendre weights", pre_p["weights"], pre_i["weights"], "as above")]
   for label, k in (("0.55 um", 0), ("channel", 1)):
      p, i = det_p[r][k], det_i[r][k]
      rows += [stage(f"{label}: Npts", p["npts"], i["npts"], "long(2*pi*(Ru-Rl)*wn/xres) > 200; wn = 1.0/wl FLOAT, xres = 0.4 FLOAT"),
               stage(f"{label}: radius abscissae R", p["r"], i["r"], "quadrature 'T' + shift_quadrature with Rl = DOUBLE(0.001 FLOAT)"),
               stage(f"{label}: size distribution W1P", p["w1p"], i["w1p"], "W1*(Nd/n)*R^alpha*exp(-b*R), S FLOAT"),
               stage(f"{label}: size parameters Dx", p["dx"], i["dx"], "2*pi*R*wn"),
               stage(f"{label}: Mie Qext", p["mie"]["qext"], i["mie"]["qext"], "mieint (gfortran vs Intel DLM) on the differing Dx"),
               stage(f"{label}: Mie F11", p["mie"]["f11"], i["mie"]["f11"], "mieint"),
               stage(f"{label}: Bext (DOUBLE)", p["bext"], i["bext"], "TOTAL(W1PA*Qext)"),
               stage(f"{label}: w (DOUBLE)", p["w"], i["w"], "Bsca/Bext"),
               stage(f"{label}: g (DOUBLE)", p["g"], i["g"], "TOTAL(W1PA*g*Qsca)/Bsca"),
               stage(f"{label}: SPM[0,*] phase function (DOUBLE)", p["spm"], i["spm"], "TOTAL(W1PA*F11*Qsca)/Bsca")]
   names = (("bext550", "Bext550 (FLOAT)"), ("bext", "Bext (FLOAT)"), ("bextrat", "BextRat = Bext/Bext550 (FLOAT/FLOAT)"),
            ("w", "w (FLOAT)"), ("g", "g (FLOAT)"), ("phs550", "Phs550 (FLOAT)"), ("phs", "Phs (FLOAT)"),
            ("amom550", "AMom550 (FLOAT)"), ("amom", "AMom = legpexp(FLOAT Phs)/(2n+1) (FLOAT)"))
   for key, label in names:
      rows.append(stage(label, out_p[key][r], out_i[key][r], "generate_scattering_properties mixing, FLTARR storage, legpexp"))
   # Mie kernels with identical (IDL) inputs
   i = det_i[r][1]
   qv = pre_i["qv"]
   intel = sr.mie_intel(i["dx"], pre_i["aerm"], qv, 0, 8)
   gfortran = sr.mie_gfortran(i["dx"], pre_i["aerm"], qv)
   for key in ("qext", "qsca", "g", "f11"):
      rows.append(stage(f"Mie kernels, identical inputs: {key} (gfortran vs Intel)", gfortran[key], intel[key], "mieint.f: gfortran -O3 vs ifort -O3 build"))
   result = {"case": case, "radius_index": r, "radius": float(sr.inputs(case)[1].efr[r]), "rows": rows}
   (sr.HERE / "results" / f"stage4_trace_{case}_r{r}.json").write_text(json.dumps(result, indent=1) + "\n")
   for row in rows:
      status = "bitwise" if row.get("bitwise") else f"{row.get('differing')} differ, max ulp {row.get('max_ulp')}, max abs {row.get('max_abs'):.3g}"
      print(f"{row['stage']:62s} {row['dtype']:10s} {str(row['shape']):12s} {status}")


if __name__ == "__main__":
   main()
