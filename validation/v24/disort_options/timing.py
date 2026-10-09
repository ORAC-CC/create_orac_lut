"""Aggregate repeated DISORT timings (validation only).

Timing runs are run_variant.py CASE EXPERIMENT REPEAT, executed in a separate,
controlled session (a fixed number of concurrent single-thread processes, no
other experiments running).  For NSTR <= 100 the experimental kernels are the
MXCMU = 100 builds, so their cost is comparable with the production library.
For NSTR > 100 only the MXCMU = 256 build can run; DISORT zeroes its
MXCMU-sized arrays on every call, a fixed per-call cost that is measured as
nstr100_k256 - nstr100 (same NSTR, same calls) and subtracted, giving
"production-equivalent" times.

   python timing.py     -> validation/v24/results/disort_options/runtime.csv
"""

import csv
import json
import statistics
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
from run_variant import CASES, EXPERIMENTS, TMP   # noqa: E402

RESULTS = ROOT / "validation" / "v24" / "results" / "disort_options"


def records(case, experiment):
   out = []
   for path in sorted((TMP / "runs" / case).glob(f"{experiment}_r*/record.json")):
      out.append(json.loads(path.read_text()))
   return out


def main():
   rows = []
   for case in CASES:
      base = records(case, "baseline")
      if not base:
         continue
      t_base = statistics.median(r["disort_cpu_total_s"] for r in base)
      k256 = records(case, "nstr100_k256")
      n100 = records(case, "nstr100")
      overhead = (statistics.median(r["disort_cpu_total_s"] for r in k256) - statistics.median(r["disort_cpu_total_s"] for r in n100)
                  if k256 and n100 else None)
      for experiment in EXPERIMENTS:
         recs = records(case, experiment)
         if not recs:
            continue
         times = [r["disort_cpu_total_s"] for r in recs]
         median = statistics.median(times)
         build = recs[0].get("kernel_build", recs[0]["kernel"])
         corrected = median - overhead if overhead is not None and build in ("single", "double") and not experiment.endswith("_k256") \
            else median
         rows.append({"case": case, "experiment": experiment, "nstreams": recs[0]["nstreams"], "kernel": recs[0]["kernel"],
                      "kernel_build": build, "repeats": len(times), "median_disort_cpu_s": corrected,
                      "raw_median_s": median, "min_s": min(times), "max_s": max(times),
                      "spread_percent": 100.0 * (max(times) - min(times)) / median,
                      "mxcmu256_overhead_subtracted_s": overhead if corrected != median else 0.0,
                      "ratio_to_baseline": corrected / t_base,
                      "direct_calls": recs[0]["disort_calls"]["direct"], "diffuse_calls": recs[0]["disort_calls"]["diffuse"],
                      "emission_calls": recs[0]["disort_calls"]["emission"],
                      "direct_share": statistics.median(r["disort_cpu_s"]["direct"] / r["disort_cpu_total_s"] for r in recs)})
   RESULTS.mkdir(parents=True, exist_ok=True)
   with open(RESULTS / "runtime.csv", "w", newline="") as handle:
      writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
      writer.writeheader()
      writer.writerows(rows)
   for row in rows:
      print(f"{row['case']:14s} {row['experiment']:22s} {row['median_disort_cpu_s']:9.1f} s  x{row['ratio_to_baseline']:6.2f}  "
            f"(n={row['repeats']}, spread {row['spread_percent']:.1f}%, build {row['kernel_build']})")


if __name__ == "__main__":
   main()
