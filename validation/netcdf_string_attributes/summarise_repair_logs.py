"""Summarise the per-file JSON records of the text-attribute repair (validation only).

    python validation/netcdf_string_attributes/summarise_repair_logs.py [LOG_DIR]

Reads every ``*.repair.json`` written by ``python -m
oraclut.repair_string_attributes --log-dir`` and prints a Markdown summary:
counts per version and status, the verification totals (variables, attributes,
data bytes compared), the conversions per file, and any file that is not
``repaired`` or ``compliant``.  Writes ``results_repair_summary.tsv`` next to
this script.
"""

from __future__ import annotations

import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
LOG_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "logs"


def main() -> int:
    records = []
    for path in sorted(LOG_DIR.glob("*.repair.json")):
        with open(path) as handle:
            records.append(json.load(handle))
    if not records:
        print(f"no *.repair.json records in {LOG_DIR}")
        return 1

    by_version_status: dict[str, Counter] = defaultdict(Counter)
    totals = Counter()
    rows = []
    for r in records:
        version = r.get("version") or "?"
        by_version_status[version][r["status"]] += 1
        v = r.get("verification") or {}
        totals["files"] += 1
        totals["conversions"] += len(r["conversions"])
        totals["required"] += r["required_conversions"]
        totals["bytes"] += v.get("bytes_compared", 0)
        totals["variables"] += v.get("variables_compared", 0)
        totals["attributes"] += v.get("attributes_compared", 0)
        totals["seconds"] += r.get("elapsed_seconds", 0.0)
        rows.append((Path(r["path"]).name, version, r["size_bytes"], r["status"], len(r["conversions"]),
                     r["required_conversions"], v.get("variables_compared", ""), v.get("attributes_compared", ""),
                     v.get("bytes_compared", ""), v.get("slabs_compared", ""), round(r.get("elapsed_seconds", 0.0), 1),
                     ",".join(f"{k}={val}" for k, val in sorted((r.get("orac_spacing_read") or {}).items())),
                     r.get("original_retained_as") or "", r["message"]))

    print("| version | " + " | ".join(sorted({s for c in by_version_status.values() for s in c})) + " | files |")
    statuses = sorted({s for c in by_version_status.values() for s in c})
    print("|---|" + "---|" * (len(statuses) + 1))
    for version in sorted(by_version_status):
        counter = by_version_status[version]
        print(f"| v{version} | " + " | ".join(str(counter.get(s, 0)) for s in statuses) + f" | {sum(counter.values())} |")
    print()
    print(f"files {totals['files']}; attributes converted {totals['conversions']} "
          f"(axis spacing {totals['required']}); verification compared {totals['variables']} variables, "
          f"{totals['attributes']} attributes and {totals['bytes']:,} data bytes; "
          f"total elapsed {totals['seconds'] / 60:.1f} min")
    problems = [r for r in records if r["status"] not in ("repaired", "compliant")]
    print(f"files not repaired/compliant: {len(problems)}")
    for r in problems:
        print(f"  {Path(r['path']).name}: {r['status']}: {r['message']}")

    out = HERE / "results_repair_summary.tsv"
    with open(out, "w") as handle:
        handle.write("\t".join(("file", "version", "size_bytes", "status", "conversions", "axis_spacing_conversions",
                                "variables_compared", "attributes_compared", "data_bytes_compared", "slabs",
                                "elapsed_s", "orac_spacing_read", "original_retained_as", "message")) + "\n")
        for row in rows:
            handle.write("\t".join(str(x) for x in row) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
