"""Command-line entry point for the LUT text-attribute audit and repair.

    python -m oraclut.repair_string_attributes [--dry-run | --check] FILE...

See ``oraclut.io.string_attributes`` for the ORAC NetCDF metadata contract,
the repair strategy and the options.  (A separate entry module avoids the
double-import warning Python emits when a module that its package already
imports is run with ``-m``.)
"""

from __future__ import annotations

import sys

from .io.string_attributes import main

if __name__ == "__main__":
    sys.exit(main())
