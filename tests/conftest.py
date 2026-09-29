"""Make the repository root and src/ importable for the tests.

The IDL-structured entry point create_orac_luts.py lives at the repository
root (next to the IDL programs it mirrors); the package lives in src/.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
for path in (ROOT, ROOT / "src"):
   if str(path) not in sys.path:
      sys.path.insert(0, str(path))
