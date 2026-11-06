#!/usr/bin/env python3
"""Run the maintained causal TOV audit without requiring an editable install."""
from pathlib import Path
import sys

PACKAGE_ROOT = Path(__file__).resolve().parents[1] / "compose_eos_pkg"
sys.path.insert(0, str(PACKAGE_ROOT))

from compose_eos.tov_cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
