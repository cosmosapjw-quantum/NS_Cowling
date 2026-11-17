#!/usr/bin/env python3
"""Generate composition-aware homogeneous DD-RMF CompOSE tables."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos.ddrmf_generator import SUPPORTED_MODELS, generate_ddrmf_core


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--models", nargs="+", choices=SUPPORTED_MODELS, default=list(SUPPORTED_MODELS))
    parser.add_argument("--n-min", type=float, required=True, help="Requested lower density in fm^-3")
    parser.add_argument("--n-max", type=float, required=True, help="Requested upper density in fm^-3")
    parser.add_argument("--points", type=int, default=129)
    args = parser.parse_args()
    if (not np.isfinite([args.n_min, args.n_max]).all()
            or not 0 < args.n_min < args.n_max or args.points < 3):
        parser.error("require 0 < n-min < n-max and at least 3 points")
    if len(set(args.models)) != len(args.models):
        parser.error("models must be unique")
    args.output.mkdir(parents=True, exist_ok=False)
    density = np.geomspace(args.n_min, args.n_max, args.points)
    for model in args.models:
        receipt = generate_ddrmf_core(model, density, args.output / model)
        print(json.dumps({"model": model, "status": receipt["status"],
                          "rows": receipt["row_count"], "composition": receipt["composition"],
                          "physical_diagnostics": receipt["physical_diagnostics"]["overall"]}),
              flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
