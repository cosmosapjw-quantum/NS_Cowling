#!/usr/bin/env python3
"""Build the twelve source-admitted E390-03B crust-only CompOSE components."""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos import build_crust_only  # noqa: E402


def main() -> int:
    output_root = Path(sys.argv[1]).expanduser().resolve() if len(sys.argv) > 1 else (
        ROOT / "build" / "e390_03b" / "crust_only"
    )
    receipts = build_crust_only(output_root)
    manifest_path = output_root / "crust_only_manifest.v1.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = {
        "SPG_M1": 365, "SPG_M2": 377, "SPG_M3": 366, "SPG_M4": 375, "SPG_M5": 363,
        "GPPVA_DD2": 283, "GPPVA_DDME2": 288, "GPPVA_FSU2": 269, "GPPVA_FSU2H": 300,
        "GPPVA_FSU2R": 296, "GPPVA_NL3wrL55": 298, "GPPVA_TM1e": 303,
    }
    observed = {item["model_id"]: item["output"]["n"] for item in receipts}
    if observed != expected:
        raise RuntimeError(f"E390-03B cohort mismatch: {observed!r}")
    if manifest.get("status") != "ALL_ADMITTED_CRUST_COMPONENTS_CONSTRUCTED_NATIVE_READBACK_VERIFIED":
        raise RuntimeError(f"E390-03B manifest is not certified: {manifest.get('status')!r}")
    if manifest.get("total_output_rows") != 3883:
        raise RuntimeError(f"E390-03B output-row total mismatch: {manifest.get('total_output_rows')!r}")
    print(json.dumps({
        "status": "E390-03B_CRUST_ONLY_ARTIFACTS_BUILT",
        "output_root": str(output_root),
        "models": observed,
        "total_output_rows": manifest["total_output_rows"],
        "manifest": str(manifest_path),
        "coverage_report": str(output_root / "crust_coverage_report.v1.json"),
        "physical_validation": {
            item["model_id"]: item["checks"]["physical_validation"]
            for item in receipts
        },
        "blocked": manifest["blocked_models"],
        "note": "Source-faithful construction/readback is separate from physical EOS admissibility.",
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
