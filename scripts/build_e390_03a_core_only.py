#!/usr/bin/env python3
"""Build the five source-admitted E390-03A core-only CompOSE components."""
from __future__ import annotations
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos import build_core_only  # noqa: E402


def main() -> int:
    output_root = Path(sys.argv[1]).expanduser().resolve() if len(sys.argv) > 1 else (
        ROOT / "build" / "e390_03a" / "core_only"
    )
    receipts = build_core_only(output_root)
    manifest_path = output_root / "core_only_manifest.v1.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = {"APR_APR": 103, "BBB_BHF2": 15, "GM_GM1": 301, "RG_KDE0v1": 274, "RG_SLY4": 54}
    observed = {item["model_id"]: item["output"]["n"] for item in receipts}
    if observed != expected:
        raise RuntimeError(f"E390-03A output cohort mismatch: {observed!r}")
    if manifest.get("status") != "ALL_REQUESTED_COMPONENTS_CONSTRUCTED_AND_NATIVE_READBACK_VERIFIED":
        raise RuntimeError(f"E390-03A manifest is not certified: {manifest.get('status')!r}")
    print(json.dumps({
        "status": "E390-03A_CORE_ONLY_ARTIFACTS_BUILT",
        "output_root": str(output_root),
        "models": observed,
        "manifest": str(manifest_path),
        "physical_validation": {
            item["model_id"]: item["checks"]["physical_validation"]
            for item in receipts
        },
        "note": "Construction/readback success is distinct from physical EOS validation; inspect each receipt.",
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
