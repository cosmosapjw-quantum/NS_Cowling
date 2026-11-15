"""Core-crust matching pipeline: maxwell + css for all (crust, core) pairs.

Runs the native Fortran maxwell-transition and css (match-core) for every
(crust, core) combination (8 crusts x 5 cores = 40 pairs x 2 methods = 80
attempts). A pair "survives" if the native backend produces a complete
combined EOS with successful CompOSE read-back.

Surviving matched EOS are written as CompOSE tables under download/matched/
and a JSON manifest of survivors is produced.

Physics:
  - maxwell: sharp first-order phase transition between the crust (low-density
    phase) and core (high-density phase) at the common P(mu_B) coexistence
    point. Requires a pressure crossing between the two phases.
  - css: smooth constant-sound-speed matching via the native match-core kernel
    with METHOD='css'. Produces a continuous thermodynamic bridge.

Both methods use drop_composition=true (the matched table is thermodynamic-only;
composition bridging is a separate subscope).
"""
from __future__ import annotations
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos.native import executable, run_native  # noqa: E402

CRUST_MODELS = [
    "SPG_M1", "SPG_M2", "SPG_M3", "SPG_M4", "SPG_M5",
    "GPPVA_DD2", "GPPVA_DDME2", "GPPVA_TW",
]
CORE_MODELS = ["APR_APR", "BBB_BHF2", "GM_GM1", "RG_SLY4", "RG_KDE0v1"]
CRUST_ROOT = ROOT / "download" / "reconstruction"
CORE_ROOT = ROOT / "eos_generator" / "eos_data"
OUT_ROOT = ROOT / "download" / "matched"
WORK_ROOT = ROOT / "download" / "matched_work"
N_BRIDGE = 8


def sha256_file(path: Path | str) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def run_pipeline(method: str, crust_dir: Path, core_dir: Path, out_dir: Path) -> dict:
    """Run one (crust, core, method) match. Returns a result dict.

    method in {"maxwell", "css"}.
    """
    work = out_dir / "work"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)
    (out_dir / "compose").mkdir(parents=True, exist_ok=True)

    crust_internal = work / "crust_internal"
    core_internal = work / "core_internal"
    matched_internal = work / "matched_internal"
    crust_internal.mkdir(exist_ok=True)
    core_internal.mkdir(exist_ok=True)
    matched_internal.mkdir(exist_ok=True)

    exe = executable()
    exe_sha = sha256_file(exe)

    try:
        # 1. read crust -> internal
        r1 = run_native("read", str(crust_dir), str(crust_internal))
        if r1.get("status") != "ok":
            return {"survived": False, "method": method, "stage": "read_crust",
                    "error": r1.get("message", "?")}
        # 2. read core -> internal
        r2 = run_native("read", str(core_dir), str(core_internal))
        if r2.get("status") != "ok":
            return {"survived": False, "method": method, "stage": "read_core",
                    "error": r2.get("message", "?")}

        if method == "maxwell":
            # maxwell-transition: low=crust, high=core, local_neutrality=true
            r3 = run_native("maxwell-transition", str(crust_internal), str(core_internal),
                            str(matched_internal), "true")
        else:  # css
            r3 = run_native("match-core", str(crust_internal), str(core_internal),
                            str(matched_internal), "css", str(N_BRIDGE), "true")
        if r3.get("status") != "ok":
            return {"survived": False, "method": method, "stage": method,
                    "error": r3.get("message", "?"),
                    "report": r3}

        combined_n = int(r3.get("n", 0))
        if combined_n < 4:
            return {"survived": False, "method": method, "stage": "match",
                    "error": f"combined n={combined_n} too small", "report": r3}

        # 4. write internal -> CompOSE output
        r4 = run_native("write", str(matched_internal), str(out_dir / "compose"))
        if r4.get("status") != "ok":
            return {"survived": False, "method": method, "stage": "write",
                    "error": r4.get("message", "?")}

        # 5. readback verification
        readback_internal = work / "readback_internal"
        readback_internal.mkdir(exist_ok=True)
        r5 = run_native("read", str(out_dir / "compose"), str(readback_internal))
        readback_n = int(r5.get("n", 0)) if r5.get("status") == "ok" else 0

        output_files = sorted(p.name for p in (out_dir / "compose").iterdir() if p.is_file())

        return {
            "survived": True,
            "method": method,
            "stage": "complete",
            "combined_n": combined_n,
            "readback_n": readback_n,
            "readback_status": r5.get("status", "?"),
            "output_files": output_files,
            "fortran_executable_sha256": exe_sha,
            "match_report": r3,
            "write_report": r4,
        }
    except Exception as exc:
        return {"survived": False, "method": method, "stage": "exception",
                "error": f"{type(exc).__name__}: {exc}"}


def main() -> int:
    if OUT_ROOT.exists():
        shutil.rmtree(OUT_ROOT)
    OUT_ROOT.mkdir(parents=True)
    if WORK_ROOT.exists():
        shutil.rmtree(WORK_ROOT)
    WORK_ROOT.mkdir(parents=True)

    results = []
    survivors = []
    for crust in CRUST_MODELS:
        crust_dir = CRUST_ROOT / crust
        if not crust_dir.is_dir():
            continue
        for core in CORE_MODELS:
            core_dir = CORE_ROOT / core
            if not core_dir.is_dir():
                continue
            for method in ("maxwell", "css"):
                pair_id = f"{crust}__{core}__{method}"
                out_dir = OUT_ROOT / pair_id
                out_dir.mkdir(parents=True, exist_ok=True)
                t0 = time.time()
                res = run_pipeline(method, crust_dir, core_dir, out_dir)
                res["crust_model"] = crust
                res["core_model"] = core
                res["pair_id"] = pair_id
                res["elapsed_s"] = round(time.time() - t0, 2)
                results.append(res)
                status = "SURVIVED" if res["survived"] else f"FAILED@{res.get('stage')}"
                print(f"  [{status:20s}] {crust:14s} + {core:12s} / {method:7s}  "
                      f"(n={res.get('combined_n','-')})  {res.get('elapsed_s','-')}s  "
                      f"{res.get('error','')[:60]}")
                if res["survived"]:
                    survivors.append({
                        "pair_id": pair_id,
                        "crust_model": crust,
                        "core_model": core,
                        "method": method,
                        "combined_n": res["combined_n"],
                        "readback_n": res["readback_n"],
                        "readback_status": res["readback_status"],
                        "output_dir": str(out_dir / "compose"),
                        "match_report": res["match_report"],
                    })
                else:
                    # clean up failed attempt dir
                    shutil.rmtree(out_dir, ignore_errors=True)

    # write manifest
    n_maxwell = sum(1 for s in survivors if s["method"] == "maxwell")
    n_css = sum(1 for s in survivors if s["method"] == "css")
    manifest = {
        "schema": "ns-cowling.core-crust-matching.v1",
        "built_date": "2026-11-14",
        "crust_models": CRUST_MODELS,
        "core_models": CORE_MODELS,
        "methods": ["maxwell", "css"],
        "total_attempts": len(results),
        "total_survivors": len(survivors),
        "survivors_maxwell": n_maxwell,
        "survivors_css": n_css,
        "survivors": survivors,
        "all_results": results,
    }
    manifest_path = OUT_ROOT / "matching_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
                             encoding="utf-8")
    print(f"\nMatching complete: {len(survivors)}/{len(results)} survived "
          f"(maxwell={n_maxwell}, css={n_css})")
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
