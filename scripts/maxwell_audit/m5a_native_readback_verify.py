"""§6 follow-up: actually run the native Fortran read/write/verify pathway on
the M5A headline barotrope so the receipt's native_fortran_compose_readback
and native_c1_c11_verification fields are TRUE with evidence, not BLOCKED.

The §6 spec says these two fields 'must be false or BLOCKED unless actual
native calls are executed successfully'. This script executes them and
records the evidence (row count match, file hashes, C1-C11 failed_checks).

The barotrope is (n_B, mu_B, P, eps). We build an EOSData from it and run:
  1. write_compose  -> CompOSE directory (native Fortran writer)
  2. read_compose   -> EOSData (native Fortran reader) -- read-back
  3. check_consistency -> C1-C11 (native Fortran verifier)

This does NOT change the analytic CSS branch formula or the barotrope values.
It only adds native-verification evidence to the receipt.
"""
from __future__ import annotations
import json
import os
import sys
import hashlib
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_ROOT", str(ROOT))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_EXE",
                      str(ROOT / "build" / "compose_eos_fortran" / "compose_eos_fortran"))

from compose_eos import read_compose, write_compose, check_consistency  # noqa: E402
from compose_eos.eos_data import EOSData  # noqa: E402
from compose_eos.native import NativeExecutionError  # noqa: E402

M_N = 939.56536
BAROTROPE = ROOT / "download" / "m5a_real_css" / "BBB_BHF2_m5a_real_css.json"
RECEIPT = ROOT / "download" / "m5a_real_css" / "E390-05A_M5A_REAL_CSS_RECEIPT_20261114.json"
NATIVE_OUT = ROOT / "download" / "m5a_real_css" / "native_readback"


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def build_eosdata_from_barotrope() -> EOSData:
    b = json.loads(BAROTROPE.read_text())
    n = np.asarray(b["n_B"], float)
    mu = np.asarray(b["mu_B"], float)
    P = np.asarray(b["P"], float)
    eps = np.asarray(b["eps"], float)
    # The barotrope has a duplicate n_B at the hadronic->CSS seam (the
    # hadronic last row and the CSS first row both carry n=n_low). The native
    # CompOSE reader requires strictly increasing n_B with >= 3 points, so we
    # drop the duplicate seam row for the native read-back ONLY. The original
    # barotrope file is NOT modified; this deduplication is a serialization
    # detail recorded in the evidence.
    keep = np.concatenate([[True], np.diff(n) > 0])
    n = n[keep]; mu = mu[keep]; P = P[keep]; eps = eps[keep]
    F = eps / n
    # real BBB_BHF2 mp (§6 finding #3)
    d_in = read_compose(str(ROOT / "eos_generator" / "eos_data" / "BBB_BHF2"))
    real_mp = float(d_in.mp)
    return EOSData(
        nb=n, mn=M_N, mp=real_mp,
        Q1=P / n, Q2=np.zeros_like(n), Q3=mu / M_N - 1.0,
        Q4=np.zeros_like(n), Q5=np.zeros_like(n),
        Q6=F / M_N - 1.0, Q7=F / M_N - 1.0,
        q_extra=np.zeros_like(n), Iphase=np.zeros_like(n, dtype=int),  # single continuous phase = Iphase 0 (CompOSE default)
        Yn=np.zeros_like(n), Yp=np.zeros_like(n), Ye=np.zeros_like(n),
        Ymu=np.zeros_like(n), Aav=np.zeros_like(n), Zav=np.zeros_like(n),
        Yh=np.zeros_like(n),
        Il=1, iT=0, iYq=1, inb_start=1, temperature_mev=0.0,
        charge_fraction_grid=0.0, q_extra_count=0,
        composition_available=False, beta_equilibrated=True)


def run_native_readback_verify() -> dict:
    """Execute the native Fortran write -> read -> verify pathway on the M5A
    headline barotrope. Returns the evidence dict."""
    if NATIVE_OUT.exists():
        import shutil
        shutil.rmtree(NATIVE_OUT)
    NATIVE_OUT.mkdir(parents=True)

    data = build_eosdata_from_barotrope()
    n_input = data.n
    # record the deduplication detail: the original barotrope had a duplicate
    # n_B at the hadronic->CSS seam; the native read-back uses the deduplicated
    # strictly-increasing version. The original file is NOT modified.
    b_raw = json.loads(BAROTROPE.read_text())
    n_raw = np.asarray(b_raw["n_B"], float)
    dedup_dropped = int(len(n_raw) - n_input)
    # 1. native write -> CompOSE directory
    compose_dir = NATIVE_OUT / "compose"
    write_compose(str(compose_dir), data)
    written_files = sorted(p.name for p in compose_dir.iterdir())
    written_hashes = {p.name: sha256_file(p) for p in compose_dir.iterdir()}

    # 2. native read-back -> EOSData
    readback = read_compose(str(compose_dir))
    n_readback = readback.n
    readback_n_match = (n_readback == n_input)

    # A matching row count is necessary but not sufficient. Compare the
    # physical barotrope arrays after the native CompOSE round-trip.
    comparison_tolerance = 5.0e-13
    reference_arrays = {
        "n_B": np.asarray(data.nb, float),
        "mu_B": np.asarray(data.mub, float),
        "P": np.asarray(data.p, float),
        "eps": np.asarray(data.eps, float),
    }
    observed_arrays = {
        "n_B": np.asarray(readback.nb, float),
        "mu_B": np.asarray(readback.mub, float),
        "P": np.asarray(readback.p, float),
        "eps": np.asarray(readback.eps, float),
    }
    readback_relative_residuals = {}
    for name, expected in reference_arrays.items():
        observed = observed_arrays[name]
        if observed.shape != expected.shape:
            readback_relative_residuals[name] = float("inf")
            continue
        scale = np.maximum(np.abs(expected), np.finfo(float).tiny)
        readback_relative_residuals[name] = float(
            np.max(np.abs(observed - expected) / scale)
        )
    readback_value_match = bool(
        readback_n_match
        and all(
            np.isfinite(value) and value <= comparison_tolerance
            for value in readback_relative_residuals.values()
        )
    )

    # 3. native C1-C11 verify
    consistency = check_consistency(readback, verbose=False)
    c1_c11_pass = bool(consistency["overall"]["pass_"])
    failed_checks = consistency["overall"]["failed_checks"]

    return {
        "native_write_executed": True,
        "native_read_executed": True,
        "native_verify_executed": True,
        "input_n_points": int(n_input),
        "readback_n_points": int(n_readback),
        "readback_n_match": bool(readback_n_match),
        "readback_value_match": bool(readback_value_match),
        "readback_comparison_relative_tolerance": comparison_tolerance,
        "readback_relative_residuals": readback_relative_residuals,
        "original_barotrope_n_points": int(len(n_raw)),
        "dedup_dropped_seam_rows": dedup_dropped,
        "dedup_note": "the barotrope has a duplicate n_B at the hadronic->CSS seam; the native CompOSE reader requires strictly increasing n_B, so the duplicate seam row is dropped for the native read-back ONLY. The original BBB_BHF2_m5a_real_css.json is NOT modified.",
        "written_compose_files": written_files,
        "written_file_hashes": written_hashes,
        "c1_c11_overall_pass": c1_c11_pass,
        "c1_c11_failed_checks": failed_checks,
        "c2_first_law_max_rel": float(consistency["C2"]["max_rel"]),
        "c6_min_dPdnb": float(consistency["C6"]["min_dPdnb"]),
        "c10_cs2_range": [float(consistency["C10"]["min_cs2"]), float(consistency["C10"]["max_cs2"])],
        "fortran_executable_sha256": sha256_file(Path(os.environ["COMPOSE_EOS_FORTRAN_EXE"])),
    }


def update_receipt(evidence: dict) -> None:
    """Patch the M5A receipt: set native_fortran_compose_readback and
    native_c1_c11_verification to TRUE with evidence, only if the native
    pathway actually succeeded. Keep output_is_compose_generated_eos as-is
    (the headline barotrope is still a barotrope; the native pathway here
    proves read-back + C1-C11, not that the barotrope is a CompOSE-generated
    multi-phase EOS)."""
    r = json.loads(RECEIPT.read_text())
    readback_ok = (
        evidence["native_write_executed"]
        and evidence["native_read_executed"]
        and evidence["readback_n_match"]
        and evidence["readback_value_match"]
    )
    verify_ok = (
        readback_ok
        and evidence["native_verify_executed"]
        and evidence["c1_c11_overall_pass"]
    )
    vbt = r["verification_backend_truthful"]
    vbt["native_fortran_compose_readback"] = bool(readback_ok)
    vbt["native_c1_c11_verification"] = bool(verify_ok)
    vbt["native_readback_scope"] = "HEADLINE_BAROTROPE_ONLY"
    # The file contains a single barotrope and the receipt contains 27 distinct
    # parameter combinations. No verified case key is derivable from this
    # artifact alone, so per-case verification remains false until each output
    # is independently generated, read back, and verified.
    vbt["native_verified_case_keys"] = []
    r["native_readback_evidence"] = evidence
    for case in r["results"].values():
        if "native_fortran_compose_readback" in case:
            case["native_fortran_compose_readback"] = False
            case["native_c1_c11_verification"] = False
            case["native_verification_status"] = "NOT_VERIFIED_INDIVIDUALLY"
    RECEIPT.write_text(json.dumps(r, indent=1))


if __name__ == "__main__":
    print("=== M5A native Fortran read/write/verify pathway ===")
    try:
        ev = run_native_readback_verify()
        print(f"  write: {ev['native_write_executed']} files={ev['written_compose_files']}")
        print(f"  read:  {ev['native_read_executed']} n_match={ev['readback_n_match']} ({ev['input_n_points']} -> {ev['readback_n_points']})")
        print(f"  verify: {ev['native_verify_executed']} c1_c11_pass={ev['c1_c11_overall_pass']} failed='{ev['c1_c11_failed_checks']}'")
        print(f"  C2 first-law max_rel={ev['c2_first_law_max_rel']:.3e}  C6 min_dPdnb={ev['c6_min_dPdnb']:.3e}  C10 cs2=[{ev['c10_cs2_range'][0]:.4e},{ev['c10_cs2_range'][1]:.4e}]")
        update_receipt(ev)
        print(f"  receipt updated: native_fortran_compose_readback={ev['readback_n_match']} native_c1_c11_verification={ev['c1_c11_overall_pass']}")
    except NativeExecutionError as e:
        print(f"  NATIVE FAILED (leaving receipt fields BLOCKED=False): {e}")
        ev = {"native_failed": str(e), "native_fortran_compose_readback": False,
              "native_c1_c11_verification": False}
    except Exception as e:
        print(f"  ERROR: {type(e).__name__}: {e}")
        ev = {"error": str(e)}
