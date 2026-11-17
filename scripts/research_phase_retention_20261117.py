"""Validate native phase transport on regenerated APR; no source-support admission.

Python assembles receipts and independent algebraic/CSS controls only. Production
matching and CompOSE I/O remain Fortran. The public APR support guard must reject.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python_frontend"))

import numpy as np

from compose_eos import EOSData, check_consistency, match_crust_core, read_compose, write_compose
from compose_eos.native import executable
from compose_cowling.css_adapter import (
    UnsupportedSourcePhaseBoundary, css_match_maxwell_eos, phase_table_from_css_match,
)

def hashes(directory):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(directory.iterdir()) if p.is_file()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apr", type=Path, required=True, help="fresh native generated CompOSE directory")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)  # preserve every prior attempt
    crust_dir = ROOT / "data/research_matching_20261117/spg/native_129_repaired/SPG_M1/generated"
    core_dir = args.apr.resolve()
    source_hashes = {"crust": hashes(crust_dir), "core": hashes(core_dir)}
    crust, core = read_compose(crust_dir), read_compose(core_dir)
    record = {
        "base_head": "4a44333377fc6252ad6bcc0c1024614277199f31",
        "scope": "native phase transport and independent invariants; full APR matching admission remains HOLD",
        "source_fidelity": "unchanged explicitly scoped APR functional/Gibbs representation",
        "CompOSE_comparison": "NOT_EVALUATED; native CompOSE serialization readback only",
        "numerical_support": "finite generated APR interval ending at 0.5 fm^-3; source boundary interpolation unapproved",
        "binary_sha256": hashlib.sha256(executable().read_bytes()).hexdigest(),
        "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in [*sorted((ROOT / "fortran_eos/src").glob("*.f90")),
                                    ROOT / "fortran_eos/app/compose_eos_fortran.f90",
                                    ROOT / "python_frontend/compose_eos/api.py",
                                    ROOT / "python_frontend/compose_cowling/css_adapter.py",
                                    ROOT / "scripts/match_report_admission.py", Path(__file__).resolve()]},
        "numpy_version": np.__version__,
        "source_hashes": source_hashes,
        "core_rows": core.n, "core_phase_ids": core.Iphase.tolist(),
        "methods": [],
    }
    # Public guard is deliberately unchanged, including immutable support seal.
    try:
        match_crust_core(crust, core, method="css", drop_composition=True)
    except ValueError as exc:
        assert "support forbids cross-phase" in str(exc)
        record["public_support_guard"] = str(exc)
    else:
        raise AssertionError("public full-APR support guard was bypassed")

    for method in ("css", "cubic_F", "cubic_P", "quintic_F", "quintic_P"):
        folder = out / method
        folder.mkdir()
        item = {"method": method, "scientific_admission": "NOT_ADMITTED_SOURCE_SUPPORT_HOLD"}
        record["methods"].append(item)
        try:
            dirs = [folder / name for name in ("crust", "core", "internal")]
            for directory in dirs:
                directory.mkdir()
            crust.write_internal_directory(dirs[0])
            core.write_internal_directory(dirs[1])
            command = [str(executable()), "match-core", *map(str, dirs), method, "40", "true"]
            run = subprocess.run(command, text=True, capture_output=True)
            (folder / "native.stdout.log").write_text(run.stdout)
            (folder / "native.stderr.log").write_text(run.stderr)
            item.update(command=command, native_exit=run.returncode)
            assert run.returncode == 0, run.stdout + run.stderr
            report = json.loads(run.stdout.splitlines()[-1])
            matched = EOSData.from_internal_directory(dirs[2])
            i1, i2 = report["i1"], report["i2"]
            start = matched.n - (core.n - i2)
            np.testing.assert_array_equal(matched.Iphase[:i1+1], crust.Iphase[:i1+1])
            np.testing.assert_array_equal(matched.Iphase[start:], core.Iphase[i2:])
            np.testing.assert_array_equal(matched.Iphase, report["output_phase_ids"])
            errors = {}
            for name, original, sl in (("crust", crust, slice(0, i1+1)), ("core", core, slice(i2, None))):
                target = slice(0, i1+1) if name == "crust" else slice(start, None)
                np.testing.assert_array_equal(matched.nb[target], original.nb[sl])
                np.testing.assert_array_equal(matched.Q1[target], original.Q1[sl])
                a, b = matched.F[target], original.F[sl]
                error = float(np.max(np.abs(a-b)))
                assert error <= 64*np.finfo(float).eps*max(1., float(np.max(np.abs(b))))
                errors[name + "_F_MeV"] = error
            euler = float(np.max(np.abs(matched.eps + matched.p - matched.nb * matched.mub)
                                 / np.abs(matched.eps + matched.p)))
            assert euler <= 2e-12
            assert not matched.composition_available
            assert not any(matched.Yn_present) and not any(matched.Yp_present)
            compose = folder / "compose"
            write_compose(compose, matched)
            back = read_compose(compose)
            for name in ("nb", "Q1", "Q3", "Q4", "Q6", "Q7", "Iphase"):
                np.testing.assert_array_equal(getattr(back, name), getattr(matched, name))
            consistency = check_consistency(back)
            assert consistency["overall"]["pass_"], consistency
            item.update(status="TRANSPORT_AND_SCOPED_CONTROLS_PASS", rows=back.n,
                        retained_energy_errors=errors, euler_relative_residual=euler,
                        output_hashes=hashes(compose), consistency=consistency,
                        native_report=report)
            if method == "css":
                hi, join = report["css_transition_high_index"], report["css_core_bridge_index"]
                s = report["css_single_sound_speed"]
                p, e, n = back.p[hi:join+1], back.eps[hi:join+1], back.nb[hi:join+1]
                predicted_e = e[0] + (p-p[0])/s
                predicted_n = n[0]*((predicted_e+p)/(e[0]+p[0]))**(1/(1+s))
                css_error = max(float(np.max(np.abs(e/predicted_e-1))),
                                float(np.max(np.abs(n/predicted_n-1))))
                assert css_error < 5e-9  # existing analytic CSS bridge contract
                item["independent_css_law_max_relative_error"] = css_error
                intake = phase_table_from_css_match(compose, report)
                np.testing.assert_array_equal(intake.Iphase, back.Iphase)
                expected = [(i-1, i) for i in range(join+1, back.n)
                            if back.Iphase[i] != back.Iphase[i-1]]
                assert len(expected) == 2
                try:
                    css_match_maxwell_eos(compose, report)
                except UnsupportedSourcePhaseBoundary as exc:
                    assert exc.boundaries == expected
                    item["stellar_intake"] = {"category": exc.category, "boundaries": exc.boundaries,
                                               "phase_labels_preserved": True}
                else:
                    raise AssertionError("undeclared APR boundary interpolation was admitted")
        except Exception as exc:
            item.update(status="FAILED", error_type=type(exc).__name__, error=str(exc))
    record["source_bytes_unchanged"] = source_hashes == {"crust": hashes(crust_dir), "core": hashes(core_dir)}
    record["status"] = ("PASS_SCOPED_TRANSPORT" if record["source_bytes_unchanged"] and
                        all(m["status"] == "TRANSPORT_AND_SCOPED_CONTROLS_PASS" for m in record["methods"]) else "FAILED")
    (out / "RECEIPT.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": record["status"], "methods": [{k: m[k] for k in
                       ("method", "status")} for m in record["methods"]]}))
    return int(record["status"] == "FAILED")


if __name__ == "__main__":
    raise SystemExit(main())
