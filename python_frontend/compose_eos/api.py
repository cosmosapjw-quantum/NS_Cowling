"""Public Python API backed by the native Fortran EOS reader/writer/verifier."""
from __future__ import annotations
import json
import tempfile
from pathlib import Path
from typing import Any
from .eos_data import EOSData
from .native import run_native, NativeFeatureUnavailable

def read_compose(directory: str | Path) -> EOSData:
    """Read a cold 1-D CompOSE slice with the Fortran reader."""
    source = Path(directory).expanduser().resolve()
    with tempfile.TemporaryDirectory(prefix="compose-eos-read-") as td:
        target = Path(td) / "internal"
        target.mkdir()
        run_native("read", source, target)
        return EOSData.from_internal_directory(target, source_dir=str(source))

def write_compose(directory: str | Path, data: EOSData, copy_from: str | Path | None = None, *,
                  composition_tolerance: float = 2e-3) -> None:
    """Write CompOSE files using the Fortran backend."""
    del copy_from, composition_tolerance
    target = Path(directory).expanduser().resolve()
    target.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="compose-eos-write-") as td:
        internal = Path(td) / "internal"
        data.write_internal_directory(internal)
        run_native("write", internal, target)

def check_consistency(data: EOSData, verbose: bool = True, *,
                      algebraic_tolerance: float = 2e-6,
                      max_log_step_for_derivatives: float = 0.25,
                      cs2_tolerance: float = 1e-6,
                      composition_tolerance: float = 2e-5) -> dict[str, Any]:
    """Run native C1-C11 checks; C7-C9 are finite-difference diagnostics."""
    del verbose
    if (algebraic_tolerance != 2e-6 or max_log_step_for_derivatives != 0.25
        or cs2_tolerance != 1e-6 or composition_tolerance != 2e-5):
        raise NativeFeatureUnavailable("custom C1-C11 tolerances are not yet exposed by the Fortran CLI")
    with tempfile.TemporaryDirectory(prefix="compose-eos-verify-") as td:
        directory = Path(td) / "internal"
        data.write_internal_directory(directory)
        output = Path(td) / "report.json"
        run_native("verify-csv", directory, output)
        raw = json.loads(output.read_text(encoding="utf-8"))
    c1 = {"temperature_mev": data.temperature_mev, "maxQ2": raw["C1_max_Q2"],
          "maxQ6mQ7": raw["C1_max_Q6_minus_Q7"], "maxmul": raw["C1_max_mu_l"],
          "pass_": raw["C1"]["pass"]}
    c2 = {"applicable": raw["C2"]["applicable"], "max_rel": raw["C2_max_relative_residual"],
          "pass_": raw["C2"]["pass"], "basis": "Fortran CompOSE thermodynamic closure"}
    def composition_check(name: str):
        key = "C3_max_absolute_residual" if name == "C3" else "C4_max_absolute_residual"
        item = raw[name]
        return {"status": item["status"], "max_abs": raw[key], "pass_": item.get("pass")}
    c5 = {"maxmul": raw["C1_max_mu_l"], "applicable": raw["C5"]["applicable"],
          "pass_": raw["C5"].get("pass")}
    c6 = {"min_dPdnb": raw["min_dP_dnb"], "min_delta_p": raw["min_delta_P"],
          "pass_": raw["C6"]["pass"]}
    def diagnostic(name: str):
        return {"max_rel": raw[name]["max_rel"], "n_points": raw[name]["n_points"],
                "status": "FINITE_DIFFERENCE_DIAGNOSTIC"}
    c10 = {"min_cs2": raw["min_cs2"], "max_cs2": raw["max_cs2"],
        "min_cs2_first_law": raw["min_cs2_first_law"], "max_cs2_first_law": raw["max_cs2_first_law"],
        "acausal_point_count": raw["acausal_point_count"],
        "acausal_point_count_first_law": raw["acausal_point_count_first_law"],
        "n_points": raw["derivative_points"], "pass_": raw["C10"]["pass"]}
    c11 = {"min_dmubdnb": raw["min_dmu_b_dnb"], "n_points": raw["derivative_points"],
           "pass_": raw["C11"]["pass"]}
    return {"C1":c1,"C2":c2,"C3":composition_check("C3"),"C4":composition_check("C4"),
        "C5":c5,"C6":c6,"C7":diagnostic("C7"),"C8":diagnostic("C8"),"C9":diagnostic("C9"),
        "C10":c10,"C11":c11,
        "overall":{"status":raw["overall"]["status"],"pass_":raw["overall"]["pass"],
            "failed_checks":[x for x in raw["failed_checks"].split(",") if x],
            "finite_primary_arrays":True,
            "note":"C7-C9 are diagnostics; the native hard gate is C1-C6/C10/C11."},
        "native_report":raw}

def normalize_compose(directory: str | Path, output_csv: str | Path) -> Path:
    output=Path(output_csv).expanduser().resolve();output.parent.mkdir(parents=True,exist_ok=True)
    run_native("normalize",Path(directory).expanduser().resolve(),output)
    return output

def find_density_gap(data: EOSData, min_ratio: float = 5.0) -> int:
    if min_ratio<=1.0:raise ValueError("min_ratio must exceed one")
    ratios=data.nb[1:]/data.nb[:-1];i=int(ratios.argmax())
    if ratios[i]<min_ratio:raise ValueError(f"no density gap found (max ratio {ratios[i]:.6g} < {min_ratio})")
    return i

def fill_density_gap(*args, **kwargs):
    raise NativeFeatureUnavailable("Fortran master-function density-gap reconstruction is not ported yet.")

def reconstruct_master(*args, **kwargs):
    raise NativeFeatureUnavailable("Full master-function smoothing/candidate selection is not yet ported to Fortran.")

def match_crust_core(*args, **kwargs):
    raise NativeFeatureUnavailable("The six-method CompOSE crust-core matching API is not yet at Fortran feature parity.")
