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

def fill_density_gap(data: EOSData, n_fill: int = 60, fit_points: int = 12, *,
                     gap_phase: int | None = None) -> EOSData:
    """Fill the largest density gap with a native First-Law master-function bridge.

    Source knots outside the gap are copied unchanged.  The inserted segment is
    built from F(ln n_B), with Q1=dF/dln(n_B), mu_B=F+Q1 and Q7=Q6.
    No legacy-Python numerical fallback is used.
    """
    if not isinstance(n_fill, int) or isinstance(n_fill, bool) or n_fill < 1:
        raise ValueError("n_fill must be >= 1")
    if not isinstance(fit_points, int) or isinstance(fit_points, bool) or fit_points < 2:
        raise ValueError("fit_points must be >= 2")
    igap = find_density_gap(data, min_ratio=5.0)
    if data.composition_available and data.Iphase[igap] != data.Iphase[igap + 1] and gap_phase is None:
        raise ValueError(
            f"density gap crosses phase IDs {int(data.Iphase[igap])}->{int(data.Iphase[igap + 1])}; "
            "pass gap_phase explicitly instead of silently assigning a phase"
        )
    with tempfile.TemporaryDirectory(prefix="compose-eos-fill-gap-") as td:
        root = Path(td)
        source = root / "source"
        target = root / "filled"
        source.mkdir()
        target.mkdir()
        data.write_internal_directory(source)
        command = ["fill-gap", source, target, str(n_fill), str(fit_points)]
        if gap_phase is not None:
            command.append(str(int(gap_phase)))
        run_native(*command)
        return EOSData.from_internal_directory(target, source_dir=data.source_dir)


def causal_completion(data: EOSData, *, cs2: float = 1.0,
                      n_max: float | None = None, n_points: int = 200,
                      return_metadata: bool = False):
    """Append a native constant-sound-speed tail at the first resolved cs2=1 crossing.

    This creates an explicit high-density model alternative, not a unique repair
    of the microscopic source EOS. The output intentionally omits composition.
    Set return_metadata=True to receive (EOSData, native_report).
    """
    import math
    if not isinstance(data, EOSData):
        raise TypeError("data must be an EOSData instance")
    if not isinstance(n_points, int) or isinstance(n_points, bool) or n_points < 4:
        raise ValueError("n_points must be an integer >= 4")
    sound_speed = float(cs2)
    if not math.isfinite(sound_speed) or not (0.0 < sound_speed <= 1.0):
        raise ValueError("cs2 must be finite and in (0, 1]")
    if n_max is None:
        n_max_arg = 0.0
    else:
        n_max_arg = float(n_max)
        if not math.isfinite(n_max_arg) or n_max_arg <= 0.0:
            raise ValueError("n_max must be finite and positive when specified")

    with tempfile.TemporaryDirectory(prefix="compose-eos-causal-completion-") as td:
        root = Path(td)
        source = root / "source"
        target = root / "completed"
        source.mkdir()
        target.mkdir()
        data.write_internal_directory(source)
        payload = run_native(
            "causal-complete", source, target,
            format(sound_speed, ".17g"), format(n_max_arg, ".17g"), str(n_points)
        )
        completed = EOSData.from_internal_directory(target, source_dir=data.source_dir)
    if return_metadata:
        return completed, payload
    return completed


def reconstruct_master(*args, **kwargs):
    raise NativeFeatureUnavailable("Full master-function smoothing/candidate selection is not yet ported to Fortran.")

def match_crust_core(d_cr: EOSData, d_co: EOSData, cr_name: str = "", co_name: str = "",
                    method: str = "quintic_P", n_bridge: int = 40, *,
                    core_compo_mode: str = "strict_thermo_only",
                    drop_composition: bool = False):
    """Join cold crust/core barotropes with the native six-method Fortran bridge.

    The native v1 implementation is deliberately strict-thermo-only. Any input
    carrying an eos.compo file or usable composition requires the caller to opt
    in to dropping that information with drop_composition=True.
    """
    import math
    if not isinstance(d_cr, EOSData) or not isinstance(d_co, EOSData):
        raise TypeError("d_cr and d_co must both be EOSData instances")
    methods = {"maxwell", "css", "cubic_F", "cubic_P", "quintic_F", "quintic_P"}
    if method not in methods:
        raise ValueError(f"method must be one of {sorted(methods)}")
    if not isinstance(n_bridge, int) or isinstance(n_bridge, bool) or n_bridge < 8:
        raise ValueError("n_bridge must be an integer >= 8")
    if core_compo_mode != "strict_thermo_only":
        raise NativeFeatureUnavailable(
            "The native matcher currently supports only core_compo_mode='strict_thermo_only'; "
            "literature_beta_eq reconstruction is not ported."
        )
    if not isinstance(drop_composition, bool):
        raise TypeError("drop_composition must be bool")
    crust_has_composition = bool(d_cr.composition_file_present or d_cr.composition_available)
    core_has_composition = bool(d_co.composition_file_present or d_co.composition_available)
    carries_composition = crust_has_composition or core_has_composition
    if carries_composition and not drop_composition:
        raise ValueError(
            "one or both EOS inputs carry composition data; pass drop_composition=True "
            "to explicitly discard it in the matched output"
        )
    if not math.isfinite(float(d_cr.temperature_mev)) or not math.isfinite(float(d_co.temperature_mev)):
        raise ValueError("input temperatures must be finite")
    if abs(float(d_cr.temperature_mev)) > 1e-12 or abs(float(d_co.temperature_mev)) > 1e-12:
        raise ValueError("native matching requires T=0 inputs")
    if not d_cr.beta_equilibrated or not d_co.beta_equilibrated:
        raise ValueError("native matching requires beta_equilibrated inputs")
    for label, data in (("crust", d_cr), ("core", d_co)):
        if not all(len(getattr(data, name)) == data.n for name in
                   ("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Iphase")):
            raise ValueError(f"{label} EOS arrays have inconsistent lengths")
        if not all(math.isfinite(float(v)) for v in data.nb):
            raise ValueError(f"{label} density grid contains nonfinite values")

    with tempfile.TemporaryDirectory(prefix="compose-eos-match-core-") as td:
        root = Path(td)
        crust_dir, core_dir, target_dir = root / "crust", root / "core", root / "matched"
        crust_dir.mkdir()
        core_dir.mkdir()
        target_dir.mkdir()
        d_cr.write_internal_directory(crust_dir)
        d_co.write_internal_directory(core_dir)
        payload = run_native(
            "match-core", crust_dir, core_dir, target_dir, method, str(n_bridge),
            "true" if drop_composition else "false",
        )
        matched = EOSData.from_internal_directory(target_dir, source_dir=d_cr.source_dir)

    # Preserve the legacy tuple return while adding the complete native report
    # and an independently generated C1-C11 consistency record.
    consistency = check_consistency(matched)
    dm_crust = float(d_co.mn - d_cr.mn)
    f1, f2 = float(payload["F1"]), float(payload["F2"])
    p1, p2 = float(payload["P1"]), float(payload["P2"])
    mu1, mu2 = float(payload["mu1"]), float(payload["mu2"])
    info = {
        "crust": cr_name, "core": co_name, "method": method,
        "i1": int(payload["i1"]), "i2": int(payload["i2"]), "n_bridge": n_bridge,
        "n1": float(payload["n1"]), "n2": float(payload["n2"]),
        "dnb": float(payload["n2"] - payload["n1"]),
        "P1": p1, "P2": p2, "dP": p2 - p1,
        "mu1": mu1, "mu2": mu2, "dmu": mu2 - mu1,
        "eps1": float(payload["n1"] * f1), "eps2": float(payload["n2"] * f2),
        "deps": float(payload["n2"] * f2 - payload["n1"] * f1),
        "seidov_ratio": None, "seidov_applicable": False,
        "is_exact_maxwell": False, "density_jump_fm3": 0.0,
        "phase_transition_order": None,
        "phase_transition_class": "not_assigned_no_density_jump",
        "construction_kind": payload["construction_kind"],
        "regularity_class": payload["construction_kind"],
        "cs2_1": float(payload["cs2_source_left"]),
        "cs2_2": float(payload["cs2_source_right"]),
        "cs2_left_br": float(payload["cs2_left"]),
        "cs2_right_br": float(payload["cs2_right"]),
        "dcs2_left": abs(float(payload["cs2_left"]) - float(payload["cs2_source_left"])),
        "dcs2_right": abs(float(payload["cs2_right"]) - float(payload["cs2_source_right"])),
        "cs2_min_bridge": float(payload["cs2_min"]),
        "cs2_max_bridge": float(payload["cs2_max"]),
        "first_law_relative_residual": float(payload["first_law_relative_residual"]),
        "pressure_span_bridge": float(payload["pressure_span"]),
        "spline_mode": payload["spline_mode"],
        "dm_crust": dm_crust, "dm_core": 0.0,
        "core_compo_mode": core_compo_mode,
        "composition_provenance": "none_strict_thermo_only",
        "core_composition_input_status": "explicitly_dropped" if core_has_composition else "input_eos_compo_absent",
        "core_composition_input_usable": False,
        "composition_dropped": bool(payload["composition_dropped"]),
        "diagnostic_only": not bool(consistency["overall"]["pass_"]),
        "consistency": consistency,
        "native_report": payload,
    }
    return matched, info
