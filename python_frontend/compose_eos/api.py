"""Public Python API backed by the native Fortran EOS reader/writer/verifier."""
from __future__ import annotations
import json
import tempfile
from pathlib import Path
from typing import Any
from .eos_data import EOSData
from .native import run_native, NativeFeatureUnavailable, NativeExecutionError

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
    """Write a complete CompOSE table using the Fortran backend."""
    if data.composition_available and not data.composition_file_present:
        raise NativeFeatureUnavailable("composition values cannot be exported without an eos.compo sidecar")
    if any(int(value) != 0 for value in data.Iphase) and not data.composition_file_present:
        raise NativeFeatureUnavailable(
            "cannot export a phase-discontinuous barotrope without phase composition; "
            "non-default Iphase labels require an eos.compo sidecar"
        )
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


def construct_maxwell_transition(phase_low: EOSData, phase_high: EOSData, *,
                                  assume_local_charge_neutrality: bool = False,
                                  return_report: bool = False):
    """Construct a genuine single-conserved-charge Maxwell transition natively.

    Inputs are two cold, beta-equilibrated, composition-free single-phase
    barotropes. All thermodynamic interpolation, stability and causality checks
    are performed by the Fortran backend. The Python layer only serializes the
    input/output tables and exposes the native report.
    """
    if not isinstance(phase_low, EOSData) or not isinstance(phase_high, EOSData):
        raise TypeError("phase_low and phase_high must both be EOSData instances")
    if not isinstance(assume_local_charge_neutrality, bool):
        raise TypeError("assume_local_charge_neutrality must be bool")
    if not assume_local_charge_neutrality:
        raise ValueError(
            "confirm that both input phases are locally charge neutral by passing "
            "assume_local_charge_neutrality=True; the thermo-only data cannot prove this from composition"
        )
    with tempfile.TemporaryDirectory(prefix="compose-eos-maxwell-") as td:
        root = Path(td)
        low_dir, high_dir, target_dir = root / "low", root / "high", root / "maxwell"
        low_dir.mkdir()
        high_dir.mkdir()
        target_dir.mkdir()
        phase_low.write_internal_directory(low_dir)
        phase_high.write_internal_directory(high_dir)
        payload = run_native("maxwell-transition", low_dir, high_dir, target_dir, "true")
        matched = EOSData.from_internal_directory(target_dir, source_dir=phase_low.source_dir)
    consistency = check_consistency(matched)
    if not consistency["overall"]["pass_"]:
        failed = ", ".join(consistency["overall"]["failed_checks"])
        raise NativeExecutionError(
            "Fortran Maxwell output failed the independent C1-C6/C10-C11 verifier: " + failed
        )
    report = dict(payload)
    report["construction"] = "single_conserved_charge_local_neutrality_maxwell"
    report["local_charge_neutrality"] = "caller_confirmed"
    report["composition_preserved"] = False
    report["compose_export_allowed"] = False
    report["consistency"] = consistency
    report["is_exact_maxwell"] = True
    report["seidov_criterion_evaluated"] = False
    report["tov_interface_support"] = "not_verified"
    return (matched, report) if return_report else matched


def reconstruct_master(data: EOSData, smooth: bool = True,
                       smooth_sigma: float = 0.05,
                       fix_dissolution: bool = True,
                       dissolution_tol: float = 0.015, *,
                       recompute_extra_cs2: bool = False,
                       return_report: bool = False):
    """Reconstruct a cold beta-equilibrated EOS from its free-energy master function.

    The numerical transformation is performed by the native Fortran backend.
    Smoothing is a local weighted-cubic regularisation in log(n_B), not a
    compatibility claim for SciPy's global UnivariateSpline. If composition
    is available, the backend enforces baryon/charge conservation and solves
    electron-muon Fermi-gas equilibrium. Composition-free tables retain their
    supplied Q4 and composition placeholders.
    """
    import math
    if not isinstance(data, EOSData):
        raise TypeError("data must be an EOSData instance")
    if not isinstance(smooth, bool) or not isinstance(fix_dissolution, bool):
        raise TypeError("smooth and fix_dissolution must be bool")
    if not isinstance(recompute_extra_cs2, bool):
        raise TypeError("recompute_extra_cs2 must be bool")
    sigma = float(smooth_sigma)
    tolerance = float(dissolution_tol)
    if not math.isfinite(sigma) or (smooth and not (0.0 < sigma <= 0.5)):
        raise ValueError("smooth_sigma must be finite and lie in (0, 0.5] when smooth=True")
    if not math.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("dissolution_tol must be finite and non-negative")
    if data.n < 5:
        raise ValueError("reconstruct_master requires at least five source points")
    with tempfile.TemporaryDirectory(prefix="compose-eos-reconstruct-master-") as td:
        root = Path(td)
        source, target = root / "source", root / "rebuilt"
        source.mkdir()
        target.mkdir()
        data.write_internal_directory(source)
        payload = run_native(
            "reconstruct-master", source, target,
            "true" if smooth else "false", format(sigma, ".17g"),
            "true" if fix_dissolution else "false", format(tolerance, ".17g"),
            "true" if recompute_extra_cs2 else "false",
        )
        rebuilt = EOSData.from_internal_directory(target, source_dir=data.source_dir)
    consistency = check_consistency(rebuilt)
    if not consistency["overall"]["pass_"]:
        failed = ", ".join(consistency["overall"]["failed_checks"])
        raise NativeExecutionError(
            "native reconstructed EOS failed independent C1-C6/C10-C11 checks: " + failed
        )
    report = dict(payload)
    report["master_function"] = "F(log(n_B))"
    report["interpolation"] = "local_weighted_cubic_then_PCHIP" if smooth else "PCHIP_of_source_F"
    report["first_law_closure"] = "Q1=dF/dln(n_B); mu_B=F+Q1; P=n_B Q1; epsilon=n_B F"
    report["composition_provenance"] = (
        "reconstructed_beta_equilibrium" if data.composition_available
        else "source_composition_absent_and_not_invented"
    )
    report["composition_preserved"] = bool(data.composition_available)
    report["q_extra_recomputed_as_cs2"] = recompute_extra_cs2
    report["consistency"] = consistency
    return (rebuilt, report) if return_report else rebuilt


def match_crust_core(d_cr: EOSData | None = None, d_co: EOSData | None = None,
                    cr_name: str = "", co_name: str = "",
                    method: str = "quintic_P", n_bridge: int = 40, *,
                    core_compo_mode: str = "strict_thermo_only",
                    drop_composition: bool = False):
    """Join cold crust/core barotropes with the native six-method Fortran bridge.

    The native v1 implementation is deliberately strict-thermo-only. Any input
    carrying an eos.compo file or usable composition requires the caller to opt
    in to dropping that information with drop_composition=True.
    """
    import math
    if d_cr is None or d_co is None:
        raise NativeFeatureUnavailable(
            "match_crust_core requires both crust and core EOSData inputs; "
            "the native smooth thermo-only matching API is available when both are supplied."
        )
    if not isinstance(d_cr, EOSData) or not isinstance(d_co, EOSData):
        raise TypeError("d_cr and d_co must both be EOSData instances")
    if method == "maxwell":
        raise NativeFeatureUnavailable(
            "method='maxwell' is disabled until the sharp composition-aware transition is integrated into match_crust_core; "
            "the legacy smooth pressure plateau is not a first-order transition"
        )
    methods = {"css", "cubic_F", "cubic_P", "quintic_F", "quintic_P"}
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
        raise ValueError("native matching requires beta-equilibrated inputs")
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


def join_crust_core_css(matched: EOSData, *,
                        p_transition: float, epsilon_low: float, nb_low: float,
                        delta_epsilon: float, cs2: float, n_css: int = 8):
    """EOS-06 native CSS parameterized high-density branch JOIN.

    Append a source-backed CSS (Alford-Han-Prakash) high-density branch to a
    matched crust+core table, producing a three-segment crust+core+CSS table.
    The matched table (the output of match_crust_core) supplies the crust+core
    base; the CSS parameters (p_transition, epsilon_low, nb_low, delta_epsilon,
    cs2) describe the constant-sound-speed high-density phase and its optional
    first-order jump at the transition. The joined table is NOT sidecar-only:
    the CSS high-density branch is concatenated into the main EOSData and
    written/read as a single CompOSE table.

    The transition jump identity Delta_epsilon = mu_t * Delta_n_B is verified
    by the Fortran kernel (compose_eos_css.f90 initialize_css_segment); the
    joined table's pressure/density grids are asserted strictly increasing;
    the assembled sound speed is checked causal (0 < cs2 <= 1) across the
    whole join; and the joined table is re-read from the Fortran internal
    representation to confirm native read-back.

    Returns (joined, info) where info carries the css_join_report and the
    native consistency record.
    """
    import math
    if not isinstance(matched, EOSData):
        raise TypeError("matched must be an EOSData instance (the output of match_crust_core)")
    if matched.n < 2:
        raise ValueError("matched crust+core table has fewer than 2 rows")
    if not math.isfinite(float(p_transition)) or p_transition < 0:
        raise ValueError("p_transition must be finite and nonnegative")
    if epsilon_low <= 0 or nb_low <= 0:
        raise ValueError("epsilon_low and nb_low must be positive")
    if not math.isfinite(float(delta_epsilon)) or delta_epsilon < 0:
        raise ValueError("delta_epsilon must be finite and nonnegative")
    if not (0 < cs2 <= 1):
        raise ValueError("cs2 must satisfy 0 < cs2 <= 1 (causal)")
    if not isinstance(n_css, int) or isinstance(n_css, bool) or n_css < 4:
        raise ValueError("n_css must be an integer >= 4")

    with tempfile.TemporaryDirectory(prefix="compose-eos-css-join-") as td:
        root = Path(td)
        matched_dir = root / "matched"
        joined_dir = root / "joined"
        matched_dir.mkdir()
        joined_dir.mkdir()
        matched.write_internal_directory(matched_dir)
        payload = run_native(
            "css-join", matched_dir, joined_dir,
            f"{p_transition:.16e}", f"{epsilon_low:.16e}", f"{nb_low:.16e}",
            f"{delta_epsilon:.16e}", f"{cs2:.16e}", str(n_css),
        )
        joined = EOSData.from_internal_directory(joined_dir, source_dir=matched.source_dir)

    consistency = check_consistency(joined)
    info = {
        "operation": "css-join",
        "n_matched": int(payload["n_matched"]),
        "n_css": int(payload["n_css"]),
        "n_total": int(payload["n_total"]),
        "interface_class": int(payload["interface_class"]),
        "interface_class_name": ("first_order" if int(payload["interface_class"]) == 1
                                 else "continuous_unclassified"),
        "p_transition": float(payload["p_transition"]),
        "epsilon_low": float(payload["epsilon_low"]),
        "nb_low": float(payload["nb_low"]),
        "delta_epsilon": float(payload["delta_epsilon"]),
        "cs2": float(payload["cs2"]),
        "mu_transition": float(payload["mu_transition"]),
        "epsilon_high": float(payload["epsilon_high"]),
        "nb_high": float(payload["nb_high"]),
        "first_law_jump_residual": float(payload["first_law_jump_residual"]),
        "jump_verified": bool(payload["jump_verified"]),
        "cs2_min": float(payload["cs2_min"]),
        "cs2_max": float(payload["cs2_max"]),
        "nb_max": float(payload["nb_max"]),
        "causal": bool(payload["causal"]),
        "pressure_strictly_increasing": bool(payload["pressure_strictly_increasing"]),
        "consistency": consistency,
        "native_report": payload,
    }
    return joined, info
