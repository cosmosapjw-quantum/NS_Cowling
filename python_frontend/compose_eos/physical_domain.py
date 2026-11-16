"""Locate sampled thermodynamic support in source-preserving cold EOS tables.

Each interval diagnostic uses adjacent source rows. Gaps and phase changes end
segments. Reported ranges require interpolation/convergence checks before use
as continuous physical domains. Model attribution requires an independent
microphysical calculation with the same conventions.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np

from .eos_data import EOSData


def _finite(value: float) -> float | None:
    return float(value) if math.isfinite(float(value)) else None


def _segments(intervals: list[dict[str, Any]], flag: str) -> list[dict[str, Any]]:
    segments: list[dict[str, Any]] = []
    for item in intervals:
        if not item[flag]:
            continue
        if segments and segments[-1]["last_row"] == item["left_row"]:
            segments[-1]["last_row"] = item["right_row"]
            segments[-1]["n_max_fm3"] = item["n_right_fm3"]
        else:
            segments.append({
                "first_row": item["left_row"], "last_row": item["right_row"],
                "n_min_fm3": item["n_left_fm3"], "n_max_fm3": item["n_right_fm3"],
            })
    return segments


def physical_domain_report(
    data: EOSData, *, max_log_step: float = 0.25,
    cs2_tolerance: float = 1e-6, euler_tolerance: float = 2e-6,
    cold_energy_tolerance: float = 2e-6,
) -> dict[str, Any]:
    """Return row and interval diagnostics without truncating an EOS.

    Sound-speed ratios are dimensionless c_s^2/c^2. Both dP/d(epsilon)
    and (dP/dn)/[(epsilon+P)/n] secant diagnostics are retained.
    C2-compatible Euler identities use the CompOSE Il and lepton conventions.
    Composition closure availability is reported independently of this cold
    thermodynamic assessment.
    """
    for name, value in (("max_log_step", max_log_step),
                        ("cs2_tolerance", cs2_tolerance),
                        ("euler_tolerance", euler_tolerance),
                        ("cold_energy_tolerance", cold_energy_tolerance)):
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"{name} must be finite and nonnegative")
    if data.temperature_mev != 0.0:
        raise ValueError("physical_domain_report currently requires T=0")
    n, pressure, energy = (np.asarray(x, dtype=float) for x in (data.nb, data.p, data.eps))
    chemical = np.asarray(data.mub, dtype=float)
    if n.ndim != 1 or n.size < 2 or any(x.shape != n.shape for x in (pressure, energy, chemical)):
        raise ValueError("at least two equal-length one-dimensional EOS arrays are required")
    if np.any(~np.isfinite(n)) or np.any(n <= 0) or np.any(np.diff(n) <= 0):
        raise ValueError("density grid must be finite, positive and strictly increasing")
    rows: list[dict[str, Any]] = []
    for i in range(n.size):
        finite = bool(np.all(np.isfinite([pressure[i], energy[i], chemical[i], data.f[i]])))
        expected = (pressure[i] + data.f[i]) / n[i]
        euler_available = True
        if data.Il == 0:
            expected -= data.charge_fraction_grid * data.muq[i]
        elif abs(data.mul[i]) <= 1e-8 * data.mn:
            pass
        elif data.Ye_present[i] and data.Ymu_present[i]:
            expected -= (data.Ye[i] + data.Ymu[i]) * data.mul[i]
        else:
            euler_available = False
        cold_residual = abs(data.f[i] - energy[i]) / (n[i] * data.mn)
        cold_pass = finite and math.isfinite(cold_residual) and cold_residual <= cold_energy_tolerance
        relative = abs(chemical[i] - expected) / max(abs(chemical[i]), abs(expected), 1e-300)
        euler_pass = finite and euler_available and math.isfinite(relative) and relative <= euler_tolerance
        rows.append({
            "row": i + 1, "source_density_row_label": data.inb_start + i,
            "n_fm3": float(n[i]), "P_MeV_fm3": _finite(pressure[i]),
            "epsilon_MeV_fm3": _finite(energy[i]), "mu_B_MeV": _finite(chemical[i]),
            "finite_positive_energy_mu": bool(finite and energy[i] > 0 and chemical[i] > 0),
            "cold_free_energy_status": "PASS" if cold_pass else "FAIL",
            "cold_Q6_minus_Q7_absolute_residual": _finite(cold_residual),
            "euler_status": "PASS" if euler_pass else ("FAIL" if euler_available else "NOT_AVAILABLE"),
            "euler_relative_residual": _finite(relative) if euler_available else None,
            "composition_status": {0: "NOT_AVAILABLE", 1: "PARTIAL", 2: "CHECKED", 3: "FAIL"}.get(int(data.row_status[i]), "UNKNOWN"),
            "baryon_closure_applicable": bool(data.C3_applicable[i]),
            "charge_closure_applicable": bool(data.C4_applicable[i]),
        })
    intervals: list[dict[str, Any]] = []
    for i in range(n.size - 1):
        dn = n[i + 1] - n[i]
        dp = pressure[i + 1] - pressure[i]
        de = energy[i + 1] - energy[i]
        enthalpy = ((energy[i] + pressure[i]) / n[i] +
                    (energy[i + 1] + pressure[i + 1]) / n[i + 1]) / 2
        reasons: list[str] = []
        if data.Iphase[i] != data.Iphase[i + 1]:
            reasons.append("PHASE_BOUNDARY")
        if math.log(n[i + 1] / n[i]) > max_log_step:
            reasons.append("COARSE_SAMPLING")
        if not all(rows[j]["finite_positive_energy_mu"] for j in (i, i + 1)):
            reasons.append("NONFINITE_OR_NONPOSITIVE_ENERGY_MU")
        sound = dp / de if de > 0 else math.nan
        sound_fl = (dp / dn) / enthalpy if enthalpy > 0 else math.nan
        pressure_roundoff = 64 * np.finfo(float).eps * max(abs(pressure[i]), abs(pressure[i + 1]), np.finfo(float).tiny)
        if de <= 0 or dp < -pressure_roundoff:
            reasons.append("MECHANICAL_INSTABILITY")
        if not math.isfinite(sound) or not math.isfinite(sound_fl):
            reasons.append("SOUND_SPEED_UNDEFINED")
        elif max(sound, sound_fl) > 1.0 + cs2_tolerance:
            reasons.append("ACAUSAL_SECANT")
        if math.isfinite(sound) and min(sound, sound_fl) < -cs2_tolerance:
            reasons.append("NEGATIVE_SOUND_SPEED_SQUARED")
        continuity_findings = [reason for reason in reasons if reason in ("PHASE_BOUNDARY", "COARSE_SAMPLING")]
        physical_findings = [reason for reason in reasons if reason not in continuity_findings]
        mechanical = not reasons
        euler_pass = all(rows[j]["euler_status"] == "PASS" for j in (i, i + 1))
        cold_pass = all(rows[j]["cold_free_energy_status"] == "PASS" for j in (i, i + 1))
        intervals.append({
            "left_row": i + 1, "right_row": i + 2,
            "n_left_fm3": float(n[i]), "n_right_fm3": float(n[i + 1]),
            "pressure_monotonic_roundoff_MeV_fm3": _finite(pressure_roundoff),
            "dP_dnb_MeV": _finite(dp / dn), "cs2_over_c2_secant": _finite(sound),
            "cs2_over_c2_first_law_secant": _finite(sound_fl),
            "mechanically_causal_sampled": mechanical,
            "thermodynamically_supported_sampled": mechanical and euler_pass and cold_pass,
            "reasons": reasons,
            "physical_findings": physical_findings,
            "continuity_findings": continuity_findings,
            "secant_diagnostic_pass": not physical_findings,
            "continuity_status": "UNRESOLVED" if continuity_findings else "ADJACENT_SAMPLED_ROWS",
            "endpoint_euler_pass": euler_pass,
            "endpoint_cold_free_energy_pass": cold_pass,
            "baryon_chemical_potential_decreases": bool(chemical[i + 1] < chemical[i] - 1e-10 * dn),
        })
    first: dict[str, Any] = {}
    for item in intervals:
        for reason in item["reasons"]:
            if reason not in first:
                first[reason] = {key: item[key] for key in (
                    "left_row", "right_row", "n_left_fm3", "n_right_fm3",
                    "cs2_over_c2_secant", "cs2_over_c2_first_law_secant")}
    failed_rows = [row["row"] for row in rows if row["euler_status"] == "FAIL"]
    unavailable_rows = [row["row"] for row in rows if row["euler_status"] == "NOT_AVAILABLE"]
    return {
        "schema": "ns-cowling.sampled-physical-domain.v1",
        "assessment": "DENSITY_RESOLVED_DIAGNOSTICS",
        "units": {"density": "fm^-3", "pressure_energy": "MeV fm^-3", "chemical_potential": "MeV"},
        "source_support_fm3": [float(n[0]), float(n[-1])],
        "tolerances": {"max_log_step": max_log_step, "cs2": cs2_tolerance, "euler": euler_tolerance, "cold_Q6_Q7": cold_energy_tolerance},
        "mechanically_causal_segments": _segments(intervals, "mechanically_causal_sampled"),
        "thermodynamically_supported_segments": _segments(intervals, "thermodynamically_supported_sampled"),
        "first_interval_findings": first,
        "pressure_roundoff_policy": "Monotonic pressure permits 64*float64 epsilon*max(abs(P_left),abs(P_right),float64 tiny); physical sound-speed tolerance is unchanged.",
        "coarse_sampling_policy": "Split connected interpolation ranges; retain endpoint values and secant diagnostics. Coarse spacing alone supplies no physical exclusion.",
        "sampled_secant_pass_count": sum(item["secant_diagnostic_pass"] for item in intervals),
        "euler_failed_rows": failed_rows, "euler_unavailable_rows": unavailable_rows,
        "composition": {
            "file_present": bool(data.composition_file_present),
            "available": bool(data.composition_available),
            "status_counts": {name: sum(row["composition_status"] == name for row in rows)
                              for name in ("NOT_AVAILABLE", "PARTIAL", "CHECKED", "FAIL", "UNKNOWN")},
        },
        "chemical_monotonicity_scope": "Baryon chemical-potential decrease is recorded independently. On charged or composition-changing paths its relation to compressibility requires the appropriate constrained chemical potential.",
        "attribution": "SOURCE_TABLE_DIAGNOSTIC; intrinsic model versus numerical or convention error requires independent reproduction",
        "continuous_domain_status": "UNRESOLVED; secants describe adjacent sampled rows and do not certify an interpolant",
        "decision": "Retain every source row and investigate supported segments before model-level exclusion",
        "rows": rows, "intervals": intervals,
    }
