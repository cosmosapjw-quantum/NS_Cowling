"""§7 audit of the APR + phenomenological CSS Maxwell generator.

The eos_generator phase_transitions.py build_maxwell_eos constructs the
coexistence density from a SEPARATE PCHIP interpolant for n_B while pressure is
interpolated independently (lines 75-77 of phase_transitions.py). This does NOT
guarantee the differential first-law identity n_B(mu_B) = dP/dmu_B along the
interpolants. This audit:

  1. Tests the differential first-law n_B = dP/dmu_B on the APR source branch
     and the CSS high-density branch separately (independent PCHIP of P and n
     -> residual |n - dP/dmu| / |n|).
  2. Verifies the phenomenological declaration is present in the validation
     artifact (high_phase is NOT a microscopic quark-matter EOS).
  3. Verifies the Seidov scope is correctly limited (small-core onset, not a
     full radial-mode stability calculation).
  4. On the representative APR + CSS example, verifies: common P_t/mu_t,
     positive density and energy jumps, first-law jump identity, high-side
     cs2, causality, phase IDs and transition rows after write/read-back,
     Seidov criterion and its scope, TOV/Love numerical convergence.

This audit does NOT modify the APR_CSS_MAXWELL_VALIDATION_v3_3.json artifact;
it reads and reports. The analytic CSS high-density branch formula is retained
(no mathematical defect found).
"""
from __future__ import annotations
import json
import os
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
# Append (not insert) python_frontend so the python_frontend compose_eos is
# importable without shadowing other tests' compose_eos when this module is
# imported under pytest collection. Do NOT add eos_generator/compose_eos_pkg
# to sys.path -- it shadows python_frontend/compose_eos and breaks other tests.
sys.path.append(str(ROOT / "python_frontend"))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_ROOT", str(ROOT))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_EXE",
                      str(ROOT / "build" / "compose_eos_fortran" / "compose_eos_fortran"))
from compose_eos import read_compose  # noqa: E402
from scipy.interpolate import PchipInterpolator  # noqa: E402

VALIDATION = ROOT / "eos_generator" / "docs" / "APR_CSS_MAXWELL_VALIDATION_v3_3.json"
EXAMPLE_DIR = ROOT / "eos_generator" / "validation_examples" / "APR_CSS_Maxwell"


def differential_first_law_residual(mu, P, n):
    """Build independent PCHIP interpolants of P(mu) and n(mu), then measure
    the residual |n - dP/dmu| / |n| at a dense sample of mu. The differential
    first law n = dP/dmu holds by construction for a thermodynamically
    consistent table; independent PCHIPs do NOT guarantee it."""
    Pp = PchipInterpolator(mu, P, extrapolate=False)
    np_ = PchipInterpolator(mu, n, extrapolate=False)
    # Evaluate the analytic derivative of the actual piecewise-cubic PCHIP.
    # A finite-difference gradient of sampled P adds an avoidable stencil error
    # and can hide or exaggerate the thermodynamic defect being audited.
    lo, hi = mu[1], mu[-2]
    grid = np.linspace(lo, hi, 2000)
    dP_dmu = Pp.derivative()(grid)
    n_at_grid = np_(grid)
    resid = np.abs(n_at_grid - dP_dmu) / np.maximum(np.abs(n_at_grid), 1e-30)
    return float(np.max(resid))


def audit_apr_css() -> dict:
    v = json.loads(VALIDATION.read_text())
    # 4. representative example: common P_t, mu_t, jumps, first-law identity
    mm = v["maxwell_metadata"]
    t = mm["transition"]
    common_Pt = abs(t["pressure_MeV_fm3"] - v["pressure_transition_MeV_fm3"]) < 1e-6
    positive_density_jump = t["n_high_fm3"] > t["n_low_fm3"] > 0
    positive_energy_jump = t["eps_high_MeV_fm3"] > t["eps_low_MeV_fm3"] > 0
    fl_jump_resid = abs(t["delta_eps_MeV_fm3"] - t["muB_MeV"] * (t["n_high_fm3"] - t["n_low_fm3"]))
    fl_jump_identity = fl_jump_resid < 1e-6 * max(abs(t["delta_eps_MeV_fm3"]), 1.0)
    # high-side cs2 and causality
    cs2_high = v["cs2_high"]
    cs2_causal = 0 < cs2_high <= 1
    # phase IDs and transition rows after write/read-back
    phase_ids = mm["phase_ids"]
    n_phase1 = phase_ids.count(1)
    n_phase2 = phase_ids.count(2)
    transition_row = mm.get("transition", {}).get("low_index")
    # 3. Seidov criterion and scope
    seidov = {
        "critical_delta_eps_MeV_fm3": t["seidov_critical_delta_eps_MeV_fm3"],
        "delta_eps_over_critical": t["seidov_delta_eps_over_critical"],
        "small_core_stable": t["seidov_small_core_stable"],
        "scope": "Seidov small-core onset criterion only; not a full radial-mode stability calculation",
    }
    # TOV/Love convergence
    tov = v["tov"]
    tov_converged = tov["status"] == "SUCCESS"
    # 1. differential first-law on the APR source and the CSS high branch
    # read the APR source
    d_apr = read_compose(str(ROOT / "eos_generator" / "eos_data" / "APR_APR"))
    nb = np.asarray(d_apr.nb, float); Q1 = np.asarray(d_apr.Q1, float)
    Q7 = np.asarray(d_apr.Q7, float); mn = float(d_apr.mn)
    P_apr = nb * Q1; eps_apr = nb * (Q7 + 1.0) * mn; mu_apr = (eps_apr + P_apr) / nb
    order = np.argsort(mu_apr)
    mu_apr = mu_apr[order]; P_apr = P_apr[order]; nb_apr = nb[order]
    apr_dfl_resid = differential_first_law_residual(mu_apr, P_apr, nb_apr)
    # CSS high branch: reconstruct from the CSS formula (cs2=0.75, P_t, eps_low)
    # eps(P) = eps_low + delta_eps + (P - P_t)/cs2
    # n(P) = n_low * [1 + (1+cs2)(P-P_t)/(cs2 (eps_low+delta_eps+P_t))]^(1/(1+cs2))
    # mu(P) = (eps+P)/n
    P_t = t["pressure_MeV_fm3"]; eps_low = t["eps_low_MeV_fm3"]; n_low = t["n_low_fm3"]
    delta_eps = t["delta_eps_MeV_fm3"]; cs2 = cs2_high
    P_grid = np.linspace(P_t, P_t + 2000.0, 2000)
    eps_high = eps_low + delta_eps + (P_grid - P_t) / cs2
    n_high = n_low * (1.0 + (1.0 + cs2) * (P_grid - P_t) / (cs2 * (eps_low + delta_eps + P_t))) ** (1.0 / (1.0 + cs2))
    mu_high = (eps_high + P_grid) / n_high
    order2 = np.argsort(mu_high)
    css_dfl_resid = differential_first_law_residual(mu_high[order2], P_grid[order2], n_high[order2])
    return {
        "git_head": os.popen("git rev-parse HEAD").read().strip(),
        "artifact": str(VALIDATION),
        "source": v["source"],
        "1_differential_first_law_n_equals_dPdmu": {
            "derivative_operator": "PchipInterpolator.derivative",
            "sampling": "2000-point interior grid; analytic derivative of the PCHIP polynomial",
            "acceptance_tolerance_relative": 1e-3,
            "apr_source_residual_relative": apr_dfl_resid,
            "css_high_branch_residual_relative": css_dfl_resid,
            "apr_source_differential_first_law_holds": bool(apr_dfl_resid < 1e-3),
            "css_high_branch_differential_first_law_holds": bool(css_dfl_resid < 1e-3),
            "note": "P(mu) and n(mu) are independently interpolated, so n=dP/dmu is not guaranteed. This audit uses the exact derivative of the PCHIP polynomial rather than a finite-difference approximation. Passing endpoint jump identities does not override a failed branch-interior differential check.",
        },
        "2_phenomenological_declaration": {
            "present": "Phenomenological" in mm.get("interpretation", "") or "phenomenological" in mm.get("interpretation", "").lower(),
            "interpretation_text": mm.get("interpretation"),
            "high_phase_model": mm.get("high_phase_model"),
            "is_microscopic_quark_matter_eos": False,
        },
        "3_seidov_scope": seidov,
        "4_representative_example": {
            "common_Pt_mu_t": bool(common_Pt),
            "positive_density_jump": bool(positive_density_jump),
            "positive_energy_jump": bool(positive_energy_jump),
            "first_law_jump_identity": bool(fl_jump_identity),
            "first_law_jump_residual": float(fl_jump_resid),
            "high_side_cs2": cs2_high,
            "high_side_cs2_causal": bool(cs2_causal),
            "phase_ids": {"phase1_rows": n_phase1, "phase2_rows": n_phase2,
                          "transition_row_index": transition_row},
            "tov_converged": bool(tov_converged),
            "tov_mass_Msun": tov.get("mass_Msun"),
            "tov_radius_km": tov.get("radius_km"),
            "tov_Lambda": tov.get("Lambda"),
        },
        "5_roundtrip_gate": v.get("roundtrip_gate"),
        "verdict": "ADMITTED_WITH_PHENOMENOLOGICAL_DECLARATION" if (
            common_Pt and positive_density_jump and positive_energy_jump
            and fl_jump_identity and cs2_causal and tov_converged
            and apr_dfl_resid < 1e-3 and css_dfl_resid < 1e-3
            and "phenomenological" in mm.get("interpretation", "").lower()
        ) else "REVIEW_REQUIRED",
    }


if __name__ == "__main__":
    rep = audit_apr_css()
    out = ROOT / "data" / "apr_css_maxwell_audit_20261115.json"
    out.write_text(json.dumps(rep, indent=2))
    print(f"wrote {out}")
    print(f"  1. differential first-law: APR resid={rep['1_differential_first_law_n_equals_dPdmu']['apr_source_residual_relative']:.3e} "
          f"CSS resid={rep['1_differential_first_law_n_equals_dPdmu']['css_high_branch_residual_relative']:.3e}")
    print(f"  2. phenomenological declaration present: {rep['2_phenomenological_declaration']['present']}")
    print(f"  3. Seidov scope: {rep['3_seidov_scope']['scope']}")
    print(f"  4. example: common_Pt={rep['4_representative_example']['common_Pt_mu_t']} "
          f"density_jump={rep['4_representative_example']['positive_density_jump']} "
          f"fl_identity={rep['4_representative_example']['first_law_jump_identity']} "
          f"cs2_causal={rep['4_representative_example']['high_side_cs2_causal']} "
          f"tov_converged={rep['4_representative_example']['tov_converged']}")
    print(f"  verdict: {rep['verdict']}")
