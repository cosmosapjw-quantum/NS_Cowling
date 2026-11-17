#!/usr/bin/env python3
"""Finite G-endpoint diagnosis and unadopted constructive H lower-EOS proposal.

This script writes only its owned thermo_construction directory. It does not
change source EOS, whole_crust.py, frozen parameters, or production contracts.
The exponential connector is phenomenological; no generated point is a TF state.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos.whole_crust import HermiteInterval

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

INPUTS = {
    "E_outer_surface_contract": ROOT / "data/stellar_extension_20261117/eos/fine/REPRESENTATION.json",
    "G_full_contract_admission": ROOT / "data/crust_refinement_20261117/stellar/FULL_CONTRACT_ADMISSION.json",
    "whole_crust_implementation": ROOT / "python_frontend/compose_eos/whole_crust.py",
    "F_frozen_parameters": ROOT / "data/new_model_maxwell_20261117/sources/PARAMETER_INPUTS.json",
    "A_outer_source_audit": ROOT / "docs/whole_crust_20261117/SOURCE_AUDIT_KO.md",
    "A_lower_derivation": ROOT / "docs/whole_crust_20261117/PHYSICS_REVIEW_KO.md",
    "G_source_physics": ROOT / "docs/crust_refinement_20261117/SOURCE_PHYSICS_KO.md",
}
before = {name: {"path": str(p.relative_to(ROOT)), "sha256": sha(p)} for name, p in INPUTS.items()}
ref = json.loads(INPUTS["E_outer_surface_contract"].read_text())
g = json.loads(INPUTS["G_full_contract_admission"].read_text())
outer = ref["outer"]
assert before["E_outer_surface_contract"]["sha256"] == g["reference_sha256"]

def state(n, p, epsilon):
    return {"n": float(n), "p": float(p), "epsilon": float(epsilon),
            "mu": float((epsilon + p) / n)}

outer_states = [state(n, p, n*f) for n, p, f in
                zip(outer["n_fm3"], outer["p_MeV_fm3"], outer["F_MeV"])]
left = outer_states[-1]
assert all(left[k] == g["left_state"][k] for k in ("n", "p", "epsilon"))

def endpoint_test(l, r):
    a, b, nl, nr = l["mu"], r["mu"], l["n"], r["n"]
    dp, d = r["p"]-l["p"], b-a
    ordered = bool(a > 0 and d > 0 and nr > nl > 0 and dp > 0
                   and r["epsilon"] > l["epsilon"] > 0)
    if not ordered:
        return {"ordered": False, "general_causal_strict_feasible": False,
                "cubic_causal_analytic": False}
    mean = dp/d
    yl, yr = nl/a, nr/b
    # Factored difference avoids cancellation in b**2-a**2.
    area = d*(a+b)/2
    dpmin, dpmax = yl*area, yr*area
    mlo = (4*nl+2*nr+nl*d/a)/6
    mhi = (2*nl+4*nr-nr*d/b)/6
    sl = (6*mean-4*nl-2*nr)/d
    sr = (-6*mean+2*nl+4*nr)/d
    bridge = HermiteInterval(nl, l["p"], l["epsilon"],
                             nr, r["p"], r["epsilon"])
    validator = bridge.validate()
    return {
        "ordered": True, "delta_mu_MeV": d, "delta_p_MeV_fm3": dp,
        "mean_n_fm3": mean, "n_over_mu_left": yl, "n_over_mu_right": yr,
        "n_over_mu_order_margin": yr-yl,
        "general_dp_lower_MeV_fm3": dpmin,
        "general_dp_upper_MeV_fm3": dpmax,
        "general_mean_lower_fm3": dpmin/d,
        "general_mean_upper_fm3": dpmax/d,
        "general_causal_strict_feasible": bool(yr > yl and dpmin < dp < dpmax),
        "cubic_mean_lower_fm3": mlo, "cubic_mean_upper_fm3": mhi,
        "cubic_lower_mean_margin_fm3": mean-mlo,
        "cubic_upper_mean_margin_fm3": mhi-mean,
        "cubic_left_dn_dmu": sl, "cubic_right_dn_dmu": sr,
        "cubic_left_causal_margin": a*sl-nl,
        "cubic_right_causal_margin": b*sr-nr,
        "cubic_causal_analytic": bool(mlo <= mean <= mhi),
        "existing_Hermite_validate": validator,
    }

def phi(x, k):
    x = np.asarray(x, dtype=float)
    if k == 0:
        return x
    if k > 0:
        return np.exp(k*(x-1))*(-np.expm1(-k*x))/(-math.expm1(-k))
    return np.expm1(k*x)/math.expm1(k)

def phi_prime(x, k):
    x = np.asarray(x, dtype=float)
    if k == 0:
        return np.ones_like(x)
    if k > 0:
        return k*np.exp(k*(x-1))/(-math.expm1(-k))
    return k*np.exp(k*x)/math.expm1(k)

def j_analytic(x, a, d, k):
    # Exact elementary primitive for the moderate nonzero parameters in this
    # finite diagnostic. k=0 has its exact polynomial limit.
    x = np.asarray(x, dtype=float)
    if k == 0:
        return a*x*x/2 + d*x*x*x/3
    return (a*(np.expm1(k*x)/k-x)
            + d*(((k*x-1)*np.exp(k*x)+1)/(k*k)-x*x/2))/math.expm1(k)

def construct(l, r):
    t = endpoint_test(l, r)
    assert t["general_causal_strict_feasible"]
    a, b = l["mu"], r["mu"]
    d = b-a
    yl, yr = l["n"]/a, r["n"]/b
    dy = yr-yl
    target = r["p"]-l["p"]
    def pressure_span(k):
        value, error = quad(lambda x: (a+d*x)*(yl+dy*float(phi(x, k))),
                            0., 1., epsabs=1e-14, epsrel=2e-13)
        return d*value
    lo, hi = -1., 1.
    while pressure_span(lo) < target:
        lo *= 2
        assert abs(lo) <= 128
    while pressure_span(hi) > target:
        hi *= 2
        assert abs(hi) <= 128
    k = float(brentq(lambda z: pressure_span(z)-target, lo, hi,
                     xtol=1e-13, rtol=1e-14))
    assert abs(k) > 1e-4  # exact k=0 separately verified below
    assert abs(k) < 100
    def evaluate(x):
        x = np.asarray(x, dtype=float)
        mu = a+d*x
        y = yl+dy*phi(x, k)
        yp = dy*phi_prime(x, k)/d
        n = mu*y
        susceptibility = y+mu*yp
        p = l["p"]+d*(yl*(a*x+d*x*x/2)+dy*j_analytic(x, a, d, k))
        epsilon = mu*n-p
        cs = y/(y+mu*yp)
        return {"mu": mu, "n": n, "p": p, "epsilon": epsilon,
                "dn_dmu": susceptibility, "cs2_over_c2": cs,
                "causal_margin": mu*mu*yp}
    x = np.linspace(0., 1., 1001)
    v = evaluate(x)
    assert np.all(v["dn_dmu"] > 0)
    assert np.all(v["causal_margin"] > 0)
    assert np.all((v["cs2_over_c2"] > 0) & (v["cs2_over_c2"] < 1))
    assert np.all(v["epsilon"] > 0)
    errors = {}
    for key in ("n", "p", "epsilon", "mu"):
        vals = v[key]
        errors[key] = max(abs(float(vals[0])-l[key])/max(abs(l[key]), 1e-300),
                          abs(float(vals[-1])-r[key])/max(abs(r[key]), 1e-300))
    quadrature_diffs = []
    for x0 in np.linspace(0, 1, 13):
        q, qerr = quad(lambda z: (a+d*z)*(yl+dy*float(phi(z, k))),
                       0., float(x0), epsabs=1e-14, epsrel=2e-13)
        quadrature_diffs.append(abs(l["p"]+d*q-float(evaluate(x0)["p"])))
    h = 1e-5
    fd_errors = []
    for x0 in (.1, .3, .7, .9):
        vp, vm, vc = evaluate(x0+h), evaluate(x0-h), evaluate(x0)
        dpdm = (float(vp["p"])-float(vm["p"]))/(2*h*d)
        dedm = (float(vp["epsilon"])-float(vm["epsilon"]))/(2*h*d)
        fd_errors.extend([abs(dpdm/float(vc["n"])-1),
                          abs(dedm/(float(vc["mu"])*float(vc["dn_dmu"]))-1)])
    assert max(errors.values()) < 2e-12
    assert max(quadrature_diffs) < 1e-12
    assert max(fd_errors) < 2e-6
    return {
        "family": "n(mu)=mu*(yL+(yR-yL)*expm1(k*x)/expm1(k))",
        "x": "(mu-muL)/(muR-muL)",
        "status": "PROPOSED_CONSTRUCTIVE_DIAGNOSTIC_NOT_PRODUCTION_CONTRACT",
        "composition": "UNSPECIFIED_PHENOMENOLOGICAL_LOWER_BRIDGE",
        "k": k, "k_determination": "unique endpoint pressure integral root; no stellar observable input",
        "root_bracket": [lo, hi],
        "endpoint_relative_errors": errors,
        "pressure_integral_absolute_residual": abs(pressure_span(k)-target),
        "analytic_vs_quad_max_absolute_pressure_error": max(quadrature_diffs),
        "first_law_central_difference_max_relative_error": max(fd_errors),
        "sample_count": len(x),
        "sample_min_cs2_over_c2": float(np.min(v["cs2_over_c2"])),
        "sample_max_cs2_over_c2": float(np.max(v["cs2_over_c2"])),
        "sample_min_epsilon": float(np.min(v["epsilon"])),
        "sample_min_dn_dmu": float(np.min(v["dn_dmu"])),
        "continuum_causality_status": "DERIVED: y>0 and dy/dmu>0 for all finite k",
    }

records = []
cut_records = []
source_hash_checks = []
for inp in g["inputs"]:
    source = ROOT/inp["compose_path"]
    for fn, expected in inp["source_files"].items():
        actual = sha(source/fn)
        assert actual == expected, (str(source/fn), actual, expected)
        source_hash_checks.append({"path": str((source/fn).relative_to(ROOT)), "sha256": actual})
    for trial in inp["lower_join_trials"]:
        right = state(trial["n"], trial["p"], trial["epsilon"])
        tests = endpoint_test(left, right)
        assert tests["existing_Hermite_validate"]["status"] == trial["bridge_validation"]["status"]
        assert tests["cubic_causal_analytic"] == (trial["bridge_validation"]["status"] == "PASS")
        construction = construct(left, right) if tests["general_causal_strict_feasible"] else None
        record = {"input_id": inp["id"], "mesh": inp["mesh"],
                  "source_row": trial["source_row"], "left": left, "right": right,
                  "endpoint_tests": tests, "constructive_example": construction}
        cuts = []
        for j, lcut in enumerate(outer_states):
            q = endpoint_test(lcut, right)
            cut = {
                "input_id": inp["id"], "right_source_row": trial["source_row"],
                "right_n_fm3": right["n"], "outer_source_row": j,
                "outer_n_fm3": lcut["n"], "outer_mu_MeV": lcut["mu"],
                "ordered": q["ordered"],
                "general_causal_strict_feasible": q["general_causal_strict_feasible"],
                "cubic_causal_analytic": q["cubic_causal_analytic"],
                "cubic_min_mean_margin": min(q["cubic_lower_mean_margin_fm3"],
                                           q["cubic_upper_mean_margin_fm3"]) if q["ordered"] else None,
            }
            if q["ordered"]:
                assert q["cubic_causal_analytic"] == (q["existing_Hermite_validate"]["status"] == "PASS")
            cuts.append(cut)
            cut_records.append(cut)
        admitted = [c for c in cuts if c["cubic_causal_analytic"]]
        best = max((c for c in cuts if c["ordered"]), key=lambda c: c["cubic_min_mean_margin"])
        record["earlier_outer_scan"] = {
            "count_including_fixed_last_endpoint": len(cuts),
            "general_causal_feasible_count": sum(c["general_causal_strict_feasible"] for c in cuts),
            "cubic_admitted_count": len(admitted),
            "last_cubic_admitted_outer_row": admitted[-1]["outer_source_row"] if admitted else None,
            "best_cubic_mean_margin": best,
            "policy": "diagnostic only; no outer cut or baseline mutation",
        }
        records.append(record)

# Independent k=0 elementary integral identity check by quadrature.
zero_checks = []
for a, d, x in ((938., 4., .2), (938., 4., 1.), (1., 3., .75)):
    q, err = quad(lambda z: (a+d*z)*z, 0, x, epsabs=1e-13)
    exact = float(j_analytic(x, a, d, 0))
    assert abs(q-exact) < 2e-12
    zero_checks.append({"a": a, "d": d, "x": x, "absolute_error": abs(q-exact)})
assert len(records) == 10
dd = [r for r in records if "DDME2" in r["input_id"]]
assert len(dd) == 6 and not any(r["endpoint_tests"]["cubic_causal_analytic"] for r in dd)
assert all(r["endpoint_tests"]["general_causal_strict_feasible"] for r in dd)
for name, p in INPUTS.items():
    assert sha(p) == before[name]["sha256"]

result = {
    "status": "FINITE_CONSTRUCTIVE_SUPPORT_VERIFIED",
    "base_G_commit": "905b845ad0f14d2f1e7f67f9eb802691144850cd",
    "scope": "G endpoint thermodynamics and unadopted mathematical lower connector; not H native TF, full EOS or stellar admission",
    "units": {"n": "fm^-3", "P_epsilon": "MeV fm^-3", "mu": "MeV", "sound_speed": "c_s^2/c^2=dP/depsilon"},
    "energy_shift": "NONE",
    "input_identities": before,
    "matched_source_files_verified": source_hash_checks,
    "fixed_E_outer_last_state": left,
    "outer_source_knot_count": len(outer_states),
    "baseline_trial_count": len(records),
    "earlier_outer_scan_count": len(cut_records),
    "general_causal_feasible_count": sum(r["endpoint_tests"]["general_causal_strict_feasible"] for r in records),
    "direct_cubic_admitted_count": sum(r["endpoint_tests"]["cubic_causal_analytic"] for r in records),
    "records": records,
    "k_zero_limit_quadrature_checks": zero_checks,
    "notes": [
        "Finite analytic endpoint tests matched existing Hermite.validate for all ordered scanned endpoint pairs.",
        "No fabricated microscopic TF state; no frozen parameter or endpoint energy mutation.",
        "Positive sound-speed sample extrema are diagnostics, not exact global extrema; the continuum inequality is proved separately.",
        "Endpoint admissibility does not establish the correct neutron-drip composition or a global free-energy minimum."
    ],
}
(OUT/"ENDPOINT_CONSTRUCTION.json").write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n")
with (OUT/"EARLIER_OUTER_SCAN.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(cut_records[0]))
    writer.writeheader()
    writer.writerows(cut_records)
print(json.dumps({
    "status": result["status"],
    "baseline_trials": len(records),
    "outer_trials": len(cut_records),
    "direct_cubic_admitted": result["direct_cubic_admitted_count"],
    "general_causal_feasible": result["general_causal_feasible_count"],
    "DDME2": [{"id": r["input_id"], "n_right": r["right"]["n"],
               "mean": r["endpoint_tests"]["mean_n_fm3"],
               "cubic_lower": r["endpoint_tests"]["cubic_mean_lower_fm3"],
               "generic_lower": r["endpoint_tests"]["general_mean_lower_fm3"],
               "generic_upper": r["endpoint_tests"]["general_mean_upper_fm3"],
               "k": r["constructive_example"]["k"],
               "cs_min": r["constructive_example"]["sample_min_cs2_over_c2"],
               "cs_max": r["constructive_example"]["sample_max_cs2_over_c2"],
               "earlier_cubic_admitted": r["earlier_outer_scan"]["cubic_admitted_count"],
               "best_outer_row": r["earlier_outer_scan"]["best_cubic_mean_margin"]["outer_source_row"]}
              for r in dd],
    "max_quad_difference": max(r["constructive_example"]["analytic_vs_quad_max_absolute_pressure_error"] for r in records),
    "max_FD_error": max(r["constructive_example"]["first_law_central_difference_max_relative_error"] for r in records),
}, indent=2))
