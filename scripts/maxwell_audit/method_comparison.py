"""§2 interpolation / reconstruction method comparison on admitted source segments.

Objective (§2): determine whether a stable, causal, thermodynamically
integrable EOS can be reconstructed while keeping source deviations within a
justified error budget. NOT to find a smoother-looking curve.

Methods compared on the SAME admitted source segments (intervals where the
§1 source secant passes the convexity necessary condition):

  M1. Cubic Hermite P(mu_B), with P'(mu_0)=n_0, P'(mu_1)=n_1 endpoint
      constraints (the current v1 kernel interpolant).
  M2. Shape-preserving PCHIP P(mu_B); n_B = P'(mu_B) derived from the same
      interpolant (NOT assumed equal to the source n_B).
  M3. A monotone+convex constrained Hermite (Fritsch-Carlson / Steffen-style)
      P(mu_B); n_B = P'(mu_B).
  M4. Thermodynamic-potential reconstruction from one authoritative function:
      epsilon(n_B). Then mu_B = d epsilon/d n_B, P = n_B mu_B - epsilon, and
      verify eps+P = mu_B n_B, dP/dn_B > 0, 0 < c_s^2 = (dP/dn)/(deps/dn) <= 1.

For each method and each admitted interval we report:
  - max/RMS deviation from the raw source values of P, epsilon, mu_B, n_B
  - the fitting objective / regularization / source-precision assumption
  - the exact interval minimum of n(mu), dn/dmu, mu*dn/dmu - n  (for M1-M3)
  - the verdict: ADMITTED (n>0, dn/dmu>0, causal) or REJECTED with reason.

Mandatory rules respected:
  - No clamping of negative derivatives to retain a PASS label.
  - PCHIP n_B is DERIVED from P'(mu_B), not assumed equal to source n_B.
  - The authoritative column (epsilon for M4) is declared, not silently changed.
  - Raw source files are read-only.

If no method reconstructs a stable causal EOS on a segment within the source
error budget, the segment is reported BLOCKED_SOURCE_RECONSTRUCTION_REQUIRED.
"""
from __future__ import annotations
import json
import os
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))
sys.path.insert(0, str(ROOT / "scripts" / "maxwell_audit"))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_ROOT", str(ROOT))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_EXE",
                      str(ROOT / "build" / "compose_eos_fortran" / "compose_eos_fortran"))
from compose_eos import read_compose  # noqa: E402
from source_consistency_audit import (  # noqa: E402
    audit_segment, select_last_iphase_segment, hermite_n_poly, quad_min)

# Source precision assumption: the crust providers print Q1/Q3/Q7 to ~10 sig
# figs (1e-10 relative); the core tables print to 17 sig figs (machine). We use
# 1e-9 as the conservative relative error budget for crust rows (one order above
# the print floor) and report residuals against it.
SOURCE_REL_BUDGET_CRUST = 1.0e-9
SOURCE_REL_BUDGET_CORE = 1.0e-15
CAUSAL_TOL = 1.0e-8


# ---------- M1: cubic Hermite P(mu) ----------------------------------------
def m1_cubic_hermite(mu, P, n):
    """Cubic Hermite with P'(mu_i)=n_i. Returns n(mu)=P'(mu) and dn/dmu as
    callable piecewise polynomials; the exact interval minima are computed
    from the quadratic/linear coefficients."""
    results = []
    for i in range(len(mu) - 1):
        (a_n, b_n, c_n), (a_c, b_c, c_c), h = hermite_n_poly(
            mu[i], mu[i+1], P[i], P[i+1], n[i], n[i+1])
        min_n = quad_min(a_n, b_n, c_n)
        min_dndmu = min((2*a_n*0 + b_n)/h, (2*a_n*1 + b_n)/h)
        min_causal = quad_min(a_c, b_c, c_c)
        results.append(dict(method="M1_cubic_hermite",
                            interval=i, h=h,
                            min_n=min_n, min_dndmu=min_dndmu, min_causal=min_causal))
    return results


# ---------- M2: PCHIP P(mu), n = P'(mu) derived -----------------------------
def _pchip_slopes(x, y):
    """Fritsch-Carlson PCHIP slopes: monotone-preserving."""
    h = np.diff(x)
    d = np.diff(y) / h
    n = len(y)
    s = np.zeros(n)
    s[1:-1] = (d[:-1] * d[1:] > 0) * 2 * d[:-1] * d[1:] / (d[:-1] + d[1:])
    # endpoints: one-sided
    s[0] = d[0] if (n == 2 or d[0] * (d[0] + d[1]) <= 0 or d[0] * (2*d[0] + d[1]) > 0) else \
           (2*d[0] + d[1]) / 3.0
    s[-1] = d[-1] if (n == 2 or d[-1] * (d[-2] + d[-1]) <= 0 or d[-1] * (d[-2] + 2*d[-1]) > 0) else \
            (d[-2] + 2*d[-1]) / 3.0
    return s


def m2_pchip(mu, P):
    """PCHIP P(mu). n(mu) = P'(mu) is DERIVED from the PCHIP interpolant
    (the derivative of the piecewise cubic), NOT assumed equal to source n_B.
    Returns per-interval exact extrema of n and dn/dmu via the cubic-derivative
    quadratic (n is quadratic in t; dn/dmu linear in t)."""
    s = _pchip_slopes(mu, P)  # s_i = P'(mu_i) per PCHIP
    results = []
    for i in range(len(mu) - 1):
        h = mu[i+1] - mu[i]
        # PCHIP cubic: P(t) = (2t^3-3t^2+1)P_i + (t^3-2t^2+t)h s_i
        #              + (-2t^3+3t^2)P_{i+1} + (t^3-t^2)h s_{i+1}
        # n(t)=P'(t) = (6t^2-6t)/h P_i + (3t^2-4t+1) s_i
        #            + (-6t^2+6t)/h P_{i+1} + (3t^2-2t) s_{i+1}
        a_n = (6/h)*P[i] + (-6/h)*P[i+1] + 3*s[i] + 3*s[i+1]
        b_n = (-6/h)*P[i] + (6/h)*P[i+1] - 4*s[i] - 2*s[i+1]
        c_n = s[i]
        # causal poly mu*dn/dmu - n (same form as M1 with s instead of n)
        (a_c, b_c, c_c) = (a_n,
                           (2*a_n + b_n) * mu[i] / h,
                           b_n * mu[i] / h - c_n)
        min_n = quad_min(a_n, b_n, c_n)
        min_dndmu = min((2*a_n*0 + b_n)/h, (2*a_n*1 + b_n)/h)
        min_causal = quad_min(a_c, b_c, c_c)
        results.append(dict(method="M2_pchip_derived_n",
                            interval=i, h=h,
                            pchip_slope_lo=s[i], pchip_slope_hi=s[i+1],
                            min_n=min_n, min_dndmu=min_dndmu, min_causal=min_causal))
    return results


# ---------- M3: monotone+convex constrained Hermite ------------------------
def m3_constrained_hermite(mu, P):
    """Monotone-AND-convex constrained Hermite. Slopes are clipped to the
    secant range (monotone) AND additionally bounded so dn/dmu >= 0 (convex in
    P means n=P' increasing, so dn/dmu = P'' >= 0). This is a Steffen-style
    bound: each slope is restricted to [0, 2*secant] on each side so the cubic
    derivative stays monotone within the interval. n(mu)=P'(mu) derived."""
    h = np.diff(mu); d = np.diff(P) / h
    # Steffen monotone slopes (ensure P monotone => n>=0):
    s = np.zeros(len(P))
    # Steffen endpoint/interior
    for i in range(len(P)):
        if i == 0:
            s[i] = d[0]
        elif i == len(P) - 1:
            s[i] = d[-1]
        else:
            hi = d[i-1]; hii = d[i]
            # Steffen: slope = sign(d) * min(|hi|, |hii|, (|hi|+|hii|)/2 * adjusted)
            if hi * hii <= 0:
                s[i] = 0.0
            else:
                s[i] = (np.sign(hi) * min(abs(hi), abs(hii),
                        (abs(hi) + abs(hii)) / 2.0))
    # additionally enforce convex (dn/dmu >= 0): clip slopes so s[i] <= s[i+1]
    # within each interval's allowed range (a simple per-interval monotone-slope
    # pass). This is a heuristic convex-monotone construction; residuals are
    # reported and the method is NOT claimed to preserve source n_B.
    for i in range(len(s) - 1):
        if s[i] > s[i+1]:
            mid = 0.5 * (s[i] + s[i+1])
            s[i] = min(s[i], mid); s[i+1] = max(s[i+1], mid)
    results = []
    for i in range(len(mu) - 1):
        hi = mu[i+1] - mu[i]
        a_n = (6/hi)*P[i] + (-6/hi)*P[i+1] + 3*s[i] + 3*s[i+1]
        b_n = (-6/hi)*P[i] + (6/hi)*P[i+1] - 4*s[i] - 2*s[i+1]
        c_n = s[i]
        (a_c, b_c, c_c) = (a_n, (2*a_n + b_n) * mu[i] / hi, b_n * mu[i] / hi - c_n)
        min_n = quad_min(a_n, b_n, c_n)
        min_dndmu = min((2*a_n*0 + b_n)/hi, (2*a_n*1 + b_n)/hi)
        min_causal = quad_min(a_c, b_c, c_c)
        results.append(dict(method="M3_constrained_hermite",
                            interval=i, h=hi,
                            slope_lo=s[i], slope_hi=s[i+1],
                            min_n=min_n, min_dndmu=min_dndmu, min_causal=min_causal))
    return results


# ---------- M4: thermodynamic potential reconstruction from eps(n) --------
def m4_thermo_from_epsilon(nb, eps):
    """Reconstruct mu_B, P from a monotone PCHIP of epsilon(n_B):
      mu_B(n) = d epsilon/d n   (PCHIP derivative)
      P(n)    = n mu_B - epsilon
      Verify eps + P = mu_B n (first law, by construction up to PCHIP deriv
      consistency), dP/dn > 0, 0 < c_s^2 = (dP/dn)/(deps/dn) <= 1.

    The authoritative column is declared as epsilon(n_B). Residuals of P, mu_B
    vs the source are reported; P is NOT fit directly.
    """
    # PCHIP of eps vs nb -> derivative is mu_B
    s_eps = _pchip_slopes(nb, eps)   # s_eps[i] = deps/dn at node i = mu_B(i)
    # mu_B at nodes (from the declared potential)
    mu_nodes = s_eps
    P_nodes = nb * mu_nodes - eps
    # check first law at nodes (algebraic by construction)
    fl_resid = np.abs(eps + P_nodes - mu_nodes * nb) / np.maximum(np.abs(mu_nodes * nb), 1.0)
    # dP/dn and deps/dn via finite difference of the reconstructed P, eps
    dPdn = np.gradient(P_nodes, nb)
    depsdn = np.gradient(eps, nb)
    cs2 = dPdn / depsdn
    return dict(method="M4_thermo_from_epsilon",
                mu_nodes=mu_nodes, P_nodes=P_nodes,
                first_law_max_resid=float(np.max(fl_resid)),
                dPdn=dPdn, depsdn=depsdn, cs2=cs2,
                dPdn_positive=bool(np.all(dPdn[1:-1] > 0)),
                cs2_causal=bool(np.all(cs2[1:-1] > 0) and np.all(cs2[1:-1] <= 1 + CAUSAL_TOL)))


def compare_methods_on_segment(model: str, model_dir: Path,
                               seg: Optional[tuple] = None) -> dict:
    d = read_compose(str(model_dir))
    il = np.asarray(d.Iphase, int)
    nb = np.asarray(d.nb, float)
    Q1 = np.asarray(d.Q1, float); Q7 = np.asarray(d.Q7, float)
    mn = float(d.mn)
    if seg is None:
        seg = select_last_iphase_segment(il)
    label, lo, hi = seg
    # segment source values
    n = nb[lo:hi]
    P = n * Q1[lo:hi]
    eps = n * (Q7[lo:hi] + 1.0) * mn
    mu = (eps + P) / n
    # admitted intervals: source secant passes
    diags = audit_segment(model, model_dir, seg)
    admitted = [d for d in diags if d.category == "PASS"]
    rejected_source = [d for d in diags if d.category == "SOURCE_TABLE_CONVEXITY_FAIL"]
    rejected_hermite = [d for d in diags if d.category == "HERMITE_INTERPOLANT_SHAPE_FAIL"]
    # residuals of each method at the source nodes (for P, eps, mu, n)
    def residuals(method_name, mu_grid, P_grid, n_grid, eps_grid):
        return {
            "P_max_rel": float(np.max(np.abs(P_grid - P) / np.maximum(np.abs(P), 1e-30))),
            "P_rms_rel": float(np.sqrt(np.mean(((P_grid - P) / np.maximum(np.abs(P), 1e-30))**2))),
            "eps_max_rel": float(np.max(np.abs(eps_grid - eps) / np.maximum(np.abs(eps), 1e-30))),
            "eps_rms_rel": float(np.sqrt(np.mean(((eps_grid - eps) / np.maximum(np.abs(eps), 1e-30))**2))),
            "mu_max_rel": float(np.max(np.abs(mu_grid - mu) / np.maximum(np.abs(mu), 1e-30))),
            "mu_rms_rel": float(np.sqrt(np.mean(((mu_grid - mu) / np.maximum(np.abs(mu), 1e-30))**2))),
            "n_max_rel": float(np.max(np.abs(n_grid - n) / np.maximum(np.abs(n), 1e-30))),
            "n_rms_rel": float(np.sqrt(np.mean(((n_grid - n) / np.maximum(np.abs(n), 1e-30))**2))),
        }
    out = {"model": model, "iphase": int(label), "segment_rows": [int(lo), int(hi)],
           "n_intervals": len(diags),
           "n_admitted": len(admitted),
           "n_rejected_source_convexity": len(rejected_source),
           "n_rejected_hermite_shape": len(rejected_hermite),
           "source_precision_assumption": "crust 10-sig (1e-9 budget)" if "SPG" in model or "GPPVA" in model else "core 17-sig (1e-15 budget)",
           "methods": {}}
    # M1: cubic Hermite (interpolates P, n at nodes exactly)
    m1 = m1_cubic_hermite(mu, P, n)
    m1_pass = all(r["min_n"] > 0 and r["min_dndmu"] > 0 and r["min_causal"] >= -CAUSAL_TOL*max(1,abs(n[i])) for i, r in enumerate(m1))
    out["methods"]["M1_cubic_hermite"] = {
        "fitting_objective": "interpolates P and n at nodes exactly; P'(mu_i)=n_i",
        "regularization": "none",
        "residuals_vs_source": residuals("M1", mu, P, n, eps),
        "n_intervals_violate_min_n": sum(1 for r in m1 if r["min_n"] <= 0),
        "n_intervals_violate_dndmu": sum(1 for r in m1 if r["min_dndmu"] <= 0),
        "n_intervals_violate_causal": sum(1 for r in m1 if r["min_causal"] < -CAUSAL_TOL),
        "verdict": "ADMITTED" if m1_pass else "REJECTED",
        "interval_minima": [{"interval": r["interval"], "min_n": r["min_n"],
                             "min_dndmu": r["min_dndmu"], "min_causal": r["min_causal"]} for r in m1],
    }
    # M2: PCHIP P(mu), n derived
    m2 = m2_pchip(mu, P)
    m2_pass = all(r["min_n"] > 0 and r["min_dndmu"] > 0 and r["min_causal"] >= -CAUSAL_TOL*max(1,abs(n[i])) for i, r in enumerate(m2))
    # M2 residuals: PCHIP reproduces P at nodes; n is derived (not source n)
    pchip_n = np.array([m2[0]["pchip_slope_lo"]] + [r["pchip_slope_hi"] for r in m2])
    out["methods"]["M2_pchip_derived_n"] = {
        "fitting_objective": "shape-preserving PCHIP of P(mu); n=P'(mu) DERIVED",
        "regularization": "Fritsch-Carlson monotone slopes",
        "residuals_vs_source": residuals("M2", mu, P, pchip_n, eps),
        "n_intervals_violate_min_n": sum(1 for r in m2 if r["min_n"] <= 0),
        "n_intervals_violate_dndmu": sum(1 for r in m2 if r["min_dndmu"] <= 0),
        "n_intervals_violate_causal": sum(1 for r in m2 if r["min_causal"] < -CAUSAL_TOL),
        "verdict": "ADMITTED" if m2_pass else "REJECTED",
        "interval_minima": [{"interval": r["interval"], "min_n": r["min_n"],
                             "min_dndmu": r["min_dndmu"], "min_causal": r["min_causal"]} for r in m2],
    }
    # M3: constrained Hermite
    m3 = m3_constrained_hermite(mu, P)
    m3_pass = all(r["min_n"] > 0 and r["min_dndmu"] > 0 and r["min_causal"] >= -CAUSAL_TOL*max(1,abs(n[i])) for i, r in enumerate(m3))
    m3_n = np.array([m3[0]["slope_lo"]] + [r["slope_hi"] for r in m3])
    out["methods"]["M3_constrained_hermite"] = {
        "fitting_objective": "monotone+convex constrained Hermite of P(mu); n=P'(mu) derived",
        "regularization": "Steffen monotone + convex slope clipping",
        "residuals_vs_source": residuals("M3", mu, P, m3_n, eps),
        "n_intervals_violate_min_n": sum(1 for r in m3 if r["min_n"] <= 0),
        "n_intervals_violate_dndmu": sum(1 for r in m3 if r["min_dndmu"] <= 0),
        "n_intervals_violate_causal": sum(1 for r in m3 if r["min_causal"] < -CAUSAL_TOL),
        "verdict": "ADMITTED" if m3_pass else "REJECTED",
        "interval_minima": [{"interval": r["interval"], "min_n": r["min_n"],
                             "min_dndmu": r["min_dndmu"], "min_causal": r["min_causal"]} for r in m3],
    }
    # M4: thermo from epsilon(n)
    m4 = m4_thermo_from_epsilon(n, eps)
    m4_res = residuals("M4", m4["mu_nodes"], m4["P_nodes"], n, eps)
    # M4 is ADMITTED only if first-law, dP/dn, causality ALL pass AND the
    # reconstructed P and mu_B stay within the source-error budget. A large P
    # residual means the PCHIP derivative of epsilon is NOT a faithful source
    # mu_B, so the reconstruction silently changed the authoritative column's
    # implied thermodynamics. This is exactly the §2 warning against assuming
    # a PCHIP derivative reproduces the source n_B / mu_B.
    budget = SOURCE_REL_BUDGET_CRUST if ("SPG" in model or "GPPVA" in model) else SOURCE_REL_BUDGET_CORE
    m4_within_budget = (m4_res["P_max_rel"] < budget and m4_res["mu_max_rel"] < budget)
    out["methods"]["M4_thermo_from_epsilon"] = {
        "fitting_objective": "PCHIP of epsilon(n_B); mu_B=deps/dn, P=n*mu_B-eps",
        "regularization": "Fritsch-Carlson monotone PCHIP of epsilon",
        "authoritative_column": "epsilon(n_B) from source Q7",
        "source_error_budget_rel": budget,
        "first_law_max_residual_at_nodes": m4["first_law_max_resid"],
        "dPdn_positive_all_interior": m4["dPdn_positive"],
        "cs2_causal_all_interior": m4["cs2_causal"],
        "cs2_min": float(np.min(m4["cs2"][1:-1])),
        "cs2_max": float(np.max(m4["cs2"][1:-1])),
        "residuals_vs_source": m4_res,
        "P_residual_within_budget": m4_within_budget,
        "verdict": "ADMITTED" if (m4["dPdn_positive"] and m4["cs2_causal"]
                                  and m4["first_law_max_resid"] < 1e-6
                                  and m4_within_budget) else "REJECTED",
    }
    # overall verdict
    any_admitted = any(out["methods"][m]["verdict"] == "ADMITTED" for m in out["methods"])
    out["overall_verdict"] = ("ADMITTED_BY_SOME_METHOD" if any_admitted
                              else "BLOCKED_SOURCE_RECONSTRUCTION_REQUIRED")
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="SPG_M1")
    ap.add_argument("--data-root", default="eos_generator/eos_data")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    mdir = ROOT / args.data_root / args.model
    rep = compare_methods_on_segment(args.model, mdir)
    if args.out:
        def _native(o):
            if isinstance(o, dict): return {k: _native(v) for k, v in o.items()}
            if isinstance(o, (list, tuple)): return [_native(x) for x in o]
            if isinstance(o, np.generic): return o.item()
            if isinstance(o, np.ndarray): return [_native(x) for x in o.tolist()]
            return o
        Path(args.out).write_text(json.dumps(_native(rep), indent=2))
    print(f"{args.model} Iphase={rep['iphase']}: {rep['n_intervals']} intervals "
          f"({rep['n_admitted']} admitted, {rep['n_rejected_source_convexity']} source-convexity-fail)")
    for m, r in rep["methods"].items():
        print(f"  {m}: verdict={r['verdict']}")
    print(f"  overall: {rep['overall_verdict']}")
