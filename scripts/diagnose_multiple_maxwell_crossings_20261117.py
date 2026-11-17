#!/usr/bin/env python3
"""Enumerate every root of the represented source P_core(mu)-P_crust(mu).

Independent diagnostic only: no production kernel, source repair, root-policy
change, or EOS output. Each union interval has one exact cubic pressure
difference. Its polynomial roots and critical points are examined explicitly.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import scipy
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "data/maxwell_three_families_20261117"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def number(value):
    return float(value.replace("D", "E").replace("d", "e"))


def read_phase(record, component):
    """Independent bounded reader of the already pinned scalar internal CSV."""
    folder = ROOT / record["internal_path"]
    for name, expected in record["internal_sha256"].items():
        if sha(folder / name) != expected:
            raise ValueError("prepared source changed after its recorded snapshot")
    text = (folder / "eos_data.nml").read_text()
    match = re.search(r"^\s*mn\s*=\s*([^,\n]+)", text, re.M | re.I)
    if match is None:
        raise ValueError("missing neutron reference energy")
    mn = number(match.group(1))
    with (folder / "eos_data.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    rows = rows[component["start_zero_based"]:component["stop_exclusive"]]
    n = np.array([number(row["nb"]) for row in rows])
    q1 = np.array([number(row["Q1"]) for row in rows])
    q7 = np.array([number(row["Q7"]) for row in rows])
    p = n*q1
    # Same arithmetic grouping as native eos_energy, without using its code.
    epsilon = n*(1+q7)*mn
    mu = (epsilon+p)/n
    if np.any(np.diff(mu) <= 0) or len(n) < 3:
        raise ValueError("diagnostic requires a monotone admitted source component")
    spline = CubicHermiteSpline(mu, p, n, extrapolate=False)
    return {"n": n, "p": p, "epsilon": epsilon, "mu": mu, "spline": spline,
            "identity": {"model": record["model"], "family": record["family"],
                         "source_path": record["source_path"],
                         "source_table_sha256": record["source_table_sha256"],
                         "internal_path": record["internal_path"],
                         "internal_sha256": record["internal_sha256"],
                         "component": component, "mn_reference_MeV": mn,
                         "density_fm3": [float(n[0]), float(n[-1])],
                         "mu_E_support_MeV": [float(mu[0]), float(mu[-1])],
                         "pressure_support_MeV_fm3": [float(p[0]), float(p[-1])]}}


def shifted_coefficients(spline, left, width):
    """Ascending coefficients in t=(mu-left)/width for an existing cubic."""
    i = min(max(int(np.searchsorted(spline.x, left, side="right"))-1, 0), len(spline.x)-2)
    c3, c2, c1, c0 = spline.c[:, i]
    d = left-spline.x[i]
    return np.array([c0+d*(c1+d*(c2+d*c3)), width*(c1+d*(2*c2+3*d*c3)),
                     width**2*(c2+3*d*c3), width**3*c3])


def unit_roots(coefficients):
    """All real roots on [0,1], including stationary/tangent roots."""
    c = np.asarray(coefficients, float)
    if np.all(c == 0):
        raise ValueError("identical pressure polynomials on one interval: continuum of coexistence, not isolated roots")
    scale = max(float(np.max(np.abs(c))), np.finfo(float).tiny)
    normalized = c/scale
    algebraic = np.polynomial.polynomial.polyroots(normalized)
    candidates = [float(z.real) for z in algebraic if abs(z.imag) < 1e-9 and -1e-10 <= z.real <= 1+1e-10]
    derivative = np.array([c[1], 2*c[2], 3*c[3]])
    critical = np.polynomial.polynomial.polyroots(derivative/scale)
    critical = sorted(float(z.real) for z in critical if abs(z.imag) < 1e-9 and 0 < z.real < 1)
    bounds = [0., *critical, 1.]
    value = lambda t: float(np.polynomial.polynomial.polyval(t, normalized))
    # A critical point whose residual is at roundoff resolves a tangential
    # root even when its two algebraic approximants acquire tiny imaginary parts.
    for t in bounds:
        if abs(value(t)) < 256*np.finfo(float).eps:
            candidates.append(t)
    for a, b in zip(bounds[:-1], bounds[1:]):
        if value(a)*value(b) < 0:
            candidates.append(brentq(value, a, b, xtol=5e-15, rtol=1e-14))
    stationary_roots = [t for t in critical if abs(value(t)) < 256*np.finfo(float).eps]
    candidates = [min(stationary_roots, key=lambda s: abs(s-t))
                  if stationary_roots and min(abs(s-t) for s in stationary_roots) < 1e-7 else t
                  for t in candidates]
    accepted = []
    for t in sorted(min(1., max(0., t)) for t in candidates):
        if not accepted or abs(t-accepted[-1]) > 1e-7:
            accepted.append(t)
    return accepted, critical, [[float(z.real), float(z.imag)] for z in algebraic]


def analytic_controls():
    cases = {"three_simple": ([.2, .5, .8], [.2, .5, .8]),
             "one_tangent_one_simple": ([.2, .5, .5], [.2, .5]),
             "endpoints": ([0., .3, 1.], [0., .3, 1.])}
    records = []
    for name, (given, unique) in cases.items():
        c = np.polynomial.polynomial.polyfromroots(given)
        found, _, _ = unit_roots(c)
        passed = len(found) == len(unique) and np.allclose(found, unique, atol=2e-7, rtol=0)
        if not passed:
            raise AssertionError(f"polynomial root control failed: {name}: {found}")
        records.append({"name": name, "expected": unique, "computed": found, "pass": bool(passed)})
    return records


def diagnose(low, high):
    lo = max(low["mu"][0], high["mu"][0])
    hi = min(low["mu"][-1], high["mu"][-1])
    if hi <= lo:
        raise ValueError("no common chemical-potential support")
    union = np.unique(np.r_[lo, hi, low["mu"][(low["mu"]>lo)&(low["mu"]<hi)],
                            high["mu"][(high["mu"]>lo)&(high["mu"]<hi)]])
    delta = lambda mu: float(high["spline"](mu)-low["spline"](mu))
    cells, candidates, extrema = [], [], [float(lo), float(hi)]
    for a, b in zip(union[:-1], union[1:]):
        coefficients = shifted_coefficients(high["spline"], a, b-a)-shifted_coefficients(low["spline"], a, b-a)
        roots, critical, algebraic = unit_roots(coefficients)
        computed = [float(a+t*(b-a)) for t in roots]
        candidates.extend(computed)
        extrema.extend(float(a+t*(b-a)) for t in critical)
        extrema.extend([float(a), float(b)])
        cells.append({"mu_interval_MeV": [float(a), float(b)],
                      "coordinate": "t=(mu-left)/(right-left)",
                      "deltaP_ascending_coefficients_MeV_fm3": coefficients.tolist(),
                      "algebraic_roots_t_real_imag": algebraic,
                      "real_roots_t_in_unit_interval": roots,
                      "critical_points_t_in_unit_interval": critical,
                      "roots_mu_MeV": computed})
    roots = []
    for mu in sorted(candidates):
        if not roots or abs(mu-roots[-1]) > 1e-10*max(1., abs(mu)):
            roots.append(mu)
    boundaries = [float(lo), *roots, float(hi)]
    ordering = []
    for a, b in zip(boundaries[:-1], boundaries[1:]):
        inside = np.array([v for v in extrema if a <= v <= b]+[a,b])
        values = np.array([delta(v) for v in inside])
        mid = delta((a+b)/2)
        ordering.append({"mu_interval_MeV": [a,b], "deltaP_midpoint_MeV_fm3": mid,
                         "preferred_phase_at_fixed_mu": "core" if mid>0 else "crust" if mid<0 else "degenerate",
                         "min_deltaP_MeV_fm3": float(values.min()),
                         "min_at_mu_MeV": float(inside[np.argmin(values)]),
                         "max_deltaP_MeV_fm3": float(values.max()),
                         "max_at_mu_MeV": float(inside[np.argmax(values)])})
    details = []
    for i, mu in enumerate(roots):
        pc, ph = float(low["spline"](mu)), float(high["spline"](mu))
        nc, nh = float(low["spline"](mu,1)), float(high["spline"](mu,1))
        ec, eh = mu*nc-pc, mu*nh-ph
        left_gap, right_gap = mu-boundaries[i], boundaries[i+2]-mu
        probe = min(left_gap,right_gap)*.01
        details.append({"index": i, "mu_MeV": mu, "P_crust_MeV_fm3": pc, "P_core_MeV_fm3": ph,
                        "pressure_residual_MeV_fm3": ph-pc,
                        "n_crust_fm3": nc, "n_core_fm3": nh, "delta_n_core_minus_crust_fm3": nh-nc,
                        "epsilon_crust_MeV_fm3": ec, "epsilon_core_MeV_fm3": eh,
                        "latent_energy_core_minus_crust_MeV_fm3": eh-ec,
                        "latent_identity_residual_MeV_fm3": (eh-ec)-mu*(nh-nc),
                        "crust_to_core_positive_density_jump": bool(nh>nc),
                        "d_deltaP_d_mu_fm3": nh-nc,
                        "d2_deltaP_d_mu2": float(high["spline"](mu,2)-low["spline"](mu,2)),
                        "left_preferred_phase": ordering[i]["preferred_phase_at_fixed_mu"],
                        "right_preferred_phase": ordering[i+1]["preferred_phase_at_fixed_mu"],
                        "probe_delta_mu_MeV": probe,
                        "deltaP_left_probe_MeV_fm3": delta(mu-probe),
                        "deltaP_right_probe_MeV_fm3": delta(mu+probe)})
    return {"common_mu_support_MeV": [float(lo),float(hi)], "union_intervals": len(cells),
            "represented_real_root_count": len(roots), "roots": details,
            "pressure_ordering_intervals": ordering, "union_cubic_polynomials": cells,
            "deltaP_definition": "P_core(mu_E)-P_crust(mu_E)",
            "equilibrium_criterion": "The phase with larger pressure at a fixed chemical potential minimizes grand-potential density.",
            "source_values_changed": False, "physical_energy_offset_applied_MeV_per_baryon": 0.,
            "classification": "REENTRANT_PRESSURE_ORDERING_IN_SUPPLIED_SUPPORT" if len(roots)>1 else "SINGLE_ROOT",
            "single_jump_policy_note": "Selecting one positive crust-to-core density-jump root does not justify retaining the core where the supplied crust branch has higher P(mu). All common-support intervals must respect the selected phase ordering. A reentrant phase chain needs an explicit different construction and remains outside the present one-jump kernel."}


def make_plot(report, low, high, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    lo,hi=report["common_mu_support_MeV"]
    grid=np.linspace(lo,hi,1601)
    delta=high["spline"](grid)-low["spline"](grid)
    fig,ax=plt.subplots(figsize=(8.5,4.8),constrained_layout=True)
    ax.axhline(0,color="0.4",lw=1)
    ax.plot(grid,delta,lw=2,color="#124f86")
    for row in report["roots"]:
        x=row["mu_MeV"]
        ax.scatter([x],[0],s=34,color="#b33d26",zorder=3)
        ax.annotate(f"$\\mu={x:.6f}$ MeV\n$\\Delta n={row['delta_n_core_minus_crust_fm3']:+.6g}$ fm$^{{-3}}$",
                    (x,0),xytext=(0 if row["index"]%2==0 else -6,35 if row["index"]%2==0 else -48),
                    textcoords="offset points",ha="center" if row["index"]%2==0 else "right",
                    fontsize=9,arrowprops={"arrowstyle":"-","lw":.6})
    ax.set(xlabel=r"Euler chemical potential $\mu_E$ [MeV]",
           ylabel=r"$P_{\rm BBB}-P_{\rm SPG\,M4}$ [MeV fm$^{-3}$]",
           title="Native SPG M4 / BBB BHF2: all common-support pressure crossings")
    ax.grid(alpha=.2)
    fig.savefig(path,dpi=180)
    plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest",type=Path,default=BASE/"prepared/MANIFEST.json")
    parser.add_argument("--cell",type=Path,default=BASE/"campaign/native_generated/SPG_M4__BBB_BHF2/CELL.json")
    parser.add_argument("--output",type=Path,default=BASE/"validation/MULTIPLE_CROSSINGS.json")
    args=parser.parse_args()
    manifest=json.loads(args.manifest.read_text());cell=json.loads(args.cell.read_text())
    records={(row["family"],row["model"]):row for row in manifest["source_records"]}
    selected=[r for r in cell["candidate_results"] if r.get("category")=="MULTIPLE_COEXISTENCE_CROSSINGS"]
    if len(selected)!=1:raise ValueError("expected exactly one rejected multiple-crossing component pair")
    candidate=selected[0]
    low=read_phase(records[cell["family"],cell["crust"]],candidate["crust_component"])
    high=read_phase(records[cell["family"],cell["core"]],candidate["core_component"])
    controls=analytic_controls()
    report={"schema":"ns_cowling.multiple_maxwell_crossings_diagnostic.v1",
            "created_at_utc":datetime.now(timezone.utc).isoformat(),
            "script_sha256":sha(__file__),"numpy_version":np.__version__,"scipy_version":scipy.__version__,
            "prepared_manifest_sha256":sha(args.manifest),"campaign_cell_sha256":sha(args.cell),
            "native_rejection_category":candidate["category"],"native_rejection_error":candidate.get("error"),
            "algorithm":"Independent SciPy CubicHermiteSpline P(mu) with slope n; analytic coefficient translations on every union interval; exhaustive cubic polynomial roots plus stationary-point/sign-bracket controls.",
            "analytic_controls":controls,"crust":low["identity"],"core":high["identity"],**diagnose(low,high)}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    plot=args.output.with_suffix(".png")
    make_plot(report,low,high,plot)
    report["plot_path"]=str(plot.relative_to(ROOT));report["plot_sha256"]=sha(plot)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"roots":report["roots"],"ordering":report["pressure_ordering_intervals"],
                      "classification":report["classification"],"union_intervals":report["union_intervals"]},indent=2))


if __name__=="__main__":main()
