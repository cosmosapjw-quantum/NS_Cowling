"""Emit CLD inner-crust reconstructions in CompOSE format.

Writes eos.nb, eos.t, eos.yq, eos.thermo (CompOSE manual v3.01 section 4.2),
a human-readable summary and a machine-readable ``reconstruction_receipt.json``.

Thermodynamics.  The CLD minimizer returns, per baryon density n_B, the cell
energy density E_cell measured from the rest mass of the *model's own*
nucleons (solve_homogeneous subtracts model.mn*n_n + model.mp*n_p; electrons
include m_e).  The total energy density is therefore

    eps(n_B) = E_cell + m_n n_n + m_p n_p ,   n_p = y_p,cell n_B .

Pressure and baryon chemical potential follow from one declared interpolant:
R(n) = eps - m_n n_B is represented by a not-a-knot cubic spline, mu_B =
m_n + R'(n_B) and P = n_B mu_B - eps (Euler exact by construction).  The
table is refused (CLDThermodynamicsError, nothing written) unless every cell
converged and n_B, eps, mu_B (convexity / mechanical stability) and P are
strictly increasing with P > 0 and row secants dP/deps <= 1.

Every Q column is encoded against the model's reference masses, which are
also written in the header; normalization never changes eps, P or mu_B.
eos.compo, eos.init and eos.mr are intentionally not written: the former
cluster encoding was not a CompOSE composition record.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from .rmf_solver import RMFModel
from .cld_cell import minimize_cld_cell, CLDCellResult

COMPOSE_COLUMNS = "iT inb iYq Q1 Q2 Q3 Q4 Q5 Q6 Q7 Nadd"
# Optimizer box used by cld_cell.minimize_cld_cell (rows on it are flagged).
R_D_BOUNDS = (2.0, 40.0)
F_BOUNDS = (1e-4, 0.95)


class CLDThermodynamicsError(ValueError):
    """A CLD table cannot be emitted as an admissible cold barotrope."""


def reconstruct_inner_crust(model: RMFModel, sigma_0: float, b: float,
                            rho_min: float = 0.001, rho_max: float = 0.08,
                            n_points: int = 20, verbose: bool = False,
                            densities=None) -> list[CLDCellResult]:
    """Run the CLD minimization on a log-spaced grid (or explicit ``densities``)."""
    if densities is None:
        if n_points == 1:
            densities = [rho_min]
        else:
            lo, hi = math.log(rho_min), math.log(rho_max)
            densities = [math.exp(lo + (hi - lo) * i / (n_points - 1)) for i in range(n_points)]
    return [minimize_cld_cell(model, float(rho), sigma_0, b, verbose=verbose)
            for rho in densities]


def thermodynamic_table(results: list[CLDCellResult], *, mn: float, mp: float) -> dict:
    """Derive (n, eps, mu_B, P) from minimized cells; raise if inadmissible."""
    if len(results) < 4:
        raise CLDThermodynamicsError("at least four CLD cells are required")
    bad = [i for i, r in enumerate(results) if not r.converged]
    if bad:
        raise CLDThermodynamicsError(f"CLD cells did not converge at rows {bad}")
    n = np.array([r.rho_B for r in results], dtype=float)
    yp = np.array([r.y_p_cell for r in results], dtype=float)
    e_cell = np.array([r.E_cell for r in results], dtype=float)
    if not (np.all(np.isfinite(n)) and np.all(np.isfinite(yp)) and np.all(np.isfinite(e_cell))):
        raise CLDThermodynamicsError("nonfinite CLD cell values")
    if np.any(n <= 0.0) or np.any(np.diff(n) <= 0.0):
        raise CLDThermodynamicsError("baryon densities must be positive and strictly increasing")
    eps = e_cell + mn * n * (1.0 - yp) + mp * n * yp
    spline = CubicSpline(n, eps - mn * n, bc_type="not-a-knot")
    mu = mn + spline(n, 1)
    p = n * mu - eps
    problems = []
    if np.any(np.diff(eps) <= 0.0):
        problems.append(f"energy density not monotone at rows {np.flatnonzero(np.diff(eps) <= 0).tolist()}")
    if np.any(np.diff(mu) <= 0.0):
        problems.append("mechanical stability violated (d mu_B/d n_B <= 0, energy not convex) at rows "
                        f"{np.flatnonzero(np.diff(mu) <= 0).tolist()}")
    if np.any(p <= 0.0):
        problems.append(f"nonpositive pressure at rows {np.flatnonzero(p <= 0).tolist()}")
    if np.any(np.diff(p) <= 0.0):
        problems.append(f"pressure not monotone at rows {np.flatnonzero(np.diff(p) <= 0).tolist()}")
    with np.errstate(divide="ignore", invalid="ignore"):
        secant = np.diff(p) / np.diff(eps)
    if np.any(~np.isfinite(secant)) or np.any(secant > 1.0):
        problems.append("row secant dP/deps exceeds 1 or is undefined")
    if problems:
        raise CLDThermodynamicsError("; ".join(problems))
    mu_solver = np.array([r.mu_n for r in results], dtype=float)
    mu_gas = np.array([getattr(r, "mu_n_gas", float("nan")) for r in results], dtype=float)
    gas_ok = np.isfinite(mu_gas)
    return {
        "gas_mu_n_minus_derived_mu_B_max_abs_MeV":
            float(np.max(np.abs(mu_gas[gas_ok] - mu[gas_ok]))) if gas_ok.any() else None,
        "gas_mu_n_rows": int(gas_ok.sum()),
        "n": n, "eps": eps, "mu": mu, "P": p, "yp": yp,
        "mu_e": np.array([r.mu_e for r in results], dtype=float),
        "solver_cluster_mu_n_minus_derived_mu_B_max_abs_MeV": float(np.max(np.abs(mu_solver - mu))),
        "interpolant": "not-a-knot cubic spline of eps - m_n n_B; mu_B = m_n + derivative; P = n_B mu_B - eps",
    }


def _box_bound_rows(results: list[CLDCellResult]) -> list[int]:
    rows = []
    for i, r in enumerate(results):
        on = any(abs(r.R_d - v) <= 1e-9 * v for v in R_D_BOUNDS) or \
             any(abs(r.f - v) <= 1e-9 * v for v in F_BOUNDS)
        if on:
            rows.append(i)
    return rows


def emit_compose_set(results: list[CLDCellResult], outdir: Path, model_key: str, *,
                     mn: float, mp: float, provenance: dict | None = None) -> dict:
    """Write the CompOSE thermodynamic file set; raises before writing if inadmissible.

    Rows whose minimum lies on the optimizer box (f or R_d on a bound) are not
    stationary points, so the envelope-theorem derivative mu_B = d eps/d n_B is
    not thermodynamic there.  Such rows are trimmed when they form contiguous
    runs at the table ends (recorded in the receipt); a bound row inside the
    table is refused.
    """
    bound = _box_bound_rows(results)
    lo, hi = 0, len(results)
    while lo < hi and lo in bound:
        lo += 1
    while hi > lo and (hi - 1) in bound:
        hi -= 1
    interior = [i for i in bound if lo <= i < hi]
    if interior:
        raise CLDThermodynamicsError(
            f"optimizer box-bound minimum inside the table at rows {interior}")
    dropped = [i for i in bound if i < lo or i >= hi]
    results = list(results[lo:hi])
    table = thermodynamic_table(results, mn=mn, mp=mp)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    n, eps, mu, p = table["n"], table["eps"], table["mu"], table["P"]
    q1 = p / n
    q3 = mu / mn - 1.0
    q4 = -table["mu_e"] / mn          # mu_q = mu_p - mu_n = -mu_e (beta equilibrium)
    q7 = eps / (n * mn) - 1.0
    (outdir / "eos.t").write_text(" 1\n 1\n   0.0000000000000000E+00\n")
    (outdir / "eos.yq").write_text(" 1\n 1\n   0.0000000000000000E+00\n")
    (outdir / "eos.nb").write_text(" 1\n " + str(len(n)) + "\n"
                                   + "".join(f"   {x:.16E}\n" for x in n))
    lines = [f"   {mn:.10E}   {mp:.10E} 1"]
    for i in range(len(n)):
        lines.append(f" 1 {i + 1:5d} 1 {q1[i]:.16E} {0.0:.16E} {q3[i]:.16E} {q4[i]:.16E} "
                     f"{0.0:.16E} {q7[i]:.16E} {q7[i]:.16E} 0")
    (outdir / "eos.thermo").write_text("\n".join(lines) + "\n")
    mu_back = (q3 + 1.0) * mn
    euler = np.abs(mu_back - ((q7 + 1.0) * mn + q1)) / np.abs(mu_back)
    summary = [f"# CLD inner-crust reconstruction: {model_key}",
               f"# reference masses m_n={mn} m_p={mp} MeV; columns {COMPOSE_COLUMNS}",
               f"{'n_B':>12} {'P':>14} {'eps':>14} {'mu_B':>12} {'E_cell/n':>10} "
               f"{'f':>7} {'R_d':>6} {'y_p':>8}"]
    for i, r in enumerate(results):
        summary.append(f"{n[i]:12.6e} {p[i]:14.7e} {eps[i]:14.7e} {mu[i]:12.5f} "
                       f"{r.E_binding_per_baryon:10.4f} {r.f:7.4f} {r.R_d:6.2f} {r.y_p_cell:8.5f}")
    (outdir / "reconstruction_summary.txt").write_text("\n".join(summary) + "\n")
    receipt = {
        "schema": "ns-cowling.cld-reconstruction-receipt.v2",
        "model": model_key,
        "columns": COMPOSE_COLUMNS,
        "reference_masses_MeV": {"neutron": mn, "proton": mp},
        "rows": int(len(n)),
        "n_B_range_fm3": [float(n[0]), float(n[-1])],
        "thermodynamics": table["interpolant"],
        "max_euler_relative_residual": float(euler.max()),
        "gas_mu_n_minus_derived_mu_B_max_abs_MeV": table["gas_mu_n_minus_derived_mu_B_max_abs_MeV"],
        "gas_mu_n_rows": table["gas_mu_n_rows"],
        "cluster_bulk_mu_n_minus_derived_mu_B_max_abs_MeV":
            table["solver_cluster_mu_n_minus_derived_mu_B_max_abs_MeV"],
        "mu_note": ("equilibrium check is gas mu_n vs derived mu_B; the cluster bulk mu_n omits "
                    "surface and Coulomb derivative terms"),
        "dropped_box_bound_rows": dropped,
        "source_row_offset": lo,
        "box_bounds": {"R_d_fm": list(R_D_BOUNDS), "f": list(F_BOUNDS)},
        "omitted_files": ["eos.compo", "eos.init", "eos.mr"],
        "provenance": provenance or {},
    }
    (outdir / "reconstruction_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt
