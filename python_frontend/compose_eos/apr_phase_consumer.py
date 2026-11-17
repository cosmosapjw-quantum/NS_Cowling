"""Bounded APR bulk-Gibbs consumer and retained-core matching assembly.

Units are n_B [fm^-3], epsilon/P [MeV fm^-3], mu_B [MeV].
The dimensionless sound_speed_sq is c_s^2/c^2, not a speed in SI.
The supplied endpoint states must be tied to the APR functional; opaque labels
alone never establish a boundary law. No Maxwell transfer is used at these
continuous endpoints. This module does not supply a complete stellar model.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq

from .eos_data import EOSData

COMPOSITION_FIELDS = (
    "Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh", "row_status",
    "Yn_present", "Yp_present", "Ye_present", "Ymu_present",
    "Abar_present", "Zbar_present", "Yh_present", "C3_applicable", "C4_applicable",
)
ROW_FIELDS = ("nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7",
              "q_extra", "Iphase") + COMPOSITION_FIELDS


def _slice(data: EOSData, first: int, last: int) -> EOSData:
    if data.q_full.size or data.q_extra_full.size:
        raise ValueError("APR phase consumer requires a cold one-dimensional table")
    result = data.copy()
    for name in ROW_FIELDS:
        value = np.asarray(getattr(data, name))
        if value.shape != (data.n,):
            raise ValueError(f"{name} must have one value per density")
        setattr(result, name, value[first:last].copy())
    return result


class APRGibbsConsumer:
    """Thermodynamically linked Hermite interpolation on each APR phase.

    epsilon(n) and epsilon'(n)=mu_B are knot data; P=n*mu_B-epsilon.
    Distinct phase interpolants share endpoint values and slopes, but preserve
    one-sided second derivatives. Inside a segment second derivatives can also
    jump at ordinary cubic knots; this is finite-grid interpolation error.
    """

    def __init__(self, table: EOSData, endpoints: dict):
        if (not table.beta_equilibrated or table.temperature_mev != 0
                or table.Il != 1):
            raise ValueError("APR consumer requires a cold beta-equilibrium lepton table")
        self.table = table.copy()
        if any(np.asarray(getattr(table, key)).shape != (table.n,)
               for key in ("nb", "Q1", "Q3", "Q7", "Iphase")):
            raise ValueError("APR arrays must have one value per density")
        if (table.n < 6 or np.any(~np.isfinite(table.nb)) or
                np.any(table.nb <= 0) or np.any(np.diff(table.nb) <= 0)):
            raise ValueError("APR density support must be finite, positive and increasing")
        if np.asarray(table.Iphase).dtype.kind not in "iu":
            raise ValueError("APR phase labels must be exact integers")
        runs = [int(table.Iphase[0])]
        runs.extend(int(table.Iphase[i]) for i in range(1, table.n)
                    if table.Iphase[i] != table.Iphase[i-1])
        if runs != [6, 8, 7]:
            raise ValueError("explicit APR 6->8->7 phase layout is required")
        lo = int(np.flatnonzero(table.Iphase == 6)[-1])
        hi = int(np.flatnonzero(table.Iphase == 7)[0])
        self.boundary_indices = (lo, hi)
        self.boundaries = (float(table.nb[lo]), float(table.nb[hi]))
        if endpoints.get("construction") != "bulk_two_charge_Gibbs" or endpoints.get("species") != "npemu":
            raise ValueError("explicit npemu bulk two-charge Gibbs endpoint law required")
        self.endpoint_checks = []
        for name, i in (("left", lo), ("right", hi)):
            ref = endpoints[name]
            got = {"n_b": table.nb[i], "pressure": table.p[i],
                   "epsilon": table.eps[i], "mu_n": table.mub[i]}
            if any(not math.isfinite(float(ref[key])) for key in got):
                raise ValueError("APR endpoint reference must be finite")
            errors = {key: abs(float(value)-float(ref[key])) /
                      max(abs(float(ref[key])), 1e-20) for key, value in got.items()}
            if max(errors.values()) > 1e-8:
                raise ValueError(f"{name} APR endpoint does not match source-bound Gibbs state")
            self.endpoint_checks.append(errors)
        if (np.any(~np.isfinite(table.eps)) or np.any(~np.isfinite(table.mub))
                or np.any(~np.isfinite(table.p))):
            raise ValueError("APR thermodynamic columns must be finite")
        euler = abs(table.eps + table.p - table.nb*table.mub)
        if np.max(euler / np.maximum(abs(table.eps), 1.)) > 1e-10:
            raise ValueError("APR source violates epsilon+P=n*mu_B")
        self.slices = ((0, lo+1), (lo, hi+1), (hi, table.n))
        self.splines = []
        for start, stop in self.slices:
            if stop-start < 2:
                raise ValueError("each APR phase needs at least two knots including endpoints")
            self.splines.append(CubicHermiteSpline(
                table.nb[start:stop], table.eps[start:stop], table.mub[start:stop],
                extrapolate=False))
        # For a cubic epsilon, epsilon'' is linear. Checking both ends of
        # every interval establishes positive epsilon'' throughout that interval.
        # mu-n*epsilon'' has derivative -n*epsilon'''; its extrema are endpoints.
        self.derivative_bounds = []
        for spline in self.splines:
            a, b, c, _ = spline.c
            h = np.diff(spline.x)
            second = np.stack((2*b, 6*a*h+2*b))
            mu = np.stack((c, 3*a*h*h+2*b*h+c))
            density = np.stack((spline.x[:-1], spline.x[1:]))
            cs2 = density*second/mu
            if np.any(mu <= 0) or np.any(second <= 0) or np.any(cs2 > 1):
                raise ValueError("APR interpolant violates positive compressibility or causality")
            self.derivative_bounds.append({
                "cs2_endpoint_min": float(np.min(cs2)),
                "cs2_endpoint_max": float(np.max(cs2)),
                "positive_compressibility_and_causal_all_intervals": True})

    def _segment(self, n: float, side: str) -> int:
        if side not in ("left", "right"):
            raise ValueError("side must be left or right")
        if not math.isfinite(n) or not self.table.nb[0] <= n <= self.table.nb[-1]:
            raise ValueError("density outside actual APR consumer support")
        return int(np.searchsorted(self.boundaries, n, side=side))

    def evaluate(self, n: float, *, side: str = "right") -> dict:
        n = float(n)
        segment = self._segment(n, side)
        spline = self.splines[segment]
        epsilon, mu, second = (float(spline(n, order)) for order in (0, 1, 2))
        return {"n_b": n, "epsilon": epsilon, "mu_n": mu,
                "pressure": n*mu-epsilon, "dP_dn": n*second,
                "sound_speed_sq": n*second/mu, "phase_id": (6, 8, 7)[segment],
                "derivative_side": side}

    def evaluate_pressure(self, pressure: float, *, side: str = "right") -> dict:
        pressure = float(pressure)
        if not math.isfinite(pressure):
            raise ValueError("pressure must be finite")
        edge_n = (float(self.table.nb[0]), *self.boundaries, float(self.table.nb[-1]))
        if side not in ("left", "right"):
            raise ValueError("side must be left or right")
        states = [self.evaluate(n) for n in edge_n]
        edge_p = [state["pressure"] for state in states]
        # P=n*mu-epsilon loses a few ulps by cancellation. Snap only within
        # its explicit floating-point error scale, preserving either side.
        for n, state in zip(edge_n, states):
            roundoff = 64*np.finfo(float).eps*max(abs(n*state["mu_n"]), 1.)
            if abs(pressure-state["pressure"]) <= roundoff:
                return self.evaluate(n, side=side)
        if not edge_p[0] <= pressure <= edge_p[-1]:
            raise ValueError("pressure outside actual APR consumer support")
        k = int(np.searchsorted(edge_p[1:-1], pressure, side=side))
        density = brentq(lambda n: self.evaluate(n)["pressure"]-pressure,
                         edge_n[k], edge_n[k+1], xtol=5e-15, rtol=1e-14)
        return self.evaluate(density, side=side)

    def boundary_report(self) -> list[dict]:
        return [{"n_b": n, "law": "continuous_bulk_Gibbs",
                 "left": self.evaluate(n, side="left"),
                 "right": self.evaluate(n, side="right"),
                 "maxwell_density_jump": 0.0,
                 "perturbation": "continuous fluid displacement and Lagrangian pressure; equilibrium closure"}
                for n in self.boundaries]


def match_retained_apr(crust: EOSData, core: EOSData, endpoints: dict, *,
                       method: str, n_bridge: int = 40):
    """Native low-branch match followed by exact APR retained-tail assembly.

    Existing SUPPORT guards apply to the native input unchanged. The explicit
    Gibbs consumer validates the retained full core separately. No global
    cross-phase spline or undeclared stencil is authorized.
    """
    from .api import match_crust_core
    consumer = APRGibbsConsumer(core, endpoints)
    stop = consumer.boundary_indices[0]+1
    low = _slice(core, 0, stop)
    matched, info = match_crust_core(crust, low, method=method,
                                    n_bridge=n_bridge, drop_composition=True)
    i2 = info["i2"]
    prefix = matched.n - (low.n-i2)
    if prefix <= 0 or not 0 <= i2 < low.n:
        raise ValueError("native retained APR source range is invalid")
    # Preserve native mass normalization while retaining physical core columns.
    tail = _slice(core, i2, core.n)
    for name in ("Q3", "Q6", "Q7"):
        setattr(tail, name, (getattr(tail, name)+1)*core.mn/matched.mn-1)
    for name in ("Q4", "Q5"):
        setattr(tail, name, getattr(tail, name)*core.mn/matched.mn)
    assembled = matched.copy()
    for name in ROW_FIELDS:
        setattr(assembled, name, np.concatenate(
            (getattr(matched, name)[:prefix], getattr(tail, name))))
    # Existing source crust rows have unchanged physical composition. Restore
    # them exactly, while all bridge abundance availability stays false.
    crust_end = info["i1"]+1
    for name in COMPOSITION_FIELDS:
        getattr(assembled, name)[:crust_end] = getattr(crust, name)[:crust_end]
    assembled.composition_file_present = True
    assembled.composition_available = bool(crust.composition_available or core.composition_available)
    assembled.source_dir = None  # composite provenance is in the returned report
    report = {"scope": "native_low_branch_match_with_retained_continuous_APR_core",
              "method": method, "native_low_branch_report": info,
              "core_source_rows": [i2, core.n-1],
              "core_output_rows": [prefix, assembled.n-1],
              "crust_output_rows": [0, crust_end-1],
              "bridge_output_rows": [crust_end, prefix-1],
              "rows": assembled.n, "output_phase_ids": assembled.Iphase.tolist(),
              "composition": "source_rows_exact_masks_bridge_unavailable",
              "APR_boundaries": consumer.boundary_report(),
              "APR_density_support_fm3": [float(core.nb[i2]), float(core.nb[-1])],
              "full_stellar_observables": "NOT_COMPUTED_BY_THIS_ASSEMBLY",
              "whole_table_smooth_interpolation": "NOT_AUTHORIZED",
              "CSS_crust_jump": "separate_native_law" if method == "css" else "not_applicable",
              "native_report_scope": "low_branch_table_before_tail_assembly"}
    return assembled, report
