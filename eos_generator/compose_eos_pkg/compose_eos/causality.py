"""Causality diagnostics and optional maximally-stiff causal completion.

A source EOS that is acausal is never silently clipped. Diagnostics sample
within every smooth tabulated interval and never interpolate across phase
boundaries or large density gaps. A completion is an explicit model alternative,
not a unique correction implied by the source EOS.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import PchipInterpolator, CubicSpline
from scipy.optimize import brentq

from .eos_data import EOSData


def _validated_arrays(data: EOSData):
    n = np.asarray(data.nb, dtype=float)
    p = np.asarray(data.p, dtype=float)
    eps = np.asarray(data.eps, dtype=float)
    mu = np.asarray(data.mub, dtype=float)
    if (n.ndim != 1 or len(n) < 4 or p.shape != n.shape or eps.shape != n.shape
            or mu.shape != n.shape):
        raise ValueError("causality analysis requires matching one-dimensional n, P, epsilon and mu_B arrays")
    if (np.any(~np.isfinite(n)) or np.any(~np.isfinite(p))
            or np.any(~np.isfinite(eps)) or np.any(~np.isfinite(mu))):
        raise ValueError("causality input arrays must be finite")
    if (np.any(n <= 0) or np.any(p <= 0) or np.any(eps <= 0) or np.any(mu <= 0)
            or np.any(np.diff(n) <= 0)):
        raise ValueError("causality input requires positive P, epsilon, mu_B and strictly increasing n_B")
    return n, p, eps, mu


def _smooth_segments(data: EOSData, n: np.ndarray, max_log_step: float):
    dln = np.diff(np.log(n))
    if np.asarray(data.Iphase).shape == n.shape:
        phase_edges = np.flatnonzero(np.diff(data.Iphase) != 0)
    else:
        phase_edges = np.empty(0, dtype=int)
    gap_edges = np.flatnonzero(dln > max_log_step)
    cuts = sorted(set(phase_edges.tolist() + gap_edges.tolist()))
    starts = [0] + [int(i) + 1 for i in cuts]
    stops = [int(i) + 1 for i in cuts] + [len(n)]
    return [np.arange(a, b, dtype=int) for a, b in zip(starts, stops) if b > a], dln, gap_edges, phase_edges


def _cs2_from_log_splines(p_spline, e_spline, xq: np.ndarray) -> np.ndarray:
    lp = np.asarray(p_spline(xq), dtype=float)
    le = np.asarray(e_spline(xq), dtype=float)
    dp = np.asarray(p_spline.derivative()(xq), dtype=float)
    de = np.asarray(e_spline.derivative()(xq), dtype=float)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        cs2 = np.exp(lp - le) * dp / de
    cs2[(~np.isfinite(cs2)) | (de <= 0)] = np.nan
    return cs2


def _quadratic_interval_min_max(poly_quad: np.ndarray, t_lo: float, t_hi: float) -> tuple[float, float, float]:
    """Exact minimum, maximum, and argmin of a quadratic polynomial on [t_lo, t_hi]."""
    cands = [t_lo, t_hi]
    a, b, _ = [float(v) for v in poly_quad]
    if abs(a) > 1e-30 * max(1.0, abs(b)):
        t_vert = -b / (2.0 * a)
        if t_lo < t_vert < t_hi:
            cands.append(float(t_vert))
    pts = np.asarray(cands, dtype=float)
    vals = np.polyval(poly_quad, pts)
    idx_min = int(np.argmin(vals))
    return float(vals[idx_min]), float(np.max(vals)), float(pts[idx_min])


def _cs2_stationary_points_in_interval(
    poly_dp: np.ndarray,
    poly_ddp: np.ndarray,
    poly_de: np.ndarray,
    poly_dde: np.ndarray,
    t_lo: float,
    t_hi: float,
) -> tuple[np.ndarray, list[float]]:
    """Construct the degree-<=6 stationary polynomial N(t)=0 and isolate real roots in (t_lo, t_hi).

    Stationary points of cs2 = exp(p - e) * p' / e' satisfy:
        N(t) = (p' - e') * p' * e' + p'' * e' - e'' * p' = 0.
    """
    poly_N = np.polyadd(
        np.polymul(np.polysub(poly_dp, poly_de), np.polymul(poly_dp, poly_de)),
        np.polysub(np.polymul(poly_ddp, poly_de), np.polymul(poly_dde, poly_dp)),
    )
    roots_in_interval: list[float] = []
    scale = float(np.max(np.abs(poly_N))) if poly_N.size else 0.0
    if scale > 0.0:
        nz = np.flatnonzero(np.abs(poly_N) > 1e-13 * scale)
        if nz.size > 0 and nz[0] < len(poly_N) - 1:
            trimmed = poly_N[nz[0]:]
            span = max(1.0, t_hi - t_lo)
            for r in np.roots(trimmed):
                if abs(r.imag) <= 1e-10 * span and np.isfinite(r.real):
                    tr = float(r.real)
                    if t_lo < tr < t_hi:
                        roots_in_interval.append(tr)
    roots_in_interval.sort()
    dedup: list[float] = []
    rtol = 1e-12 * max(1.0, t_hi - t_lo)
    for r in roots_in_interval:
        if not dedup or abs(r - dedup[-1]) > rtol:
            dedup.append(r)
    return poly_N, dedup


def _analyze_log_spline_intervals(
    ps,
    es,
    xs: np.ndarray,
    *,
    n_min: float = 0.0,
    max_log_step: float = 0.25,
    tolerance: float = 1e-6,
    contact_tol: float = 1e-9,
) -> dict:
    """Exact interval-extrema analysis of cs2(x) over one smooth log-spline segment."""
    min_cs2 = float("inf")
    max_cs2 = float("-inf")
    max_cs2_nb = None
    min_cs2_nb = None
    min_dp = float("inf")
    min_de = float("inf")
    dp_positive = True
    de_positive = True
    nonpositive_deriv_nb: list[float] = []

    crossing_events: list[dict] = []
    tangential_contacts: list[dict] = []
    acausal_intervals: list[dict] = []
    interval_reports: list[dict] = []

    def _eval_interval_cs2(a_poly, b_poly, dp_poly, de_poly, t_val: float) -> float:
        p_val = float(np.polyval(a_poly, t_val))
        e_val = float(np.polyval(b_poly, t_val))
        dp_v = float(np.polyval(dp_poly, t_val))
        de_v = float(np.polyval(de_poly, t_val))
        if de_v <= 0.0 or not np.isfinite(p_val + e_val + dp_v + de_v):
            return float("nan")
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            val = float(np.exp(p_val - e_val) * dp_v / de_v)
        return val if np.isfinite(val) else float("nan")

    x_min_filter = float(np.log(n_min)) if n_min > 0.0 else float("-inf")
    stationary_root_count = 0

    for j in range(len(xs) - 1):
        h = float(xs[j + 1] - xs[j])
        if h > max_log_step:
            continue
        x0 = float(xs[j])
        if x0 + h < x_min_filter:
            continue
        t_lo = max(0.0, x_min_filter - x0)
        if t_lo >= h - 1e-14:
            continue

        a3, a2, a1, a0 = [float(v) for v in ps.c[:, j]]
        b3, b2, b1, b0 = [float(v) for v in es.c[:, j]]
        poly_p = np.array([a3, a2, a1, a0], dtype=float)
        poly_dp = np.array([3.0 * a3, 2.0 * a2, a1], dtype=float)
        poly_ddp = np.array([6.0 * a3, 2.0 * a2], dtype=float)
        poly_e = np.array([b3, b2, b1, b0], dtype=float)
        poly_de = np.array([3.0 * b3, 2.0 * b2, b1], dtype=float)
        poly_dde = np.array([6.0 * b3, 2.0 * b2], dtype=float)

        dp_min, _, dp_arg = _quadratic_interval_min_max(poly_dp, t_lo, h)
        de_min, _, de_arg = _quadratic_interval_min_max(poly_de, t_lo, h)
        min_dp = min(min_dp, dp_min)
        min_de = min(min_de, de_min)
        if dp_min <= 0.0:
            dp_positive = False
            nonpositive_deriv_nb.append(float(np.exp(x0 + dp_arg)))
        if de_min <= 0.0:
            de_positive = False
            nonpositive_deriv_nb.append(float(np.exp(x0 + de_arg)))

        _, st_roots = _cs2_stationary_points_in_interval(
            poly_dp, poly_ddp, poly_de, poly_dde, t_lo, h
        )
        stationary_root_count += len(st_roots)

        partition = [t_lo, *st_roots, h]
        vals = [_eval_interval_cs2(poly_p, poly_e, poly_dp, poly_de, t_v) for t_v in partition]
        finite_pairs = [(t_v, c_v) for t_v, c_v in zip(partition, vals) if np.isfinite(c_v)]
        if not finite_pairs:
            continue

        iv_min = float(min(c_v for _, c_v in finite_pairs))
        iv_max = float(max(c_v for _, c_v in finite_pairs))
        for t_v, c_v in finite_pairs:
            nb_pt = float(np.exp(x0 + t_v))
            if c_v < min_cs2:
                min_cs2 = c_v
                min_cs2_nb = nb_pt
            if c_v > max_cs2:
                max_cs2 = c_v
                max_cs2_nb = nb_pt

        # Because every interior stationary point splits [t_lo, h] into strictly monotone
        # subintervals, any cs2 = 1 crossing in (t_r, t_{r+1}) is isolated uniquely.
        for r in range(len(partition) - 1):
            t_a, t_b = partition[r], partition[r + 1]
            c_a, c_b = vals[r], vals[r + 1]
            if not (np.isfinite(c_a) and np.isfinite(c_b)):
                continue
            f_a, f_b = c_a - 1.0, c_b - 1.0
            if (f_a < -contact_tol and f_b > tolerance) or (f_a <= tolerance and f_b > tolerance and c_a <= 1.0):
                root_t = t_a if abs(f_a) <= 1e-14 else float(
                    brentq(
                        lambda z: _eval_interval_cs2(poly_p, poly_e, poly_dp, poly_de, z) - 1.0,
                        t_a, t_b, xtol=1e-14, rtol=1e-14,
                    )
                )
                crossing_events.append({
                    "kind": "true_crossing",
                    "direction": "upward",
                    "interval_index": int(j),
                    "x": float(x0 + root_t),
                    "t": float(root_t),
                    "nb": float(np.exp(x0 + root_t)),
                    "n_left": float(np.exp(x0 + t_a)),
                    "n_right": float(np.exp(x0 + t_b)),
                })
            elif f_a > tolerance and f_b < -contact_tol:
                root_t = float(
                    brentq(
                        lambda z: _eval_interval_cs2(poly_p, poly_e, poly_dp, poly_de, z) - 1.0,
                        t_a, t_b, xtol=1e-14, rtol=1e-14,
                    )
                )
                crossing_events.append({
                    "kind": "true_crossing",
                    "direction": "downward",
                    "interval_index": int(j),
                    "x": float(x0 + root_t),
                    "t": float(root_t),
                    "nb": float(np.exp(x0 + root_t)),
                    "n_left": float(np.exp(x0 + t_a)),
                    "n_right": float(np.exp(x0 + t_b)),
                })

        # Detect smooth interior tangential contacts at cs2 = 1 where cs2 does not cross 1.
        for r in range(1, len(partition) - 1):
            c_mid = vals[r]
            c_prev, c_next = vals[r - 1], vals[r + 1]
            if not (np.isfinite(c_mid) and np.isfinite(c_prev) and np.isfinite(c_next)):
                continue
            if abs(c_mid - 1.0) <= contact_tol:
                if (c_prev < 1.0 - 1e-12 and c_next < 1.0 - 1e-12) or (c_prev > 1.0 + 1e-12 and c_next > 1.0 + 1e-12):
                    tangential_contacts.append({
                        "kind": "tangential_contact",
                        "interval_index": int(j),
                        "x": float(x0 + partition[r]),
                        "nb": float(np.exp(x0 + partition[r])),
                        "cs2": float(c_mid),
                        "side": "from_below" if c_prev < 1.0 else "from_above",
                    })

        if iv_max > 1.0 + tolerance:
            lo_ac = float(np.exp(x0 + finite_pairs[0][0]))
            hi_ac = float(np.exp(x0 + finite_pairs[-1][0]))
            for ev in crossing_events:
                if ev["interval_index"] == j and ev["direction"] == "upward":
                    lo_ac = ev["nb"]
                    break
            acausal_intervals.append({
                "interval_index": int(j),
                "nb_lo": lo_ac,
                "nb_hi": hi_ac,
                "max_cs2": iv_max,
            })

        interval_reports.append({
            "interval_index": int(j),
            "nb_left": float(np.exp(x0 + t_lo)),
            "nb_right": float(np.exp(x0 + h)),
            "min_cs2": iv_min,
            "max_cs2": iv_max,
            "n_stationary": int(len(st_roots)),
            "min_dp_dlnn": dp_min,
            "min_de_dlnn": de_min,
        })

    upward = [ev for ev in crossing_events if ev["direction"] == "upward"]
    first_acausal_nb = None
    acausal_nb_range = None
    if acausal_intervals:
        first_acausal_nb = float(upward[0]["nb"]) if upward else float(acausal_intervals[0]["nb_lo"])
        acausal_nb_range = [float(first_acausal_nb), float(acausal_intervals[-1]["nb_hi"])]
    elif np.isfinite(max_cs2) and max_cs2 > 1.0 + tolerance:
        first_acausal_nb = float(max_cs2_nb)
        acausal_nb_range = [float(max_cs2_nb), float(max_cs2_nb)]

    if crossing_events:
        crossing_kind = "true_crossing"
    elif tangential_contacts:
        crossing_kind = "tangential_contact"
    elif np.isfinite(min_cs2) and min_cs2 >= 1.0 + tolerance:
        crossing_kind = "always_superluminal"
    else:
        crossing_kind = "none"

    return {
        "min_cs2": float(min_cs2) if np.isfinite(min_cs2) else None,
        "max_cs2": float(max_cs2) if np.isfinite(max_cs2) else None,
        "min_cs2_nb": min_cs2_nb,
        "max_cs2_nb": max_cs2_nb,
        "dp_positive": bool(dp_positive),
        "de_positive": bool(de_positive),
        "min_dp_dlnn": float(min_dp) if np.isfinite(min_dp) else None,
        "min_de_dlnn": float(min_de) if np.isfinite(min_de) else None,
        "nonpositive_derivative_nb": nonpositive_deriv_nb,
        "stationary_point_count": int(stationary_root_count),
        "interval_count": int(len(interval_reports)),
        "crossing_kind": crossing_kind,
        "crossing_events": crossing_events,
        "tangential_contacts": tangential_contacts,
        "acausal": bool(len(acausal_intervals) > 0 or (np.isfinite(max_cs2) and max_cs2 > 1.0 + tolerance)),
        "acausal_interval_count": int(len(acausal_intervals)),
        "first_acausal_nb": first_acausal_nb,
        "acausal_nb_range": acausal_nb_range,
        "intervals": interval_reports,
    }


def causality_profile(
    data: EOSData, *, n_min=0.0, max_log_step=0.25,
    tolerance=1e-6, samples_per_interval=16
):
    """Compare PCHIP, natural-cubic and centered-FD c_s^2 diagnostics.

    For piecewise-cubic splines (PCHIP and natural cubic), exact interval
    extrema, derivative positivity (p'>0, e'>0), and true vs tangential cs2=1
    crossings are evaluated analytically across every admissible interval via
    the degree-<=6 stationary polynomial (p'-e')*p'*e' + p''*e' - e''*p' = 0.
    The uniform sub-grid controlled by ``samples_per_interval`` is computed only
    for plotting/diagnostics and is never the hard acceptance criterion.
    """
    n, p, eps, _ = _validated_arrays(data)
    if (not np.isfinite(n_min) or n_min < 0 or not np.isfinite(max_log_step)
            or max_log_step <= 0 or not np.isfinite(tolerance) or tolerance < 0
            or int(samples_per_interval) < 2):
        raise ValueError("invalid causality-profile sampling configuration")
    samples_per_interval = int(samples_per_interval)
    x = np.log(n)
    lp = np.log(p)
    le = np.log(eps)
    segments, dln, gap_edges, phase_edges = _smooth_segments(data, n, max_log_step)
    out = {}

    for name in ("pchip", "cubic_spline", "centered_fd"):
        if name == "centered_fd":
            cs_values = []
            n_values = []
            for ids in segments:
                if len(ids) < 3:
                    continue
                xs = x[ids]
                dp = np.gradient(lp[ids], xs, edge_order=2)
                de = np.gradient(le[ids], xs, edge_order=2)
                with np.errstate(divide="ignore", invalid="ignore"):
                    cs = (p[ids] / eps[ids]) * dp / de
                valid = (n[ids] >= n_min) & np.isfinite(cs) & (de > 0)
                cs_values.extend(cs[valid].tolist())
                n_values.extend(n[ids][valid].tolist())
            cs_all = np.asarray(cs_values, dtype=float)
            n_all = np.asarray(n_values, dtype=float)
            ids_ac = np.flatnonzero(cs_all > 1.0 + tolerance)
            out[name] = {
                "max_cs2": float(np.max(cs_all)) if cs_all.size else None,
                "min_cs2": float(np.min(cs_all)) if cs_all.size else None,
                "acausal_count": int(ids_ac.size),
                "acausal_nb_range": ([float(n_all[ids_ac[0]]), float(n_all[ids_ac[-1]])]
                                      if ids_ac.size else None),
                "first_acausal_nb": float(n_all[ids_ac[0]]) if ids_ac.size else None,
                "n_valid": int(cs_all.size),
                "acausal": bool(ids_ac.size),
                "sampling": "original knots only, per smooth segment",
            }
            continue

        grid_cs_values = []
        grid_n_values = []
        seg_analyses = []
        for ids in segments:
            if len(ids) < 3:
                continue
            xs = x[ids]
            if name == "pchip":
                ps = PchipInterpolator(xs, lp[ids], extrapolate=False)
                es = PchipInterpolator(xs, le[ids], extrapolate=False)
            else:
                ps = CubicSpline(xs, lp[ids], bc_type="natural", extrapolate=False)
                es = CubicSpline(xs, le[ids], bc_type="natural", extrapolate=False)

            seg_analyses.append(
                _analyze_log_spline_intervals(
                    ps, es, xs,
                    n_min=n_min,
                    max_log_step=max_log_step,
                    tolerance=tolerance,
                )
            )
            chunks = [np.linspace(xs[j], xs[j + 1], samples_per_interval + 2)
                      for j in range(len(xs) - 1)]
            xq = np.unique(np.concatenate(chunks))
            cs = _cs2_from_log_splines(ps, es, xq)
            nq = np.exp(xq)
            valid = (nq >= n_min) & np.isfinite(cs)
            grid_cs_values.extend(cs[valid].tolist())
            grid_n_values.extend(nq[valid].tolist())

        valid_segs = [sa for sa in seg_analyses if sa["interval_count"] > 0 and sa["max_cs2"] is not None]
        if valid_segs:
            max_cs2 = float(max(sa["max_cs2"] for sa in valid_segs))
            min_cs2 = float(min(sa["min_cs2"] for sa in valid_segs))
            dp_pos = all(sa["dp_positive"] for sa in valid_segs)
            de_pos = all(sa["de_positive"] for sa in valid_segs)
            min_dp_seg = float(min(sa["min_dp_dlnn"] for sa in valid_segs if sa["min_dp_dlnn"] is not None))
            min_de_seg = float(min(sa["min_de_dlnn"] for sa in valid_segs if sa["min_de_dlnn"] is not None))
            all_crossings = [ev for sa in valid_segs for ev in sa["crossing_events"]]
            all_tangential = [ev for sa in valid_segs for ev in sa["tangential_contacts"]]
            acausal_segs = [sa for sa in valid_segs if sa["acausal"]]
            first_acausal_nb = acausal_segs[0]["first_acausal_nb"] if acausal_segs else None
            acausal_nb_range = (
                [acausal_segs[0]["acausal_nb_range"][0], acausal_segs[-1]["acausal_nb_range"][1]]
                if acausal_segs else None
            )
            acausal_count = sum(sa["acausal_interval_count"] for sa in valid_segs)
            interval_count = sum(sa["interval_count"] for sa in valid_segs)
            stationary_count = sum(sa["stationary_point_count"] for sa in valid_segs)
            nonpos_deriv_nb = [v for sa in valid_segs for v in sa["nonpositive_derivative_nb"]]
            if any(sa["crossing_kind"] == "true_crossing" for sa in valid_segs):
                crossing_kind = "true_crossing"
            elif any(sa["crossing_kind"] == "tangential_contact" for sa in valid_segs):
                crossing_kind = "tangential_contact"
            elif all(sa["crossing_kind"] == "always_superluminal" for sa in valid_segs):
                crossing_kind = "always_superluminal"
            else:
                crossing_kind = "none"
        else:
            max_cs2 = min_cs2 = min_dp_seg = min_de_seg = None
            dp_pos = de_pos = True
            all_crossings = all_tangential = nonpos_deriv_nb = []
            first_acausal_nb = acausal_nb_range = None
            acausal_count = interval_count = stationary_count = 0
            crossing_kind = "none"

        cs_grid = np.asarray(grid_cs_values, dtype=float)
        out[name] = {
            "max_cs2": max_cs2,
            "min_cs2": min_cs2,
            "dp_positive": bool(dp_pos),
            "de_positive": bool(de_pos),
            "min_dp_dlnn": min_dp_seg,
            "min_de_dlnn": min_de_seg,
            "nonpositive_derivative_nb": nonpos_deriv_nb,
            "acausal_count": int(acausal_count),
            "acausal_interval_count": int(acausal_count),
            "acausal_nb_range": acausal_nb_range,
            "first_acausal_nb": first_acausal_nb,
            "crossing_kind": crossing_kind,
            "crossing_events": all_crossings,
            "tangential_contacts": all_tangential,
            "stationary_point_count": int(stationary_count),
            "interval_count": int(interval_count),
            "n_valid": int(interval_count),
            "acausal": bool(first_acausal_nb is not None or (max_cs2 is not None and max_cs2 > 1.0 + tolerance)),
            "grid_diagnostic_max_cs2": float(np.max(cs_grid)) if cs_grid.size else None,
            "grid_diagnostic_min_cs2": float(np.min(cs_grid)) if cs_grid.size else None,
            "grid_diagnostic_n_points": int(cs_grid.size),
            "sampling": "exact degree<=6 interval extrema (uniform subgrid diagnostic-only)",
        }

    robust = all(out[k]["acausal"] and out[k]["n_valid"] > 0
                 for k in ("pchip", "cubic_spline", "centered_fd"))
    onset = [out[k]["first_acausal_nb"] for k in out
             if out[k]["first_acausal_nb"] is not None]
    out["classification"] = {
        "robust_acausality": robust,
        "onset_spread_fm3": float(max(onset) - min(onset)) if onset else 0.0,
        "root_cause": "source_EOS_acausality" if robust else
                      ("derivative_scheme_sensitive" if onset else "no_acausality_detected"),
        "phase_boundary_indices": phase_edges.astype(int).tolist(),
        "density_gap_indices": gap_edges.astype(int).tolist(),
        "max_log_step": float(max_log_step),
        "samples_per_interval": samples_per_interval,
        "acceptance_basis": "exact_interval_extrema",
    }
    return out


def causal_completion(
    data: EOSData, *, cs2=1.0, n_max=None, n_points=200,
    max_log_step=0.25, samples_per_interval=32
):
    """Create a separate thermodynamically matched constant-sound-speed tail.

    Uses exact interval extrema and the degree-<=6 stationary polynomial on the
    PCHIP log-spline interpolant to locate the first true upward c_s^2=1
    crossing and distinguish it from tangential contacts or non-crossing tables.
    Accepts one smooth phase only: refuses to interpolate across Iphase edges.
    """
    n, P, eps, mu = _validated_arrays(data)
    s = float(cs2)
    if not np.isfinite(s) or not 0 < s <= 1:
        raise ValueError("cs2 must be in (0,1]")
    if (not np.isfinite(max_log_step) or max_log_step <= 0
            or int(samples_per_interval) < 4 or int(n_points) < 4):
        raise ValueError("invalid causal-completion sampling configuration")
    if np.asarray(data.Iphase).shape == n.shape and np.any(np.diff(data.Iphase) != 0):
        raise ValueError("causal completion is not phase-aware; split the source at its Iphase phase boundary first")

    x = np.log(n)
    lp = np.log(P)
    le = np.log(eps)
    p_spline = PchipInterpolator(x, lp, extrapolate=False)
    e_spline = PchipInterpolator(x, le, extrapolate=False)

    analysis = _analyze_log_spline_intervals(
        p_spline, e_spline, x,
        n_min=0.0,
        max_log_step=float(max_log_step),
        tolerance=1e-12,
        contact_tol=1e-9,
    )
    if analysis["interval_count"] == 0:
        raise ValueError("no valid smooth intervals for causal-completion analysis")
    if not analysis["dp_positive"] or not analysis["de_positive"]:
        raise ValueError(
            "source log-spline has non-positive derivative (p'<=0 or e'<=0) inside an interval; "
            "cannot define sound-speed completion"
        )
    cs0 = float(_cs2_from_log_splines(p_spline, e_spline, np.asarray([x[0]]))[0])
    if cs0 >= 1.0:
        raise ValueError("source begins at or above the causal limit; no causal low-density branch exists")

    upward = [ev for ev in analysis["crossing_events"] if ev["direction"] == "upward"]
    if not upward:
        if analysis["tangential_contacts"]:
            raise ValueError(
                "tangential c_s^2=1 contact detected without crossing into superluminal regime; "
                "causal completion not applied"
            )
        raise ValueError("no resolved first c_s^2=1 crossing found; causal completion not applied")

    first = upward[0]
    i = int(first["interval_index"])
    root = float(first["x"])
    nt = float(np.exp(root))
    pt = float(np.exp(p_spline(root)))
    et = float(np.exp(e_spline(root)))
    mu_t = (et + pt) / nt
    mu_source = float(PchipInterpolator(x, mu, extrapolate=False)(root))
    mu_match_rel = abs(mu_t - mu_source) / max(abs(mu_t), 1.0)
    if mu_match_rel > 5e-3:
        raise ValueError(f"first-law chemical potential mismatch at causal completion point: {mu_match_rel:.3e}")

    nmax = float(n[-1] if n_max is None else n_max)
    if not np.isfinite(nmax) or nmax <= nt:
        raise ValueError("n_max must be finite and exceed causal-limit density")
    q = np.geomspace(1.0, nmax / nt, int(n_points))
    tail = nt * q
    w = et + pt
    p_tail = pt + s / (1.0 + s) * w * (q ** (1.0 + s) - 1.0)
    mu_tail = mu_t * q ** s
    e_tail = mu_tail * tail - p_tail
    keep = n < nt
    n_out = np.r_[n[keep], tail]
    p_out = np.r_[P[keep], p_tail]
    e_out = np.r_[eps[keep], e_tail]
    mu_out = (e_out + p_out) / n_out
    if (np.any(np.diff(n_out) <= 0) or np.any(np.diff(p_out) <= 0)
            or np.any(np.diff(e_out) <= 0) or np.any(~np.isfinite(n_out + p_out + e_out + mu_out))):
        raise ValueError("causal completion is not finite and strictly monotone")

    mn = float(data.mn)
    F = e_out / n_out
    eos = EOSData(
        nb=n_out, mn=mn, mp=float(data.mp), Q1=p_out / n_out,
        Q2=np.zeros_like(n_out), Q3=mu_out / mn - 1.0,
        Q4=np.zeros_like(n_out), Q5=np.zeros_like(n_out),
        Q6=F / mn - 1.0, Q7=F / mn - 1.0,
        q_extra=np.zeros_like(n_out),
        # The source/completion join can have a finite c_s^2 kink. Preserve
        # source labels and give the phenomenological tail a distinct phase ID
        # so derivative stencils never interpolate across that model boundary.
        Iphase=np.r_[np.asarray(data.Iphase, dtype=int)[keep],
                     np.full(len(tail), (int(np.max(data.Iphase)) + 1)
                             if np.asarray(data.Iphase).shape == n.shape else 1, dtype=int)],
        Yn=np.zeros_like(n_out), Yp=np.zeros_like(n_out), Ye=np.zeros_like(n_out),
        Ymu=np.zeros_like(n_out), Aav=np.zeros_like(n_out), Zav=np.zeros_like(n_out),
        Yh=np.zeros_like(n_out), Il=int(data.Il), iT=int(data.iT),
        iYq=int(data.iYq), inb_start=1, temperature_mev=0,
        charge_fraction_grid=0, q_extra_count=0,
        composition_available=False, beta_equilibrated=True,
    )
    meta = {
        "schema": "compose-eos-causal-completion-v2",
        "method": "constant_sound_speed_extension",
        "cs2_extension": s, "n_match_fm3": nt, "P_match_MeV_fm3": pt,
        "eps_match_MeV_fm3": et, "n_max_fm3": nmax,
        "source_composition_omitted": bool(data.composition_available),
        "source_max_cs2": float(analysis["max_cs2"]),
        "source_min_cs2": float(analysis["min_cs2"]),
        "source_crossing_kind": analysis["crossing_kind"],
        "source_first_limiter_index": i,
        "chemical_potential_match_relative_residual": float(mu_match_rel),
        "sound_speed_continuity_at_join": bool(abs(s - 1.0) < 1e-12),
        "source_sampling": {
            "method": "exact_interval_extrema",
            "max_log_step": float(max_log_step),
            "diagnostic_samples_per_interval": int(samples_per_interval),
        },
        "warning": "This is an explicit constant-cs alternative, not a unique repair of the microscopic source EOS. If cs2_extension differs from 1, c_s^2 has a finite kink at the first causal-limit point even though P, epsilon and mu_B remain continuous.",
    }
    return eos, meta
