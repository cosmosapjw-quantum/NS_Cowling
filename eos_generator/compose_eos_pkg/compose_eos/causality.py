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


def causality_profile(
    data: EOSData, *, n_min=0.0, max_log_step=0.25,
    tolerance=1e-6, samples_per_interval=16
):
    """Compare PCHIP, natural-cubic and centered-FD c_s^2 diagnostics.

    PCHIP/cubic estimates are evaluated at knots and at interior points of
    every admissible interval. Each smooth phase and density-gap segment is
    interpolated independently. The centered finite-difference diagnostic is
    evaluated only at original tabulated knots within each segment.
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
        cs_values = []
        n_values = []
        for ids in segments:
            if len(ids) < 3:
                continue
            xs = x[ids]
            if name == "centered_fd":
                dp = np.gradient(lp[ids], xs, edge_order=2)
                de = np.gradient(le[ids], xs, edge_order=2)
                with np.errstate(divide="ignore", invalid="ignore"):
                    cs = (p[ids] / eps[ids]) * dp / de
                valid = (n[ids] >= n_min) & np.isfinite(cs) & (de > 0)
                cs_values.extend(cs[valid].tolist())
                n_values.extend(n[ids][valid].tolist())
                continue

            if name == "pchip":
                ps = PchipInterpolator(xs, lp[ids], extrapolate=False)
                es = PchipInterpolator(xs, le[ids], extrapolate=False)
            else:
                ps = CubicSpline(xs, lp[ids], bc_type="natural", extrapolate=False)
                es = CubicSpline(xs, le[ids], bc_type="natural", extrapolate=False)

            chunks = [np.linspace(xs[j], xs[j + 1], samples_per_interval + 2)
                      for j in range(len(xs) - 1)]
            xq = np.unique(np.concatenate(chunks))
            cs = _cs2_from_log_splines(ps, es, xq)
            nq = np.exp(xq)
            valid = (nq >= n_min) & np.isfinite(cs)
            cs_values.extend(cs[valid].tolist())
            n_values.extend(nq[valid].tolist())

        cs_all = np.asarray(cs_values, dtype=float)
        n_all = np.asarray(n_values, dtype=float)
        ids = np.flatnonzero(cs_all > 1.0 + tolerance)
        out[name] = {
            "max_cs2": float(np.max(cs_all)) if cs_all.size else None,
            "min_cs2": float(np.min(cs_all)) if cs_all.size else None,
            "acausal_count": int(ids.size),
            "acausal_nb_range": ([float(n_all[ids[0]]), float(n_all[ids[-1]])]
                                  if ids.size else None),
            "first_acausal_nb": float(n_all[ids[0]]) if ids.size else None,
            "n_valid": int(cs_all.size),
            "acausal": bool(ids.size),
            "sampling": "between-knot per smooth segment" if name != "centered_fd"
                        else "original knots only, per smooth segment",
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
    }
    return out


def causal_completion(
    data: EOSData, *, cs2=1.0, n_max=None, n_points=200,
    max_log_step=0.25, samples_per_interval=32
):
    """Create a separate thermodynamically matched constant-sound-speed tail.

    The source EOS is sampled between knots to locate the first c_s^2=1
    crossing. This implementation accepts one smooth phase only: it refuses
    to interpolate across an Iphase discontinuity. The original EOS is never
    modified, and a completion remains a model-dependent alternative.
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
    dln = np.diff(x)
    samples_per_interval = int(samples_per_interval)
    candidate_x = []
    candidate_cs = []
    candidate_i = []
    for i in range(len(n) - 1):
        if dln[i] > max_log_step:
            continue
        xq = np.linspace(x[i], x[i + 1], samples_per_interval + 2)
        csq = _cs2_from_log_splines(p_spline, e_spline, xq)
        candidate_x.extend(xq.tolist())
        candidate_cs.extend(csq.tolist())
        candidate_i.extend([i] * len(xq))

    xq = np.asarray(candidate_x, dtype=float)
    csq = np.asarray(candidate_cs, dtype=float)
    iq = np.asarray(candidate_i, dtype=int)
    valid = np.isfinite(csq)
    xq, csq, iq = xq[valid], csq[valid], iq[valid]
    if not xq.size:
        raise ValueError("no valid smooth intervals for causal-completion analysis")
    order = np.argsort(xq)
    xq, csq, iq = xq[order], csq[order], iq[order]
    if csq[0] >= 1.0:
        raise ValueError("source begins at or above the causal limit; no causal low-density branch exists")

    crossing = np.flatnonzero((csq[:-1] <= 1.0) & (csq[1:] > 1.0))
    if not crossing.size:
        raise ValueError("no resolved first c_s^2=1 crossing found; causal completion not applied")
    k = int(crossing[0])
    i = int(iq[k])
    def cs2_minus_one(z):
        val = float(_cs2_from_log_splines(p_spline, e_spline, np.asarray([z]))[0])
        if not np.isfinite(val):
            raise FloatingPointError("non-finite sound-speed interpolation at causal crossing")
        return val - 1.0
    root = brentq(cs2_minus_one, float(xq[k]), float(xq[k + 1]), xtol=1e-13, rtol=1e-13)
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
        q_extra=np.zeros_like(n_out), Iphase=np.zeros(len(n_out), dtype=int),
        Yn=np.zeros_like(n_out), Yp=np.zeros_like(n_out), Ye=np.zeros_like(n_out),
        Ymu=np.zeros_like(n_out), Aav=np.zeros_like(n_out), Zav=np.zeros_like(n_out),
        Yh=np.zeros_like(n_out), Il=int(data.Il), iT=int(data.iT),
        iYq=int(data.iYq), inb_start=1, temperature_mev=0,
        charge_fraction_grid=0, q_extra_count=0,
        composition_available=False, beta_equilibrated=True,
    )
    meta = {
        "schema": "compose-eos-causal-completion-v1",
        "method": "constant_sound_speed_extension",
        "cs2_extension": s, "n_match_fm3": nt, "P_match_MeV_fm3": pt,
        "eps_match_MeV_fm3": et, "n_max_fm3": nmax,
        "source_composition_omitted": bool(data.composition_available),
        "source_max_cs2": float(np.max(csq)),
        "source_first_limiter_index": i,
        "chemical_potential_match_relative_residual": float(mu_match_rel),
        "sound_speed_continuity_at_join": bool(abs(s - 1.0) < 1e-12),
        "source_sampling": {"max_log_step": float(max_log_step),
                            "samples_per_interval": samples_per_interval},
        "warning": "This is an explicit constant-cs alternative, not a unique repair of the microscopic source EOS. If cs2_extension differs from 1, c_s^2 has a finite kink at the first causal-limit point even though P, epsilon and mu_B remain continuous.",
    }
    return eos, meta
