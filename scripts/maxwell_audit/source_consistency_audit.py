"""Source-table convexity and thermodynamic-consistency pre-interpolation audit.

Implements §1 of the Maxwell Stability audit. Computes, directly from the raw
CompOSE source values (no interpolation), the per-interval diagnostics that
distinguish a SOURCE_TABLE_CONVEXITY_FAIL from a HERMITE_INTERPOLANT_SHAPE_FAIL
and the other §1 categories.

CompOSE / repository conventions (compose_eos_types.f90):
  P        = n_b * Q1
  epsilon  = n_b * (Q7 + 1) * m_n
  mu_Euler = (epsilon + P) / n_b
  mu_Q3    = (Q3 + 1) * m_n             (Il=1: mu_Q3 = mu_Euler when Q5~0)
           (Il=0: mu_Q3 = mu_Euler - Yq*mu_q, recovered before comparison)

For each interval i -> i+1 of a selected single-Iphase phase segment:
  s_i = (P_{i+1} - P_i) / (mu_{i+1} - mu_i)        (source secant slope)

Mechanical-stability necessary condition for a single phase:
  n_B = dP/dmu_B,  dn_B/dmu_B >= 0   ==>
  min(n_i, n_{i+1}) <= s_i <= max(n_i, n_{i+1})

This is a PRE-INTERPOLATION test on the raw source values. It does not depend
on Hermite, PCHIP, or any reconstruction. A failure here is a
SOURCE_TABLE_CONVEXITY_FAIL: the source table itself is not convex in (mu, P)
and no interpolant that reproduces the source endpoints can be mechanically
stable. A convex source that the cubic Hermite then distorts is a
HERMITE_INTERPOLANT_SHAPE_FAIL (computed separately).

Categories (per §1):
  SOURCE_TABLE_CONVEXITY_FAIL    secant necessary condition fails
  HERMITE_INTERPOLANT_SHAPE_FAIL source secant passes; cubic Hermite violates
                                 n>0 or dn/dmu_B>0 (exact polynomial extrema)
  SOURCE_THERMO_IDENTITY_FAIL    Q3/Euler or first-law source residual > tol
  SOURCE_CAUSALITY_FAIL          source or reconstructed phase violates c_s^2<=1
  SOURCE_PHASE_LAYOUT_UNSUPPORTED phase/provider layout cannot be read safely
"""
from __future__ import annotations
import json
import hashlib
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
import os
import sys
sys.path.insert(0, str(ROOT / "python_frontend"))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_ROOT", str(ROOT))
os.environ.setdefault("COMPOSE_EOS_FORTRAN_EXE",
                      str(ROOT / "build" / "compose_eos_fortran" / "compose_eos_fortran"))
from compose_eos import read_compose  # noqa: E402

# §1 tolerances. The Q3/Euler residual tolerance is the same as the v1 kernel
# (euler_tolerance = 5e-4, empirically justified in compose_eos_phase_transition.f90
# lines 9-33: worst observed 1.15e-4 across 18-model corpus; 5e-4 = 4.3x margin).
EULER_REL_TOL = 5.0e-4
FIRST_LAW_REL_TOL = 1.0e-10   # |eps+P - mu*n|/|mu*n|; source printed to ~10 sig figs
CAUSAL_TOL = 1.0e-8
SECANT_MARGIN = 0.0           # exact: s_i must lie in [min, max] of endpoint n


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def source_hashes(model_dir: Path) -> dict:
    return {name: sha256_file(model_dir / name)
            for name in ("eos.nb", "eos.thermo", "eos.compo", "eos.t", "eos.yq", "eos.mr", "eos.init")
            if (model_dir / name).exists()}


@dataclass
class IntervalDiag:
    model: str
    iphase: int
    interval_index: int            # 0-based within the segment
    source_row_lo: int             # absolute row index in the source table
    source_row_hi: int
    n0: float
    n1: float
    mu0: float
    mu1: float
    P0: float
    P1: float
    eps0: float
    eps1: float
    mu_q3_0: float                  # source Q3 chemical potential (recovered)
    mu_q3_1: float
    euler_residual_0: float         # |mu_Q3 - mu_Euler|/|mu_Euler| at lo
    euler_residual_1: float
    first_law_residual_0: float      # |eps+P - mu*n|/|mu*n| at lo
    first_law_residual_1: float
    dP: float
    dmu: float
    secant: float                   # dP/dmu
    n_min: float
    n_max: float
    secant_in_range: bool           # min(n0,n1) - margin <= s <= max(n0,n1) + margin
    # Hermite polynomial (cubic, P'=n_B endpoint constraints) exact extrema
    hermite_min_n: float            # exact min of n(mu) over the interval
    hermite_min_dndmu: float         # exact min of dn/dmu (linear in mu => endpoint min)
    hermite_min_causal: float       # exact min of mu*dn/dmu - n (quadratic)
    # diagnoses (independent)
    category: str                   # one of the §1 categories or PASS
    reason: str                     # human-readable detail


def hermite_n_poly(mu0, mu1, p0, p1, n0, n1):
    """Cubic Hermite of P(mu) with P'(mu0)=n0, P'(mu1)=n1.

    Returns the polynomial coefficients of n(mu)=P'(mu) as a quadratic in t
    and of dn/dmu as linear in t, plus the causal polynomial mu*dn/dmu - n(mu)
    evaluated analytically. t in [0,1] parametrizes [mu0,mu1], h=mu1-mu0.
    """
    h = mu1 - mu0
    # n(t) = ((6t^2-6t)/h) p0 + (3t^2-4t+1) n0 + ((-6t^2+6t)/h) p1 + (3t^2-2t) n1
    # collect as a*t^2 + b*t + c
    a_n = (6.0/h)*p0 + (-6.0/h)*p1 + 3.0*n0 + 3.0*n1
    b_n = (-6.0/h)*p0 + (6.0/h)*p1 - 4.0*n0 - 2.0*n1
    c_n = n0
    # dn/dmu = d(n)/dmu = (1/h) d/dt. n is quadratic in t => dn/dmu linear in t.
    # dn/dmu = (1/h)(2 a_n t + b_n)
    # causal poly: mu*dn/dmu - n. mu = mu0 + h t.
    # = (mu0 + h t)*(2 a_n t + b_n)/h - (a_n t^2 + b_n t + c_n)
    # collect in t:
    # t^2: (2 a_n) + ... let's expand:
    #   (mu0 + h t)(2 a_n t + b_n)/h = mu0(2 a_n t + b_n)/h + t(2 a_n t + b_n)
    #     = (2 a_n mu0/h) t + (b_n mu0/h) + 2 a_n t^2 + b_n t
    #   minus (a_n t^2 + b_n t + c_n)
    # t^2: 2 a_n - a_n = a_n
    # t^1: 2 a_n mu0/h + b_n mu0/h + b_n - b_n = (2 a_n + b_n) mu0/h
    #      wait: b_n mu0/h + b_n - b_n = b_n mu0/h; plus 2 a_n mu0/h => (2 a_n + b_n) mu0/h
    # t^0: b_n mu0/h - c_n
    a_c = a_n
    b_c = (2.0*a_n + b_n) * mu0 / h
    c_c = b_n * mu0 / h - c_n
    return (a_n, b_n, c_n), (a_c, b_c, c_c), h


def quad_min(a, b, c, t_lo=0.0, t_hi=1.0):
    """Exact minimum of a t^2 + b t + c on [t_lo, t_hi]."""
    # endpoints
    ymin = min(a*t_lo*t_lo + b*t_lo + c, a*t_hi*t_hi + b*t_hi + c)
    # vertex t* = -b/(2a), if a != 0 and inside
    if abs(a) > 1e-15 * max(1.0, abs(b), abs(c)):
        tv = -b / (2.0 * a)
        if t_lo < tv < t_hi:
            yv = a*tv*tv + b*tv + c
            if yv < ymin:
                ymin = yv
    return ymin


def linear_min(a, b, t_lo=0.0, t_hi=1.0):
    """Minimum of a t + b on [t_lo, t_hi] (linear => endpoint min)."""
    return min(a*t_lo + b, a*t_hi + b)


def select_last_iphase_segment(il: np.ndarray) -> tuple[int, int, int]:
    """Return (label, start_idx, end_idx_exclusive) of the last/highest-density
    contiguous Iphase segment (rows ordered by increasing nb)."""
    segs = []
    cur = int(il[0]); start = 0
    for i in range(1, len(il)):
        if int(il[i]) != cur:
            segs.append((cur, start, i)); cur = int(il[i]); start = i
    segs.append((cur, start, len(il)))
    return segs[-1]


def audit_segment(model: str, model_dir: Path, seg: Optional[tuple] = None) -> list[IntervalDiag]:
    """Audit every interval of one Iphase segment. If seg is None, the last
    highest-density contiguous segment is selected (the §1 convention)."""
    d = read_compose(str(model_dir))
    il = np.asarray(d.Iphase, int)
    nb = np.asarray(d.nb, float)
    Q = np.vstack([np.asarray(d.Q1, float), np.asarray(d.Q2, float),
                   np.asarray(d.Q3, float), np.asarray(d.Q4, float),
                   np.asarray(d.Q5, float), np.asarray(d.Q6, float),
                   np.asarray(d.Q7, float)])  # shape (7, n)
    q = Q  # alias: q[k] is Q_{k+1}
    mn = float(d.mn)
    il_flag = int(getattr(d, "Il", 1))
    yq = float(getattr(d, "charge_fraction_grid", 0.0))
    if seg is None:
        seg = select_last_iphase_segment(il)
    label, lo, hi = seg
    diags: list[IntervalDiag] = []
    for i in range(lo, hi - 1):
        n0, n1 = nb[i], nb[i+1]
        P0 = n0 * q[0, i]; P1 = n1 * q[0, i+1]
        eps0 = n0 * (q[6, i] + 1.0) * mn
        eps1 = n1 * (q[6, i+1] + 1.0) * mn
        mu0 = (eps0 + P0) / n0
        mu1 = (eps1 + P1) / n1
        # source Q3 chemical potential, recovered per Il convention
        mu_q3_0 = (q[2, i] + 1.0) * mn
        mu_q3_1 = (q[2, i+1] + 1.0) * mn
        if il_flag == 0:
            mu_q3_0 = mu_q3_0 + yq * q[3, i] * mn
            mu_q3_1 = mu_q3_1 + yq * q[3, i+1] * mn
        euler_res_0 = abs(mu_q3_0 - mu0) / max(abs(mu0), 1.0)
        euler_res_1 = abs(mu_q3_1 - mu1) / max(abs(mu1), 1.0)
        fl_res_0 = abs(eps0 + P0 - mu0 * n0) / max(abs(mu0 * n0), 1.0)
        fl_res_1 = abs(eps1 + P1 - mu1 * n1) / max(abs(mu1 * n1), 1.0)
        dP = P1 - P0; dmu = mu1 - mu0
        secant = dP / dmu if dmu != 0 else float("inf")
        n_min = min(n0, n1); n_max = max(n0, n1)
        secant_ok = (n_min - SECANT_MARGIN) <= secant <= (n_max + SECANT_MARGIN)
        # Hermite exact extrema
        (a_n, b_n, c_n), (a_c, b_c, c_c), h = hermite_n_poly(mu0, mu1, P0, P1, n0, n1)
        herm_min_n = quad_min(a_n, b_n, c_n)
        # dn/dmu = (1/h)(2 a_n t + b_n) -> linear in t, min at endpoints
        herm_min_dndmu = min((2*a_n*0 + b_n)/h, (2*a_n*1 + b_n)/h)
        herm_min_causal = quad_min(a_c, b_c, c_c)
        # classify
        category = "PASS"; reason = ""
        if not secant_ok:
            category = "SOURCE_TABLE_CONVEXITY_FAIL"
            where = "above" if secant > n_max else "below"
            reason = (f"source secant dP/dmu={secant:.9e} lies {where} endpoint "
                     f"density range [{n_min:.9e}, {n_max:.9e}]")
        elif herm_min_n <= 0 or herm_min_dndmu <= 0:
            category = "HERMITE_INTERPOLANT_SHAPE_FAIL"
            if herm_min_n <= 0:
                reason = f"cubic Hermite n(mu) min={herm_min_n:.6e} <= 0"
            else:
                reason = f"cubic Hermite dn/dmu min={herm_min_dndmu:.6e} <= 0"
        elif herm_min_causal < -CAUSAL_TOL * max(1.0, abs(n0), abs(n1), abs(mu0), abs(mu1)):
            category = "SOURCE_CAUSALITY_FAIL"
            reason = f"cubic Hermite causal poly min={herm_min_causal:.6e} < 0"
        elif max(euler_res_0, euler_res_1) > EULER_REL_TOL:
            category = "SOURCE_THERMO_IDENTITY_FAIL"
            reason = (f"Q3/Euler residual {max(euler_res_0, euler_res_1):.3e} > "
                     f"tolerance {EULER_REL_TOL:.0e}")
        elif max(fl_res_0, fl_res_1) > FIRST_LAW_REL_TOL:
            category = "SOURCE_THERMO_IDENTITY_FAIL"
            reason = (f"first-law residual {max(fl_res_0, fl_res_1):.3e} > "
                     f"tolerance {FIRST_LAW_REL_TOL:.0e}")
        diags.append(IntervalDiag(
            model=model, iphase=int(label), interval_index=i - lo,
            source_row_lo=int(i), source_row_hi=int(i+1),
            n0=n0, n1=n1, mu0=mu0, mu1=mu1, P0=P0, P1=P1, eps0=eps0, eps1=eps1,
            mu_q3_0=mu_q3_0, mu_q3_1=mu_q3_1,
            euler_residual_0=euler_res_0, euler_residual_1=euler_res_1,
            first_law_residual_0=fl_res_0, first_law_residual_1=fl_res_1,
            dP=dP, dmu=dmu, secant=secant, n_min=n_min, n_max=n_max,
            secant_in_range=secant_ok,
            hermite_min_n=herm_min_n, hermite_min_dndmu=herm_min_dndmu,
            hermite_min_causal=herm_min_causal,
            category=category, reason=reason))
    return diags


def _to_native(o):
    """Recursively convert numpy scalars/arrays to native python for JSON."""
    import numpy as np
    if isinstance(o, dict):
        return {k: _to_native(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_to_native(x) for x in o]
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return [_to_native(x) for x in o.tolist()]
    return o


def audit_model(model: str, model_dir: Path, seg: Optional[tuple] = None) -> dict:
    diags = audit_segment(model, model_dir, seg)
    from collections import Counter
    cats = Counter(d.category for d in diags)
    seg_used = seg if seg is not None else select_last_iphase_segment(
        np.asarray(read_compose(str(model_dir)).Iphase, int))
    return {
        "model": model,
        "source_dir": str(model_dir),
        "source_hashes": source_hashes(model_dir),
        "segment": list(seg_used),
        "n_intervals": len(diags),
        "category_counts": dict(cats),
        "intervals": [asdict(d) for d in diags],
    }


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="SPG_M1")
    ap.add_argument("--data-root", default="eos_generator/eos_data")
    ap.add_argument("--out", default=None)
    ap.add_argument("--seg", default=None, help="label,lo,hi (1-based inclusive) or 'last'")
    args = ap.parse_args()
    mdir = ROOT / args.data_root / args.model
    seg = None
    if args.seg and args.seg != "last":
        lab, lo, hi = [int(x) for x in args.seg.split(",")]
        seg = (lab, lo - 1, hi)
    rep = _to_native(audit_model(args.model, mdir, seg))
    if args.out:
        Path(args.out).write_text(json.dumps(rep, indent=2))
    print(f"{args.model}: {rep['n_intervals']} intervals, categories={rep['category_counts']}")
    # show first few failing intervals
    fails = [iv for iv in rep["intervals"] if iv["category"] != "PASS"]
    print(f"  failing intervals: {len(fails)}")
    for iv in fails[:5]:
        print(f"  int{iv['interval_index']} (rows {iv['source_row_lo']}-{iv['source_row_hi']}): "
              f"{iv['category']} -- {iv['reason']}")
        print(f"    n0={iv['n0']:.9e} n1={iv['n1']:.9e} secant={iv['secant']:.9e} "
              f"range=[{iv['n_min']:.3e},{iv['n_max']:.3e}]")
