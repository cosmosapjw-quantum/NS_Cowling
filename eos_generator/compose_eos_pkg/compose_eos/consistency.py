"""Thermodynamic and composition consistency checks for cold barotropic EOSs.

The algebraic identities are checked pointwise. Derivative checks are explicitly
reported as *finite-grid diagnostics*: points adjacent to large density gaps or
phase-index changes are excluded from smooth-branch derivatives. This avoids
mistaking a finite difference across a discontinuity / 100x interpolation gap
for a local derivative.

Convention: EOSData.Yh is an internal nuclear-group baryon fraction, i.e.
Yh=Abar*Y_I with CompOSE's standard nuclear-group number abundance Y_I.
"""
from __future__ import annotations

from typing import Dict, List
import numpy as np

from .eos_data import EOSData


def _relerr(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Symmetric relative residual with a scale floor to avoid 0/0."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    scale = np.maximum(np.maximum(np.abs(a), np.abs(b)), 1e-30)
    return (a - b) / scale


def _maxabs(values: np.ndarray, mask: np.ndarray = None) -> float:
    v = np.asarray(values, dtype=float)
    if mask is not None:
        v = v[np.asarray(mask, dtype=bool)]
    v = v[np.isfinite(v)]
    return float(np.max(np.abs(v))) if v.size else float("nan")


def _finite_difference_mask(data: EOSData, max_log_step: float) -> tuple:
    """Mask stencils that cross a large density gap or a phase boundary."""
    n = data.n
    mask = np.ones(n, dtype=bool)
    dln = np.diff(np.log(data.nb))
    gap_edges = np.where(dln > max_log_step)[0]
    for i in gap_edges:
        mask[i:i + 2] = False
    phase_edges = np.array([], dtype=int)
    if data.composition_available and data.Iphase.size == n:
        phase_edges = np.where(np.diff(data.Iphase) != 0)[0]
        for i in phase_edges:
            mask[i:i + 2] = False
    return mask, dln, gap_edges, phase_edges


def check_consistency(data: EOSData, verbose: bool = True,
                      *, algebraic_tolerance: float = 2e-6,
                      max_log_step_for_derivatives: float = 0.25,
                      cs2_tolerance: float = 1e-6,
                      composition_tolerance: float = 2e-5) -> Dict[str, object]:
    """Run consistency checks on a cold, one-dimensional CompOSE slice.

    Returns the historical C1--C11 keys plus an explicit `overall` verdict,
    applicability statuses, and derivative stencil masks. C3/C4 are NOT
    AVAILABLE when eos.compo was absent, which is a valid CompOSE case.
    """
    nb = np.asarray(data.nb, dtype=float)
    n = len(nb)
    if n < 3 or np.any(~np.isfinite(nb)) or np.any(nb <= 0) or np.any(np.diff(nb) <= 0):
        raise ValueError("EOS density grid must be finite, positive and strictly increasing")
    p, f, eps = data.p, data.f, data.eps
    mub, mul = data.mub, data.mul
    Q2, Q5, Q6, Q7 = data.Q2, data.Q5, data.Q6, data.Q7
    Fb = f / nb

    # Ensure all primary table arrays agree in length and are finite.
    primaries = [p, f, eps, mub, mul, Q2, Q5, Q6, Q7]
    if any(np.asarray(x).shape != (n,) for x in primaries):
        raise ValueError("EOS primary arrays have inconsistent shapes")
    finite_primary = all(np.all(np.isfinite(x)) for x in primaries)

    # C1: coldness and zero-temperature free/internal energy equality.
    c1 = dict(
        temperature_mev=float(data.temperature_mev),
        maxQ2=_maxabs(Q2),
        maxQ6mQ7=_maxabs(Q6 - Q7),
        maxmul=_maxabs(mul),
        pass_=bool(
            abs(float(data.temperature_mev)) <= 1e-12
            and _maxabs(Q2) <= algebraic_tolerance
            and _maxabs(Q6 - Q7) <= algebraic_tolerance
            and (not data.beta_equilibrated or _maxabs(mul) <= algebraic_tolerance * data.mn)
        ),
    )

    # C2: Euler/Gibbs-Duhem relation with the correct CompOSE thermodynamic
    # potential terms. For no-lepton tables, f=-p+nb*(mu_b+Yq*mu_q).
    # For charge-neutral lepton tables, f=-p+nb*(mu_b+Y_l*mu_l); in beta
    # equilibrium mu_l=0. Y_l is read from eos.compo when needed.
    if int(data.Il) == 0:
        mu_expected = (p + f) / nb - float(data.charge_fraction_grid) * data.muq
        c2_applicable = True
        c2_basis = "no-lepton: mu_b + Yq mu_q"
    else:
        if np.max(np.abs(mul)) <= 1e-8 * data.mn:
            mu_expected = (p + f) / nb
            c2_applicable = True
            c2_basis = "charge-neutral beta equilibrium: mu_l=0"
        elif data.composition_available:
            yl = data.Ye + data.Ymu
            mu_expected = (p + f) / nb - yl * mul
            c2_applicable = True
            c2_basis = "charge-neutral leptons: mu_b + Y_l mu_l"
        else:
            mu_expected = np.full(n, np.nan)
            c2_applicable = False
            c2_basis = "unresolved: nonzero mu_l but no lepton fractions"
    r2 = _relerr(mub, mu_expected)
    c2 = dict(
        basis=c2_basis,
        applicable=bool(c2_applicable),
        max_rel=_maxabs(r2) if c2_applicable else float("nan"),
        pass_=bool(c2_applicable and _maxabs(r2) <= algebraic_tolerance),
    )

    # C3/C4: standard CompOSE number-abundance closure expressed through the
    # internal baryonic fraction Yh=Abar*YI.
    if data.composition_available:
        if any(np.asarray(getattr(data, name)).shape != (n,) for name in
               ("Iphase", "Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh")):
            raise ValueError("composition arrays have inconsistent lengths")
        crust = (data.Aav > 0)
        core = ~crust
        ZA = np.divide(data.Zav, data.Aav, out=np.zeros(n), where=data.Aav > 0)
        r3 = data.Yp + ZA * data.Yh - data.Ye - data.Ymu
        r4 = data.Yn + data.Yp + data.Yh - 1.0
        c3 = dict(
            status="CHECKED",
            crust=_maxabs(r3, crust), core=_maxabs(r3, core),
            max_abs=_maxabs(r3), pass_=bool(_maxabs(r3) <= composition_tolerance),
            formula="Yp + (Zbar/Abar) Yh = Ye + Ymu; Yh=Abar*YI",
        )
        c4 = dict(
            status="CHECKED",
            crust=_maxabs(r4, crust), core=_maxabs(r4, core),
            max_abs=_maxabs(r4), pass_=bool(_maxabs(r4) <= composition_tolerance),
            formula="Yn + Yp + Yh = 1; Yh=Abar*YI",
        )
    else:
        crust = core = np.zeros(n, dtype=bool)
        r3 = r4 = np.full(n, np.nan)
        c3 = dict(status="NOT_AVAILABLE", crust=float("nan"), core=float("nan"),
                  max_abs=float("nan"), pass_=None,
                  reason="eos.compo is optional in CompOSE and absent from this table")
        c4 = dict(status="NOT_AVAILABLE", crust=float("nan"), core=float("nan"),
                  max_abs=float("nan"), pass_=None,
                  reason="eos.compo is optional in CompOSE and absent from this table")

    # C5: beta equilibrium, when it is part of the table contract.
    c5 = dict(maxmul=_maxabs(mul), applicable=bool(data.beta_equilibrated),
              pass_=(bool(_maxabs(mul) <= algebraic_tolerance * data.mn)
                     if data.beta_equilibrated else None))

    # C6: global discrete pressure monotonicity (not a local derivative).
    dp = np.diff(p)
    pscale = max(float(np.max(np.abs(p))), 1.0)
    dp_tol = 1e-12 * pscale
    min_secant = float(np.min(dp / np.diff(nb)))
    c6 = dict(min_dPdnb=min_secant, min_delta_p=float(np.min(dp)),
              pass_=bool(np.all(dp >= -dp_tol)))

    # Smooth-branch derivatives. Do not evaluate local derivative identities
    # across a large logarithmic density gap or across a change of Iphase.
    smooth_mask, dln, gap_edges, phase_edges = _finite_difference_mask(
        data, max_log_step_for_derivatives
    )
    dPdnb = np.gradient(p, nb, edge_order=2)
    dmubdnb = np.gradient(mub, nb, edge_order=2)
    dFbdnb = np.gradient(Fb, nb, edge_order=2)
    depsdnb = np.gradient(eps, nb, edge_order=2)
    derivative_mask = smooth_mask & np.isfinite(dPdnb) & np.isfinite(dmubdnb) & np.isfinite(dFbdnb) & np.isfinite(depsdnb)
    derivative_mask &= depsdnb > 0

    r7 = _relerr(dPdnb, nb * dmubdnb)
    r8 = _relerr(p, nb ** 2 * dFbdnb)
    r9 = _relerr(depsdnb, mub)
    c7 = dict(max_rel=_maxabs(r7, derivative_mask), n_points=int(derivative_mask.sum()),
              status="FINITE_DIFFERENCE_DIAGNOSTIC")
    c8 = dict(max_rel=_maxabs(r8, derivative_mask), n_points=int(derivative_mask.sum()),
              status="FINITE_DIFFERENCE_DIAGNOSTIC")
    c9 = dict(max_rel=_maxabs(r9, derivative_mask), n_points=int(derivative_mask.sum()),
              status="FINITE_DIFFERENCE_DIAGNOSTIC")

    # Primary discrete definition c_s^2 = (dP/dn)/(d epsilon/dn).  A second,
    # independent thermodynamic form uses d epsilon/dn = mu_b at T=0:
    # c_s^2 = (dP/dn)/mu_b.  Agreement is a useful cross-check on C9.
    cs2 = np.divide(dPdnb, depsdnb, out=np.full(n, np.nan), where=depsdnb > 0)
    cs2_first_law = np.divide(dPdnb, mub, out=np.full(n, np.nan), where=mub > 0)
    cs2_valid = cs2[derivative_mask]
    cs2_fl_valid = cs2_first_law[derivative_mask]
    cs2_min = float(np.min(cs2_valid)) if cs2_valid.size else float("nan")
    cs2_max = float(np.max(cs2_valid)) if cs2_valid.size else float("nan")
    cs2_fl_min = float(np.min(cs2_fl_valid)) if cs2_fl_valid.size else float("nan")
    cs2_fl_max = float(np.max(cs2_fl_valid)) if cs2_fl_valid.size else float("nan")
    cs2_cross = np.where(derivative_mask & (cs2 > 1.0 + cs2_tolerance))[0]
    cs2_fl_cross = np.where(derivative_mask & (cs2_first_law > 1.0 + cs2_tolerance))[0]
    both_valid = bool(cs2_valid.size and cs2_fl_valid.size)
    c10 = dict(
        min_cs2=cs2_min, max_cs2=cs2_max,
        min_cs2_first_law=cs2_fl_min, max_cs2_first_law=cs2_fl_max,
        max_abs_difference_between_forms=(float(np.max(np.abs(cs2_valid-cs2_fl_valid)))
                                           if both_valid else float("nan")),
        acausal_point_count=int(cs2_cross.size),
        acausal_nb_range=([float(nb[cs2_cross[0]]), float(nb[cs2_cross[-1]])]
                          if cs2_cross.size else None),
        acausal_point_count_first_law=int(cs2_fl_cross.size),
        n_points=int(cs2_valid.size),
        pass_=bool(
            both_valid
            and cs2_min >= -cs2_tolerance and cs2_max < 1.0 + cs2_tolerance
            and cs2_fl_min >= -cs2_tolerance and cs2_fl_max < 1.0 + cs2_tolerance
        ),
        status="SMOOTH_BRANCH_FINITE_DIFFERENCE; checked dP/deps and dP/(mu_b dn)",
    )

    dmu_valid = dmubdnb[derivative_mask]
    min_dmu = float(np.min(dmu_valid)) if dmu_valid.size else float("nan")
    c11 = dict(min_dmubdnb=min_dmu, n_points=int(dmu_valid.size),
               pass_=bool(dmu_valid.size and min_dmu >= -1e-10),
               status="SMOOTH_BRANCH_FINITE_DIFFERENCE")

    grid = dict(
        max_logstep=float(np.max(dln)), median_logstep=float(np.median(dln)),
        max_log10step=float(np.max(dln) / np.log(10.0)),
        median_log10step=float(np.median(dln) / np.log(10.0)),
        gap_indices=gap_edges.astype(int).tolist(),
        phase_boundary_indices=phase_edges.astype(int).tolist(),
        derivative_points=int(derivative_mask.sum()),
        derivative_masked_points=int(n - derivative_mask.sum()),
        derivative_gap_threshold=float(max_log_step_for_derivatives),
    )

    hard_checks = {
        "C1": c1["pass_"], "C2": c2["pass_"], "C6": c6["pass_"],
        "C10": c10["pass_"],
    }
    if data.composition_available:
        hard_checks.update({"C3": c3["pass_"], "C4": c4["pass_" ]})
    if data.beta_equilibrated:
        hard_checks["C5"] = c5["pass_"]
    # Thermodynamic stability is independently expressed by dmu_b/dnb >= 0.
    # A finite-difference sample set is required; no-sample is not silently PASS.
    hard_checks["C11"] = c11["pass_"]
    failures = [k for k, v in hard_checks.items() if not v]
    overall = dict(status="PASS" if not failures and finite_primary else "FAIL",
                   pass_=bool(not failures and finite_primary),
                   failed_checks=failures,
                   finite_primary_arrays=bool(finite_primary),
                   note="C7-C9/C11 are reported as finite-difference diagnostics, not standalone proofs")

    report = dict(
        C1=c1, C2=c2, C3=c3, C4=c4, C5=c5, C6=c6, C7=c7, C8=c8, C9=c9,
        C10=c10, C11=c11, grid=grid, overall=overall,
        cs2=cs2, cs2_first_law=cs2_first_law, derivative_mask=derivative_mask,
        dPdnb=dPdnb, dmubdnb=dmubdnb, depsdnb=depsdnb,
        gd_residual=r2, charge_residual=r3, baryon_residual=r4,
    )
    if verbose:
        _print_report(report)
    return report


def _print_report(r: dict) -> None:
    c1, c2, c3, c4, c5 = (r[k] for k in ("C1", "C2", "C3", "C4", "C5"))
    print("=" * 76)
    print(" COLD BAROTROPIC EOS CONSISTENCY AUDIT")
    print("=" * 76)
    print(f"[C1] T=0, Q2=0, Q6=Q7" +
          (f", beta equilibrium |mu_l|={c1['maxmul']:.3e} MeV" if c5['applicable'] else ""))
    print(f"     PASS={c1['pass_']}  max|Q2|={c1['maxQ2']:.3e}  max|Q6-Q7|={c1['maxQ6mQ7']:.3e}")
    print(f"[C2] Euler/Gibbs-Duhem ({c2['basis']})")
    print(f"     status={'CHECKED' if c2['applicable'] else 'NOT_AVAILABLE'}  max rel={c2['max_rel']:.3e}  PASS={c2['pass_']}")
    for key, c in (("C3", c3), ("C4", c4)):
        if c["status"] == "CHECKED":
            print(f"[{key}] {c['formula']} max abs={c['max_abs']:.3e} PASS={c['pass_']}")
        else:
            print(f"[{key}] NOT_AVAILABLE: {c['reason']}")
    if c5["applicable"]:
        print(f"[C5] beta equilibrium max|mu_l|={c5['maxmul']:.3e} MeV PASS={c5['pass_']}")
    print(f"[C6] monotone pressure min secant dP/dn={r['C6']['min_dPdnb']:.3e} PASS={r['C6']['pass_']}")
    for k in ("C7", "C8", "C9"):
        print(f"[{k}] finite-difference max rel={r[k]['max_rel']:.3e} ({r[k]['n_points']} smooth points)")
    print(f"[C10] smooth-branch c_s^2 range=[{r['C10']['min_cs2']:.5g},{r['C10']['max_cs2']:.5g}] PASS={r['C10']['pass_']}")
    print(f"[C11] smooth-branch min dmu_b/dn={r['C11']['min_dmubdnb']:.3e} PASS={r['C11']['pass_']}")
    print(f"[grid] excluded large-gap edges={r['grid']['gap_indices']}; phase edges={r['grid']['phase_boundary_indices']}; derivative points={r['grid']['derivative_points']}")
    print(f"[OVERALL] {r['overall']['status']} failed={r['overall']['failed_checks']}")


__all__ = ["check_consistency"]
