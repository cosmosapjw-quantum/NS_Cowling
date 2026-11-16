"""Thermodynamically explicit Maxwell/Gibbs constructions and phase-aware TOV.

The legacy matcher remains a smooth core--crust interpolation. This module adds
separate, physically distinct sharp Maxwell and global-neutrality Gibbs paths.
Maxwell data carry an explicit pressure-continuous, density-discontinuous event;
Gibbs data require two-charge pressure surfaces P(mu_B,mu_Q), never a 1-D beta-
equilibrium barotrope. The l=2 tidal jump uses the pressure-corrected
Takatsy--Kovacs (2020) condition.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional
import hashlib, json, os, shutil, tempfile
import numpy as np
from scipy.interpolate import PchipInterpolator, RectBivariateSpline
from scipy.optimize import brentq, root_scalar
from scipy.integrate import solve_ivp
from .eos_data import EOSData
from .io import read_compose, write_compose, verify_compose_readback, _compose_dir_hashes, _sha256_file
from .matching import MAXWELL_METHOD_ID

@dataclass
class MaxwellEOS:
    data: EOSData
    transition: dict
    metadata: dict = field(default_factory=dict)

@dataclass
class GibbsResult:
    eos: EOSData
    metadata: dict

@dataclass
class CausalCompletion:
    eos: EOSData
    metadata: dict


def _thermodynamic_pressure_from_density(mu, n, p):
    """Return the pressure primitive P(mu)=P(mu0)+integral n_B dmu_B.

    Interpolating n_B as the derivative of the thermodynamic potential and
    integrating it enforces dP/dmu_B=n_B identically. The deviation from the
    supplied tabulated pressure is retained as an independent source-table
    consistency diagnostic.
    """
    ns = PchipInterpolator(mu, n, extrapolate=False)
    primitive = ns.antiderivative()
    offset = float(p[0]) - float(primitive(mu[0]))
    return ns, primitive, offset


def _polynomial_roots_on_unit_interval(poly):
    """Isolate real polynomial roots by recursively partitioning at extrema.

    Used only for degree <= 3 density differences.  The normalized coordinate
    avoids scaling the root test by a chemical potential of order 1 GeV.
    This is floating-point isolation of the declared interpolant, not a claim
    to recover unrepresented structure in the microscopic EOS.
    """
    poly = poly.trim(tol=0.0)
    if poly.degree() <= 0:
        return []
    if poly.degree() == 1:
        root = -float(poly.coef[0]) / float(poly.coef[1])
        return [root] if 0.0 < root < 1.0 else []
    stationary = _polynomial_roots_on_unit_interval(poly.deriv())
    partition = [0.0, *stationary, 1.0]
    tol = 64.0*np.finfo(float).eps*float(np.sum(np.abs(poly.coef)))
    roots = [t for t in stationary if abs(float(poly(t))) <= tol]
    for a, b in zip(partition[:-1], partition[1:]):
        fa, fb = float(poly(a)), float(poly(b))
        if abs(fa) > tol and abs(fb) > tol and np.signbit(fa) != np.signbit(fb):
            roots.append(float(brentq(poly, a, b, xtol=4*np.finfo(float).eps,
                                     rtol=4*np.finfo(float).eps)))
    return sorted(roots)


def _maxwell_pressure_equality_roots(n_low, n_high, pressure_low, pressure_high, lo, hi):
    """Return all isolated (mu, tangential) pressure contacts on the common domain.

    P is the quartic primitive of each cubic PCHIP density.  On the union of
    their knot intervals, Delta n = d(Delta P)/dmu is cubic.  Its roots split
    Delta P into monotone pieces, so close root pairs cannot hide between
    samples. A non-isolated pressure-equality interval is explicitly rejected.
    """
    from numpy.polynomial import Polynomial
    knots = np.unique(np.r_[lo, n_low.x, n_high.x, hi])
    knots = knots[(knots >= lo) & (knots <= hi)]
    contacts = []
    epsmach = np.finfo(float).eps
    for left, right in zip(knots[:-1], knots[1:]):
        width = float(right-left)
        def density_poly(spline):
            index = min(max(int(np.searchsorted(spline.x, left, side='right'))-1, 0),
                        spline.c.shape[1]-1)
            return Polynomial(spline.c[:, index][::-1])(
                Polynomial([float(left-spline.x[index]), width]))
        delta_n = density_poly(n_low)-density_poly(n_high)
        delta_p = (width*delta_n).integ()
        pl, ph = float(pressure_low(left)), float(pressure_high(left))
        delta_p.coef[0] = pl-ph
        pressure_scale = max(abs(pl), abs(ph), abs(float(pressure_low(right))),
                             abs(float(pressure_high(right))), np.finfo(float).tiny)
        ptol = 64.0*epsmach*pressure_scale
        if np.all(delta_n.coef == 0.0):
            if abs(float(delta_p.coef[0])) <= ptol:
                raise ValueError('Maxwell pressure equality contains a non-isolated root interval')
            continue
        stationary = _polynomial_roots_on_unit_interval(delta_n)
        partition = [0.0, *stationary, 1.0]
        values = [float(delta_p(t)) for t in partition]
        for i, (t, value) in enumerate(zip(partition, values)):
            if abs(value) <= ptol:
                density_scale = max(abs(float(n_low(left+width*t))),
                                    abs(float(n_high(left+width*t))))
                tangential = (0 < i < len(partition)-1 or
                              abs(float(delta_n(t))) <= 64.0*epsmach*density_scale)
                contacts.append((float(left+width*t), bool(tangential)))
        for a, b, fa, fb in zip(partition[:-1], partition[1:], values[:-1], values[1:]):
            if abs(fa) > ptol and abs(fb) > ptol and np.signbit(fa) != np.signbit(fb):
                root = float(brentq(delta_p, a, b, xtol=4*epsmach, rtol=4*epsmach))
                contacts.append((float(left+width*root), False))
    contacts.sort()
    unique = []
    mutol = 8.0*epsmach*max(abs(lo), abs(hi), 1.0)
    for mu, tangential in contacts:
        if unique and abs(mu-unique[-1][0]) <= mutol:
            unique[-1] = (unique[-1][0], unique[-1][1] or tangential)
        else:
            unique.append((mu, tangential))
    return unique


MAXWELL_CS2_TOLERANCE = 1e-6   # fixed causal gate c_s^2 <= 1 + 1e-6 (unchanged)


def _real_roots_in_unit_interval(poly):
    """Real roots of a numpy Polynomial strictly inside (0, 1)."""
    poly = poly.trim(tol=0.0)
    if poly.degree() <= 0:
        return []
    out = []
    for r in poly.roots():
        if abs(r.imag) <= 1e-9 * max(1.0, abs(r.real)) and 0.0 < r.real < 1.0:
            out.append(float(r.real))
    return out


def _branchwise_first_law_report(mu, p, n, eps, *, samples_per_interval=8,
                                 cs2_tolerance=MAXWELL_CS2_TOLERANCE):
    """Exact interval-extrema audit of the declared branch interpolant.

    Representation: n_B(mu_B) is one monotone PCHIP (cubic per interval) and
    P(mu_B) = P(mu_0) + integral n_B dmu_B.  Hence n_B = dP/dmu_B holds *by
    construction* on the whole branch; it is not reported as a residual.  The
    independent source-table checks are the pressure and energy reconstruction
    residuals at the knots.

    Stability and causality are decided on the continuous interpolant, per
    interval s in [0, 1] (mu = mu_j + h s), without sampling:
      * dn/dmu is quadratic: extrema at s=0, s=1 and its vertex;
      * c_s^2 <= 1 + tol  <=>  g(s) = (1+tol) mu dn/dmu - n >= 0 when
        dn/dmu > 0; g is cubic, so its minimum is at the endpoints or at the
        real roots of the quadratic g'(s) (closed form);
      * the reported c_s^2 extrema use the roots of the at-most-quartic stationarity
        polynomial of n/(mu dn/dmu) plus every candidate above.
    Root solving and evaluation use float64; exact refers to the algebraic
    candidate set, not a directed-rounding interval proof.
    ``samples_per_interval`` is accepted for configuration compatibility only
    and is not used for any gate.
    """
    from numpy.polynomial import Polynomial
    mu = np.asarray(mu, dtype=float); p = np.asarray(p, dtype=float)
    n = np.asarray(n, dtype=float); eps = np.asarray(eps, dtype=float)
    if (mu.ndim != 1 or mu.size < 2 or n.shape != mu.shape or p.shape != mu.shape
            or eps.shape != mu.shape):
        raise ValueError("branch audit requires matching 1-D arrays with >= 2 rows")
    if not np.all(np.isfinite(mu)) or np.any(mu <= 0.0) or np.any(np.diff(mu) <= 0.0):
        raise ValueError("branch audit requires strictly increasing finite positive mu_B; "
                         "a phase jump must be split into one-sided branches")
    if not (np.isfinite(cs2_tolerance) and cs2_tolerance >= 0.0):
        raise ValueError("invalid causal tolerance")
    ns, primitive, offset = _thermodynamic_pressure_from_density(mu, n, p)
    p_reconstructed = np.asarray(primitive(mu), dtype=float) + offset
    p_scale = np.maximum(np.maximum(np.abs(p), abs(float(p[0]))), 1e-8)
    p_rel = np.abs(p_reconstructed - p) / p_scale
    eps_reconstructed = mu * n - p_reconstructed
    eps_scale = np.maximum(np.abs(eps), 1.0)
    eps_rel = np.abs(eps_reconstructed - eps) / eps_scale

    epsmach = np.finfo(float).eps
    min_dn = np.inf; arg_dn = float(mu[0])
    min_n = np.inf
    max_cs2 = -np.inf; arg_max = float(mu[0])
    min_cs2 = np.inf
    min_margin_rel = np.inf; arg_margin = float(mu[0])
    causal = True
    finite = bool(np.all(np.isfinite(ns.c)))
    one_plus = 1.0 + float(cs2_tolerance)
    for j in range(len(mu) - 1):
        h = float(mu[j + 1] - mu[j]); m0 = float(mu[j])
        c3, c2, c1, c0 = (float(v) for v in ns.c[:, j])
        # n(s) on s in [0, 1]
        n_s = Polynomial([c0, c1 * h, c2 * h * h, c3 * h ** 3])
        dn_ds = n_s.deriv()                               # = h dn/dmu
        mu_s = Polynomial([m0, h])
        cand = {0.0, 1.0}
        cand.update(_real_roots_in_unit_interval(dn_ds.deriv()))   # dn/dmu extrema
        cand.update(_real_roots_in_unit_interval(dn_ds))           # n extrema
        g = one_plus * mu_s * dn_ds - h * n_s                      # h*g
        g_cand = {0.0, 1.0, *_real_roots_in_unit_interval(g.deriv())}
        num, den = h * n_s, mu_s * dn_ds                           # cs2 = num/den
        stat = num.deriv() * den - num * den.deriv()
        cand.update(_real_roots_in_unit_interval(stat))
        cand.update(g_cand)
        for s in sorted(cand):
            nv = float(n_s(s)); dnv = float(dn_ds(s)) / h; muv = m0 + h * s
            if not (np.isfinite(nv) and np.isfinite(dnv)):
                finite = False; continue
            if dnv < min_dn:
                min_dn, arg_dn = dnv, muv
            min_n = min(min_n, nv)
            cs2 = nv / (muv * dnv) if dnv > 0.0 else np.inf
            if cs2 > max_cs2:
                max_cs2, arg_max = cs2, muv
            min_cs2 = min(min_cs2, cs2)
        for s in g_cand:
            gv = float(g(s)) / h
            scale = max(abs(one_plus * (m0 + h * s) * float(dn_ds(s)) / h), abs(float(n_s(s))))
            rel = gv / scale if scale > 0.0 else -np.inf
            if rel < min_margin_rel:
                min_margin_rel, arg_margin = rel, m0 + h * s
            if gv < -64.0 * epsmach * scale:
                causal = False
    stable = finite and min_dn > 0.0 and min_n > 0.0
    causal = causal and stable and max_cs2 <= one_plus
    if not finite or not (np.all(np.isfinite(p_rel)) and np.all(np.isfinite(eps_rel))):
        return {"method": "exact_interval_extrema", "stability_pass": False,
                "causality_pass": False,
                "pressure_reconstruction_max_relative_residual": float("inf"),
                "energy_reconstruction_max_relative_residual": float("inf"),
                "min_dP_dmu_fm3": float("nan"), "min_dn_dmu_fm3_per_MeV": float("nan"),
                "min_cs2_thermodynamic": float("nan"), "max_cs2_thermodynamic": float("inf"),
                "interval_count": int(len(mu) - 1),
                "identity": "n_B = dP/dmu_B by construction (P is the antiderivative of n_B)"}
    return {
        "method": "exact_interval_extrema",
        "stability_pass": bool(stable),
        "causality_pass": bool(causal),
        "causal_tolerance": float(cs2_tolerance),
        "pressure_reconstruction_max_relative_residual": float(np.max(p_rel)),
        "pressure_reconstruction_rms_relative_residual": float(np.sqrt(np.mean(p_rel ** 2))),
        "energy_reconstruction_max_relative_residual": float(np.max(eps_rel)),
        "energy_reconstruction_rms_relative_residual": float(np.sqrt(np.mean(eps_rel ** 2))),
        "min_dP_dmu_fm3": float(min_n),
        "min_dn_dmu_fm3_per_MeV": float(min_dn),
        "argmin_dn_dmu_mu_MeV": float(arg_dn),
        "min_cs2_thermodynamic": float(min_cs2),
        "max_cs2_thermodynamic": float(max_cs2),
        "argmax_cs2_mu_MeV": float(arg_max),
        "min_causal_margin_relative": float(min_margin_rel),
        "argmin_causal_margin_mu_MeV": float(arg_margin),
        "interval_count": int(len(mu) - 1),
        "samples_per_interval_ignored": int(samples_per_interval),
        "identity": "n_B = dP/dmu_B by construction (P is the antiderivative of n_B)",
        "interpolant": "P(mu_B)=P(mu_0)+integral PchipInterpolator(mu_B,n_B) dmu_B",
    }

def _phase_arrays(data: EOSData, *, differential_tolerance=5e-3, differential_samples_per_interval=8):
    n=np.asarray(data.nb,float); mu=np.asarray(data.mub,float); p=np.asarray(data.p,float); eps=np.asarray(data.eps,float)
    if n.ndim!=1 or len(n)<5 or eps.shape!=n.shape or np.any(~np.isfinite(n+mu+p+eps)) or np.any(n<=0) or np.any(mu<=0) or np.any(p<=0) or np.any(eps<=0):
        raise ValueError('phase EOS needs finite positive 1-D n, P, epsilon and mu arrays')
    if not np.isfinite(data.temperature_mev) or abs(float(data.temperature_mev))>1e-12:
        raise ValueError('Maxwell phase input must be a cold T=0 slice')
    if int(data.Il)!=1:
        raise ValueError('Maxwell beta-equilibrated phase input must include leptons (CompOSE Il=1)')
    if not np.isfinite(data.charge_fraction_grid) or abs(float(data.charge_fraction_grid))>1e-10:
        raise ValueError('Maxwell requires a locally charge-neutral Yq=0 source slice')
    if np.any(np.diff(n)<=0) or np.any(np.diff(mu)<=0) or np.any(np.diff(p)<=0) or np.any(np.diff(eps)<=0):
        raise ValueError('each phase must have strictly increasing n_b, mu_B, P and epsilon')
    if np.asarray(data.Iphase).shape == n.shape and np.any(np.diff(data.Iphase) != 0):
        raise ValueError('Maxwell input phase contains an internal Iphase phase boundary; split it into smooth branches first')
    if (not np.isfinite(differential_tolerance) or differential_tolerance <= 0
            or int(differential_samples_per_interval) < 2):
        raise ValueError('invalid branchwise first-law validation configuration')
    differential = _branchwise_first_law_report(
        mu, p, n, eps, samples_per_interval=int(differential_samples_per_interval)
    )
    if (not np.isfinite(differential['pressure_reconstruction_max_relative_residual'])
            or not np.isfinite(differential['energy_reconstruction_max_relative_residual'])
            or differential['pressure_reconstruction_max_relative_residual'] > differential_tolerance
            or differential['energy_reconstruction_max_relative_residual'] > differential_tolerance
            or differential['min_dP_dmu_fm3'] <= 0):
        raise ValueError(
            'phase thermodynamic-potential consistency failed (dP/dmu_B = n_B by construction; '
            'pressure reconstruction residual '
            f"P residual {differential['pressure_reconstruction_max_relative_residual']:.3e}, "
            f"epsilon residual {differential['energy_reconstruction_max_relative_residual']:.3e}; "
            f"tolerance {differential_tolerance:.3e})"
        )
    # The actual Maxwell interpolant is n_B(mu_B) with P obtained by
    # integration. Check mechanical stability and causality on that same
    # interpolant between knots, not only on a finite-difference knot stencil.
    if (not differential['stability_pass'] or not differential['causality_pass']
            or differential['min_dn_dmu_fm3_per_MeV'] <= 0
            or differential['min_cs2_thermodynamic'] <= 0
            or differential['max_cs2_thermodynamic'] > 1.0 + MAXWELL_CS2_TOLERANCE):
        raise ValueError(
            'thermodynamic n_B(mu_B) interpolant is unstable or acausal '
            f"(min dn/dmu={differential['min_dn_dmu_fm3_per_MeV']:.3e}, "
            f"c_s^2 range=[{differential['min_cs2_thermodynamic']:.6g},"
            f"{differential['max_cs2_thermodynamic']:.6g}])"
        )
    # A sharp coexistence construction is only meaningful when each input is
    # already a stable, causal, cold barotrope with its first-law identity.
    euler=mu*n-p
    euler_rel=float(np.max(np.abs(euler-eps)/np.maximum(np.abs(eps),1.0)))
    if euler_rel>5e-3:
        raise ValueError(f'phase EOS violates epsilon=mu_B*n_B-P (max relative residual {euler_rel:.3e})')
    if not bool(data.beta_equilibrated):
        raise ValueError('1-D Maxwell input must be a beta-equilibrated, locally neutral phase barotrope')
    from .consistency import check_consistency
    rep=check_consistency(data,verbose=False)
    rep['C12'] = {'pass_': True, **differential}
    rep['C13'] = {
        'pass_': True,
        'identity': 'c_s^2 = n_B/(mu_B dn_B/dmu_B)',
        'max_cs2': differential['max_cs2_thermodynamic'],
        'min_cs2': differential['min_cs2_thermodynamic'],
        'min_dn_dmu': differential['min_dn_dmu_fm3_per_MeV'],
        'method': 'exact interval extrema of the declared PCHIP n_B(mu_B) on every interval',
        'causality_pass': differential['causality_pass'],
    }
    for key in ('C1','C2','C6','C10','C11'):
        if not bool(rep[key].get('pass_',False)):
            raise ValueError(f'phase {key} hard gate failed; Maxwell construction requires stable causal source phases')
    if data.composition_available and (not rep['C3']['pass_'] or not rep['C4']['pass_']):
        raise ValueError('phase composition fails local baryon/charge closure')
    return n,mu,p,eps,rep


def build_maxwell_eos(low_phase: EOSData, high_phase: EOSData, *, low_name='low', high_name='high',
                      mu_samples=800, crossing_index=None, differential_tolerance=5e-3,
                      differential_samples_per_interval=8) -> MaxwellEOS:
    """Construct one exact local-neutrality Maxwell transition from two 1-D phases.

    Coexistence solves P_low(mu_B)=P_high(mu_B). The returned table contains
    two adjacent rows at the same P_t and mu_t but n_high>n_low and eps_high>eps_low.
    The two input EOS are understood as locally charge-neutral phase barotropes.
    Multiple roots require explicit ``crossing_index`` selection. Such output
    represents that interface only; it does not certify global pressure
    dominance or include the other phase transitions on the retained domain.
    """
    if int(mu_samples) < 2:
        raise ValueError('mu_samples must be at least 2')
    nl,mul,pl,el,repl=_phase_arrays(low_phase, differential_tolerance=differential_tolerance,
        differential_samples_per_interval=differential_samples_per_interval)
    nh,muh,ph,eh,reph=_phase_arrays(high_phase, differential_tolerance=differential_tolerance,
        differential_samples_per_interval=differential_samples_per_interval)
    lo=max(mul[0],muh[0]); hi=min(mul[-1],muh[-1])
    if not lo<hi: raise ValueError('the two phases have no common baryon-chemical-potential interval')
    nli, pli_primitive, pli_offset = _thermodynamic_pressure_from_density(mul, nl, pl)
    nhi, phi_primitive, phi_offset = _thermodynamic_pressure_from_density(muh, nh, ph)
    Pl=lambda x: pli_primitive(x)+pli_offset
    Ph=lambda x: phi_primitive(x)+phi_offset
    eli=PchipInterpolator(mul,el,extrapolate=False); ehi_sp=PchipInterpolator(muh,eh,extrapolate=False)
    # mu_samples is retained for API compatibility; acceptance no longer
    # depends on a uniform scan resolution.
    roots = _maxwell_pressure_equality_roots(nli, nhi, Pl, Ph, lo, hi)
    if not roots: raise ValueError('no P_low(mu_B)=P_high(mu_B) coexistence root in common interval')
    if crossing_index is None:
        if len(roots) != 1:
            raise ValueError(f'{len(roots)} coexistence roots found; choose crossing_index explicitly; '
                             'a single selected interface does not certify global pressure dominance')
        crossing_index = 0
    if isinstance(crossing_index, (bool, np.bool_)) or not isinstance(crossing_index, (int, np.integer)):
        raise ValueError('crossing_index must be an integer or None')
    if crossing_index<0 or crossing_index>=len(roots): raise ValueError(f'{len(roots)} coexistence roots found; choose crossing_index explicitly')
    mut, tangential = roots[crossing_index]
    if tangential:
        raise ValueError('tangential pressure equality has no resolved positive density jump')
    pt=float(0.5*(Pl(mut)+Ph(mut))); nlo=float(nli(mut)); nhi_t=float(nhi(mut))
    if not (pt>0 and mut>0 and nlo>0 and nhi_t>nlo):
        raise ValueError(f'coexistence root is not a physical density jump: P={pt}, mu={mut}, n_low={nlo}, n_high={nhi_t}')
    elo=float(mut*nlo-pt); epshi_t=float(mut*nhi_t-pt); de=epshi_t-elo
    euler_low_rel=abs(elo-float(eli(mut)))/max(abs(elo),1.0)
    euler_high_rel=abs(epshi_t-float(ehi_sp(mut)))/max(abs(epshi_t),1.0)
    if max(euler_low_rel,euler_high_rel)>5e-3:
        raise ValueError('interpolated coexistence point violates the phase first-law identity beyond 5e-3')
    if not (elo>0 and epshi_t>elo and de>0): raise ValueError('Maxwell energy-density jump is non-positive')
    # Keep the untouched tabulated branches and add exact interpolated endpoints.
    ml=mul<mut; mh=muh>mut
    if ml.sum()<3 or mh.sum()<3: raise ValueError('not enough tabulated points on both sides of the coexistence point')
    # Export the same thermodynamic pressure primitive used to locate the
    # coexistence root. Keep the original input immutable; quantify how much
    # P and epsilon were reconstructed from n_B(mu_B).
    p_low_consistent = np.asarray(Pl(mul), dtype=float)
    p_high_consistent = np.asarray(Ph(muh), dtype=float)
    eps_low_consistent = mul * nl - p_low_consistent
    eps_high_consistent = muh * nh - p_high_consistent
    if (np.any(~np.isfinite(p_low_consistent + eps_low_consistent))
            or np.any(~np.isfinite(p_high_consistent + eps_high_consistent))
            or np.any(eps_low_consistent <= 0) or np.any(eps_high_consistent <= 0)):
        raise ValueError('thermodynamic-potential reconstruction produced non-finite/non-positive branch energy')
    n_out=np.r_[nl[ml],nlo,nhi_t,nh[mh]]
    mu_out=np.r_[mul[ml],mut,mut,muh[mh]]
    p_out=np.r_[p_low_consistent[ml],pt,pt,p_high_consistent[mh]]
    phase_ids=np.r_[np.ones(ml.sum()+1,dtype=int),np.full(1+mh.sum(),2,dtype=int)]
    eps_out=mu_out*n_out-p_out
    if np.any(np.diff(n_out)<=0) or np.any(np.diff(p_out)<-1e-10) or np.any(np.diff(eps_out)<=0):
        raise ValueError('constructed Maxwell table is not monotone in n_b, P, and epsilon')
    mn=float(low_phase.mn); F=eps_out/n_out
    data=EOSData(nb=n_out,mn=mn,mp=float(low_phase.mp),
        Q1=p_out/n_out,Q2=np.zeros_like(n_out),Q3=mu_out/mn-1.0,
        Q4=np.zeros_like(n_out),Q5=np.zeros_like(n_out),Q6=F/mn-1.0,Q7=F/mn-1.0,
        q_extra=np.zeros_like(n_out),Iphase=phase_ids,
        Yn=np.zeros_like(n_out),Yp=np.zeros_like(n_out),Ye=np.zeros_like(n_out),Ymu=np.zeros_like(n_out),
        Aav=np.zeros_like(n_out),Zav=np.zeros_like(n_out),Yh=np.zeros_like(n_out),
        Il=1,iT=0,iYq=1,inb_start=1,temperature_mev=0.0,charge_fraction_grid=0.0,
        q_extra_count=0,composition_available=False,beta_equilibrated=True)
    il=int(ml.sum()); ih=il+1
    critical=0.5*elo+1.5*pt
    transition={'kind':'maxwell','method_id':MAXWELL_METHOD_ID,'pressure_MeV_fm3':pt,'muB_MeV':mut,
       'n_low_fm3':nlo,'n_high_fm3':nhi_t,'eps_low_MeV_fm3':elo,'eps_high_MeV_fm3':epshi_t,
       'delta_eps_MeV_fm3':de,'low_index':il,'high_index':ih,
       'seidov_critical_delta_eps_MeV_fm3':critical,
       'seidov_delta_eps_over_critical':de/critical,
       'seidov_small_core_stable':bool(de<critical),
       'tidal_jump_formula':'Takatsy-Kovacs-2020: Delta y=-4 pi r_t^3 Delta epsilon/(m_t+4 pi r_t^3 P_t)',
       'charge_neutrality':'local in each phase (Maxwell)'}
    # Store provenance and interpolation consistency at the coexistence point.
    source_e_low=float(eli(mut)); source_e_high=float(ehi_sp(mut))
    meta={'schema':'compose-eos-maxwell-v1','method_id':MAXWELL_METHOD_ID,'low_phase':low_name,'high_phase':high_name,
          'coexistence_root_count':len(roots),'selected_root_index':int(crossing_index),
          'coexistence_root_method':'union-knot quartic pressure primitives; cubic density stationary partition',
          'coexistence_roots_MeV':[float(root[0]) for root in roots],
          'coexistence_scope':'single_isolated_interface' if len(roots)==1 else 'explicitly_selected_interface',
          'global_pressure_dominance_validated':False,
          'differential_tolerance':float(differential_tolerance),
          'differential_samples_per_interval':int(differential_samples_per_interval),
          'low_energy_identity_residual':elo-source_e_low,'high_energy_identity_residual':epshi_t-source_e_high,
          'low_source_hard_gates':{k:bool(repl[k].get('pass_',False)) for k in ('C1','C2','C6','C10','C11','C12','C13')},
          'high_source_hard_gates':{k:bool(reph[k].get('pass_',False)) for k in ('C1','C2','C6','C10','C11','C12','C13')},
          'low_source_differential_first_law':repl['C12'],
          'high_source_differential_first_law':reph['C12'],
          'low_source_thermodynamic_causality':repl['C13'],
          'high_source_thermodynamic_causality':reph['C13'],
          'transition':transition,'phase_ids':phase_ids.tolist(),
          'output_pressure_interpolation':'P(mu_B)=P(mu_0)+integral PchipInterpolator(mu_B,n_B) dmu_B',
          'output_energy_identity':'epsilon=mu_B*n_B-P, enforced at every exported row',
          'low_branch_source_pressure_reconstruction_max_relative_residual':repl['C12']['pressure_reconstruction_max_relative_residual'],
          'high_branch_source_pressure_reconstruction_max_relative_residual':reph['C12']['pressure_reconstruction_max_relative_residual'],
          'low_branch_source_energy_reconstruction_max_relative_residual':repl['C12']['energy_reconstruction_max_relative_residual'],
          'high_branch_source_energy_reconstruction_max_relative_residual':reph['C12']['energy_reconstruction_max_relative_residual'],
          'interpretation':'Sharp Maxwell interface between locally-neutral 1-D barotropes using an n_B(mu_B) thermodynamic pressure primitive; source P/epsilon reconstruction residuals are explicitly recorded.'}
    return MaxwellEOS(data=data,transition=transition,metadata=meta)


def _stage_and_promote_compose_with_sidecar(
    directory: str,
    data: EOSData,
    sidecar_name: str,
    sidecar_payload: dict,
    native_revalidator: Callable[[str], object],
    *,
    kind: str,
    source_hashes: Optional[dict] = None,
) -> dict:
    """Write CompOSE files + sidecar in staging, verify native readback, then promote."""
    os.makedirs(directory, exist_ok=True)
    stage_dir = tempfile.mkdtemp(prefix=f".{kind}_stage_", dir=directory)
    backup_dir = tempfile.mkdtemp(prefix=f".{kind}_backup_", dir=directory)
    failure_receipt = os.path.join(directory, "eos.readback_failure_receipt.json")
    committed: list[str] = []
    existed: dict[str, bool] = {}
    try:
        receipt = write_compose(stage_dir, data)
        sidecar_stage = Path(stage_dir) / sidecar_name
        sidecar_stage.write_text(
            json.dumps(sidecar_payload, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        try:
            native_revalidator(stage_dir)
        except Exception as exc:
            receipt["status"] = "FAIL"
            receipt.setdefault("errors", []).append(f"{kind} native readback failed: {exc}")
            Path(failure_receipt).write_text(
                json.dumps(receipt, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            raise

        files_to_promote = [
            fn for fn in ("eos.nb", "eos.t", "eos.yq", "eos.thermo", "eos.compo", sidecar_name)
            if os.path.isfile(os.path.join(stage_dir, fn))
        ]
        all_targets = list(files_to_promote)
        if not data.composition_available and os.path.isfile(os.path.join(directory, "eos.compo")):
            all_targets.append("eos.compo")

        for fn in all_targets:
            dst = os.path.join(directory, fn)
            existed[fn] = os.path.isfile(dst)
            if existed[fn]:
                shutil.copy2(dst, os.path.join(backup_dir, fn))

        for fn in files_to_promote:
            os.replace(os.path.join(stage_dir, fn), os.path.join(directory, fn))
            committed.append(fn)

        if not data.composition_available:
            cp_dst = os.path.join(directory, "eos.compo")
            if os.path.exists(cp_dst):
                os.remove(cp_dst)
                committed.append("eos.compo")

        final_hashes = _compose_dir_hashes(directory)
        final_hashes[sidecar_name] = _sha256_file(os.path.join(directory, sidecar_name))
        receipt["directory"] = os.path.abspath(directory)
        receipt["operation"] = kind
        receipt["sidecar"] = sidecar_name
        receipt["output_sha256"] = final_hashes
        if source_hashes:
            receipt["source_sha256"] = source_hashes
        receipt["native_sidecar_readback"] = "PASS"
        Path(directory, "eos.readback_receipt.json").write_text(
            json.dumps(receipt, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return receipt
    except Exception:
        for fn in committed:
            dst = os.path.join(directory, fn)
            if existed.get(fn, False):
                shutil.copy2(os.path.join(backup_dir, fn), dst)
            elif os.path.exists(dst):
                os.remove(dst)
        raise
    finally:
        shutil.rmtree(stage_dir, ignore_errors=True)
        shutil.rmtree(backup_dir, ignore_errors=True)


def write_maxwell_eos(directory: str, eos: MaxwellEOS) -> dict:
    """Write CompOSE arrays + Maxwell sidecar via staged native readback verification."""
    payload = {
        **eos.metadata,
        "transition": eos.transition,
        "phase_ids": eos.data.Iphase.astype(int).tolist(),
    }
    return _stage_and_promote_compose_with_sidecar(
        directory,
        eos.data,
        "eos.phase_transitions.json",
        payload,
        read_maxwell_eos,
        kind="maxwell",
    )


def read_maxwell_eos(directory: str) -> MaxwellEOS:
    """Read and revalidate a sharp Maxwell EOS and its phase-event sidecar.

    The sidecar is not trusted as a validation receipt: each smooth phase is
    checked again against the pressure primitive generated from n_B(mu_B).
    This catches post-write corruption away from the transition rows.
    """
    data = read_compose(directory)
    path = Path(directory) / "eos.phase_transitions.json"
    if not path.is_file():
        raise FileNotFoundError(
            "EOS directory has no eos.phase_transitions.json sharp-interface metadata"
        )
    meta = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(meta, dict):
        raise ValueError("Maxwell sidecar must contain a JSON object")
    transition = meta.get("transition")
    if (meta.get("schema") != "compose-eos-maxwell-v1"
            or not isinstance(transition, dict)
            or transition.get("kind") != "maxwell"
            or meta.get("method_id", MAXWELL_METHOD_ID) != MAXWELL_METHOD_ID):
        raise ValueError("unsupported or malformed Maxwell sidecar schema")

    # The current Maxwell exporter is deliberately thermo-only. An eos.compo
    # file appearing beside this output is not silently accepted as validated.
    if data.composition_available:
        raise ValueError("Maxwell readback is thermo-only; unexpected eos.compo is present")

    numeric_keys = (
        "pressure_MeV_fm3", "muB_MeV", "n_low_fm3", "n_high_fm3",
        "eps_low_MeV_fm3", "eps_high_MeV_fm3", "delta_eps_MeV_fm3",
    )
    try:
        numeric = {key: float(transition[key]) for key in numeric_keys}
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Maxwell sidecar is missing a numeric transition field") from exc
    if not all(np.isfinite(value) for value in numeric.values()):
        raise ValueError("Maxwell sidecar transition fields must all be finite")

    def _exact_index(key):
        try:
            value = float(transition[key])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Maxwell sidecar {key} must be an integer index") from exc
        if not np.isfinite(value) or value != np.rint(value):
            raise ValueError(f"Maxwell sidecar {key} must be an integer index")
        return int(value)

    il = _exact_index("low_index")
    ih = _exact_index("high_index")
    nrows = data.n
    if ih != il + 1 or not (0 <= il < ih < nrows):
        raise ValueError("Maxwell sidecar indices do not identify adjacent in-range transition rows")

    try:
        raw_ids = np.asarray(meta.get("phase_ids", []), dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("Maxwell sidecar phase_ids must be an integer array") from exc
    if (raw_ids.shape != (nrows,) or np.any(~np.isfinite(raw_ids))
            or np.any(raw_ids != np.rint(raw_ids))):
        raise ValueError("Maxwell sidecar phase_ids must be finite integers with one ID per EOS row")
    ids = raw_ids.astype(int)
    expected_ids = np.r_[np.ones(ih, dtype=int), np.full(nrows - ih, 2, dtype=int)]
    if not np.array_equal(ids, expected_ids):
        raise ValueError("Maxwell sidecar phase_ids must be contiguous low-phase 1 then high-phase 2")
    data.Iphase = ids

    n = np.asarray(data.nb, dtype=float)
    mu = np.asarray(data.mub, dtype=float)
    pressure = np.asarray(data.p, dtype=float)
    eps = np.asarray(data.eps, dtype=float)
    if (any(array.shape != (nrows,) for array in (n, mu, pressure, eps))
            or any(np.any(~np.isfinite(array)) for array in (n, mu, pressure, eps))
            or any(np.any(array <= 0.0) for array in (n, mu, pressure, eps))):
        raise ValueError("Maxwell readback requires finite, positive n_B, mu_B, P and epsilon arrays")
    if np.any(np.diff(n) <= 0.0) or np.any(np.diff(eps) <= 0.0):
        raise ValueError("Maxwell readback density and energy density must increase across the full table")

    pt = numeric["pressure_MeV_fm3"]
    mut = numeric["muB_MeV"]
    comparisons = (
        ("P_low", pressure[il], pt),
        ("P_high", pressure[ih], pt),
        ("mu_low", mu[il], mut),
        ("mu_high", mu[ih], mut),
        ("n_low", n[il], numeric["n_low_fm3"]),
        ("n_high", n[ih], numeric["n_high_fm3"]),
        ("eps_low", eps[il], numeric["eps_low_MeV_fm3"]),
        ("eps_high", eps[ih], numeric["eps_high_MeV_fm3"]),
    )
    for label, actual, expected in comparisons:
        if abs(float(actual) - float(expected)) > 2e-8 * max(abs(float(expected)), 1.0):
            raise ValueError(f"Maxwell sidecar {label} mismatch: table={actual}, metadata={expected}")

    if not (pressure[il] > 0 and mut > 0 and n[ih] > n[il] and eps[ih] > eps[il]):
        raise ValueError("Maxwell transition must have positive P, mu_B and upward density/energy jumps")

    delta_n = float(n[ih] - n[il])
    delta_eps = float(eps[ih] - eps[il])
    identity_value = mut * delta_n
    identity_scale = max(abs(delta_eps), abs(identity_value), 1.0)
    declared_delta_eps = numeric["delta_eps_MeV_fm3"]
    if (declared_delta_eps <= 0.0
            or abs(declared_delta_eps - delta_eps) > 2e-8 * identity_scale
            or abs(delta_eps - identity_value) > 2e-8 * identity_scale):
        raise ValueError("Maxwell readback violates Delta epsilon = mu_B Delta n_B")

    try:
        differential_tolerance = float(meta.get("differential_tolerance", 5e-3))
        samples_per_interval = int(meta.get("differential_samples_per_interval", 8))
    except (TypeError, ValueError) as exc:
        raise ValueError("Maxwell readback has invalid branchwise validation settings") from exc
    if (not np.isfinite(differential_tolerance) or differential_tolerance <= 0.0
            or samples_per_interval < 2):
        raise ValueError("Maxwell readback has invalid branchwise validation settings")

    branch_reports = {}
    for label, lo, hi in (("low", 0, ih), ("high", ih, nrows)):
        if hi - lo < 4:
            raise ValueError(f"Maxwell readback {label} branch needs at least four rows")
        if (np.any(np.diff(mu[lo:hi]) <= 0.0)
                or np.any(np.diff(pressure[lo:hi]) <= 0.0)
                or np.any(np.diff(eps[lo:hi]) <= 0.0)):
            raise ValueError(f"Maxwell readback {label} branch is not strictly monotone in mu_B, P and epsilon")
        report = _branchwise_first_law_report(
            mu[lo:hi], pressure[lo:hi], n[lo:hi], eps[lo:hi],
            samples_per_interval=samples_per_interval,
        )
        pressure_residual = report["pressure_reconstruction_max_relative_residual"]
        energy_residual = report["energy_reconstruction_max_relative_residual"]
        if (not np.isfinite(pressure_residual) or not np.isfinite(energy_residual)
                or pressure_residual > differential_tolerance
                or energy_residual > differential_tolerance):
            raise ValueError(
                f"Maxwell readback {label} branch pressure reconstruction failed "
                f"(P residual={pressure_residual:.3e}, epsilon residual={energy_residual:.3e}, "
                f"tolerance={differential_tolerance:.3e})"
            )
        if (not report["stability_pass"] or not report["causality_pass"]
                or not np.isfinite(report["min_dn_dmu_fm3_per_MeV"])
                or report["min_dn_dmu_fm3_per_MeV"] <= 0.0
                or not np.isfinite(report["min_cs2_thermodynamic"])
                or report["min_cs2_thermodynamic"] <= 0.0
                or not np.isfinite(report["max_cs2_thermodynamic"])
                or report["max_cs2_thermodynamic"] > 1.0 + MAXWELL_CS2_TOLERANCE):
            raise ValueError(
                f"Maxwell readback {label} branch fails mechanical-stability/causality checks "
                f"(dn/dmu min={report['min_dn_dmu_fm3_per_MeV']:.3e}, "
                f"cs2=[{report['min_cs2_thermodynamic']:.6g},"
                f"{report['max_cs2_thermodynamic']:.6g}])"
            )
        branch_reports[label] = {
            "pressure_reconstruction_max_relative_residual": float(pressure_residual),
            "energy_reconstruction_max_relative_residual": float(energy_residual),
            "min_dn_dmu_fm3_per_MeV": float(report["min_dn_dmu_fm3_per_MeV"]),
            "min_cs2": float(report["min_cs2_thermodynamic"]),
            "max_cs2": float(report["max_cs2_thermodynamic"]),
            "intervals": int(report["interval_count"]),
            "argmax_cs2_mu_MeV": float(report["argmax_cs2_mu_MeV"]),
        }

    meta = dict(meta)
    meta["readback_validation"] = {
        "status": "PASS",
        "method": "branchwise n_B(mu_B) PCHIP pressure primitive; exact interval extrema",
        "differential_tolerance": differential_tolerance,
        "branches": branch_reports,
        "transition_identity_relative_residual": float(
            abs(delta_eps - identity_value) / identity_scale
        ),
    }
    return MaxwellEOS(data, transition, meta)


def tidal_y_jump(y_inner: float, radius_km: float, mass_geometric_km: float,
                 pressure_geometric_km2: float, delta_epsilon_geometric_km2: float) -> float:
    """Apply the pressure-corrected l=2 tidal jump in G=c=1 units.

    ``mass_geometric_km`` is m(r_t) in km; pressure and energy density are km^-2.
    The returned value is y just outside the interface.
    """
    vals=np.asarray([y_inner,radius_km,mass_geometric_km,pressure_geometric_km2,delta_epsilon_geometric_km2],float)
    if np.any(~np.isfinite(vals)) or radius_km<=0 or mass_geometric_km<0 or delta_epsilon_geometric_km2<=0:
        raise ValueError('invalid tidal jump arguments')
    den=mass_geometric_km+4*np.pi*radius_km**3*pressure_geometric_km2
    if den<=0: raise ValueError('non-positive tidal-jump denominator')
    return float(y_inner-4*np.pi*radius_km**3*delta_epsilon_geometric_km2/den)


# Exact unit conversion: (MeV fm^-3) -> km^-2 in G=c=1.
_G=6.67430e-11; _C=299792458.0; _MEV=1.602176634e-13; _MSUN=1.988409870698051e30
_MEVFM3_TO_KM2=_G*_MEV/(1e-45*_C**4)*1e6
_MSUN_GEOM_KM=_G*_MSUN/_C**2/1000.0


def solve_tov_with_phase_transition(eos: MaxwellEOS, central_pressure_MeV_fm3: float,
                                    *, rtol=2e-8, max_step_km=0.12) -> dict:
    """Solve TOV plus l=2 Love equation with one sharp Maxwell interface.

    Integrates the high-density phase to P_t, applies the corrected y-jump,
    then integrates the low-density phase to the finite tabulated surface.
    All stellar structure equations are in G=c=1 geometric units.
    """
    from .tov import BarotropicEOS, find_causal_domain
    if not np.isfinite(rtol) or not 0.0 < rtol < 1.0:
        raise ValueError('rtol must be finite and in (0, 1)')
    if not np.isfinite(max_step_km) or max_step_km <= 0.0:
        raise ValueError('max_step_km must be finite and positive')
    d=eos.data; t=eos.transition; pt=float(t['pressure_MeV_fm3'])
    pc=float(central_pressure_MeV_fm3); pgrid=d.p; eps=d.eps
    if not np.isfinite(pc) or pc<=pgrid[0] or pc>pgrid[-1]: raise ValueError('central pressure outside tabulated EOS range')
    il=int(t['low_index']); ih=int(t['high_index'])
    def make_branch(start,stop):
        pp=np.asarray(pgrid[start:stop],float); ee=np.asarray(eps[start:stop],float)
        if len(pp)<4 or np.any(np.diff(pp)<=0) or np.any(np.diff(ee)<=0): raise ValueError('each smooth phase branch needs >=4 strictly increasing P/epsilon points')
        barotrope = BarotropicEOS(np.asarray(d.nb[start:stop], float), pp, ee)
        sp = barotrope._lneps_of_lnp
        def state(P):
            values = np.asarray(P, dtype=float)
            if np.any(~np.isfinite(values)):
                raise FloatingPointError('non-finite pressure query in phase branch')
            # solve_ivp can return P_t +/- a few ulps at the phase interface.
            # Permit only roundoff-scale endpoint excursions; never extrapolate
            # the EOS interpolant across a physically finite pressure interval.
            endpoint_tol = 64.0 * np.finfo(float).eps * max(1.0, abs(float(pp[0])), abs(float(pp[-1])))
            if np.any(values < pp[0] - endpoint_tol) or np.any(values > pp[-1] + endpoint_tol):
                raise FloatingPointError('pressure query lies outside selected phase branch')
            clipped = np.clip(values, pp[0], pp[-1])
            lp=np.log(clipped); ev=np.exp(sp(lp)); dedp=ev/clipped*sp.derivative()(lp)
            if np.any(~np.isfinite(ev+dedp)) or np.any(dedp<=0): raise FloatingPointError('invalid d epsilon/dP in phase branch')
            return ev,dedp
        return (pp,ee,state,barotrope)
    low=make_branch(0,il+1); high=make_branch(ih,d.n)
    active_high=pc>pt
    # A causal center does not imply a causal shell. Certify each traversed
    # pressure range on the actual epsilon(P) RHS interpolant. The independent
    # Maxwell n(mu) potential has a different interpolation contract.
    traversed = [('low', low, pt if active_high else pc)]
    if active_high:
        traversed.append(('high', high, pc))
    causal_domains = {}
    for label, phase, pressure_end in traversed:
        domain = find_causal_domain(phase[3], tolerance=1e-5)
        if pressure_end > domain.pressure_max * (1.0 + 2e-12):
            raise ValueError(f'{label} phase pressure path exceeds its causal domain '
                             f'ending at {domain.pressure_max:g} MeV fm^-3')
        causal_domains[label] = domain.to_dict()
    branch=high if active_high else low
    ec,de_dp=branch[2](pc); ec=float(ec); de_dp=float(de_dp); cs2c=1.0/de_dp
    if not (0<cs2c<=1+1e-5): raise ValueError(f'central sound speed is noncausal: {cs2c:g}')
    pcg=pc*_MEVFM3_TO_KM2; ecg=ec*_MEVFM3_TO_KM2
    r0=1e-4
    p0=pcg-(2*np.pi/3)*(ecg+pcg)*(ecg+3*pcg)*r0*r0
    m0=4*np.pi/3*ecg*r0**3
    y0=2-(4*np.pi/7)*(5*ecg+9*pcg+(ecg+pcg)/cs2c)*r0*r0
    if p0<=0 or (active_high and p0<=pt*_MEVFM3_TO_KM2): raise ValueError('central expansion step is too large or starts below the Maxwell interface')
    transitions=[]
    def integrate(beg_state, r_start, target_P, phase):
        pmin=float(phase[0][0]); pmax=float(phase[0][-1])
        if target_P<pmin*(1-1e-10) or target_P>pmax*(1+1e-10): raise ValueError('TOV event lies outside selected phase table')
        target_P_geom = target_P * _MEVFM3_TO_KM2
        def rhs(r,s):
            # The integrated pressure state is geometric (km^-2). Convert back
            # to CompOSE MeV/fm^3 before querying the phase interpolator.
            P_geom=max(float(s[0]),target_P_geom)
            P_mev=P_geom/_MEVFM3_TO_KM2
            eps_v,dedp_v=phase[2](P_mev)
            pg=P_geom; eg=float(eps_v)*_MEVFM3_TO_KM2
            m,y=s[1],s[2]; fac=1-2*m/r
            if fac<=0: raise FloatingPointError('trapped-surface condition reached')
            F=(1-4*np.pi*r*r*(eg-pg))/fac
            nu_prime=2*(m+4*np.pi*r**3*pg)/(r*(r-2*m))
            Q=4*np.pi*(5*eg+9*pg+(eg+pg)*float(dedp_v))/fac-6/(r*r*fac)-nu_prime**2
            dpdr=-(eg+pg)*(m+4*np.pi*r**3*pg)/(r*(r-2*m))
            dmdr=4*np.pi*r*r*eg
            dydr=-(y*y+y*F+r*r*Q)/r
            return [dpdr,dmdr,dydr]
        def event(r,s): return s[0]-target_P_geom
        event.terminal=True; event.direction=-1
        sol=solve_ivp(rhs,(r_start,40.0),beg_state,method='DOP853',rtol=rtol,atol=[1e-13,1e-12,1e-10],max_step=max_step_km,events=event)
        if not sol.success or not sol.t_events or len(sol.t_events[0])!=1: raise RuntimeError(f'TOV phase integration missed pressure event {target_P:g}: {sol.message}')
        return float(sol.t_events[0][0]),[float(v) for v in sol.y_events[0][0]]
    r=r0; state=[p0,m0,y0]
    if active_high:
        r,state=integrate(state,r,pt,high)
        de=float(t['delta_eps_MeV_fm3'])*_MEVFM3_TO_KM2
        pgeom=pt*_MEVFM3_TO_KM2
        y_before=state[2]; state[2]=tidal_y_jump(y_before,r,state[1],pgeom,de)
        transitions.append({'radius_km':r,'mass_geometric_km':float(state[1]),'pressure_geometric_km2':pgeom,
          'delta_epsilon_geometric_km2':de,'y_inside':y_before,'y_outside':state[2],'delta_y':state[2]-y_before,
          'delta_y_formula':'-4*pi*r_t^3*Delta_epsilon/(m_t+4*pi*r_t^3*P_t)'})
        r,state=integrate(state,r,float(pgrid[0]),low)
    else:
        r,state=integrate(state,r,float(pgrid[0]),low)
    R=r; Mgeom=state[1]; yR=state[2]; C=Mgeom/R
    if not (0<C<0.5): raise RuntimeError(f'unphysical compactness {C:g}')
    t1=(8/5)*C**5*(1-2*C)**2*(2+2*C*(yR-1)-yR)
    t2=2*C*(6-3*yR+3*C*(5*yR-8))
    t3=4*C**3*(13-11*yR+C*(3*yR-2)+2*C**2*(1+yR))
    t4=3*(1-2*C)**2*(2-yR+2*C*(yR-1))*np.log1p(-2*C)
    den=t2+t3+t4
    if abs(den)<1e-30: raise RuntimeError('Love-number denominator is singular')
    k2=float(t1/den); lam=float((2/3)*k2/C**5)
    return {'status':'SUCCESS','mass_Msun':Mgeom/_MSUN_GEOM_KM,'radius_km':R,'compactness':C,
            'k2':k2,'Lambda':lam,'central_pressure_MeV_fm3':pc,'surface_y':yR,
            'phase_interfaces':transitions,'tidal_jump_formula':'Takatsy-Kovacs-2020 Eq. (11)',
            'rhs_causal_domains':causal_domains,
            'seidov_small_core_stable':bool(eos.transition['seidov_small_core_stable']),
            'radial_stability_scope':'Seidov small-core onset criterion only; no full radial-mode eigenfrequency calculation is performed.'}


@dataclass
class ChargeResolvedPhase:
    """Tabulated two-charge phase surface P(mu_B,mu_Q), in MeV fm^-3."""
    mu_b: np.ndarray
    mu_q: np.ndarray
    pressure: np.ndarray
    label: str = 'phase'
    source_id: str = 'in_memory_table'
    source_file_path: Optional[str] = None
    source_sha256: Optional[str] = None
    model_id: Optional[str] = None
    phase_identity: str = ''
    mu_b_units: str = 'MeV'
    mu_q_units: str = 'MeV'
    pressure_units: str = 'MeV fm^-3'
    normalization_convention: str = 'grand-potential pressure P(mu_B, mu_Q) relative to zero-pressure vacuum'
    interpolation_method: str = 'RectBivariateSpline'
    interpolation_order: tuple[int, int] = (3, 3)
    lepton_contributions: str = 'excluded_baryonic_only'
    _spline: object = field(init=False, repr=False)

    def __post_init__(self):
        self.mu_b = np.asarray(self.mu_b, float)
        self.mu_q = np.asarray(self.mu_q, float)
        self.pressure = np.asarray(self.pressure, float)
        if self.mu_b.ndim != 1 or self.mu_q.ndim != 1 or len(self.mu_b) < 4 or len(self.mu_q) < 4 or self.pressure.shape != (len(self.mu_b), len(self.mu_q)):
            raise ValueError('two-charge surface needs >=4x4 pressure grid with shape (mu_b,mu_q)')
        if (np.any(~np.isfinite(self.mu_b)) or np.any(~np.isfinite(self.mu_q))
                or np.any(~np.isfinite(self.pressure))):
            raise ValueError('phase surface grids and pressure values must be finite')
        if np.any(np.diff(self.mu_b) <= 0) or np.any(np.diff(self.mu_q) <= 0):
            raise ValueError('phase surface grids must be strictly increasing')
        if not self.phase_identity:
            self.phase_identity = str(self.label)
        if self.model_id is None:
            self.model_id = str(self.label)
        if self.source_file_path:
            fpath = os.path.abspath(str(self.source_file_path))
            self.source_file_path = fpath
            if self.source_sha256 is None and os.path.isfile(fpath):
                self.source_sha256 = _sha256_file(fpath)
        order = tuple(int(k) for k in self.interpolation_order)
        if len(order) != 2:
            raise ValueError('interpolation_order must be a 2-tuple (kx, ky)')
        self.interpolation_order = order
        self._spline = RectBivariateSpline(self.mu_b, self.mu_q, self.pressure, kx=order[0], ky=order[1], s=0)

    def provenance_dict(self) -> dict:
        return {
            'source_identifier': str(self.source_id),
            'source_file_path': self.source_file_path,
            'source_sha256': self.source_sha256,
            'model_parameter_set_identifier': str(self.model_id),
            'phase_identity': str(self.phase_identity),
            'mu_b_units': str(self.mu_b_units),
            'mu_q_units': str(self.mu_q_units),
            'pressure_units': str(self.pressure_units),
            'normalization_convention': str(self.normalization_convention),
            'interpolation_method': str(self.interpolation_method),
            'interpolation_order': list(self.interpolation_order),
            'lepton_contributions': str(self.lepton_contributions),
            'mu_b_range': [float(self.mu_b[0]), float(self.mu_b[-1])],
            'mu_q_range': [float(self.mu_q[0]), float(self.mu_q[-1])],
            'grid_shape': [int(len(self.mu_b)), int(len(self.mu_q))],
        }

    @classmethod
    def from_npz(cls, path, label='phase', **metadata):
        path_str = os.path.abspath(str(path))
        file_sha = _sha256_file(path_str)
        z = np.load(path_str, allow_pickle=False)
        def _opt_str(key: str, default):
            if key in metadata and metadata[key] is not None:
                return metadata[key]
            if key in z:
                val = z[key]
                return str(val.item() if np.ndim(val) == 0 else val)
            return default
        return cls(
            z['mu_b'],
            z['mu_q'],
            z['pressure'],
            label=label,
            source_id=_opt_str('source_id', Path(path_str).stem),
            source_file_path=path_str,
            source_sha256=file_sha,
            model_id=_opt_str('model_id', label),
            phase_identity=_opt_str('phase_identity', label),
            mu_b_units=_opt_str('mu_b_units', 'MeV'),
            mu_q_units=_opt_str('mu_q_units', 'MeV'),
            pressure_units=_opt_str('pressure_units', 'MeV fm^-3'),
            normalization_convention=_opt_str(
                'normalization_convention',
                'grand-potential pressure P(mu_B, mu_Q) relative to zero-pressure vacuum',
            ),
            lepton_contributions=_opt_str('lepton_contributions', 'excluded_baryonic_only'),
        )
    def evaluate(self,mu_b,mu_q):
        if not(self.mu_b[0]<=mu_b<=self.mu_b[-1] and self.mu_q[0]<=mu_q<=self.mu_q[-1]): raise ValueError('chemical potentials outside tabulated phase surface')
        P=float(self._spline.ev(mu_b,mu_q)); n=float(self._spline.ev(mu_b,mu_q,dx=1,dy=0)); q=float(self._spline.ev(mu_b,mu_q,dx=0,dy=1))
        hbb=float(self._spline.ev(mu_b,mu_q,dx=2,dy=0)); hbq=float(self._spline.ev(mu_b,mu_q,dx=1,dy=1)); hqq=float(self._spline.ev(mu_b,mu_q,dx=0,dy=2))
        hess=np.array([[hbb,hbq],[hbq,hqq]],float); eig=np.linalg.eigvalsh(hess)
        scale=max(float(np.max(np.abs(hess))),1e-12)
        if np.min(eig)<-1e-7*scale: raise ValueError(f'{self.label} pressure surface is not convex in (mu_B,mu_Q): Hessian eigenvalues {eig}')
        if not np.isfinite(P+n+q) or n<=0: raise ValueError(f'{self.label} has non-positive/non-finite baryon density')
        return P,n,q

@dataclass
class LeptonSector:
    pressure: Callable[[float],float]
    charge_density: Callable[[float],float]
    label: str='custom-lepton-sector'
    @staticmethod
    def electron_muon():
        hbarc=197.3269804; masses=(0.51099895,105.6583755)
        def p_of_mu_q(mu_q):
            mu_e=-float(mu_q); out=0.0
            for mass in masses:
                if mu_e>mass:
                    k=np.sqrt(mu_e**2-mass**2)
                    out+=(mu_e*k*(2*mu_e**2-5*mass**2)+3*mass**4*np.log((mu_e+k)/mass))/(24*np.pi**2*hbarc**3)
            return float(out)
        def q_of_mu_q(mu_q):
            mu_e=-float(mu_q); out=0.0
            for mass in masses:
                if mu_e>mass: out-=(mu_e**2-mass**2)**1.5/(3*np.pi**2*hbarc**3)
            return float(out)
        return LeptonSector(p_of_mu_q,q_of_mu_q,'zero-temperature electron+muon Fermi gases')



def _gibbs_pressure_equality_roots(hadronic, high_density, mu_b, q_lo, q_hi):
    """Isolate all roots of Delta P(mu_Q) for the declared bicubic surfaces.

    At fixed mu_B, Delta P is piecewise cubic in mu_Q and its derivative
    q_H-q_Q is piecewise quadratic. The union of spline knots partitions the
    domain, allowing all stationary points and hence all roots to be isolated.
    """
    if not q_lo < q_hi:
        return []

    def dp(q):
        return float(hadronic.evaluate(mu_b, q)[0] - high_density.evaluate(mu_b, q)[0])

    def dq(q):
        return float(hadronic.evaluate(mu_b, q)[2] - high_density.evaluate(mu_b, q)[2])

    parts=[np.asarray([q_lo,q_hi],float)]
    for phase in (hadronic,high_density):
        _,knots=phase._spline.get_knots()
        knots=np.asarray(knots,float)
        parts.append(knots[(knots>=q_lo)&(knots<=q_hi)])
    breaks=np.unique(np.concatenate(parts))
    stationary=[]; eps=np.finfo(float).eps
    for left,right in zip(breaks[:-1],breaks[1:]):
        width=float(right-left)
        if width<=0: continue
        mid=0.5*(left+right)
        d0,dm,d1=dq(left),dq(mid),dq(right)
        a=2.0*(d1+d0-2.0*dm); b=d1-d0-a; c=d0
        scale=max(abs(d0),abs(dm),abs(d1),1e-300); tol=128.0*eps*scale
        if max(abs(a),abs(b),abs(c))<=tol:
            vals=[dp(left),dp(mid),dp(right)]
            pscale=max(1.0,abs(hadronic.evaluate(mu_b,mid)[0]),abs(high_density.evaluate(mu_b,mid)[0]))
            if max(abs(v) for v in vals)<=128.0*eps*pscale:
                raise ValueError(
                    "Gibbs pressure equality contains a non-isolated root interval "
                    f"at mu_B={mu_b:g}, mu_Q=[{left:g},{right:g}]"
                )
            continue
        if abs(a)<=tol:
            tr=[] if abs(b)<=tol else [-c/b]
        else:
            tr=[float(r.real) for r in np.roots([a,b,c])
                if abs(float(r.imag))<=1e-10*max(1.0,abs(float(r.real)))]
        candidates=sorted(left+width*t for t in tr if np.isfinite(t) and 0.0<t<1.0)
        for i, guess in enumerate(candidates):
            lo = left if i == 0 else 0.5 * (candidates[i - 1] + guess)
            hi = right if i + 1 == len(candidates) else 0.5 * (guess + candidates[i + 1])
            if np.signbit(dq(lo)) != np.signbit(dq(hi)):
                guess = float(brentq(dq, lo, hi, xtol=1e-15, rtol=4.0 * eps))
            stationary.append(guess)

    points = np.unique(np.asarray([q_lo, q_hi, *stationary], float))
    roots = []
    for q in points:
        ph = hadronic.evaluate(mu_b, float(q))[0]
        pq = high_density.evaluate(mu_b, float(q))[0]
        if abs(float(ph - pq)) <= 128.0 * eps * max(1.0, abs(float(ph)), abs(float(pq))):
            roots.append(float(q))
    for lo, hi in zip(points[:-1], points[1:]):
        flo, fhi = dp(float(lo)), dp(float(hi))
        if np.signbit(flo) != np.signbit(fhi):
            roots.append(float(brentq(dp, float(lo), float(hi), xtol=1e-15, rtol=4.0 * eps)))
    roots.sort()
    qtol = 128.0 * eps * max(1.0, abs(q_lo), abs(q_hi))
    unique = []
    for root in roots:
        if not unique or abs(root - unique[-1]) > qtol:
            unique.append(root)
    return unique


def validate_gibbs_barotrope(
    mu_b_mev,
    nb_fm3,
    p_mev_fm3,
    eps_mev_fm3,
    *,
    stage: str = "post_construction",
    euler_rtol: float = 1e-6,
    euler_atol: float = 1e-10,
    pressure_rtol: float = 5e-4,
    pressure_atol: float = 5e-4,
    causality_margin_tol: float = 1e-11,
    cs2_tol: float = 1e-10,
) -> dict:
    """Validate thermodynamic identity, positivity, and exact interval causality
    for a Gibbs mixed-phase barotrope using one PCHIP representation of n_B(mu_B).

    Checks over the full piecewise-cubic spline domain (never endpoints alone):
      1. dP/dmu_B = n_B via exact piecewise integration of PCHIP n_B(mu_B);
      2. Exact interval extrema of dn_B/dmu_B > 0 (thermodynamic stability);
      3. Exact interval extrema of causality margin g = mu_B * dn_B/dmu_B - n_B >= 0
         and 0 < c_s^2 = n_B / (mu_B * dn_B/dmu_B) <= 1;
      4. Euler relation eps + P = mu_B * n_B.
    """
    mu_b = np.asarray(mu_b_mev, dtype=float)
    nb = np.asarray(nb_fm3, dtype=float)
    p = np.asarray(p_mev_fm3, dtype=float)
    eps = np.asarray(eps_mev_fm3, dtype=float)
    if (mu_b.ndim != 1 or nb.shape != mu_b.shape or p.shape != mu_b.shape
            or eps.shape != mu_b.shape or mu_b.size < 4):
        raise ValueError(f"Gibbs {stage} validation requires >= 4 finite points of matching 1D shape")
    if any(np.any(~np.isfinite(arr)) for arr in (mu_b, nb, p, eps)):
        raise ValueError(f"Gibbs {stage} validation found non-finite thermodynamic values")
    if np.any(mu_b <= 0.0) or np.any(np.diff(mu_b) <= 0.0):
        raise ValueError(f"Gibbs {stage} requires strictly increasing positive mu_B grid")
    if np.any(nb <= 0.0) or np.any(p <= 0.0) or np.any(eps <= 0.0):
        raise ValueError(f"Gibbs {stage} requires positive n_B, P, and energy density")
    if np.any(np.diff(nb) <= 0.0) or np.any(np.diff(p) <= 0.0) or np.any(np.diff(eps) <= 0.0):
        raise ValueError(f"Gibbs {stage} requires strictly increasing (n_B, P, eps)")

    # Euler identity: eps + P = mu_B * n_B
    euler_lhs = eps + p
    euler_rhs = mu_b * nb
    euler_abs = np.abs(euler_lhs - euler_rhs)
    euler_scale = np.maximum(np.maximum(np.abs(euler_lhs), np.abs(euler_rhs)), 1e-30)
    euler_rel = euler_abs / euler_scale
    euler_max_abs = float(np.max(euler_abs))
    euler_max_rel = float(np.max(euler_rel))
    euler_rms_rel = float(np.sqrt(np.mean(euler_rel ** 2)))
    if not np.allclose(euler_lhs, euler_rhs, rtol=euler_rtol, atol=euler_atol):
        raise ValueError(
            f"Gibbs {stage} violates Euler identity eps + P = mu_B * n_B "
            f"(max_rel={euler_max_rel:.3e}, max_abs={euler_max_abs:.3e})"
        )

    # Single declared interpolation representation: PCHIP of n_B(mu_B)
    nb_spline = PchipInterpolator(mu_b, nb, extrapolate=False)
    antideriv = nb_spline.antiderivative()
    p_reconstructed = p[0] + (antideriv(mu_b) - antideriv(mu_b[0]))
    p_abs = np.abs(p_reconstructed - p)
    p_scale = np.maximum(np.maximum(np.abs(p_reconstructed), np.abs(p)), 1e-30)
    p_rel = p_abs / p_scale
    p_max_abs = float(np.max(p_abs))
    p_rms_abs = float(np.sqrt(np.mean(p_abs ** 2)))
    p_max_rel = float(np.max(p_rel))
    p_rms_rel = float(np.sqrt(np.mean(p_rel ** 2)))
    if not np.allclose(p_reconstructed, p, rtol=pressure_rtol, atol=pressure_atol):
        raise ValueError(
            f"Gibbs {stage} fails dP/dmu_B = n_B integral reconstruction "
            f"(max_abs={p_max_abs:.3e}, max_rel={p_max_rel:.3e})"
        )

    min_dn = float("inf")
    max_dn = float("-inf")
    min_margin = float("inf")
    max_margin = float("-inf")
    min_cs2 = float("inf")
    max_cs2 = float("-inf")
    worst_cs2_mub = float(mu_b[0])
    worst_margin_mub = float(mu_b[0])
    worst_dn_mub = float(mu_b[0])

    n_intervals = len(mu_b) - 1
    for j in range(n_intervals):
        h = float(mu_b[j + 1] - mu_b[j])
        mu0 = float(mu_b[j])
        c3, c2, c1, c0 = [float(v) for v in nb_spline.c[:, j]]
        poly_n = np.array([c3, c2, c1, c0], dtype=float)
        poly_dn = np.array([3.0 * c3, 2.0 * c2, c1], dtype=float)
        poly_ddn = np.array([6.0 * c3, 2.0 * c2], dtype=float)
        poly_mu = np.array([1.0, mu0], dtype=float)

        # Stationary points of dn/dmu_B AND of g = mu_B * dn/dmu_B - n_B:
        # g'(t) = (mu0 + t) * n''(t), and mu0 + t > 0 => only inflection point t* = -c2/(3*c3).
        cand_dn = [0.0, h]
        if abs(c3) > 1e-30 * max(1.0, abs(c2), abs(c1)):
            t_flex = -c2 / (3.0 * c3)
            if 0.0 < t_flex < h:
                cand_dn.append(float(t_flex))
        ts_dn = np.asarray(sorted(set(cand_dn)), dtype=float)
        n_vals = np.polyval(poly_n, ts_dn)
        dn_vals = np.polyval(poly_dn, ts_dn)
        g_vals = (mu0 + ts_dn) * dn_vals - n_vals

        for t_val, dnv, gv in zip(ts_dn, dn_vals, g_vals):
            if dnv < min_dn:
                min_dn = float(dnv)
                worst_dn_mub = float(mu0 + t_val)
            if dnv > max_dn:
                max_dn = float(dnv)
            if gv < min_margin:
                min_margin = float(gv)
                worst_margin_mub = float(mu0 + t_val)
            if gv > max_margin:
                max_margin = float(gv)

        if np.min(dn_vals) <= 0.0:
            raise ValueError(
                f"Gibbs {stage} violates thermodynamic stability dn_B/dmu_B > 0 "
                f"(min dn_B/dmu_B={min_dn:.3e} at mu_B={worst_dn_mub:.6g} MeV)"
            )

        # Stationary points of cs2 = n(t) / ((mu0 + t) * n'(t)):
        # M(t) = (mu0+t)*(n')^2 - n*(n' + (mu0+t)*n'') = 0 (degree <= 5).
        poly_M = np.polysub(
            np.polymul(poly_mu, np.polymul(poly_dn, poly_dn)),
            np.polymul(poly_n, np.polyadd(poly_dn, np.polymul(poly_mu, poly_ddn))),
        )
        cand_cs2 = list(ts_dn)
        m_scale = float(np.max(np.abs(poly_M))) if poly_M.size else 0.0
        if m_scale > 0.0:
            nonzero_idx = np.flatnonzero(np.abs(poly_M) > 1e-14 * m_scale)
            if nonzero_idx.size > 0 and nonzero_idx[0] < len(poly_M) - 1:
                trimmed = poly_M[nonzero_idx[0]:]
                for r in np.roots(trimmed):
                    if abs(r.imag) <= 1e-10 * max(1.0, h) and np.isfinite(r.real):
                        tr = float(r.real)
                        if 0.0 < tr < h:
                            cand_cs2.append(tr)
        ts_cs2 = np.asarray(sorted(set(cand_cs2)), dtype=float)
        n_c = np.polyval(poly_n, ts_cs2)
        dn_c = np.polyval(poly_dn, ts_cs2)
        if np.any(dn_c <= 0.0):
            raise ValueError(
                f"Gibbs {stage} encounters non-positive dn_B/dmu_B inside interval [{mu0:.6g}, {mu0+h:.6g}]"
            )
        cs2_c = n_c / ((mu0 + ts_cs2) * dn_c)
        for t_val, cval in zip(ts_cs2, cs2_c):
            if cval < min_cs2:
                min_cs2 = float(cval)
            if cval > max_cs2:
                max_cs2 = float(cval)
                worst_cs2_mub = float(mu0 + t_val)

    margin_scale = causality_margin_tol * max(1.0, float(np.max(nb)))
    if min_dn <= 0.0:
        raise ValueError(
            f"Gibbs {stage} violates dn_B/dmu_B > 0 (min={min_dn:.3e} at mu_B={worst_dn_mub:.6g} MeV)"
        )
    if min_margin < -margin_scale:
        raise ValueError(
            f"Gibbs {stage} violates causality margin g = mu_B*dn_B/dmu_B - n_B >= 0 "
            f"(min_g={min_margin:.3e} at mu_B={worst_margin_mub:.6g} MeV, max_cs2={max_cs2:.6g})"
        )
    if min_cs2 <= 0.0 or max_cs2 > 1.0 + cs2_tol:
        raise ValueError(
            f"Gibbs {stage} violates causal sound speed 0 < cs2 <= 1 "
            f"(range=[{min_cs2:.6g}, {max_cs2:.6g}], max at mu_B={worst_cs2_mub:.6g} MeV)"
        )

    return {
        "stage": stage,
        "representation": "PchipInterpolator(mu_b_MeV, nb_fm3)",
        "interval_count": int(n_intervals),
        "euler_residual_max_abs_MeV_fm3": euler_max_abs,
        "euler_residual_max_rel": euler_max_rel,
        "euler_residual_rms_rel": euler_rms_rel,
        "pressure_reconstruction_max_abs_MeV_fm3": p_max_abs,
        "pressure_reconstruction_rms_abs_MeV_fm3": p_rms_abs,
        "pressure_reconstruction_max_rel": p_max_rel,
        "pressure_reconstruction_rms_rel": p_rms_rel,
        "min_dnb_dmu_b_fm3_per_MeV": min_dn,
        "max_dnb_dmu_b_fm3_per_MeV": max_dn,
        "worst_dnb_dmu_b_mu_b_MeV": worst_dn_mub,
        "min_causality_margin_fm3": min_margin,
        "max_causality_margin_fm3": max_margin,
        "worst_causality_margin_mu_b_MeV": worst_margin_mub,
        "min_cs2": min_cs2,
        "max_cs2": max_cs2,
        "worst_cs2_mu_b_MeV": worst_cs2_mub,
        "tolerances": {
            "euler_rtol": float(euler_rtol),
            "euler_atol": float(euler_atol),
            "pressure_rtol": float(pressure_rtol),
            "pressure_atol": float(pressure_atol),
            "causality_margin_tol_scaled": float(margin_scale),
            "cs2_tol": float(cs2_tol),
        },
        "pass_": True,
    }


def gibbs_mixed_phase(hadronic: ChargeResolvedPhase, high_density: ChargeResolvedPhase,
                      mu_b_values, *, lepton_sector: Optional[LeptonSector]=None,
                      mu_q_bounds=None) -> GibbsResult:
    """Construct a globally neutral Gibbs mixed phase from two-charge surfaces.

    Requires P_i(mu_B,mu_Q), including charge-resolved baryon/charge densities.
    A one-dimensional beta-equilibrated CompOSE barotrope is insufficient by design.
    """
    if not isinstance(hadronic, ChargeResolvedPhase) or not isinstance(high_density, ChargeResolvedPhase):
        raise TypeError('Gibbs construction requires two ChargeResolvedPhase P(mu_B,mu_Q) surfaces; 1-D EOSData is insufficient')
    for ph in (hadronic, high_density):
        if (ph.mu_b_units, ph.mu_q_units, ph.pressure_units) != ('MeV', 'MeV', 'MeV fm^-3'):
            raise ValueError(
                f'{ph.label}: unsupported units ({ph.mu_b_units}, {ph.mu_q_units}, {ph.pressure_units}); '
                'expected (MeV, MeV, MeV fm^-3)'
            )
    lb = max(hadronic.mu_b[0], high_density.mu_b[0])
    ub = min(hadronic.mu_b[-1], high_density.mu_b[-1])
    mus = np.asarray(mu_b_values, float)
    if mus.ndim != 1 or len(mus) < 4 or np.any(np.diff(mus) <= 0) or mus[0] < lb or mus[-1] > ub:
        raise ValueError('mu_b_values must be an increasing grid inside both phase surfaces')
    qlo = max(hadronic.mu_q[0], high_density.mu_q[0])
    qhi = min(hadronic.mu_q[-1], high_density.mu_q[-1])
    if mu_q_bounds is not None:
        qlo = max(qlo, float(mu_q_bounds[0]))
        qhi = min(qhi, float(mu_q_bounds[1]))
    if qlo >= qhi:
        raise ValueError('two phases have no common charge-chemical-potential interval')
    leptons = lepton_sector or LeptonSector.electron_muon()
    out = []
    pressure_equal_residuals = []
    for mub in mus:
        roots = _gibbs_pressure_equality_roots(hadronic, high_density, float(mub), qlo, qhi)
        if len(roots) != 1:
            raise ValueError(
                f'Gibbs pressure equality has {len(roots)} roots at mu_B={mub:g}; '
                'a unique coexistence branch is required'
            )
        qmu = float(roots[0])
        ph, nh, qh = hadronic.evaluate(mub, qmu)
        pq, nq, qq = high_density.evaluate(mub, qmu)
        pressure_equal_residuals.append(abs(ph - pq))
        ql = float(leptons.charge_density(qmu))
        pl = float(leptons.pressure(qmu))
        for phase_obj in (hadronic, high_density):
            if str(phase_obj.lepton_contributions).lower().startswith('included') and (abs(ql) > 1e-12 or abs(pl) > 1e-12):
                raise ValueError(
                    f'{phase_obj.label} already declares lepton_contributions={phase_obj.lepton_contributions!r}; '
                    'refusing to add an external LeptonSector on top'
                )
        denom = qq - qh
        denom_tol = max(1e-10, 1e-10 * max(1.0, abs(qh), abs(qq), abs(ql)))
        if abs(denom) <= denom_tol:
            raise ValueError('phase charge densities are indistinguishable; mixed fraction is undefined')
        chi = float(-(qh + ql) / denom)
        if chi < -1e-8 or chi > 1.0 + 1e-8:
            raise ValueError(f'global neutrality requires mixed fraction chi={chi:g} outside [0,1] at mu_B={mub:g}')
        chi = float(np.clip(chi, 0.0, 1.0))
        n = (1.0 - chi) * nh + chi * nq
        P = 0.5 * (ph + pq) + pl
        charge = (1.0 - chi) * qh + chi * qq + ql
        eps = mub * n - P
        if n <= 0.0 or P <= 0.0 or eps <= 0.0 or abs(charge) > 2e-7:
            raise ValueError('Gibbs point violates positivity or global charge neutrality')
        out.append((n, P, eps, mub, qmu, chi, charge))
    a = np.asarray(out)
    n, P, eps, mub, qmu, chi, qres = a.T
    if np.any(np.diff(n) <= 0) or np.any(np.diff(P) <= 0) or np.any(np.diff(eps) <= 0):
        raise ValueError('Gibbs mixed phase is not monotone/stable on requested mu_B grid')

    post_construction_val = validate_gibbs_barotrope(
        mub, n, P, eps, stage="post_construction"
    )

    mn = float(939.0)
    F = eps / n
    eos = EOSData(
        nb=n, mn=mn, mp=mn, Q1=P / n, Q2=np.zeros_like(n), Q3=mub / mn - 1.0,
        Q4=qmu / mn, Q5=np.zeros_like(n), Q6=F / mn - 1.0, Q7=F / mn - 1.0,
        q_extra=np.zeros_like(n), Iphase=np.full(len(n), 3, dtype=int),
        Yn=np.zeros_like(n), Yp=np.zeros_like(n), Ye=np.zeros_like(n), Ymu=np.zeros_like(n),
        Aav=np.zeros_like(n), Zav=np.zeros_like(n), Yh=np.zeros_like(n),
        Il=1, iT=0, iYq=1, inb_start=1,
        temperature_mev=0.0, charge_fraction_grid=0.0, q_extra_count=0,
        composition_available=False, beta_equilibrated=True,
    )
    meta = {
        'schema': 'compose-eos-gibbs-v2',
        'hadronic_phase': hadronic.label,
        'high_density_phase': high_density.label,
        'lepton_sector': leptons.label,
        'composition_status': 'thermo_only',
        'composition_note': (
            'Gibbs output provides equilibrium thermodynamics only (composition_available=False, '
            'Iphase=[3] in sidecar) until explicit species-resolved mixed-phase composition is implemented.'
        ),
        'normalization_masses_mev': {
            'mn': mn,
            'mp': mn,
            'role': 'CompOSE normalization reference only; does not represent constituent masses of either underlying phase',
        },
        'provenance': {
            'hadronic': hadronic.provenance_dict(),
            'high_density': high_density.provenance_dict(),
            'lepton_sector': {'label': leptons.label},
        },
        'mu_b_MeV': mub.tolist(),
        'mu_q_MeV': qmu.tolist(),
        'mixed_fraction_chi': chi.tolist(),
        'global_charge_residual': qres.tolist(),
        'pressure_equality_residual_max_MeV_fm3': float(max(pressure_equal_residuals)),
        'validation': {
            'post_construction': post_construction_val,
            'post_readback': None,
        },
        'construction': 'two conserved charges; pressure equality plus global charge neutrality',
        'input_requirement': 'P_H(mu_B,mu_Q), P_Q(mu_B,mu_Q), derivatives n_B=partial P/partial mu_B and q=partial P/partial mu_Q',
    }
    return GibbsResult(eos, meta)


def write_gibbs_eos(directory: str, result: GibbsResult) -> dict:
    """Write CompOSE arrays + Gibbs provenance sidecar via staged native readback."""
    src_hashes = {}
    prov = result.metadata.get('provenance', {})
    for role in ('hadronic', 'high_density'):
        info = prov.get(role, {})
        if info.get('source_sha256'):
            src_hashes[role] = info['source_sha256']
    payload = {
        **result.metadata,
        'phase_ids': result.eos.Iphase.astype(int).tolist(),
    }
    return _stage_and_promote_compose_with_sidecar(
        directory,
        result.eos,
        'eos.gibbs_construction.json',
        payload,
        read_gibbs_eos,
        kind='gibbs',
        source_hashes=src_hashes or None,
    )


def read_gibbs_eos(directory: str) -> GibbsResult:
    data = read_compose(directory)
    path = Path(directory) / 'eos.gibbs_construction.json'
    if not path.is_file():
        raise FileNotFoundError('EOS directory has no eos.gibbs_construction.json provenance sidecar')
    meta = json.loads(path.read_text(encoding='utf-8'))
    if meta.get('schema') not in ('compose-eos-gibbs-v1', 'compose-eos-gibbs-v2'):
        raise ValueError('unsupported Gibbs sidecar schema')
    fields = ('mu_b_MeV', 'mu_q_MeV', 'mixed_fraction_chi', 'global_charge_residual')
    if any(len(meta.get(k, [])) != data.n for k in fields):
        raise ValueError('Gibbs chemical-potential/mixture arrays do not match EOS row count')
    arrays = {k: np.asarray(meta[k], dtype=float) for k in fields}
    if any(v.shape != (data.n,) or np.any(~np.isfinite(v)) for v in arrays.values()):
        raise ValueError('Gibbs sidecar arrays must be finite and match EOS row count')
    if np.max(np.abs(arrays['global_charge_residual'])) > 2e-7:
        raise ValueError('Gibbs sidecar violates global charge neutrality')
    if np.any(arrays['mixed_fraction_chi'] < -1e-8) or np.any(arrays['mixed_fraction_chi'] > 1.0 + 1e-8):
        raise ValueError('Gibbs mixed fraction lies outside [0,1]')
    for label, actual, expected in (
        ('mu_B', arrays['mu_b_MeV'], data.mub),
        ('mu_Q', arrays['mu_q_MeV'], data.muq),
    ):
        residual = np.abs(actual - np.asarray(expected, dtype=float))
        scale = np.maximum(np.maximum(np.abs(actual), np.abs(expected)), 1.0)
        if np.any(residual > 2e-8 * scale):
            raise ValueError(f'Gibbs sidecar {label} does not match CompOSE table chemical potential')
    if 'phase_ids' in meta:
        ids = np.asarray(meta['phase_ids'], dtype=int)
        if ids.shape != (data.n,):
            raise ValueError('Gibbs sidecar phase_ids length mismatch')
        data.Iphase = ids
    post_readback_val = validate_gibbs_barotrope(
        data.mub, data.nb, data.p, data.eps, stage="post_readback"
    )
    meta = {
        **meta,
        'validation': {
            **(meta.get('validation') or {}),
            'post_readback': post_readback_val,
        },
    }
    return GibbsResult(data, meta)


def maxwell_tidal_jump_report(eos: MaxwellEOS, radius_km, mass_geometric_km):
    t=eos.transition; return {'transition':dict(t),'radius_km':float(radius_km),'mass_geometric_km':float(mass_geometric_km),
        'delta_y_for_y_inside_1':tidal_y_jump(1.0,radius_km,mass_geometric_km,
          t['pressure_MeV_fm3']*_MEVFM3_TO_KM2,t['delta_eps_MeV_fm3']*_MEVFM3_TO_KM2)-1.0}


def build_maxwell_css_eos(hadronic: EOSData, pressure_transition_MeV_fm3: float,
                          delta_epsilon_MeV_fm3: float, cs2_high: float, *,
                          n_max_fm3: float, n_points: int=240,
                          overlap_fraction: float=0.01, hadronic_name='hadronic',
                          high_phase_name='CSS-high-density') -> MaxwellEOS:
    """Create an exact Maxwell transition with a phenomenological CSS phase.

    The hadronic phase is taken from the supplied cold beta-equilibrated table
    below P_t. The high-density phase obeys dP/dε=cs2_high and is fixed by the
    specified energy-density jump Δε at equal P_t and μ_B. A tiny metastable
    overlap brackets the chemical-potential crossing and is removed from the
    final table. This is a generic CSS parameterization, not a microscopic QCD EOS.
    """
    from scipy.interpolate import PchipInterpolator
    n=np.asarray(hadronic.nb,float); p=np.asarray(hadronic.p,float); eps=np.asarray(hadronic.eps,float)
    _phase_arrays(hadronic)
    pt=float(pressure_transition_MeV_fm3); de=float(delta_epsilon_MeV_fm3); cs=float(cs2_high); nmax=float(n_max_fm3)
    if not hadronic.beta_equilibrated: raise ValueError('CSS Maxwell builder requires a beta-equilibrated one-dimensional hadronic EOS')
    if not np.all(np.isfinite([pt,de,cs,nmax,overlap_fraction])) or not p[0]<pt<p[-1] or de<=0 or not 0<cs<=1:
        raise ValueError('invalid transition pressure, energy jump or high-phase sound speed')
    if n_points<40 or not 0<overlap_fraction<0.1: raise ValueError('n_points>=40 and overlap_fraction in (0,0.1) required')
    x=np.log(n); lp=PchipInterpolator(x,np.log(p)); le=PchipInterpolator(x,np.log(eps))
    lognt=brentq(lambda z: float(lp(z)-np.log(pt)),float(x[0]),float(x[-1]),xtol=1e-13)
    nt=float(np.exp(lognt)); et=float(np.exp(le(lognt))); mut=float((et+pt)/nt)
    mu_source=float(PchipInterpolator(x,np.asarray(hadronic.mub,float))(lognt))
    mu_residual=abs(mut-mu_source)/max(abs(mut),1.0)
    if mu_residual>5e-3: raise ValueError(f'hadronic phase violates mu_B=(epsilon+P)/n_B at P_t by {mu_residual:.3e}')
    # Limit the metastable overlap to keep the CSS pressure positive below P_t.
    epsq=et+de; nq=(epsq+pt)/mut
    alpha=min(float(overlap_fraction),0.2*pt/(cs*(epsq+pt)))
    if alpha<1e-7: raise ValueError('P_t is too small to resolve a positive-pressure CSS metastable overlap')
    ncut=nt*(1+alpha)
    if ncut>=n[-1] or nmax<=nq: raise ValueError('hadronic table does not extend beyond P_t overlap or n_max does not exceed n_Q')
    nlow0=n[n<ncut]
    if len(nlow0)<5: raise ValueError('insufficient hadronic points below the transition')
    plow0=p[n<ncut]; elow0=eps[n<ncut]
    nlow=np.r_[nlow0,ncut]; plow=np.r_[plow0,float(np.exp(lp(np.log(ncut))))]; elow=np.r_[elow0,float(np.exp(le(np.log(ncut))))]
    mulow=(elow+plow)/nlow
    q4sp=PchipInterpolator(x,np.asarray(hadronic.muq,float)); q5sp=PchipInterpolator(x,np.asarray(hadronic.mul,float))
    q4low=q4sp(np.log(nlow))/float(hadronic.mn); q5low=q5sp(np.log(nlow))/float(hadronic.mn)
    def make_phase(nn,pp,ee,mm,q4,q5):
        z=np.zeros_like(nn); F=ee/nn; mn=float(hadronic.mn)
        return EOSData(nb=nn,mn=mn,mp=float(hadronic.mp),Q1=pp/nn,Q2=z,Q3=mm/mn-1,
          Q4=q4,Q5=q5,Q6=F/mn-1,Q7=F/mn-1,q_extra=z,Iphase=np.zeros(len(nn),int),
          Yn=z,Yp=z,Ye=z,Ymu=z,Aav=z,Zav=z,Yh=z,Il=int(hadronic.Il),iT=0,iYq=1,inb_start=1,
          temperature_mev=0,charge_fraction_grid=0,q_extra_count=0,composition_available=False,beta_equilibrated=True)
    low=make_phase(nlow,plow,elow,mulow,q4low,q5low)
    nmin=nq/(1+alpha); npre=np.geomspace(nmin,nq,max(8,n_points//3)); npost=np.geomspace(nq,nmax,max(8,n_points-len(npre)+1))[1:]
    nh=np.r_[npre,npost]; ratio=nh/nq; muh=mut*ratio**cs
    ph=pt+cs/(1+cs)*(epsq+pt)*(ratio**(1+cs)-1); eh=muh*nh-ph
    if np.any(ph<=0) or np.any(eh<=0): raise ValueError('CSS high phase became non-positive')
    q4h=np.full(len(nh),float(q4sp(lognt))/float(hadronic.mn)); q5h=np.full(len(nh),float(q5sp(lognt))/float(hadronic.mn))
    high=make_phase(nh,ph,eh,muh,q4h,q5h)
    result=build_maxwell_eos(low,high,low_name=hadronic_name,high_name=high_phase_name)
    result.metadata.update({'high_phase_model':'constant_sound_speed','cs2_high':cs,
      'input_transition_pressure_MeV_fm3':pt,'input_delta_epsilon_MeV_fm3':de,
      'input_epsilon_hadronic_transition_MeV_fm3':et,'input_n_hadronic_transition_fm3':nt,
      'input_muB_transition_MeV':mut,'n_high_phase_transition_fm3':nq,
      'metastable_overlap_fraction':alpha,'hadronic_mu_match_relative_residual':mu_residual,
      'charge_potential_extension':'mu_Q and effective lepton chemical potential are held at the transition values as a one-fluid placeholder; no high-phase microscopic composition is claimed.',
      'interpretation':'Phenomenological CSS high-density phase. Its transition pressure, energy jump and sound speed are model parameters, not a microscopic quark-matter calculation.'})
    return result
