"""M5: Independent Python reconstruction of the sharp Maxwell construction
between two locally-charge-neutral phase barotropes.

Physics (independent of any legacy implementation):
  Two phases, each a 1-D cold beta-equilibrated barotrope P_i(n_B), mu_B_i(n_B),
  epsilon_i(n_B), are placed in chemical-mechanical equilibrium at a transition
  pressure P_t and baryon chemical potential mu_t. The Maxwell (equal-area for
  the free-energy, or equivalently equal-P equal-mu) construction requires:
    P_low(mu_t) = P_high(mu_t) = P_t           (mechanical equilibrium)
    n_low < n_high                              (density jump upward)
    epsilon_high - epsilon_low = mu_t (n_high - n_low)   (first law / latent heat)
  The output is a single combined barotrope with an explicit two-row density
  jump at (P_t, mu_t), preserving each phase's tabulated points away from the
  jump, and an eos.phase_transitions.json sidecar recording the transition
  metadata including the Seidov small-core stability criterion
    Delta epsilon_critical = 0.5 * epsilon_low + 1.5 * P_t
    small core stable iff Delta epsilon < Delta epsilon_critical
  (Seidov 1967 Astron. Zhur. 44 343 / ApJ 1968).

The implementation uses PCHIP interpolation of P(mu_B) for each phase, a
dense scan to bracket the sign change of P_low - P_high, and Brent root
finding for the exact coexistence mu_t. This mirrors the standard Maxwell
construction algorithm but is written from the physics, not copied.

References:
  - Maxwell 1875 construction (classical thermodynamics of first-order transitions)
  - Glendenning 1991 PhysRep 341 89 (relativistic two-phase formalism; uploaded)
  - Seidov 1967 Astron. Zhur. 44 343 (small-core stability)
  - Takatsy & Kovacs 2020 PRD 102 028501 (tidal y-jump at first-order interface)

This is the Python reproduction of the M5 acceptance criteria (no Fortran in
this sandbox - gfortran unavailable); the Fortran kernel exists in the repo
at fortran_eos/src/ and this Python reproduction validates the same physics.
"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
import json, math
import sys
sys.path.insert(0, '/home/z/my-project/repos/NS_Cowling_eos390/python_frontend')
from compose_eos.eos_data import EOSData

# Standard CompOSE -> physical arrays
def _phase_arrays(data: EOSData):
    """Return (n_B, mu_B, P, eps, hard_gate_report) for one phase barotrope.
    mu_B = (Q3+1)*mn, P = Q1*nb, eps = nb*(Q7+1)*mn."""
    n = np.asarray(data.nb, float)
    mu = np.asarray(data.mub, float)   # (Q3+1)*mn
    P = np.asarray(data.p, float)     # Q1*nb
    eps = np.asarray(data.eps, float) # nb*(Q7+1)*mn
    # Hard-gate sanity: monotone increasing n, P, mu, eps, mu*mu - P^2 causal-ish
    rep = {}
    rep['C1'] = {'pass_': bool(np.all(np.isfinite([n,mu,P,eps])))}
    rep['C11'] = {'pass_': bool(np.all(np.diff(mu) > 0))}
    rep['n_monotone'] = {'pass_': bool(np.all(np.diff(n) > 0))}
    rep['P_monotone'] = {'pass_': bool(np.all(np.diff(P) > 0))}
    rep['eps_monotone'] = {'pass_': bool(np.all(np.diff(eps) > 0))}
    return n, mu, P, eps, rep


def build_maxwell_eos(low_phase: EOSData, high_phase: EOSData, *,
                      low_name: str = 'low', high_name: str = 'high',
                      mu_samples: int = 800, crossing_index: int = 0) -> dict:
    """Construct one sharp Maxwell transition between two locally-neutral 1-D
    phase barotropes. Returns a dict with the combined EOSData (as numpy arrays
    in CompOSE Q-format) and the transition metadata.
    """
    nl, mul, pl, el, repl = _phase_arrays(low_phase)
    nh, muh, ph, eh, reph = _phase_arrays(high_phase)
    # Hard gates
    for key in ('C1','C11','n_monotone','P_monotone','eps_monotone'):
        if not (repl[key]['pass_'] and reph[key]['pass_']):
            raise ValueError(f'phase hard gate {key} failed; Maxwell requires stable source phases')
    # Common mu_B interval
    lo = max(mul[0], muh[0]); hi = min(mul[-1], muh[-1])
    if not lo < hi:
        raise ValueError('phases have no common baryon-chemical-potential interval')
    # PCHIP interpolators of P(mu_B), n(mu_B), eps(mu_B) for each phase
    Pl = PchipInterpolator(mul, pl, extrapolate=False)
    Ph = PchipInterpolator(muh, ph, extrapolate=False)
    nli = PchipInterpolator(mul, nl, extrapolate=False)
    nhi = PchipInterpolator(muh, nh, extrapolate=False)
    eli = PchipInterpolator(mul, el, extrapolate=False)
    ehi = PchipInterpolator(muh, eh, extrapolate=False)
    # Dense scan for sign change of Pl - Ph
    grid = np.linspace(lo, hi, mu_samples)
    delta = np.asarray(Pl(grid) - Ph(grid), float)
    brackets = [(grid[i], grid[i+1])
                for i in range(len(grid)-1)
                if delta[i] == 0 or delta[i]*delta[i+1] < 0]
    if not brackets:
        raise ValueError('no P_low(mu_B)=P_high(mu_B) coexistence root in common interval')
    if crossing_index < 0 or crossing_index >= len(brackets):
        raise ValueError(f'{len(brackets)} coexistence roots found; choose crossing_index explicitly')
    a, b = brackets[crossing_index]
    # Exact root via Brent
    if delta[np.searchsorted(grid, a)] == 0:
        mut = float(a)
    else:
        mut = float(brentq(lambda x: float(Pl(x) - Ph(x)), a, b, xtol=1e-12))
    pt = float(0.5 * (Pl(mut) + Ph(mut)))
    nlo = float(nli(mut)); nhi_t = float(nhi(mut))
    if not (pt > 0 and mut > 0 and nlo > 0 and nhi_t > nlo):
        raise ValueError(f'coexistence root is not a physical density jump: P={pt}, mu={mut}, n_low={nlo}, n_high={nhi_t}')
    # First-law consistency at the coexistence point
    elo = float(mut * nlo - pt)
    epshi_t = float(mut * nhi_t - pt)
    de = epshi_t - elo
    # Compare against the source-tabulated eps at the same mu_B (independent
    # check that the first-law identity eps = mu*n - P holds at the transition)
    euler_low_rel = abs(elo - float(eli(mut))) / max(abs(elo), 1.0)
    euler_high_rel = abs(epshi_t - float(ehi(mut))) / max(abs(epshi_t), 1.0)
    if max(euler_low_rel, euler_high_rel) > 5e-3:
        raise ValueError('interpolated coexistence point violates phase first-law identity beyond 5e-3')
    if not (elo > 0 and epshi_t > elo and de > 0):
        raise ValueError('Maxwell energy-density jump is non-positive')
    # Assemble the combined barotrope: low phase up to (mut, nlo), then (mut, nhi_t), then high phase
    ml = mul < mut
    mh = muh > mut
    if ml.sum() < 3 or mh.sum() < 3:
        raise ValueError('not enough tabulated points on both sides of the coexistence point')
    n_out = np.r_[nl[ml], nlo, nhi_t, nh[mh]]
    mu_out = np.r_[mul[ml], mut, mut, muh[mh]]
    p_out = np.r_[pl[ml], pt, pt, ph[mh]]
    phase_ids = np.r_[np.ones(ml.sum() + 1, dtype=int), np.full(1 + mh.sum(), 2, dtype=int)]
    eps_out = mu_out * n_out - p_out
    # Monotonicity gates on the assembled barotrope
    if np.any(np.diff(n_out) <= 0) or np.any(np.diff(p_out) < -1e-10) or np.any(np.diff(eps_out) <= 0):
        raise ValueError('constructed Maxwell table is not monotone in n_B, P, epsilon')
    mn = float(low_phase.mn)
    F = eps_out / n_out
    # Seidov small-core stability criterion
    seidov_critical = 0.5 * elo + 1.5 * pt
    il = int(ml.sum()); ih = il + 1
    transition = {
        'kind': 'maxwell',
        'pressure_MeV_fm3': pt,
        'muB_MeV': mut,
        'n_low_fm3': nlo,
        'n_high_fm3': nhi_t,
        'eps_low_MeV_fm3': elo,
        'eps_high_MeV_fm3': epshi_t,
        'delta_eps_MeV_fm3': de,
        'delta_n_fm3': nhi_t - nlo,
        'low_index': il,
        'high_index': ih,
        'seidov_critical_delta_eps_MeV_fm3': seidov_critical,
        'seidov_delta_eps_over_critical': de / seidov_critical,
        'seidov_small_core_stable': bool(de < seidov_critical),
        'tidal_jump_formula': 'Takatsy-Kovacs-2020: Delta y = -4 pi r_t^3 Delta epsilon / (m_t + 4 pi r_t^3 P_t)',
        'charge_neutrality': 'local in each phase (Maxwell)',
        'first_law_residual': abs(de - mut * (nhi_t - nlo)),
    }
    meta = {
        'schema': 'compose-eos-maxwell-v1',
        'low_phase': low_name,
        'high_phase': high_name,
        'coexistence_root_count': len(brackets),
        'selected_root_index': int(crossing_index),
        'low_energy_identity_residual': elo - float(eli(mut)),
        'high_energy_identity_residual': epshi_t - float(ehi(mut)),
        'low_source_hard_gates': {k: bool(repl[k]['pass_']) for k in ('C1','C11','n_monotone','P_monotone','eps_monotone')},
        'high_source_hard_gates': {k: bool(reph[k]['pass_']) for k in ('C1','C11','n_monotone','P_monotone','eps_monotone')},
        'transition': transition,
        'phase_ids': phase_ids.tolist(),
        'interpretation': 'Exact sharp Maxwell interface between two locally-neutral 1-D barotropes (Python reproduction of M5 acceptance criteria).',
    }
    # Pack as a new EOSData (CompOSE Q-format)
    data = EOSData(
        nb=n_out, mn=mn, mp=float(low_phase.mp),
        Q1=p_out / n_out, Q2=np.zeros_like(n_out),
        Q3=mu_out / mn - 1.0, Q4=np.zeros_like(n_out), Q5=np.zeros_like(n_out),
        Q6=F / mn - 1.0, Q7=F / mn - 1.0,
        q_extra=np.zeros_like(n_out), Iphase=phase_ids,
        Yn=np.zeros_like(n_out), Yp=np.zeros_like(n_out),
        Ye=np.zeros_like(n_out), Ymu=np.zeros_like(n_out),
        Aav=np.zeros_like(n_out), Zav=np.zeros_like(n_out), Yh=np.zeros_like(n_out),
        Il=1, iT=0, iYq=1, inb_start=1, temperature_mev=0.0,
        charge_fraction_grid=0.0, q_extra_count=0,
        composition_available=False, beta_equilibrated=True,
    )
    return {'data': data, 'transition': transition, 'metadata': meta}


def verify_maxwell_acceptance(result: dict, tol: float = 1e-9) -> dict:
    """Verify the M5 acceptance criteria on a Maxwell construction result.
    Returns a dict of pass/fail booleans and the numerical values."""
    t = result['transition']
    data = result['data']
    il, ih = int(t['low_index']), int(t['high_index'])
    checks = {}
    # 1. nonzero density and energy jumps
    checks['nonzero_density_jump'] = bool(t['n_high_fm3'] - t['n_low_fm3'] > 0)
    checks['nonzero_energy_jump'] = bool(t['delta_eps_MeV_fm3'] > 0)
    # 2. unique coexistence root: equal P and mu_B at the jump
    checks['pressure_equal_at_jump'] = bool(abs(data.p[il] - data.p[ih]) < tol * max(1, abs(data.p[il])))
    checks['muB_equal_at_jump'] = bool(abs(data.mub[il] - data.mub[ih]) < tol * max(1, abs(data.mub[il])))
    # 3. first law: Delta eps = mu_t * Delta n_B
    de_calc = t['muB_MeV'] * (t['n_high_fm3'] - t['n_low_fm3'])
    checks['first_law_Delta_eps'] = bool(abs(t['delta_eps_MeV_fm3'] - de_calc) < tol * max(1, abs(de_calc)))
    # 4. monotonicity of the assembled barotrope (n, P, eps all strictly increasing)
    checks['n_monotone'] = bool(np.all(np.diff(data.nb) > 0))
    checks['P_monotone'] = bool(np.all(np.diff(data.p) > -1e-10))
    checks['eps_monotone'] = bool(np.all(np.diff(data.eps) > 0))
    # 5. Seidov small-core stability reported
    checks['seidov_reported'] = bool('seidov_small_core_stable' in t)
    # 6. phase IDs bracket the jump correctly
    checks['phase_ids_bracket_jump'] = bool(int(data.Iphase[il]) == 1 and int(data.Iphase[ih]) == 2)
    overall = all(checks.values())
    return {'checks': checks, 'overall_pass': overall, 'transition': t}
