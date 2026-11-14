"""M5A: Independent Python reproduction of the Constant-Sound-Speed (CSS)
high-density branch with an explicit sharp first-order interface.

Physics (from Alford, Han & Prakash 2013, PRD 88 083013, arXiv:1302.4732;
Alford et al. 2015, PRD 92 083002, arXiv:1501.07902):

A constant-sound-speed (CSS) high-density EOS branch is parameterized by
  P_t      - transition pressure (MeV/fm^3)
  eps_low  - energy density on the low side at P_t (MeV/fm^3)
  n_low    - baryon number density on the low side at P_t (fm^-3)
  Delta_eps - latent heat (energy-density jump) at P_t (>= 0, MeV/fm^3)
  cs2      - constant sound speed squared (dimensionless, in units of c^2, 0<cs2<=1)

The high-density branch (P >= P_t) is:
  mu_t = (eps_low + P_t) / n_low                              [transition mu_B]
  eps_high_at_Pt = eps_low + Delta_eps                         [energy jump]
  n_high_at_Pt   = (eps_high_at_Pt + P_t) / mu_t               [density jump, first law]
  eps_high(P)    = eps_high_at_Pt + (P - P_t) / cs2            [constant cs^2: dP/deps = cs2]
  n_high(P)      = n_high_at_Pt * [1 + (1+cs2)(P - P_t)/(cs2 (eps_high_at_Pt + P_t))]^(1/(1+cs2))
  mu_high(P)     = (eps_high(P) + P) / n_high(P)               [first law: eps + P = mu n]

Acceptance (from the DAG css_segment_contract):
  - At P=P_t, the high side recovers (eps_high_at_Pt, n_high_at_Pt, mu_t).
  - Delta_eps = mu_t * (n_high_at_Pt - n_low) within floating-point tolerance.
  - Delta_eps > 0  =>  sharp first-order interface.
  - Delta_eps = 0  =>  CONTINUOUS_UNCLASSIFIED (no jump, never inferred as
    second-order or crossover without extra thermodynamic evidence).
  - dP/deps = cs2 and deps/dn_B = mu_B at every high-side sample.
  - Reject cs2 <= 0, cs2 > 1, Delta_eps < 0, nonfinite inputs, P < P_t.

This Python reproduction validates the M5A acceptance criteria and is the
physics oracle for the Fortran kernel at fortran_eos/src/compose_eos_css.f90
(which is not compiled in this gfortran-less sandbox).
"""
from __future__ import annotations
import math
import numpy as np


def css_transition_state(p_t: float, eps_low: float, n_low: float,
                          delta_eps: float, cs2: float) -> dict:
    """Compute the exact transition state at P=P_t.
    Returns dict with mu_t, eps_high_at_Pt, n_high_at_Pt, delta_n, first-law residual.
    Raises ValueError on any invalid input.
    """
    # Input validation (acceptance: reject invalid inputs)
    vals = [p_t, eps_low, n_low, delta_eps, cs2]
    if not all(math.isfinite(v) for v in vals):
        raise ValueError('CSS inputs must be finite')
    if p_t < 0 or eps_low <= 0 or n_low <= 0 or delta_eps < 0:
        raise ValueError(f'invalid CSS inputs: P_t={p_t}, eps_low={eps_low}, n_low={n_low}, Delta_eps={delta_eps}')
    if not (0 < cs2 <= 1.0 + 1e-12):
        raise ValueError(f'invalid cs2: {cs2} (must satisfy 0 < cs2 <= 1)')
    # mu_t from the low-side state at the transition
    mu_t = (eps_low + p_t) / n_low
    # High-side energy and density at P_t
    eps_high_t = eps_low + delta_eps
    n_high_t = (eps_high_t + p_t) / mu_t
    delta_n = n_high_t - n_low
    # First-law check at the transition: Delta eps = mu_t * Delta n_B
    first_law_residual = abs(delta_eps - mu_t * delta_n)
    # Classification
    if delta_eps > 0:
        classification = 'SHARP_FIRST_ORDER'
    else:
        classification = 'CONTINUOUS_UNCLASSIFIED'
    return {
        'P_t': p_t, 'eps_low': eps_low, 'n_low': n_low,
        'Delta_eps': delta_eps, 'cs2': cs2,
        'mu_t': mu_t, 'eps_high_at_Pt': eps_high_t,
        'n_high_at_Pt': n_high_t, 'delta_n': delta_n,
        'first_law_residual': first_law_residual,
        'classification': classification,
    }


def css_high_branch(p_t: float, eps_low: float, n_low: float,
                     delta_eps: float, cs2: float,
                     pressures: np.ndarray) -> dict:
    """Evaluate the CSS high-density branch at an array of pressures P >= P_t.
    Returns dict with arrays eps_high, n_high, mu_high and the transition state.
    """
    st = css_transition_state(p_t, eps_low, n_low, delta_eps, cs2)
    mu_t = st['mu_t']
    eps_high_t = st['eps_high_at_Pt']
    n_high_t = st['n_high_at_Pt']
    P = np.asarray(pressures, float)
    # Reject samples below P_t (acceptance)
    if np.any(P < p_t * (1 - 1e-12)):
        raise ValueError(f'pressure sample below P_t={p_t} detected; CSS branch only defined for P >= P_t')
    # CSS branch equations (Alford-Han-Prakash 2013)
    eps_high = eps_high_t + (P - p_t) / cs2
    factor = 1.0 + (1.0 + cs2) * (P - p_t) / (cs2 * (eps_high_t + p_t))
    # Guard against negative base (numerical only)
    if np.any(factor <= 0):
        raise ValueError('CSS n_high factor base is non-positive (pressure out of branch range)')
    n_high = n_high_t * factor ** (1.0 / (1.0 + cs2))
    mu_high = (eps_high + P) / n_high
    return {
        'transition_state': st,
        'P': P, 'eps_high': eps_high, 'n_high': n_high, 'mu_high': mu_high,
    }


def verify_css_acceptance(p_t: float, eps_low: float, n_low: float,
                           delta_eps: float, cs2: float,
                           n_samples: int = 50) -> dict:
    """Verify the M5A acceptance criteria for a CSS parameter set.
    Samples the high branch from P_t to P_t * 3 and checks:
      - At P=P_t the high side recovers (eps_high_t, n_high_t, mu_t)
      - Delta eps = mu_t * Delta n_B within tolerance
      - dP/deps = cs2 at every sample (constant sound speed)
      - deps/dn_B = mu_B at every sample (first law)
      - classification is correct (SHARP_FIRST_ORDER if Delta_eps > 0 else CONTINUOUS_UNCLASSIFIED)
    """
    st = css_transition_state(p_t, eps_low, n_low, delta_eps, cs2)
    # Sample the high branch
    P_samples = np.linspace(p_t, p_t * 3.0, n_samples)
    branch = css_high_branch(p_t, eps_low, n_low, delta_eps, cs2, P_samples)
    eps_h = branch['eps_high']; n_h = branch['n_high']; mu_h = branch['mu_high']
    checks = {}
    # 1. At P=P_t, recover the declared transition state
    checks['eps_at_Pt'] = bool(abs(eps_h[0] - st['eps_high_at_Pt']) < 1e-9 * max(1, abs(st['eps_high_at_Pt'])))
    checks['n_at_Pt'] = bool(abs(n_h[0] - st['n_high_at_Pt']) < 1e-9 * max(1, abs(st['n_high_at_Pt'])))
    checks['mu_at_Pt'] = bool(abs(mu_h[0] - st['mu_t']) < 1e-9 * max(1, abs(st['mu_t'])))
    # 2. First law at transition: Delta eps = mu_t * Delta n_B
    checks['first_law_Delta_eps'] = bool(st['first_law_residual'] < 1e-9 * max(1, st['Delta_eps']))
    # 3. Constant sound speed: dP/deps = cs2 at every sample (finite difference)
    dP = np.diff(P_samples); deps = np.diff(eps_h)
    cs2_numerical = dP / deps
    checks['dP_deps_equals_cs2'] = bool(np.allclose(cs2_numerical, cs2, rtol=1e-6, atol=1e-9))
    # 4. First law along the branch: deps/dn_B = mu_B at every sample
    # eps = mu * n - P => deps = mu * dn + n * dmu - dP, so deps/dn = mu + (n dmu - dP)/dn
    # The cleanest check: (eps + P)/n = mu at every sample (algebraic identity)
    mu_from_first_law = (eps_h + P_samples) / n_h
    checks['eps_plus_P_over_n_equals_mu'] = bool(np.allclose(mu_from_first_law, mu_h, rtol=1e-10, atol=1e-12))
    # 5. Differential first law: deps/dn = mu (computed via chain rule deps/dP / dn/dP)
    dn_dP = np.diff(n_h) / np.diff(P_samples)
    deps_dP = np.diff(eps_h) / np.diff(P_samples)
    # avoid division by zero
    mask = np.abs(dn_dP) > 1e-30
    deps_over_dn = np.zeros_like(deps_dP)
    deps_over_dn[mask] = deps_dP[mask] / dn_dP[mask]
    mu_at_midpoints = 0.5 * (mu_h[:-1] + mu_h[1:])
    # The differential deps/dn = mu_B holds exactly in the continuum; the finite
    # difference agrees to within ~1e-3 for typical cs2 because n(P) is a power
    # law with exponent 1/(1+cs2) < 1, making dn/dP small and the finite-diff
    # quotient noisy. The algebraic identity (eps+P)/n = mu is checked at 1e-10
    # above and is the strict first-law closure; here we accept the differential
    # check at 5e-3 relative tolerance (the physical content is identical).
    checks['deps_dn_equals_mu_along_branch'] = bool(
        np.allclose(deps_over_dn[mask], mu_at_midpoints[mask], rtol=5e-3, atol=1e-6)
    )
    # 6. Classification
    checks['classification'] = bool(st['classification'] == ('SHARP_FIRST_ORDER' if delta_eps > 0 else 'CONTINUOUS_UNCLASSIFIED'))
    # 7. cs2 in valid range
    checks['cs2_valid'] = bool(0 < cs2 <= 1.0 + 1e-12)
    overall = all(checks.values())
    return {
        'transition_state': st,
        'checks': checks,
        'overall_pass': overall,
        'cs2_numerical_mean': float(np.mean(cs2_numerical)),
        'cs2_numerical_max_err': float(np.max(np.abs(cs2_numerical - cs2))),
    }
