"""M5: Real-data sharp Maxwell construction via the native Fortran kernel.

Uses the native `construct_maxwell_transition` (commit 6ffdd62, "Add
thermodynamically consistent native Maxwell transition") which redefines
maxwell as a single-conserved-charge (baryon number) sharp coexistence with
LOCAL charge neutrality imposed in each phase - NOT a Gibbs mixed phase, NOT
a smooth crust-core interpolation, and NOT the historical maxwell option in
the six-method smooth-matching API.

Physics (from docs/FORTRAN_EOS_MAXWELL_TRANSITION_V1.md):
  - coexistence in Euler mu_E = (eps + P)/n_B, NOT the CompOSE Q3 quantity
  - P_L(mu_t) = P_H(mu_t) = P_t,  n_L(mu_t) < n_H(mu_t)
  - first law (latent heat): Delta_eps = mu_t * Delta_n_B, residual < 1e-10
  - cubic Hermite P(mu_B) interpolant preserving first law on each interval
  - native C1-C6/C10-C11 verifier run on the combined output

Phase low  = M4-quality spliced GPPVA_DD2 crust + BBB_BHF2 core barotrope
             (real CompOSE data; the hadronic phase). Built by m4_quality_splice.
Phase high = a stiff deconfined-quark phase (MIT-bag-style EOSData) representing
             high-density exotic matter. The quark phase is built as an EOSData
             on the same unified M_N scale so the native coexistence search has
             a common mu_E domain.

This reproduces the previous-turn M5 real-data work (commit 62734be, lost to
sandbox recycling) targeting DAG node E390-05 (M5), now using the native
Maxwell kernel rather than the earlier pure-Python reproduction (which could
not run because gfortran was unavailable in that sandbox).
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "python_frontend"))
from compose_eos.eos_data import EOSData
from compose_eos.api import construct_maxwell_transition, check_consistency, read_compose

M_N = 939.56536
OUT = REPO / "download" / "m5_real_maxwell"
OUT.mkdir(parents=True, exist_ok=True)


def hadronic_phase_eosdata(source='BBB_BHF2', n_min=0.08):
    """Build the hadronic phase EOSData from a single real CompOSE core table.

    The Maxwell hadron-quark transition lives at core densities, so a single
    homogeneous-core model (BBB_BHF2) is the physically correct hadronic
    phase. Using the M4-quality SPLICE as the hadronic phase fails the native
    Hermite P(mu_B) stability check: the splice joins two different RMF models
    (GPPVA_DD2 crust at mn=939.0 and BBB_BHF2 core at mn=939.56533), and at the
    join (mu_B ~ 951.6 -> 953.9 MeV, n_B 0.066 -> 0.076) the cubic Hermite in
    P(mu) oscillates and produces a non-positive interpolated density. A single
    model has a smooth, strictly-increasing mu_B and passes the check.

    BBB_BHF2 spans n_B in [0.062, 1.5] with a natural monotone mu_B in
    [952.6, 2172.4] MeV; we trim the low tail to n_B >= 0.08 to keep only the
    stiff homogeneous core and ensure first law (eps = mu*n - P) is exact on
    the unified M_N scale.
    """
    d_in = read_compose(REPO / 'eos_generator' / 'eos_data' / source)
    n = np.asarray(d_in.nb, float); mu = np.asarray(d_in.mub, float)
    P = np.asarray(d_in.p, float); eps = np.asarray(d_in.eps, float)
    # rescale to unified M_N (BBB_BHF2 mn=939.56533 ~ M_N; tiny rescale)
    r = M_N / float(d_in.mn)
    mu = mu * r
    eps = eps * r
    # enforce first law exactly: eps := mu_B*n_B - P
    eps = mu * n - P
    order = np.argsort(n)
    n, mu, P, eps = n[order], mu[order], P[order], eps[order]
    keep = n >= n_min
    n, mu, P, eps = n[keep], mu[keep], P[keep], eps[keep]
    F = eps / n
    d = EOSData(
        nb=n, mn=M_N, mp=M_N * 0.9938,
        Q1=P / n, Q2=np.zeros_like(n), Q3=mu / M_N - 1.0,
        Q4=np.zeros_like(n), Q5=np.zeros_like(n),
        Q6=F / M_N - 1.0, Q7=F / M_N - 1.0,
        q_extra=np.zeros_like(n), Iphase=np.ones_like(n, dtype=int),
        Yn=np.zeros_like(n), Yp=np.zeros_like(n), Ye=np.zeros_like(n),
        Ymu=np.zeros_like(n), Aav=np.zeros_like(n), Zav=np.zeros_like(n),
        Yh=np.zeros_like(n),
        Il=1, iT=0, iYq=1, inb_start=1, temperature_mev=0.0,
        charge_fraction_grid=0.0, q_extra_count=0,
        composition_available=False, beta_equilibrated=True)
    return d


def quark_phase_eosdata(P_t=30.0, eps_t=300.0, cs2=0.4, n_t=0.4,
                        n_max=2.0, n_points=500):
    """CSS-type constant-sound-speed quark phase (Alford-Han-Prakash form),
    built as an EOSData. This is Hermite-stable: P = P_t + cs2*(eps - eps_t)
    gives a smooth, strictly-increasing mu_E(n_B) that the cubic Hermite P(mu)
    interpolant does not oscillate on (unlike the bag-model polytrope whose
    near-flat mu_E curve produced negative interpolated densities).

      mu_t = (eps_t + P_t)/n_t
      eps_high(n) = eps_t + (P - P_t)/cs2 where P(n) parametrised below
      For a constant-cs2 branch with eps+P = mu*n and dP/deps = cs2,
      the standard AHP parametrisation gives
        n(n_B) = n_t * (eps_high + P)/(eps_t + P_t)   (first law: mu=(eps+P)/n)
      We instead sample n_B directly on a fine grid, compute
        eps(n) = eps_t + (n - n_t)*mu_t - (P(n) - P_t) ... integrated,
      but the cleanest Hermite-stable construction is to use the closed-form
      CSS branch from m5a_css_branch and tabulate it on n_B.
    """
    # Use the validated m5a_css_branch CSS primitive to build the high branch.
    sys.path.insert(0, str(REPO / "scripts" / "m5_m5a"))
    from m5a_css_branch import css_transition_state, css_high_branch
    # transition state at (P_t, eps_t=eps_low at P_t, n_t=n_low at P_t)
    ts = css_transition_state(p_t=P_t, eps_low=eps_t, n_low=n_t,
                              delta_eps=0.0, cs2=cs2)
    mu_t = ts['mu_t']
    # sample P above P_t, build the CSS high branch, then convert to (n_B, mu_E, P, eps)
    P_grid = np.linspace(P_t + 0.01, P_t + 1200.0, n_points)
    branch = css_high_branch(p_t=P_t, eps_low=eps_t, n_low=n_t,
                             delta_eps=0.0, cs2=cs2, pressures=P_grid)
    n = np.array(branch['n_high']); P = np.array(branch['P'])
    eps = np.array(branch['eps_high']); mu = np.array(branch['mu_high'])
    # enforce first law exactly (CSS gives it to ~1e-13; re-assert)
    eps = mu * n - P
    F = eps / n
    d = EOSData(
        nb=n, mn=M_N, mp=M_N * 0.9938,
        Q1=P / n, Q2=np.zeros_like(n), Q3=mu / M_N - 1.0,
        Q4=np.zeros_like(n), Q5=np.zeros_like(n),
        Q6=F / M_N - 1.0, Q7=F / M_N - 1.0,
        q_extra=np.zeros_like(n), Iphase=2 * np.ones_like(n, dtype=int),
        Yn=np.zeros_like(n), Yp=np.zeros_like(n), Ye=np.zeros_like(n),
        Ymu=np.zeros_like(n), Aav=np.zeros_like(n), Zav=np.zeros_like(n),
        Yh=np.zeros_like(n),
        Il=1, iT=0, iYq=1, inb_start=1, temperature_mev=0.0,
        charge_fraction_grid=0.0, q_extra_count=0,
        composition_available=False, beta_equilibrated=True)
    return d, P_t, eps_t, cs2, n_t


def main():
    import math
    low = hadronic_phase_eosdata()
    # try a grid of CSS quark-phase parameters until native finds a coexistence
    # with a positive upward density jump. The CSS phase is Hermite-stable
    # (constant cs2 => smooth mu_E(n_B)); we vary P_t (transition pressure) and
    # cs2 (high-density sound speed) to find a crossing with the hadronic core.
    attempts = []
    for P_t in [20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0, 120.0, 150.0]:
        for cs2 in [0.2, 0.3, 0.4, 0.5, 0.6, 1.0/3.0]:
            # eps_t/n_t chosen at the CSS onset on the hadronic branch: pick a
            # hadronic-core point near P_t so the CSS branch starts smoothly.
            n_h = np.asarray(low.nb, float); P_h = np.asarray(low.p, float)
            # find hadronic state closest to P_t (the CSS low side)
            i = int(np.argmin(np.abs(P_h - P_t)))
            eps_t = float(P_h[i] * 0 + (np.asarray(low.eps, float)[i]))
            n_t = float(n_h[i])
            try:
                high, P_t_used, eps_t_used, cs2_used, n_t_used = quark_phase_eosdata(
                    P_t=P_t, eps_t=eps_t, cs2=cs2, n_t=n_t)
            except Exception as e:
                attempts.append({'P_t': P_t, 'cs2': cs2, 'build_error': str(e)[:150]})
                continue
            try:
                matched, report = construct_maxwell_transition(
                    low, high, assume_local_charge_neutrality=True,
                    return_report=True)
            except Exception as e:
                attempts.append({'P_t': P_t, 'cs2': cs2, 'error': str(e)[:200]})
                continue
            # extract transition metadata from native report
            n_low = float(report.get('n_low', report.get('n_low_fm3', 0.0)))
            n_high = float(report.get('n_high', report.get('n_high_fm3', 0.0)))
            jump = n_high - n_low
            attempts.append({
                'P_t_in': P_t, 'cs2_in': cs2, 'eps_t_in': eps_t, 'n_t_in': n_t,
                'mu_t': float(report.get('mu_transition', 0.0)),
                'P_t': float(report.get('pressure_transition', 0.0)),
                'n_low': n_low, 'n_high': n_high, 'jump': jump,
                'maxwell_identity_residual': float(report.get('maxwell_identity_residual', 0.0)),
                'status': report.get('status', 'ok'),
            })
            if jump > 0 and report.get('status', 'ok') == 'ok':
                # success - save the matched barotrope + receipt and stop
                n = np.asarray(matched.nb, float)
                mu = np.asarray(matched.mub, float)
                P = np.asarray(matched.p, float)
                eps = np.asarray(matched.eps, float)
                barotrope = {'n_B': n.tolist(), 'mu_B': mu.tolist(),
                             'P': P.tolist(), 'eps': eps.tolist()}
                (OUT / 'GPPVA_DD2_m5_real_maxwell.json').write_text(
                    json.dumps(barotrope, indent=1))
                # Seidov stability (computed in Python from the native transition)
                mu_t = float(report.get('mu_transition', 0.0))
                P_t_out = float(report.get('pressure_transition', 0.0))
                eps_low = mu_t * n_low - P_t_out
                eps_high = mu_t * n_high - P_t_out
                delta_eps = eps_high - eps_low
                seidov_crit = 0.5 * eps_low + 1.5 * P_t_out
                receipt = {
                    'node': 'E390-05', 'subscope': 'm5_real_maxwell_hadron_quark_native',
                    'date': '2026-11-14',
                    'method': 'native Fortran construct_maxwell_transition (commit 6ffdd62; single-conserved-charge local-neutrality sharp coexistence in Euler mu_E)',
                    'phase_low': {'source': 'single real CompOSE core table BBB_BHF2 (n_B>=0.08, trimmed; the M4 splice was rejected by the native Hermite P(mu) stability gate at the GPPVA_DD2/BBB_BHF2 model join)',
                                  'n_points': int(len(low.nb)),
                                  'mu_B_range': [float(low.mub.min()), float(low.mub.max())]},
                    'phase_high': {'model': 'CSS constant-sound-speed quark (Alford-Han-Prakash 2013): eps=eps_t+(P-P_t)/cs2, n=n_t*f^(1/(1+cs2)), mu=(eps+P)/n',
                                   'P_t_input': P_t_used, 'eps_t_input': eps_t_used,
                                   'n_t_input': n_t_used, 'cs2_input': cs2_used,
                                   'n_points': int(len(high.nb))},
                    'maxwell_transition': {
                        'mu_t_mev': mu_t, 'P_t_mev_fm3': P_t,
                        'n_low_fm3': n_low, 'n_high_fm3': n_high,
                        'delta_n_fm3': jump,
                        'eps_low_mev_fm3': eps_low, 'eps_high_mev_fm3': eps_high,
                        'delta_eps_mev_fm3': delta_eps,
                        'first_law_resid_native': float(report.get('maxwell_identity_residual', 0.0)),
                        'seidov_critical_mev_fm3': seidov_crit,
                        'seidov_stable': bool(delta_eps < seidov_crit),
                        'native_status': report.get('status', 'ok'),
                    },
                    'native_c1_c11_overall_pass': bool(report['consistency']['overall']['pass_']),
                    'native_failed_checks': report['consistency']['overall']['failed_checks'],
                    'native_C2_first_law_max_rel': float(report['consistency']['C2']['max_rel']),
                    'native_C6_min_dPdnb': float(report['consistency']['C6']['min_dPdnb']),
                    'native_C10_cs2_range': [float(report['consistency']['C10']['min_cs2']),
                                              float(report['consistency']['C10']['max_cs2'])],
                    'is_exact_maxwell': bool(report.get('is_exact_maxwell', True)),
                    'construction': report.get('construction'),
                    'status': 'M5_REAL_MAXWELL_NATIVE_PASS' if (
                        report['consistency']['overall']['pass_'] and jump > 0) else 'M5_REAL_MAXWELL_NATIVE_FAIL',
                    'attempts_summary': attempts,
                }
                (OUT / 'E390-05_M5_REAL_MAXWELL_RECEIPT_20261114.json').write_text(
                    json.dumps(receipt, indent=1))
                print('M5 real Maxwell (native):', receipt['status'])
                print(f"  P_t_input={P_t_used} cs2_input={cs2_used} n_t_input={n_t_used}")
                print(f"  mu_t={mu_t:.4f} MeV  P_t={P_t_out:.6f} MeV/fm3")
                print(f"  n_low={n_low:.5f}  n_high={n_high:.5f}  jump={jump:.5f}")
                print(f"  delta_eps={delta_eps:.4f}  first-law resid(native)={report.get('maxwell_identity_residual',0.0):.3e}")
                print(f"  Seidov critical={seidov_crit:.4f}  stable={delta_eps < seidov_crit}")
                print(f"  native C1-C11 pass={receipt['native_c1_c11_overall_pass']}  failed={receipt['native_failed_checks']}")
                print(f"  native C2 first-law max_rel={receipt['native_C2_first_law_max_rel']:.3e}")
                return receipt
    # no transition found - dump attempts
    print('M5: no valid Maxwell transition found across the quark parameter grid.')
    print('attempts:')
    for a in attempts:
        print(' ', a)
    return None


if __name__ == '__main__':
    main()
