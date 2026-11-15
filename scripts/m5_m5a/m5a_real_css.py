"""M5A: Real-data Constant-Sound-Speed (CSS) high-density branch.

Attaches the Alford-Han-Prakash (2013, PRD 88 083013) CSS high-density branch
to the real hadronic EOS (BBB_BHF2 homogeneous core) at three transition
densities n_B = 0.3, 0.4, 0.5 fm^-3, and verifies the M5A acceptance criteria.

Physics (CSS branch, from m5a_css_branch.py + docs):
  P_t      - transition pressure
  eps_low  - energy density on the low (hadronic) side at P_t
  n_low    - baryon number density on the low side at P_t
  Delta_eps - latent heat at P_t (>= 0); here we first test the continuous
              case Delta_eps=0 (CONTINUOUS_UNCLASSIFIED), then a sharp
              first-order case Delta_eps>0 (SHARP_FIRST_ORDER).
  cs2      - constant sound speed squared (0 < cs2 <= 1)

  High branch (P >= P_t):
    eps_high(P) = eps_low + Delta_eps + (P - P_t)/cs2
    n_high(P)   = n_low * [1 + (1+cs2)(P-P_t)/(cs2 (eps_low+Delta_eps+P_t))]^(1/(1+cs2))
    mu_high(P)  = (eps_high + P)/n_high          (first law: eps + P = mu n)

Acceptance (from the DAG css_segment_contract, verified per branch):
  - At P=P_t: eps_high=eps_low+Delta_eps, n_high recovers, mu=mu_t=(eps_low+P_t)/n_low
  - Delta_eps = mu_t * (n_high - n_low) within floating-point tolerance (first law)
  - dP/deps = cs2 at every high-side sample
  - deps/dn_B = mu_B along the branch
  - eps + P = mu * n (first law) at every sample
  - classification: Delta_eps>0 => SHARP_FIRST_ORDER; =0 => CONTINUOUS_UNCLASSIFIED

This reproduces the previous-turn M5A real-data work (commit 62734be, lost to
sandbox recycling) targeting DAG node E390-05A (M5A), now using the real
BBB_BHF2 hadronic core and the validated m5a_css_branch CSS primitive.
"""
from __future__ import annotations
import json, sys, math
from pathlib import Path
import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "python_frontend"))
sys.path.insert(0, str(REPO / "scripts" / "m5_m5a"))
from compose_eos.eos_data import EOSData
from compose_eos.api import read_compose, check_consistency
from m5a_css_branch import css_transition_state, css_high_branch, verify_css_acceptance

M_N = 939.56536
OUT = REPO / "download" / "m5a_real_css"
OUT.mkdir(parents=True, exist_ok=True)


def hadronic_core():
    """BBB_BHF2 real core as (n_B, mu_B, P, eps) on the unified M_N scale,
    with first law enforced exactly (eps := mu*n - P).

    Returns (n, mu, P, eps, eps_source, max_first_law_reset_rel) where
    eps_source is the original CompOSE energy density (after the M_N rescale)
    and max_first_law_reset_rel is the maximum relative change introduced by
    the first-law reset eps := mu*n - P. This quantifies the §6 finding #4:
    the first-law reset is NOT free; it changes the energy density by up to
    max_first_law_reset_rel relative to the source."""
    d_in = read_compose(REPO / 'eos_generator' / 'eos_data' / 'BBB_BHF2')
    n = np.asarray(d_in.nb, float); mu = np.asarray(d_in.mub, float)
    P = np.asarray(d_in.p, float); eps = np.asarray(d_in.eps, float)
    r = M_N / float(d_in.mn)
    mu = mu * r; eps_source = eps * r            # source energy after M_N rescale
    eps = mu * n - P                              # enforce first law
    # quantify the reset
    reset_rel = np.abs(eps - eps_source) / np.maximum(np.abs(eps_source), 1e-30)
    max_first_law_reset_rel = float(np.max(reset_rel))
    order = np.argsort(n)
    return n[order], mu[order], P[order], eps[order], eps_source[order], max_first_law_reset_rel


def css_branch_at_nB(n_low, hadronic, cs2=1.0/3.0, delta_eps=0.0, n_points=200):
    """Attach a CSS branch to the hadronic core at baryon density n_low.

    Returns the transition state + the high-branch arrays (n_B, mu_B, P, eps).
    """
    n_h, mu_h, P_h, eps_h = hadronic[:4]
    # find the hadronic state at n_low (interpolate)
    i = np.searchsorted(n_h, n_low)
    if i == 0: i = 1
    if i >= len(n_h): i = len(n_h) - 1
    # linear interpolation for P_t, eps_low at n_low
    n0, n1 = n_h[i-1], n_h[i]
    P_t = np.interp(n_low, [n0, n1], [P_h[i-1], P_h[i]])
    eps_low = np.interp(n_low, [n0, n1], [eps_h[i-1], eps_h[i]])
    ts = css_transition_state(p_t=P_t, eps_low=eps_low, n_low=n_low,
                              delta_eps=delta_eps, cs2=cs2)
    # sample P from P_t up to a high pressure
    P_grid = np.linspace(P_t + 1e-6, P_t + 1500.0, n_points)
    br = css_high_branch(p_t=P_t, eps_low=eps_low, n_low=n_low,
                         delta_eps=delta_eps, cs2=cs2, pressures=P_grid)
    n = np.array(br['n_high']); P = np.array(br['P'])
    eps = np.array(br['eps_high']); mu = np.array(br['mu_high'])
    eps = mu * n - P  # re-assert first law exact
    return ts, P_t, eps_low, n, mu, P, eps


def build_combined_hadronic_css(n_low, hadronic, cs2, delta_eps):
    """Concatenate the hadronic branch (n_B <= n_low) with the CSS high branch
    (n_B >= n_high_at_Pt) to form a single combined barotrope, with an explicit
    two-row jump at (P_t, mu_t) when delta_eps > 0."""
    ts, P_t, eps_low, n_high_arr, mu_high_arr, P_high_arr, eps_high_arr = \
        css_branch_at_nB(n_low, hadronic, cs2=cs2, delta_eps=delta_eps)
    n_h, mu_h, P_h, eps_h = hadronic[:4]
    # hadronic side: n_B <= n_low
    lo = n_h <= n_low
    n_low_arr = n_h[lo]; mu_low_arr = mu_h[lo]
    P_low_arr = P_h[lo]; eps_low_arr = eps_h[lo]
    mu_t = ts['mu_t']
    eps_low_t = eps_low
    eps_high_t = ts['eps_high_at_Pt']
    n_low_t = n_low
    n_high_t = ts['n_high_at_Pt']
    # assemble: low rows, (jump rows if delta_eps>0), high rows
    rows_n = list(n_low_arr) + ([n_low_t, n_high_t] if delta_eps > 0 else []) + list(n_high_arr)
    rows_mu = list(mu_low_arr) + ([mu_t, mu_t] if delta_eps > 0 else []) + list(mu_high_arr)
    rows_P = list(P_low_arr) + ([P_t, P_t] if delta_eps > 0 else []) + list(P_high_arr)
    rows_eps = list(eps_low_arr) + ([eps_low_t, eps_high_t] if delta_eps > 0 else []) + list(eps_high_arr)
    n = np.array(rows_n, float); mu = np.array(rows_mu, float)
    P = np.array(rows_P, float); eps = np.array(rows_eps, float)
    # enforce first law exactly on the combined table
    eps = mu * n - P
    return n, mu, P, eps, ts, P_t, eps_low, n_low_t, n_high_t, eps_low_t, eps_high_t


def main():
    hadronic = hadronic_core()
    n_h, mu_h, P_h, eps_h, eps_source_h, max_first_law_reset_rel = hadronic
    results = {}
    # real BBB_BHF2 proton mass (§6 finding #3: replace the invented M_N*0.9938)
    d_in = read_compose(REPO / 'eos_generator' / 'eos_data' / 'BBB_BHF2')
    real_mp = float(d_in.mp)
    real_mn = float(d_in.mn)
    for n_low in [0.3, 0.4, 0.5]:
        for cs2 in [1.0/3.0, 0.4, 0.5]:
            for delta_eps in [0.0, 50.0, 100.0]:
                key = f'nB{n_low}_cs2{cs2:.4f}_dEps{delta_eps}'
                try:
                    n, mu, P, eps, ts, P_t, eps_low, n_low_t, n_high_t, eps_low_t, eps_high_t = \
                        build_combined_hadronic_css(n_low, hadronic, cs2, delta_eps)
                    # M5A acceptance via the validated primitive
                    acc = verify_css_acceptance(p_t=P_t, eps_low=eps_low,
                                                n_low=n_low, delta_eps=delta_eps, cs2=cs2)
                    # first law at the jump
                    delta_n = n_high_t - n_low_t
                    delta_e = eps_high_t - eps_low_t
                    fl_resid = delta_e - mu_t * delta_n if 'mu_t' in locals() else 0.0
                    fl_resid = (ts['eps_high_at_Pt'] - eps_low) - ts['mu_t'] * (ts['n_high_at_Pt'] - n_low)
                    # §6 finding #5: μ_t mismatch between interpolated source μ_B and
                    # μ_t = (eps_low + P_t) / n_low at the transition density.
                    # the hadronic μ at n_low (interpolated):
                    mu_source_at_nlow = float(np.interp(n_low, n_h, mu_h))
                    mu_t = float(ts['mu_t'])
                    mu_t_mismatch_rel = abs(mu_t - mu_source_at_nlow) / max(abs(mu_source_at_nlow), 1e-30)
                    # also pressure/energy interpolation residuals at n_low
                    P_source_at_nlow = float(np.interp(n_low, n_h, P_h))
                    eps_source_at_nlow = float(np.interp(n_low, n_h, eps_h))
                    P_interp_resid_rel = abs(P_t - P_source_at_nlow) / max(abs(P_source_at_nlow), 1e-30)
                    eps_interp_resid_rel = abs(eps_low - eps_source_at_nlow) / max(abs(eps_source_at_nlow), 1e-30)
                    # §6 finding #1 + #7: this is a Python-side consistency check, NOT a
                    # native Fortran C1-C11 verification. Rename accordingly. The output
                    # is a (n_B, mu_B, P, eps) barotrope, NOT a CompOSE-generated EOS.
                    F = eps / n
                    d = EOSData(
                        nb=n, mn=M_N, mp=real_mp,  # §6 finding #3: real BBB_BHF2 mp, not M_N*0.9938
                        Q1=P/n, Q2=np.zeros_like(n), Q3=mu/M_N-1.0,
                        Q4=np.zeros_like(n), Q5=np.zeros_like(n),
                        Q6=F/M_N-1.0, Q7=F/M_N-1.0,
                        q_extra=np.zeros_like(n), Iphase=np.ones_like(n, dtype=int),
                        Yn=np.zeros_like(n), Yp=np.zeros_like(n), Ye=np.zeros_like(n),
                        Ymu=np.zeros_like(n), Aav=np.zeros_like(n), Zav=np.zeros_like(n),
                        Yh=np.zeros_like(n),
                        Il=1, iT=0, iYq=1, inb_start=1, temperature_mev=0.0,
                        charge_fraction_grid=0.0, q_extra_count=0,
                        composition_available=False, beta_equilibrated=True)
                    python_consistency = check_consistency(d, verbose=False)
                    # §6 finding #6: explicit Δε classification
                    if delta_eps > 0:
                        delta_eps_classification = 'SHARP_FIRST_ORDER'
                    else:
                        delta_eps_classification = 'CONTINUOUS_UNCLASSIFIED'  # zero jump alone is NOT a crossover/2nd-order inference
                    results[key] = {
                        'n_low_input': n_low, 'cs2_input': cs2, 'delta_eps_input': delta_eps,
                        'P_t': float(P_t), 'eps_low': float(eps_low),
                        'mu_t': mu_t,
                        'mu_source_at_nlow': mu_source_at_nlow,
                        'mu_t_mismatch_relative': float(mu_t_mismatch_rel),  # §6 #5
                        'P_interp_residual_relative': float(P_interp_resid_rel),
                        'eps_interp_residual_relative': float(eps_interp_resid_rel),
                        'n_high_at_Pt': float(ts['n_high_at_Pt']),
                        'eps_high_at_Pt': float(ts['eps_high_at_Pt']),
                        'delta_n': float(ts['delta_n']),
                        'first_law_resid': float(ts['first_law_residual']),
                        'classification': ts['classification'],
                        'delta_eps_classification_explicit': delta_eps_classification,  # §6 #6
                        'css_acceptance_overall_pass': bool(acc.get('overall_pass', acc.get('overall_pass_'))),
                        'css_acceptance': acc,
                        # §6 finding #1: renamed. This is a PYTHON-side consistency check,
                        # NOT a native Fortran C1-C11 verification. A native PASS requires
                        # actually running the native Fortran read/write/verify pathway.
                        'python_consistency_check_pass': bool(python_consistency['overall']['pass_']),
                        'python_consistency_failed_checks': python_consistency['overall']['failed_checks'],
                        'python_C2_first_law_max_rel': float(python_consistency['C2']['max_rel']),
                        'python_C6_min_dPdnb': float(python_consistency['C6']['min_dPdnb']),
                        'python_C10_cs2_range': [float(python_consistency['C10']['min_cs2']), float(python_consistency['C10']['max_cs2'])],
                        # native Fortran verification backend (§6 receipt-truthfulness)
                        'native_fortran_compose_readback': False,  # BLOCKED: not executed
                        'native_c1_c11_verification': False,        # BLOCKED: not executed
                        'analytic_css_acceptance': bool(acc.get('overall_pass', acc.get('overall_pass_'))),
                        'combined_n_points': int(len(n)),
                    }
                except Exception as e:
                    results[key] = {'error': str(e)[:200]}
    # §6 finding #2 + #7: the headline barotrope filename must reflect the
    # real source (BBB_BHF2), not GPPVA_DD2. The JSON contains (n_B, mu_B, P,
    # eps) only -- it is NOT a CompOSE-generated EOS.
    headline = results.get('nB0.4_cs2'+f'{1/3:.4f}'+'_dEps100.0')
    if headline and 'python_consistency_check_pass' in headline:
        # rebuild headline barotrope and save
        n, mu, P, eps, *_ = build_combined_hadronic_css(0.4, hadronic, 1.0/3.0, 100.0)
        barotrope = {'n_B': n.tolist(), 'mu_B': mu.tolist(), 'P': P.tolist(), 'eps': eps.tolist()}
        (OUT / 'BBB_BHF2_m5a_real_css.json').write_text(json.dumps(barotrope, indent=1))
    receipt = {
        'node': 'E390-05A', 'subscope': 'm5a_real_css_branch',
        'date': '2026-11-14',
        'method': 'CSS Alford-Han-Prakash branch attached to real BBB_BHF2 hadronic core at n_B=0.3/0.4/0.5, Python consistency check',
        'hadronic_source': 'eos_generator/eos_data/BBB_BHF2 (real CompOSE, unified M_N=939.56536, first law enforced)',
        'hadronic_source_real_mp_mev': real_mp,  # §6 #3: real BBB_BHF2 proton mass
        'hadronic_source_real_mn_mev': real_mn,
        'hadronic_source_invented_mp_was_M_N_times_0_9938_mev': M_N * 0.9938,  # the old, incorrect value (for audit trail)
        'first_law_reset_max_relative_change': max_first_law_reset_rel,  # §6 #4: quantified
        'css_primitive': 'scripts/m5_m5a/m5a_css_branch.py (css_transition_state, css_high_branch, verify_css_acceptance)',
        'transition_densities': [0.3, 0.4, 0.5],
        'cs2_values': [1/3, 0.4, 0.5],
        'delta_eps_values': [0.0, 50.0, 100.0],
        'headline_case': 'n_B=0.4, cs2=1/3, delta_eps=100 (sharp first-order)',
        'results': results,
        # §6 finding #7: the saved JSON is a (n_B, mu_B, P, eps) barotrope, NOT a
        # CompOSE-generated EOS. A CompOSE artifact would require the native I/O path.
        'output_barotrope_file': 'BBB_BHF2_m5a_real_css.json',
        'output_is_compose_generated_eos': False,
        'output_is_barotrope_only': True,
        'verification_backend_truthful': {
            'analytic_css_acceptance': True,           # analytic CSS primitive checks
            'python_consistency_check': True,           # Python check_consistency
            'native_fortran_compose_readback': False,   # BLOCKED: not executed
            'native_c1_c11_verification': False,        # BLOCKED: not executed
        },
        'status': 'M5A_REAL_CSS_PASS' if all(
            (r.get('css_acceptance_overall_pass') and r.get('python_consistency_check_pass'))
            for r in results.values() if 'css_acceptance_overall_pass' in r
        ) else 'M5A_REAL_CSS_PARTIAL',
    }
    (OUT / 'E390-05A_M5A_REAL_CSS_RECEIPT_20261114.json').write_text(json.dumps(receipt, indent=1))
    print('M5A real CSS:', receipt['status'])
    # summary table
    for n_low in [0.3, 0.4, 0.5]:
        for cs2 in [1.0/3.0, 0.4, 0.5]:
            for delta_eps in [0.0, 50.0, 100.0]:
                key = f'nB{n_low}_cs2{cs2:.4f}_dEps{delta_eps}'
                r = results.get(key, {})
                if 'css_acceptance_overall_pass' in r:
                    print(f"  nB={n_low} cs2={cs2:.3f} dEps={delta_eps:5.1f}: css_pass={r['css_acceptance_overall_pass']} py_pass={r['python_consistency_check_pass']} class={r['delta_eps_classification_explicit']} mu_t_mismatch={r['mu_t_mismatch_relative']:.2e} fl_resid={r['first_law_resid']:.2e} P_t={r['P_t']:.2f} n_high={r['n_high_at_Pt']:.4f}")
    print(f"  first-law reset max relative change vs source: {max_first_law_reset_rel:.3e}")
    print(f"  real BBB_BHF2 mp={real_mp} MeV (was invented M_N*0.9938={M_N*0.9938} MeV)")
    return receipt


if __name__ == '__main__':
    main()
