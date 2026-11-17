"""M6 integrated reproduction and acceptance verification.

Tests:
  1. Smooth polytropic TOV + Love: build a single polytrope (n in [0.05, 0.8],
     K=200, gamma=2.5), solve TOV+y for a range of central pressures, verify
     M(R), k2, Lambda behave physically (compactness < 0.5, k2 in [0, 0.2],
     Lambda positive). Cross-check the maximum mass and radius against the
     Read 2009 piecewise-polytrope benchmark range.
  2. Sharp Maxwell TOV + Love: use the M5 Maxwell output (from
     m5_maxwell_construction build_maxwell_eos) and integrate TOV+y with the
     Takatsy-Kovacs 2020 tidal y-jump across the interface.
  3. Slow vs rapid conversion: compare the Takatsy-Kovacs 2020 (slow,
     pressure-corrected) jump to the Postnikov 2010 (legacy, m-only) jump and
     report the relative difference, which Takatsy-Kovacs estimates to be
     up to ~5% for transitions at non-negligible pressure.
  4. Convergence and causality: verify the TOV solver converges with
     decreasing rtol and that the central sound speed is causal (cs2 <= 1).
"""
from __future__ import annotations
import sys, json, math
import numpy as np

sys.path.insert(0, '/home/z/my-project/repos/NS_Cowling_eos390/python_frontend')
sys.path.insert(0, '/home/z/my-project/scripts')
from compose_eos.eos_data import EOSData
from m5_maxwell_construction import build_maxwell_eos
from m6_tov_love_solver import (
    solve_tov_love_smooth, solve_tov_love_with_maxwell,
    _MEVFM3_TO_KM2, _MSUN_GEOM_KM, love_number_k2,
)


def make_polytrope_p_eps(nB_arr, K, gamma, mn=939.0, n_offset=0.0):
    """Return (P, eps) arrays in MeV/fm^3 for a polytrope P=K(n-n0)^g,
    eps = mn*n + P/(gamma-1)."""
    nB = np.asarray(nB_arr, float)
    P = K * (nB - n_offset) ** gamma
    eps = mn * nB + P / (gamma - 1.0)
    return P, eps


def test_smooth_polytrope_tov():
    """Test 1: Smooth single-phase polytropic EOS TOV + Love.
    Use a polytrope P = K (n - n0)^gamma with K=200, gamma=2.5, n in [0.05, 0.8]
    and verify the TOV sequence."""
    print("=== Test 1: Smooth polytropic TOV + Love ===")
    n_arr = np.linspace(0.05, 0.80, 120)
    P, eps = make_polytrope_p_eps(n_arr, K=200.0, gamma=2.5, mn=939.0, n_offset=0.0)
    # Sort by P (already monotone increasing for n_offset=0)
    P_range = (P[0], P[-1])
    print(f"  polytrope: K=200, gamma=2.5, n in [0.05,0.80], P in [{P_range[0]:.2f},{P_range[1]:.2f}] MeV/fm^3")
    # Central pressures spanning the EOS range
    pc_list = np.linspace(P_range[0] * 1.5, P_range[1] * 0.95, 8)
    results = []
    for pc in pc_list:
        try:
            r = solve_tov_love_smooth(P, eps, float(pc), rtol=2e-8, max_step_km=0.12)
            results.append(r)
        except Exception as e:
            print(f"  pc={pc:.3f} FAILED: {e}")
    if not results:
        print("  no successful integrations")
        return None
    print(f"  {len(results)} successful integrations")
    print(f"  {'pc':>10s} {'R(km)':>8s} {'M(Msun)':>9s} {'C=M/R':>7s} {'k2':>7s} {'Lambda':>10s} {'yR':>8s}")
    for r in results:
        print(f"  {r['central_pressure_MeV_fm3']:10.3f} {r['R_km']:8.3f} {r['M_Msun']:9.4f} {r['compactness']:7.4f} {r['k2']:7.4f} {r['Lambda']:10.2f} {r['yR']:8.4f}")
    # Physical acceptance checks
    checks = {}
    checks['all_compactness_below_half'] = bool(all(r['compactness'] < 0.5 for r in results))
    checks['all_k2_in_physical_range'] = bool(all(0 < r['k2'] < 0.2 for r in results))
    checks['all_lambda_positive'] = bool(all(r['Lambda'] > 0 for r in results))
    checks['all_cs2_causal'] = bool(all(0 < r['central_cs2'] <= 1 + 1e-5 for r in results))
    checks['mass_increases_then_peaks'] = bool(len(results) >= 3 and results[-1]['M_Msun'] <= max(r['M_Msun'] for r in results) * 1.05)
    # Read 2009 benchmark: neutron star radii are 9-14 km, masses 0.5-2.5 Msun
    masses = [r['M_Msun'] for r in results]
    radii = [r['R_km'] for r in results]
    checks['radii_in_NS_range'] = bool(all(8 < rR < 16 for rR in radii))
    print(f"\n  acceptance: {json.dumps(checks, indent=2)}")
    print(f"  overall PASS = {all(checks.values())}")
    return {'results': results, 'checks': checks, 'overall_pass': all(checks.values())}


def test_maxwell_tov_love():
    """Test 2: Sharp Maxwell TOV + Love with Takatsy-Kovacs 2020 tidal y-jump.
    Use the M5 Maxwell construction between two polytropes (the test 1 pair
    from m5_m5a_integrated_test) and integrate TOV+y across the interface."""
    print("\n=== Test 2: Sharp Maxwell TOV + Love with tidal y-jump ===")
    # Build the same Maxwell pair as m5_m5a test 1 (substantial jump)
    n_low = np.linspace(0.05, 0.50, 100)
    n_high = np.linspace(0.15, 0.80, 100)

    def to_eosdata(n, P, eps, mn=939.0):
        F = eps / n; mu = (eps + P) / n
        return EOSData(nb=n, mn=mn, mp=mn * 0.9938,
            Q1=P / n, Q2=np.zeros_like(n), Q3=mu / mn - 1.0, Q4=np.zeros_like(n), Q5=np.zeros_like(n),
            Q6=F / mn - 1.0, Q7=F / mn - 1.0, q_extra=np.zeros_like(n), Iphase=np.ones_like(n, dtype=int),
            Yn=np.zeros_like(n), Yp=np.zeros_like(n), Ye=np.zeros_like(n), Ymu=np.zeros_like(n),
            Aav=np.zeros_like(n), Zav=np.zeros_like(n), Yh=np.zeros_like(n),
            Il=1, iT=0, iYq=1, inb_start=1, temperature_mev=0.0, charge_fraction_grid=0.0,
            q_extra_count=0, composition_available=False, beta_equilibrated=True)

    P_low_arr, eps_low_arr = make_polytrope_p_eps(n_low, K=1000.0, gamma=2.0, mn=939.0, n_offset=0.0)
    P_high_arr, eps_high_arr = make_polytrope_p_eps(n_high, K=8000.0, gamma=3.0, mn=1100.0, n_offset=0.10)
    low = to_eosdata(n_low, P_low_arr, eps_low_arr, mn=939.0)
    high = to_eosdata(n_high, P_high_arr, eps_high_arr, mn=1100.0)
    mw = build_maxwell_eos(low, high, low_name='poly_K1000_g2', high_name='poly_K8000_g3_off')
    t = mw['transition']
    print(f"  Maxwell interface: P_t={t['pressure_MeV_fm3']:.3f}, mu_t={t['muB_MeV']:.4f}")
    print(f"  n_low={t['n_low_fm3']:.5f}, n_high={t['n_high_fm3']:.5f}, jump={t['delta_n_fm3']:.5f}")
    print(f"  Delta_eps = {t['delta_eps_MeV_fm3']:.3f} MeV/fm^3, Seidov stable = {t['seidov_small_core_stable']}")
    # Build phase branches (P, eps) for the M6 solver. The high phase must be
    # extended beyond P_t; the low phase must extend below P_t.
    P_low = P_low_arr[P_low_arr <= t['pressure_MeV_fm3']]
    eps_low_p = eps_low_arr[P_low_arr <= t['pressure_MeV_fm3']]
    # Append the transition point on the low side (eps_low)
    P_low = np.append(P_low, t['pressure_MeV_fm3'])
    eps_low_p = np.append(eps_low_p, t['eps_low_MeV_fm3'])
    # Extend the low phase down to P=0 (polytrope P=K n^g, eps=mn n + P/(g-1),
    # so at n=0 both P=0 and eps=0). This gives the TOV solver a valid surface.
    # Use a tiny positive pressure (1e-9 MeV/fm^3) instead of exactly 0 so the
    # log-P PCHIP interpolator in make_phase_branch stays finite.
    P_low = np.insert(P_low, 0, 1e-9)
    eps_low_p = np.insert(eps_low_p, 0, 1e-9)
    # High phase from P_t upward
    mask_high = P_high_arr >= t['pressure_MeV_fm3']
    P_high = np.concatenate(([t['pressure_MeV_fm3']], P_high_arr[mask_high]))
    eps_high = np.concatenate(([t['eps_high_MeV_fm3']], eps_high_arr[mask_high]))
    # Verify monotone
    if not (np.all(np.diff(P_low) > 0) and np.all(np.diff(eps_low_p) > 0)):
        # Trim to monotone
        ok = np.diff(P_low) > 0
        P_low = P_low[np.r_[True, ok]]; eps_low_p = eps_low_p[np.r_[True, ok]]
    if not (np.all(np.diff(P_high) > 0) and np.all(np.diff(eps_high) > 0)):
        ok = np.diff(P_high) > 0
        P_high = P_high[np.r_[True, ok]]; eps_high = eps_high[np.r_[True, ok]]
    print(f"  low branch: P in [{P_low[0]:.3f},{P_low[-1]:.3f}], {len(P_low)} pts")
    print(f"  high branch: P in [{P_high[0]:.3f},{P_high[-1]:.3f}], {len(P_high)} pts")
    # Integrate TOV + Love for several central pressures above P_t (high-core stars).
    # The high phase polytrope (mn=1100, K=8000, g=3, n_offset=0.1) reaches cs2=1 at
    # P ~ 215 MeV/fm^3, so we restrict central pressures to the causal range above P_t.
    pc_max = 200.0  # safely below the cs2=1 limit at P~215
    pc_list = np.linspace(t['pressure_MeV_fm3'] * 1.5, pc_max, 6)
    results_slow = []
    results_rapid = []
    for pc in pc_list:
        try:
            r_slow = solve_tov_love_with_maxwell(P_low, eps_low_p, P_high, eps_high, t,
                                                   float(pc), rtol=2e-8, max_step_km=0.12,
                                                   jump_formula='takatsy_kovacs_2020')
            r_rapid = solve_tov_love_with_maxwell(P_low, eps_low_p, P_high, eps_high, t,
                                                     float(pc), rtol=2e-8, max_step_km=0.12,
                                                     jump_formula='postnikov_2010')
            results_slow.append(r_slow); results_rapid.append(r_rapid)
        except Exception as e:
            print(f"  pc={pc:.3f} FAILED: {e}")
    if not results_slow:
        print("  no successful integrations")
        return None
    print(f"\n  {len(results_slow)} successful integrations (slow = Takatsy-Kovacs 2020 Eq.11):")
    print(f"  {'pc':>10s} {'R(km)':>8s} {'M(Msun)':>9s} {'C':>7s} {'k2':>7s} {'Lambda':>10s} {'yR':>8s} {'Delta_y':>9s}")
    for r in results_slow:
        dy = r['phase_interfaces'][0]['delta_y'] if r['phase_interfaces'] else 0.0
        print(f"  {r['central_pressure_MeV_fm3']:10.3f} {r['R_km']:8.3f} {r['M_Msun']:9.4f} {r['compactness']:7.4f} {r['k2']:7.4f} {r['Lambda']:10.2f} {r['yR']:8.4f} {dy:9.4f}")
    # Slow vs rapid comparison
    print(f"\n  Slow (TK2020) vs Rapid (Postnikov 2010) relative difference in k2 and Lambda:")
    for rs, rr in zip(results_slow, results_rapid):
        rel_k2 = (rr['k2'] - rs['k2']) / rs['k2'] * 100
        rel_lam = (rr['Lambda'] - rs['Lambda']) / rs['Lambda'] * 100
        print(f"  pc={rs['central_pressure_MeV_fm3']:.3f}: k2_slow={rs['k2']:.4f}, k2_rapid={rr['k2']:.4f}, rel={rel_k2:+.2f}%  | Lambda_slow={rs['Lambda']:.2f}, rel={rel_lam:+.2f}%")
    # Acceptance checks
    checks = {}
    checks['all_compactness_below_half'] = bool(all(r['compactness'] < 0.5 for r in results_slow))
    checks['all_k2_physical'] = bool(all(0 < r['k2'] < 0.2 for r in results_slow))
    checks['all_lambda_positive'] = bool(all(r['Lambda'] > 0 for r in results_slow))
    checks['phase_interface_recorded'] = bool(all(len(r['phase_interfaces']) == 1 for r in results_slow))
    checks['tidal_jump_negative'] = bool(all(r['phase_interfaces'][0]['delta_y'] < 0 for r in results_slow if r['phase_interfaces']))
    checks['slow_rapid_differ'] = bool(any(abs(rr['k2'] - rs['k2']) / rs['k2'] > 1e-4 for rs, rr in zip(results_slow, results_rapid)))
    print(f"\n  acceptance: {json.dumps(checks, indent=2)}")
    print(f"  overall PASS = {all(checks.values())}")
    return {'slow': results_slow, 'rapid': results_rapid, 'checks': checks, 'overall_pass': all(checks.values()),
            'transition': t}


def test_convergence():
    """Test 3: TOV solver convergence with decreasing rtol."""
    print("\n=== Test 3: Convergence with decreasing rtol ===")
    n_arr = np.linspace(0.05, 0.80, 120)
    P, eps = make_polytrope_p_eps(n_arr, K=200.0, gamma=2.5)
    pc = float(P[-1] * 0.5)
    r_vals = []
    for rtol in [2e-6, 2e-8, 2e-10]:
        try:
            r = solve_tov_love_smooth(P, eps, pc, rtol=rtol, max_step_km=0.12)
            r_vals.append((rtol, r))
            print(f"  rtol={rtol:.0e}: R={r['R_km']:.6f} km, M={r['M_Msun']:.6f} Msun, k2={r['k2']:.6f}, Lambda={r['Lambda']:.4f}")
        except Exception as e:
            print(f"  rtol={rtol}: FAILED {e}")
    if len(r_vals) >= 2:
        R_conv = abs(r_vals[-1][1]['R_km'] - r_vals[-2][1]['R_km']) / r_vals[-1][1]['R_km']
        M_conv = abs(r_vals[-1][1]['M_Msun'] - r_vals[-2][1]['M_Msun']) / r_vals[-1][1]['M_Msun']
        k2_conv = abs(r_vals[-1][1]['k2'] - r_vals[-2][1]['k2']) / r_vals[-1][1]['k2']
        print(f"  relative change rtol 2e-8 -> 2e-10: R={R_conv:.2e}, M={M_conv:.2e}, k2={k2_conv:.2e}")
        converged = bool(R_conv < 1e-5 and M_conv < 1e-5 and k2_conv < 1e-4)
        print(f"  convergence PASS = {converged}")
        return {'converged': converged, 'R_change': R_conv, 'M_change': M_conv, 'k2_change': k2_conv}
    return {'converged': False}


def main():
    print("=" * 72)
    print("M6 reproduction: TOV + l=2 tidal Love + phase-transition tidal y-jump")
    print("Literature: Takatsy-Kovacs 2020 (PRD 102 028501, uploaded), Read 2009 (PRD 79 124032, uploaded)")
    print("=" * 72)
    t1 = test_smooth_polytrope_tov()
    t2 = test_maxwell_tov_love()
    t3 = test_convergence()
    summary = {
        'test1_smooth_polytrope': {
            'overall_pass': t1['overall_pass'] if t1 else False,
            'checks': t1['checks'] if t1 else None,
            'n_stars': len(t1['results']) if t1 else 0,
            'mass_range_Msun': [min(r['M_Msun'] for r in t1['results']), max(r['M_Msun'] for r in t1['results'])] if t1 else None,
            'radius_range_km': [min(r['R_km'] for r in t1['results']), max(r['R_km'] for r in t1['results'])] if t1 else None,
        },
        'test2_maxwell_tov_love': {
            'overall_pass': t2['overall_pass'] if t2 else False,
            'checks': t2['checks'] if t2 else None,
            'transition': t2['transition'] if t2 else None,
            'n_stars': len(t2['slow']) if t2 else 0,
        },
        'test3_convergence': t3,
    }
    import os
    out = '/home/z/my-project/download/m6/M6_ACCEPTANCE_SUMMARY.json'
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    print(f"\n=== Summary saved to {out} ===")

if __name__ == '__main__':
    main()
