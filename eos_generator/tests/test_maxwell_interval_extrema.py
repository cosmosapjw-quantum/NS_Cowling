"""Exact interval-extrema admission of the declared Maxwell branch interpolant.

The Maxwell branch representation is n_B(mu_B) = PCHIP (piecewise cubic) and
P(mu_B) its antiderivative.  Mechanical stability and causality must hold on
that whole continuous representation.  A finite set of samples per interval
can miss a narrow acausal pocket; the counterexample below was admitted by the
former 8-samples-per-interval gate although the interpolant reaches
c_s^2 = 1 + 6.1e-6 (gate: 1 + 1e-6).
"""
import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator

from compose_eos import EOSData
from compose_eos.phase_transitions import (
    MAXWELL_CS2_TOLERANCE,
    _branchwise_first_law_report,
    _phase_arrays,
)

# Former sampled gate passed this perturbation at max c_s^2 = 1 + 5.0e-7.
POCKET_A = -0.0034416369635021592


def _branch(a=0.0, *, cs2=0.9, k=9, bump_index=4, mu0=1000.0, mu1=1400.0):
    mu = np.linspace(mu0, mu1, k)
    n = 0.3 * (mu / mu0) ** (1.0 / cs2)
    n[bump_index] *= 1.0 + a
    spline = PchipInterpolator(mu, n)
    prim = spline.antiderivative()
    p = 10.0 + prim(mu) - prim(mu[0])     # exact potential of the declared n_B
    return mu, n, p, mu * n - p


def _eos(mu, n, p, eps):
    f = eps / n
    z = np.zeros_like(n)
    return EOSData(nb=n, mn=939.0, mp=939.0, Q1=p / n, Q2=z, Q3=mu / 939.0 - 1.0,
                   Q4=z, Q5=z, Q6=f / 939.0 - 1.0, Q7=f / 939.0 - 1.0, q_extra=z,
                   Iphase=np.zeros(len(n), int), Yn=z, Yp=z, Ye=z, Ymu=z, Aav=z,
                   Zav=z, Yh=z, composition_available=False, beta_equilibrated=True)


def _dense(mu, n, points=400_001):
    s = PchipInterpolator(mu, n)
    q = np.linspace(mu[0], mu[-1], points)
    dn = s.derivative()(q)
    return float(np.max(s(q) / (q * dn))), float(np.min(dn))


def test_narrow_acausal_pocket_between_former_samples_is_rejected():
    mu, n, p, eps = _branch(POCKET_A)
    report = _branchwise_first_law_report(mu, p, n, eps)
    dense_max, _ = _dense(mu, n)
    assert dense_max > 1.0 + 6.0e-6
    assert report["method"] == "exact_interval_extrema"
    assert report["max_cs2_thermodynamic"] >= dense_max - 1e-12
    assert report["causality_pass"] is False
    with pytest.raises(ValueError, match="acausal"):
        _phase_arrays(_eos(mu, n, p, eps))


def test_causal_neighbour_is_admitted_with_exact_maximum():
    mu, n, p, eps = _branch(-0.0034)
    report = _branchwise_first_law_report(mu, p, n, eps)
    dense_max, _ = _dense(mu, n)
    assert report["causality_pass"] is True and report["stability_pass"] is True
    assert dense_max <= report["max_cs2_thermodynamic"] <= dense_max + 1e-9
    _phase_arrays(_eos(mu, n, p, eps))


@pytest.mark.parametrize("a", [-0.004, -0.002, 0.0, 0.002, 0.01])
def test_exact_extrema_bound_every_dense_sample(a):
    mu, n, p, eps = _branch(a)
    report = _branchwise_first_law_report(mu, p, n, eps)
    dense_max, dense_min_dn = _dense(mu, n)
    assert report["max_cs2_thermodynamic"] >= dense_max - 1e-12
    assert report["min_dn_dmu_fm3_per_MeV"] <= dense_min_dn + 1e-15
    assert report["max_cs2_thermodynamic"] - dense_max < 1e-8


def test_tangency_at_luminal_bound_is_admitted():
    # n_B proportional to mu_B gives c_s^2 = 1 identically; PCHIP reproduces
    # linear data exactly, so the bound is touched on the whole branch.
    mu = np.linspace(1000.0, 1400.0, 9)
    n = 3.0e-4 * mu
    prim = PchipInterpolator(mu, n).antiderivative()
    p = 10.0 + prim(mu) - prim(mu[0])
    report = _branchwise_first_law_report(mu, p, n, mu * n - p)
    assert abs(report["max_cs2_thermodynamic"] - 1.0) < 1e-12
    assert report["causality_pass"] is True


@pytest.mark.parametrize("excess, admitted", [(0.5e-6, True), (3.0e-6, False)])
def test_constant_speed_on_either_side_of_the_fixed_gate(excess, admitted):
    mu, n, p, eps = _branch(cs2=1.0 + excess, k=33)
    report = _branchwise_first_law_report(mu, p, n, eps)
    assert abs(report["max_cs2_thermodynamic"] - (1.0 + excess)) < 2e-7
    assert report["causality_pass"] is admitted
    assert MAXWELL_CS2_TOLERANCE == 1e-6


def test_compressibility_minimum_inside_interval_is_found():
    # Strongly unequal neighbouring secants push the PCHIP dn/dmu minimum into
    # the interval interior, away from both knots.
    mu = np.array([1000.0, 1010.0, 1020.0, 1030.0, 1200.0, 1210.0])
    n = np.array([0.30, 0.31, 0.32, 0.33, 0.331, 0.50])
    prim = PchipInterpolator(mu, n).antiderivative()
    p = 10.0 + prim(mu) - prim(mu[0])
    report = _branchwise_first_law_report(mu, p, n, mu * n - p)
    _, dense_min_dn = _dense(mu, n)
    knot_dn = PchipInterpolator(mu, n).derivative()(mu)
    assert report["min_dn_dmu_fm3_per_MeV"] <= dense_min_dn + 1e-15
    assert report["min_dn_dmu_fm3_per_MeV"] < np.min(knot_dn)
    assert mu[3] < report["argmin_dn_dmu_mu_MeV"] < mu[4]


def test_pocket_in_interval_ending_at_one_sided_jump_boundary_is_rejected():
    # The final interval terminates at the branch end (the coexistence row).
    mu, n, p, eps = _branch(POCKET_A, bump_index=7)
    report = _branchwise_first_law_report(mu, p, n, eps)
    dense_max, _ = _dense(mu, n)
    assert report["max_cs2_thermodynamic"] >= dense_max - 1e-12
    if dense_max > 1.0 + MAXWELL_CS2_TOLERANCE:
        assert report["causality_pass"] is False


def test_report_refuses_to_span_a_density_jump():
    mu, n, p, eps = _branch()
    mu2 = np.r_[mu[:5], mu[4], mu[5:] + 1e-9]
    n2 = np.r_[n[:5], n[4] * 1.1, n[5:] * 1.1]
    p2 = np.r_[p[:5], p[4], p[5:]]
    with pytest.raises(ValueError, match="strictly increasing"):
        _branchwise_first_law_report(mu2, p2, n2, mu2 * n2 - p2)


def test_vacuous_derivative_residual_is_not_reported_as_evidence():
    mu, n, p, eps = _branch()
    report = _branchwise_first_law_report(mu, p, n, eps)
    assert "max_relative_residual" not in report
    assert "by construction" in report["identity"]


@pytest.mark.parametrize("mu", [
    np.array([0.0, 1.0, 2.0, 3.0]),
    np.array([-2.0, -1.0, 0.0, 1.0]),
])
def test_branch_audit_rejects_nonpositive_chemical_potential(mu):
    n = np.array([0.1, 0.2, 0.3, 0.4])
    p = np.array([0.01, 0.02, 0.03, 0.04])
    with pytest.raises(ValueError, match="positive"):
        _branchwise_first_law_report(mu, p, n, mu * n - p)
