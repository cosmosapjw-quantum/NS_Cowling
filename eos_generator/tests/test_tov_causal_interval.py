"""Certified causal domain of the interpolant the TOV right-hand side integrates.

The TOV RHS evaluates epsilon(P) through the log-log PCHIP ln eps(ln P).  The
former ``find_causal_domain`` sampled c_s^2 from a *different* PCHIP,
ln P(ln eps), at 16 points per interval.  Two counterexamples were declared
FULL_TABLE_CAUSAL:

* SAMPLE_MISS: sampled max c_s^2 = 1 + 5e-9 (gate 1 + 1e-8), but the sampled
  interpolant itself reaches 1 + 9.2e-5 and the RHS interpolant 1 + 9.4e-4;
* MISMATCH: the sampled interpolant is causal (max 0.9995) while the RHS
  interpolant reaches 1 + 4.6e-4.
"""
import numpy as np
import pytest

from compose_eos.tov import BarotropicEOS, find_causal_domain

SAMPLE_MISS_A = 0.01973334388796515
MISMATCH_A = -0.02234662920626299
TOL = 1e-8


def _table(a, k=12, cs2=0.9, bump=6):
    eps = np.geomspace(200.0, 2000.0, k)
    p = cs2 * (eps - 150.0)
    p[bump] *= 1.0 + a
    h = (eps + p)
    nb = 0.16 * np.exp(np.cumsum(np.r_[0.0, np.diff(eps) / ((h[1:] + h[:-1]) / 2.0)]))
    return BarotropicEOS(nb, p, eps)


def _dense_rhs_cs2(eos, lo=None, hi=None, points=400_001):
    lo = np.log(eos.pressure[0]) if lo is None else lo
    hi = np.log(eos.pressure[-1]) if hi is None else hi
    x = np.linspace(lo, hi, points)
    z = eos._lneps_of_lnp
    return x, np.exp(x - z(x)) / z.derivative()(x)


@pytest.mark.parametrize("a", [SAMPLE_MISS_A, MISMATCH_A])
def test_former_full_table_counterexamples_are_cut(a):
    eos = _table(a)
    _, cs2 = _dense_rhs_cs2(eos)
    assert cs2.max() > 1.0 + 1e-4
    domain = find_causal_domain(eos)
    assert domain.status == "ACASUAL_CUTOFF"
    assert domain.pressure_max < eos.pressure[-1]


@pytest.mark.parametrize("a", [SAMPLE_MISS_A, MISMATCH_A, 0.03, -0.03])
def test_rhs_interpolant_is_causal_below_the_cutoff(a):
    eos = _table(a)
    domain = find_causal_domain(eos)
    x, cs2 = _dense_rhs_cs2(eos, hi=np.log(domain.pressure_max))
    assert cs2.max() <= 1.0 + TOL
    assert abs(float(eos.cs2_of_pressure(domain.pressure_max)) - 1.0) <= 1e-9
    assert abs(domain.cs2_at_limit - 1.0) <= 1e-9
    # epsilon at the limit comes from the same RHS interpolant
    assert abs(domain.eps_max - float(eos.energy_of_pressure(domain.pressure_max))) <= 1e-9 * domain.eps_max


def test_certified_full_table_bounds_every_dense_sample():
    eos = _table(0.0)
    domain = find_causal_domain(eos)
    _, cs2 = _dense_rhs_cs2(eos)
    assert domain.status == "FULL_TABLE_CAUSAL"
    assert "exact interval bounds" in domain.reason
    assert cs2.max() <= 1.0 + TOL
    assert cs2.max() <= domain.cs2_certified_max <= 1.0 + TOL


def test_cs2_of_pressure_matches_rhs_derivative():
    eos = _table(0.01)
    p = np.geomspace(eos.pressure[0], eos.pressure[-1], 37)[1:-1]
    eps = eos.energy_of_pressure(p)
    dp = 1e-6 * p
    fd = 2.0 * dp / (eos.energy_of_pressure(p + dp) - eos.energy_of_pressure(p - dp))
    assert np.allclose(eos.cs2_of_pressure(p), fd, rtol=1e-6)
    assert np.all(eps > 0)
