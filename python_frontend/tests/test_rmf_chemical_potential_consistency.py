"""mu_n and mu_p of the homogeneous RMF solver must equal d eps/d n_n, d eps/d n_p.

Counterexample (SPG_GOGELEIN): the energy density used the no-1/2 isovector
convention (b = g_rho rho_3/m_rho^2, which reproduces J = 32.1 MeV), while the
chemical potentials and rearrangement term used an explicit 1/2; mu_n differed
from d eps/d n by 4-44 MeV in pure neutron matter, so the CLD cell
beta-equilibrium and the gas chemical potential were wrong for every SPG crust.
"""
from pathlib import Path

import pytest

from compose_eos.rmf_solver import load_models, solve_homogeneous

ROOT = Path(__file__).resolve().parents[2]
MODELS = load_models(ROOT / "data" / "eos390_rmf_parameter_sets.v1.json")
NAMES = [k for k in ("SPG_M1", "SPG_M3", "SPG_M5", "GPPVA_DD2", "GPPVA_DDME2", "GPPVA_TW")
         if k in MODELS]


def _eps(model, n_n, n_p):
    return solve_homogeneous(model, n_n + n_p, n_p).energy_density


@pytest.mark.parametrize("name", NAMES)
@pytest.mark.parametrize("n, yp", [(0.02, 0.0), (0.08, 0.05), (0.16, 0.5), (0.3, 0.1)])
def test_chemical_potentials_are_energy_derivatives(name, n, yp):
    model = MODELS[name]
    n_p, n_n = yp * n, (1.0 - yp) * n
    r = solve_homogeneous(model, n, n_p)
    h = 1e-6 * n
    mu_n_fd = (_eps(model, n_n + h, n_p) - _eps(model, n_n - h, n_p)) / (2 * h)
    assert r.mu_n == pytest.approx(mu_n_fd, abs=2e-4)
    if yp > 0:
        mu_p_fd = (_eps(model, n_n, n_p + h) - _eps(model, n_n, n_p - h)) / (2 * h)
        assert r.mu_p == pytest.approx(mu_p_fd, abs=2e-4)
