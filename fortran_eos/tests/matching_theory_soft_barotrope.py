"""Regression: decreasing P/n need not imply mechanical instability.

n is in fm^-3; F, q=P/n and mu are in MeV; epsilon and P are MeV/fm^3.
The reported s is c_s^2/c^2, not c_s^2 in dimensional units.
Run explicitly with pytest (the filename records ownership of this audit).
"""
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos import EOSData, match_crust_core


def soft_branch(n):
    """Exact cold barotrope: F=940+10*x-x^2/2, q=10-x, x=ln(n/.04)."""
    x = np.log(n / 0.04)
    q = 10.0 - x
    f = 940.0 + 10.0 * x - 0.5 * x**2
    mu = f + q
    z = np.zeros_like(n)
    return EOSData(
        nb=n, mn=939.0, mp=939.0,
        Q1=q, Q2=z.copy(), Q3=mu / 939.0 - 1.0,
        Q4=z.copy(), Q5=z.copy(), Q6=f / 939.0 - 1.0,
        Q7=f / 939.0 - 1.0, q_extra=z.copy(),
        Iphase=np.zeros(n.size, dtype=int),
        Yn=z.copy(), Yp=z.copy(), Ye=z.copy(), Ymu=z.copy(),
        Aav=z.copy(), Zav=z.copy(), Yh=z.copy(),
        Il=1, iT=1, iYq=1, inb_start=1,
        temperature_mev=0.0, charge_fraction_grid=0.0, q_extra_count=0,
        composition_file_present=False, composition_available=False,
        beta_equilibrated=True,
    )


@pytest.mark.parametrize("method", ["cubic_F", "quintic_F", "cubic_P"])
def test_decreasing_pressure_per_baryon_has_causal_bridge(method):
    crust = soft_branch(np.geomspace(0.005, 0.04, 41))
    core = soft_branch(np.linspace(0.045, 0.12, 61))
    # Independently known analytic stiffness is positive despite dq/dx=-1 MeV.
    for source in (crust, core):
        s = (source.Q1 - 1.0) / ((source.Q6 + 1.0) * source.mn + source.Q1)
        assert np.all((0.0 < s) & (s < 1.0))
        assert np.all(np.diff(source.Q1) < 0.0)
    matched, info = match_crust_core(
        crust, core, method=method, n_bridge=40, drop_composition=True
    )
    assert (info["i1"], info["i2"]) == (40, 0)
    r = info["native_report"]
    assert 0.0 < r["cs2_min"] <= r["cs2_max"] <= 1.0
    assert info["first_law_relative_residual"] < 5e-9
    assert np.all(np.diff(matched.nb) > 0.0)
    assert np.all(np.diff(matched.p) > 0.0)
    assert np.all(np.diff(matched.eps) > 0.0)
    # Both source endpoint states must survive the join.
    np.testing.assert_allclose(matched.p[[40, 81]], [crust.p[-1], core.p[0]], rtol=2e-12)
    np.testing.assert_allclose(matched.eps[[40, 81]], [crust.eps[-1], core.eps[0]], rtol=2e-12)
