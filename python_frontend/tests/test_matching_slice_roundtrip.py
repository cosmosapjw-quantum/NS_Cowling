"""Native matching preserves physical charge slices and their CompOSE encoding."""
from __future__ import annotations

import numpy as np
import pytest

from compose_eos import (
    EOSData, check_consistency, construct_maxwell_transition, match_crust_core,
    read_compose, write_compose,
)


def quadratic_source(nlo, nhi, count, *, rest_energy=939.0, stiffness=100.0):
    """Independent analytic reference: epsilon=m*n+K*n^2 and P=K*n^2."""
    n = np.linspace(nlo, nhi, count)
    pressure = stiffness * n**2
    epsilon = rest_energy * n + pressure
    mu = (epsilon + pressure) / n
    z = np.zeros_like(n)
    free = epsilon / (n * 939.0) - 1.0
    return EOSData(
        nb=n, mn=939.0, mp=938.0, Q1=pressure / n, Q2=z.copy(),
        Q3=mu / 939.0 - 1.0, Q4=z.copy(), Q5=z.copy(),
        Q6=free.copy(), Q7=free.copy(), q_extra=z.copy(),
        Iphase=np.zeros(n.size, dtype=int), Yn=z.copy(), Yp=z.copy(),
        Ye=z.copy(), Ymu=z.copy(), Aav=z.copy(), Zav=z.copy(), Yh=z.copy(),
        Il=1, iT=1, iYq=1, inb_start=1, temperature_mev=0.0,
        charge_fraction_grid=0.0, q_extra_count=0,
        composition_file_present=False, composition_available=False,
        beta_equilibrated=True,
    )


@pytest.mark.parametrize("method", ["cubic_F", "cubic_P", "quintic_F", "quintic_P"])
def test_il0_nonzero_charge_slice_and_baryon_potential_roundtrip(method, tmp_path):
    crust = quadratic_source(.01, .04, 20)
    core = quadratic_source(.06, .15, 30)
    reference, _ = match_crust_core(crust, core, method=method, drop_composition=True)
    for data, index in ((crust, 3), (core, 9)):
        data.Il = 0
        data.iYq = index
        data.charge_fraction_grid = .1
        data.Q4[:] = -.01
        data.Q3 -= data.charge_fraction_grid * data.Q4
    matched, report = match_crust_core(crust, core, method=method, drop_composition=True)
    assert matched.Il == 0
    assert matched.iYq == 9
    assert matched.charge_fraction_grid == pytest.approx(.1)
    assert report["diagnostic_only"] is False
    assert report["native_report"]["source_euler_residual_meV"] < 2e-10
    np.testing.assert_allclose(matched.eps, reference.eps, rtol=2e-13, atol=2e-12)
    np.testing.assert_allclose(matched.p, reference.p, rtol=2e-13, atol=2e-12)
    mu_euler = (matched.eps + matched.p) / matched.nb
    np.testing.assert_allclose(
        matched.mub, mu_euler - .1 * matched.Q4 * matched.mn,
        rtol=2e-13, atol=2e-12,
    )
    left, right = report["i1"] + 1, core.n - report["i2"]
    np.testing.assert_allclose(matched.Q3[:left], crust.Q3[:left], rtol=2e-13, atol=2e-14)
    np.testing.assert_allclose(matched.Q3[-right:], core.Q3[report["i2"]:], rtol=2e-13, atol=2e-14)
    target = tmp_path / "matched"
    write_compose(target, matched)
    restored = read_compose(target)
    assert restored.Il == matched.Il
    assert restored.iYq == matched.iYq
    assert restored.charge_fraction_grid == pytest.approx(.1)
    np.testing.assert_allclose(restored.Q3, matched.Q3, rtol=2e-13, atol=2e-14)
    np.testing.assert_allclose(restored.eps, matched.eps, rtol=2e-13, atol=2e-12)
    assert check_consistency(restored)["overall"]["pass_"] is True


def test_maxwell_same_physical_charge_slice_accepts_different_index_labels(tmp_path):
    low = quadratic_source(.3, 1.4, 81, stiffness=1.0)
    high = quadratic_source(1.3, 2.8, 81, rest_energy=940.0, stiffness=.25)
    reference, reference_report = construct_maxwell_transition(
        low, high, assume_local_charge_neutrality=True, return_report=True,
    )
    low.iYq, high.iYq = 3, 9
    matched, report = construct_maxwell_transition(
        low, high, assume_local_charge_neutrality=True, return_report=True,
    )
    assert report["status"] == "ok"
    assert matched.iYq == 3
    assert matched.charge_fraction_grid == 0.0
    np.testing.assert_allclose(matched.nb, reference.nb, rtol=2e-13, atol=2e-14)
    np.testing.assert_allclose(matched.eps, reference.eps, rtol=2e-13, atol=2e-12)
    assert report["pressure_transition"] == pytest.approx(reference_report["pressure_transition"], rel=2e-13)
    target = tmp_path / "maxwell"
    write_compose(target, matched)
    restored = read_compose(target)
    assert restored.iYq == 3
    assert restored.charge_fraction_grid == 0.0
    assert check_consistency(restored)["overall"]["pass_"] is True
