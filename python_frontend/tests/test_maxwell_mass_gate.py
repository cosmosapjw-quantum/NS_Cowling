"""Maxwell reference encoding must preserve physical EOS data.

Header mn/mp are not a license to fit an additive physical energy offset.
These analytic controls replace the former metadata-equality gate tests:
different encodings of the same P, epsilon and chemical potentials must give
the same coexistence; a physical shift remains in the source curve.
"""
from __future__ import annotations

import sys
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos import EOSData, NativeExecutionError, construct_maxwell_transition


def physical_pair(il=1):
    """Exactly quadratic P(mu), one upward density jump at mu=990 MeV."""
    mu = np.array([950., 960., 980., 1000., 1030., 1040.])
    x = mu - 990.
    phases = []
    for upper in (False, True):
        p = 8. + .12*x + .001*x*x + (0.01*x if upper else 0.)
        n = .12 + .002*x + (0.01 if upper else 0.)
        eps = mu*n-p
        zeros = np.zeros_like(n)
        q4 = np.full_like(n, -40./939.)
        q3 = mu/939.-1. - (.1*q4 if il == 0 else 0.)
        q7 = eps/(n*939.)-1.
        phases.append(EOSData(
            nb=n, mn=939., mp=939., Q1=p/n, Q2=zeros.copy(), Q3=q3,
            Q4=q4, Q5=zeros.copy(), Q6=q7.copy(), Q7=q7.copy(),
            q_extra=zeros.copy(), Iphase=np.zeros(n.size, dtype=int),
            Yn=zeros.copy(), Yp=zeros.copy(), Ye=zeros.copy(), Ymu=zeros.copy(),
            Aav=zeros.copy(), Zav=zeros.copy(), Yh=zeros.copy(),
            Il=il, iT=1, iYq=1, inb_start=1, temperature_mev=0.,
            charge_fraction_grid=.1, q_extra_count=0,
            composition_file_present=False, composition_available=False,
            beta_equilibrated=True))
    return phases


def reencode(data, mn, mp):
    """Analytic test-input transformation; not a production EOS calculation."""
    result = data.copy()
    ratio = data.mn/mn
    for name in ("Q3", "Q6", "Q7"):
        setattr(result, name, (getattr(data, name)+1.)*ratio-1.)
    for name in ("Q4", "Q5"):
        setattr(result, name, getattr(data, name)*ratio)
    result.mn, result.mp = mn, mp
    return result


@pytest.mark.parametrize("il", [0, 1])
@pytest.mark.parametrize("mn,mp", [
    (939., 939.), (939.56536, 939.), (939., 938.27203),
    (939.56536, 938.27203), (939.5654, 938.272), (950., 900.)])
def test_reference_encoding_does_not_change_physical_maxwell(il, mn, mp):
    low, high = physical_pair(il)
    original = {name:getattr(high, name).copy() for name in ("Q1","Q3","Q4","Q6","Q7")}
    reference, _ = construct_maxwell_transition(
        low, high, assume_local_charge_neutrality=True, return_report=True)
    encoded = reencode(high, mn, mp)
    out, report = construct_maxwell_transition(
        low, encoded, assume_local_charge_neutrality=True, return_report=True)
    for name in ("nb", "p", "eps", "mub", "muq"):
        np.testing.assert_allclose(getattr(out, name), getattr(reference, name),
                                   rtol=0., atol=2e-10)
    np.testing.assert_array_equal(out.Iphase, reference.Iphase)
    np.testing.assert_allclose(encoded.eps, high.eps, rtol=0., atol=2e-11)
    for name, before in original.items():
        np.testing.assert_array_equal(getattr(high, name), before)
    assert report["mu_transition"] == pytest.approx(990., abs=2e-10)
    assert report["physical_energy_shift_mev"] == 0.
    assert report["mass_normalization_semantics"] == "physical_energy_preserving_Q_reencoding"
    assert report["rest_energy_reference_mev"] == {
        "low":[low.mn, low.mp], "high":[mn, mp], "output":[out.mn, out.mp]}


@pytest.mark.parametrize("il", [0, 1])
def test_physical_energy_offset_is_not_removed_by_normalization(il):
    low, high = physical_pair(il)
    shifted = high.copy()
    for name in ("Q3", "Q6", "Q7"):
        setattr(shifted, name, getattr(high, name)+100./high.mn)
    np.testing.assert_allclose(shifted.eps-high.eps, 100.*high.nb, atol=2e-11, rtol=0.)
    with pytest.raises(NativeExecutionError, match="no overlapping baryon-chemical-potential"):
        construct_maxwell_transition(
            low, shifted, assume_local_charge_neutrality=True, return_report=True)


def test_header_only_neutron_change_is_not_equivalent_reencoding():
    _, high = physical_pair()
    changed = high.copy()
    changed.mn = 950.
    assert np.max(np.abs(changed.eps-high.eps)) > 1.
    encoded = reencode(high, 950., 938.27203)
    np.testing.assert_allclose(encoded.eps, high.eps, rtol=0., atol=2e-11)
