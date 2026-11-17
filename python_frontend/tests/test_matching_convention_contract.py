"""Matching cannot silently change the thermodynamic content of a source slice."""
from __future__ import annotations

import numpy as np
import pytest

from compose_eos import EOSData, api


def _branch(lo: float, hi: float) -> EOSData:
    n = np.linspace(lo, hi, 20)
    z = np.zeros_like(n)
    # Exact cold polytrope: epsilon = m*n + K*n**2, P = K*n**2.
    q = 100.0 * n
    f = 939.0 + q
    return EOSData(
        nb=n, mn=939.0, mp=938.0, Q1=q, Q2=z.copy(),
        Q3=(f + q) / 939.0 - 1.0, Q4=z.copy(), Q5=z.copy(),
        Q6=f / 939.0 - 1.0, Q7=f / 939.0 - 1.0,
        q_extra=z.copy(), Iphase=np.zeros(n.size, dtype=int),
        **{name: z.copy() for name in ("Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh")},
        beta_equilibrated=True,
    )


@pytest.mark.parametrize("method", ["css", "cubic_F", "cubic_P", "quintic_F", "quintic_P"])
@pytest.mark.parametrize("invalid", ["lepton_content", "charge_fraction", "nan_charge", "infinite_charge"])
def test_public_matching_rejects_incompatible_thermodynamic_slices(method, invalid, monkeypatch):
    crust, core = _branch(.01, .04), _branch(.06, .15)
    if invalid == "lepton_content":
        core.Il = 0
    elif invalid == "charge_fraction":
        core.charge_fraction_grid = .1
    elif invalid == "nan_charge":
        core.charge_fraction_grid = float("nan")
    else:
        core.charge_fraction_grid = float("inf")

    def unexpected_native_call(*args):
        raise AssertionError("incompatible source slice entered native matching")

    monkeypatch.setattr(api, "run_native", unexpected_native_call)
    with pytest.raises(ValueError, match="lepton-content|charge-fraction"):
        api.match_crust_core(crust, core, method=method, n_bridge=8)
