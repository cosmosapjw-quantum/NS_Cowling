"""Frozen-fixture parity of the native smooth crust-core bridges (EOS-05a).

Scope (characterization only; this is NOT a production PASS):
  * Native Fortran ``match-core`` (composition dropped, ``drop_composition=True``)
    versus the frozen historical outputs in ``data/eos_original_normalized``.
  * Tokens cubic_F, cubic_P, quintic_F, quintic_P. ``css`` is excluded (its
    identity is governed by the cloud DAG) and ``maxwell`` is covered separately.
  * Cores BBB_BHF2 and GM_GM1 (26 pairs, 104 cases) are compared point by point.
  * RG_SLY4 cores are rejected by the native guard ``beta_equilibrated`` (the
    source carries a constant non-zero Q5). The exact rejection is asserted so
    that any future change is a conscious one; nothing is relaxed here.

No legacy Python is executed: the frozen CSV files are the reference.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest

from compose_eos import read_compose, match_crust_core

ROOT = Path(__file__).resolve().parents[2]
MODELS = ROOT / "eos_generator" / "eos_data_normalized"
FROZEN = ROOT / "data" / "eos_original_normalized"

CRUSTS = (
    "SPG_M1", "SPG_M2", "SPG_M3", "SPG_M4", "SPG_M5",
    "GPPVA_DD2", "GPPVA_DDME2", "GPPVA_FSU2", "GPPVA_FSU2H",
    "GPPVA_FSU2R", "GPPVA_NL3wrL55", "GPPVA_TM1e", "GPPVA_TW",
)
COMPARED_CORES = ("BBB_BHF2", "GM_GM1")
REJECTED_CORE = "RG_SLY4"
TOKENS = ("cubic_F", "cubic_P", "quintic_F", "quintic_P")
N_BRIDGE = 40
RTOL, ATOL = 8e-7, 2e-9  # same tolerances as the live differential test

COLUMNS = {
    "nb": "n_b_fm3",
    "Q1": "q1_p_over_nb_mev",
    "Q3": "q3_mu_b_over_mn_minus_1",
    "Q4": "q4_mu_q_over_mn",
    "Q6": "q6_f_over_nb_mn_minus_1",
    "Q7": "q7_e_over_nb_mn_minus_1",
}
SPLINE_MODES = {
    "exact_single_cubic", "monotone_1knot_cubic",
    "exact_single_quintic", "adaptive_2knot_c2_hermite",
}
RG_REJECTION = "native matching requires beta-equilibrated inputs"

PAIRS = [(cr, co) for co in COMPARED_CORES for cr in CRUSTS]
CASES = [(cr, co, m) for (cr, co) in PAIRS for m in TOKENS]


@lru_cache(maxsize=None)
def _model(name: str):
    return read_compose(MODELS / name)


@lru_cache(maxsize=None)
def _native(crust: str, core: str, method: str):
    return match_crust_core(
        _model(crust), _model(core), cr_name=crust, co_name=core,
        method=method, n_bridge=N_BRIDGE, drop_composition=True,
    )


@lru_cache(maxsize=None)
def _frozen(crust: str, core: str, method: str):
    path = FROZEN / f"{crust}__{core}__{method}.csv"
    assert path.is_file(), f"frozen fixture missing: {path}"
    return np.genfromtxt(path, delimiter=",", names=True)


def test_fixture_inventory_is_complete():
    """13 crusts x 3 cores x 6 tokens are tracked; absence is a failure, not a skip."""
    for core in COMPARED_CORES + (REJECTED_CORE,):
        for crust in CRUSTS:
            for token in TOKENS + ("css", "maxwell"):
                assert (FROZEN / f"{crust}__{core}__{token}.csv").is_file()
    for name in CRUSTS + COMPARED_CORES + (REJECTED_CORE,):
        assert (MODELS / name / "eos.thermo").is_file()


@pytest.mark.parametrize("crust,core,method", CASES)
def test_native_bridge_matches_frozen_fixture(crust, core, method):
    matched, info = _native(crust, core, method)
    frozen = _frozen(crust, core, method)

    assert info["spline_mode"] in SPLINE_MODES
    assert info["method"] == method
    assert info["composition_dropped"] is (
        bool(_model(crust).composition_file_present or _model(crust).composition_available)
    )
    assert info["is_exact_maxwell"] is False
    assert matched.n == len(frozen)
    for key, column in COLUMNS.items():
        np.testing.assert_allclose(
            getattr(matched, key), frozen[column], rtol=RTOL, atol=ATOL,
            err_msg=f"{crust}x{core} {method}: native/frozen mismatch in {key}",
        )
    assert np.all(np.diff(matched.nb) > 0.0)
    assert np.all(np.diff(matched.p) > 0.0)


@pytest.mark.parametrize("crust,core", PAIRS)
def test_negative_control_fixture_comparison_discriminates_methods(crust, core):
    """Native cubic_F must NOT reproduce the frozen quintic_F bridge rows."""
    matched, info = _native(crust, core, "cubic_F")
    frozen_other = _frozen(crust, core, "quintic_F")
    i1 = info["i1"]
    bridge = slice(i1 + 1, i1 + 1 + N_BRIDGE)
    assert not np.allclose(
        matched.Q1[bridge], frozen_other[COLUMNS["Q1"]][bridge], rtol=RTOL, atol=ATOL
    )


def test_all_spline_branches_are_exercised():
    seen = {_native(cr, co, m)[1]["spline_mode"] for (cr, co, m) in CASES}
    assert seen == SPLINE_MODES


@pytest.mark.parametrize("crust", CRUSTS)
def test_rg_sly4_core_is_rejected_by_native_guard_with_exact_message(crust):
    """Documented limitation, not parity: the RG source has a constant Q5 != 0."""
    core = _model(REJECTED_CORE)
    assert core.beta_equilibrated is False
    assert float(np.min(core.Q5)) == float(np.max(core.Q5)) != 0.0
    for method in TOKENS:
        with pytest.raises(ValueError, match=RG_REJECTION):
            match_crust_core(
                _model(crust), core, cr_name=crust, co_name=REJECTED_CORE,
                method=method, n_bridge=N_BRIDGE, drop_composition=True,
            )
