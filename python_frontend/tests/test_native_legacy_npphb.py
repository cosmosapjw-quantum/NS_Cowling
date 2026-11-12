"""Native route for the historical ``maxwell`` key = LEGACY_NPPHB (EOS-05b).

Contract (data/ns_cowling_cloud_agent_dag.v1.json, method_ontology.historical_maxwell):
the historical ``maxwell`` token is the project-defined smooth near-pressure-plateau
Hermite bridge (NPPHB): continuous density, no latent heat, no common-P/common-mu
coexistence. It is NOT a sharp first-order Maxwell transition.

The native backend exposes this construction under the internal kernel name
``legacy_npphb``. It is not a matrix token. The bare API/CLI token ``maxwell``
stays rejected (see test_native_crust_core_matching.py) and ``sharp_maxwell`` is
a separate operation that is not touched here.

Parity reference: frozen ``*__maxwell.csv`` in data/eos_original_normalized, with the
same tolerances as the other native/legacy differential tests. Composition is dropped
(drop_composition=True), so this is characterization, not a production PASS.
"""
from __future__ import annotations

import subprocess
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest

from compose_eos import NativeFeatureUnavailable, executable, match_crust_core, read_compose

ROOT = Path(__file__).resolve().parents[2]
MODELS = ROOT / "eos_generator" / "eos_data_normalized"
FROZEN = ROOT / "data" / "eos_original_normalized"

KERNEL = "legacy_npphb"
LABEL = "legacy_npphb_smooth_plateau"  # <= 32 characters (Fortran len=32 fields)
CRUSTS = (
    "SPG_M1", "SPG_M2", "SPG_M3", "SPG_M4", "SPG_M5",
    "GPPVA_DD2", "GPPVA_DDME2", "GPPVA_FSU2", "GPPVA_FSU2H",
    "GPPVA_FSU2R", "GPPVA_NL3wrL55", "GPPVA_TM1e", "GPPVA_TW",
)
COMPARED_CORES = ("BBB_BHF2", "GM_GM1")
PAIRS = [(cr, co) for co in COMPARED_CORES for cr in CRUSTS]
RTOL, ATOL = 8e-7, 2e-9
N_BRIDGE = 40
COLUMNS = {
    "nb": "n_b_fm3", "Q1": "q1_p_over_nb_mev", "Q3": "q3_mu_b_over_mn_minus_1",
    "Q4": "q4_mu_q_over_mn", "Q6": "q6_f_over_nb_mn_minus_1", "Q7": "q7_e_over_nb_mn_minus_1",
}


@lru_cache(maxsize=None)
def _model(name: str):
    return read_compose(MODELS / name)


@lru_cache(maxsize=None)
def _native(crust: str, core: str):
    return match_crust_core(
        _model(crust), _model(core), cr_name=crust, co_name=core,
        method=KERNEL, n_bridge=N_BRIDGE, drop_composition=True,
    )


@pytest.mark.parametrize("crust,core", PAIRS)
def test_legacy_npphb_matches_frozen_historical_maxwell_output(crust, core):
    matched, info = _native(crust, core)
    path = FROZEN / f"{crust}__{core}__maxwell.csv"
    assert path.is_file(), f"frozen fixture missing: {path}"
    frozen = np.genfromtxt(path, delimiter=",", names=True)

    assert matched.n == len(frozen)
    for key, column in COLUMNS.items():
        np.testing.assert_allclose(
            getattr(matched, key), frozen[column], rtol=RTOL, atol=ATOL,
            err_msg=f"{crust}x{core} legacy_npphb: native/frozen mismatch in {key}",
        )


@pytest.mark.parametrize("crust,core", PAIRS)
def test_legacy_npphb_is_labelled_smooth_and_never_a_sharp_transition(crust, core):
    matched, info = _native(crust, core)
    assert info["method"] == KERNEL
    assert info["construction_kind"] == LABEL
    assert info["spline_mode"] == LABEL
    assert info["is_exact_maxwell"] is False
    assert info["seidov_applicable"] is False
    assert info["density_jump_fm3"] == 0.0
    assert info["phase_transition_order"] is None
    assert info["composition_dropped"] is bool(
        _model(crust).composition_file_present or _model(crust).composition_available
    )
    # Continuous density: no latent heat, strictly increasing monotone barotrope.
    assert np.all(np.diff(matched.nb) > 0.0)
    assert np.all(np.diff(matched.p) > 0.0)
    assert np.all(np.diff(matched.eps) > 0.0)


def test_label_fits_the_32_character_fortran_fields_without_truncation():
    assert len(LABEL) <= 32


def test_bare_maxwell_token_remains_rejected_at_the_api():
    with pytest.raises(NativeFeatureUnavailable, match="not a first-order transition"):
        match_crust_core(_model("GPPVA_FSU2H"), _model("GM_GM1"), method="maxwell",
                         drop_composition=True)


def test_native_cli_accepts_the_internal_kernel_name_but_not_the_bare_token(tmp_path):
    cdir, kdir = tmp_path / "crust", tmp_path / "core"
    for d in (cdir, kdir):
        d.mkdir()
    _model("GPPVA_FSU2H").write_internal_directory(cdir)
    _model("GM_GM1").write_internal_directory(kdir)

    ok_out = tmp_path / "ok"
    ok_out.mkdir()
    accepted = subprocess.run(
        [str(executable()), "match-core", str(cdir), str(kdir), str(ok_out), KERNEL, "40", "true"],
        capture_output=True, text=True, check=False,
    )
    assert accepted.returncode == 0, accepted.stdout
    assert f'"construction_kind":"{LABEL}"' in accepted.stdout
    assert f'"method":"{KERNEL}"' in accepted.stdout

    bad_out = tmp_path / "bad"
    bad_out.mkdir()
    rejected = subprocess.run(
        [str(executable()), "match-core", str(cdir), str(kdir), str(bad_out), "maxwell", "40", "true"],
        capture_output=True, text=True, check=False,
    )
    assert rejected.returncode != 0
    assert "maxwell requires a sharp coexistence construction" in rejected.stdout


@pytest.mark.parametrize("crust", CRUSTS)
def test_rg_sly4_core_still_rejected_by_native_guard_for_legacy_npphb(crust):
    """Same documented limitation as the smooth tokens (EOS-05a F-EOS05A-1)."""
    with pytest.raises(ValueError, match="native matching requires beta-equilibrated inputs"):
        match_crust_core(_model(crust), _model("RG_SLY4"), cr_name=crust, co_name="RG_SLY4",
                         method=KERNEL, n_bridge=N_BRIDGE, drop_composition=True)
