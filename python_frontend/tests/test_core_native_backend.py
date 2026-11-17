"""Native core generation: independent physics checks and research-oracle parity.

The Python research implementations are used only in tests, never called by
the Fortran production generator.
"""
from pathlib import Path
from functools import lru_cache
import shutil
import subprocess
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
MODELS = ("APR_APR", "BBB_BHF2", "GM_GM1", "RG_KDE0v1", "RG_SLY4")


@pytest.fixture(scope="module")
def native_core(tmp_path_factory):
    if shutil.which("gfortran") is None:
        pytest.skip("gfortran required for native core contract")
    build = tmp_path_factory.mktemp("native_core")
    common = ["gfortran", "-std=f2018", "-O2", "-ffree-line-length-none",
              "-fimplicit-none", "-fcheck=all", "-J" + str(build), "-I" + str(build),
              str(ROOT / "fortran_eos/src/compose_eos_types.f90"),
              str(ROOT / "fortran_eos/src/compose_eos_core_native.f90")]
    probe = build / "probe"
    subprocess.run(common + [str(ROOT / "fortran_eos/tests/core_native_probe.f90"),
                             "-o", str(probe)], check=True, capture_output=True, text=True)
    unit = build / "unit"
    subprocess.run(common + [str(ROOT / "fortran_eos/tests/test_core_native.f90"),
                             "-o", str(unit)], check=True, capture_output=True, text=True)
    subsaturation = build / "subsaturation"
    subprocess.run(common + [str(ROOT / "fortran_eos/tests/test_core_subsaturation.f90"),
                             "-o", str(subsaturation)], check=True, capture_output=True, text=True)
    return probe, unit, subsaturation


def native_rows(native_core, model, *args):
    run = subprocess.run([str(native_core[0]), model, *map(str, args)],
                         check=True, capture_output=True, text=True)
    return np.atleast_2d(np.loadtxt(run.stdout.splitlines()))


@lru_cache(maxsize=None)
def oracle(model):
    from compose_eos.core_regeneration import state_solver
    return state_solver(model)


def test_native_independent_first_law_conservation_and_gibbs(native_core):
    run = subprocess.run([str(native_core[1])], check=True, capture_output=True, text=True)
    assert "PASS native core functional, Gibbs, composition and generation tests" in run.stdout


def test_native_subsaturation_domain_and_independent_first_law(native_core):
    run = subprocess.run([str(native_core[2])], check=True, capture_output=True, text=True)
    assert "PASS native subsaturation domain and independent first-law checks" in run.stdout


def test_computational_support_is_distinct_from_default_sampling():
    from compose_eos.generation_dispatch import generation_request
    from compose_eos.native_model_contracts import MODEL_CONTRACTS
    for model in ("GM_GM1", "RG_KDE0v1", "RG_SLY4"):
        assert generation_request(model)["n_min"] == .08
        assert MODEL_CONTRACTS[model]["computational_interval_fm3"][0] == .02
        assert generation_request(model, n_min=.02)["n_min"] == .02
    assert MODEL_CONTRACTS["APR_APR"]["default_interval_fm3"][0] == .085
    assert MODEL_CONTRACTS["BBB_BHF2"]["default_interval_fm3"][0] == .08


@pytest.mark.parametrize("model", MODELS)
def test_native_states_match_declared_functionals(native_core, model):
    solve, _ = oracle(model)
    for density in (.16, .213, .47, .8):
        native = native_rows(native_core, model, density)[0]
        reference = solve(density)
        keys = ("n_b", "epsilon", "pressure", "mu_n", "mu_p", "mu_e",
                "Y_n", "Y_p", "Y_e", "Y_mu")
        expected = np.array([reference[key] for key in keys])
        np.testing.assert_allclose(native[:10], expected, rtol=5e-10, atol=2e-8)
        assert abs(native[7]-native[8]-native[9]) < 1e-10
        assert abs(native[1]+native[2]-native[0]*native[3]) < 2e-8


def test_native_apr_inserts_only_in_window_gibbs_endpoints(native_core):
    full = native_rows(native_core, "APR_APR", 3, .085, 1.34)
    _, context = oracle("APR_APR")
    assert len(full) == 5
    for side in ("left", "right"):
        assert min(abs(full[:, 0]-context[side]["n_b"])) < 1e-10
    narrow = native_rows(native_core, "APR_APR", 4, .5, .6)
    assert len(narrow) == 4
    assert np.all(narrow[:, 16] == 7)
    np.testing.assert_allclose(narrow[:, 0], np.geomspace(.5, .6, 4), rtol=2e-15)


@pytest.mark.parametrize("model", MODELS)
def test_native_requested_default_support_is_regenerable(native_core, model):
    from compose_eos.core_regeneration import MODEL_CONTRACTS
    lo, hi = MODEL_CONTRACTS[model]["default_interval_fm3"]
    rows = native_rows(native_core, model, 161, lo, hi)
    assert len(rows) == (163 if model == "APR_APR" else 161)
    assert np.all(np.isfinite(rows))
    assert np.all(np.diff(rows[:, 0]) > 0)
    assert np.max(abs(rows[:, 7]-rows[:, 8]-rows[:, 9])) < 1e-10
    assert np.max(abs(rows[:, 1]+rows[:, 2]-rows[:, 0]*rows[:, 3])) < 2e-7


def test_bbb_embedded_inputs_are_published_energies(native_core):
    import json
    source = json.loads((ROOT / "data/research_generator_20261116/"
                        "bbb2_published_tables.json").read_text())
    for density, symmetric, neutron in zip(
            source["n"], source["symmetric_EA"], source["neutron_EA"]):
        snm = native_rows(native_core, "BBB_BHF2", "nuclear", density/2, density/2, 1)[0]
        pnm = native_rows(native_core, "BBB_BHF2", "nuclear", density, 0, 1)[0]
        assert abs(snm[0]/density-939-symmetric) < 5e-13
        assert abs(pnm[0]/density-939-neutron) < 5e-13
