from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from compose_eos import (
    NativeExecutionError,
    NativeFeatureUnavailable,
    executable,
    match_crust_core,
    read_compose,
)

ROOT = Path(__file__).resolve().parents[2]
CRUST = ROOT / "eos_generator" / "eos_data_normalized" / "GPPVA_FSU2H"
CORE = ROOT / "eos_generator" / "eos_data_normalized" / "GM_GM1"
LEGACY_ROOT = ROOT / "eos_generator" / "compose_eos_pkg"
METHODS = ("css", "cubic_F", "cubic_P", "quintic_F", "quintic_P")


def _require_reference_inputs():
    required = (CRUST / "eos.nb", CRUST / "eos.thermo", CRUST / "eos.compo",
                CORE / "eos.nb", CORE / "eos.thermo")
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        pytest.skip("repository benchmark EOS inputs are unavailable: " + ", ".join(missing))


def _legacy_match(method: str, tmp_path: Path):
    _require_reference_inputs()
    npz = tmp_path / f"legacy_{method}.npz"
    metadata = tmp_path / f"legacy_{method}.json"
    script = r'''
import json, os, sys, types
from pathlib import Path
import numpy as np

# The legacy package __init__ eagerly imports optional plotting backends.
# Install a minimal package namespace so the reference reader/matcher can be
# imported without matplotlib; their numerical dependencies are NumPy/SciPy.
legacy_pkg = types.ModuleType("compose_eos")
legacy_pkg.__path__ = [str(Path(os.environ["LEGACY_ROOT"]) / "compose_eos")]
sys.modules["compose_eos"] = legacy_pkg
from compose_eos.io import read_compose
from compose_eos.matching import match_crust_core

cr = read_compose(sys.argv[1])
co = read_compose(sys.argv[2])
matched, info = match_crust_core(
    cr, co, cr_name="GPPVA_FSU2H", co_name="GM_GM1",
    method=sys.argv[3], n_bridge=40, core_compo_mode="strict_thermo_only")
np.savez(sys.argv[4], nb=matched.nb, Q1=matched.Q1, Q3=matched.Q3,
         Q4=matched.Q4, Q6=matched.Q6, Q7=matched.Q7, Iphase=matched.Iphase)
with open(sys.argv[5], "w", encoding="utf-8") as stream:
    json.dump(info, stream, allow_nan=False)
'''
    env = os.environ.copy()
    env["PYTHONPATH"] = str(LEGACY_ROOT)
    env["LEGACY_ROOT"] = str(LEGACY_ROOT)
    completed = subprocess.run(
        [sys.executable, "-c", script, str(CRUST), str(CORE), method, str(npz), str(metadata)],
        cwd=ROOT, env=env, capture_output=True, text=True, check=False,
    )
    assert completed.returncode == 0, (
        f"legacy reference failed for {method}: stdout={completed.stdout}\n"
        f"stderr={completed.stderr}"
    )
    return np.load(npz), json.loads(metadata.read_text(encoding="utf-8"))


@pytest.mark.parametrize("method", METHODS)
def test_native_smooth_matching_differential_against_legacy(method, tmp_path):
    _require_reference_inputs()
    legacy, reference_info = _legacy_match(method, tmp_path)
    crust = read_compose(CRUST)
    core = read_compose(CORE)
    matched, info = match_crust_core(
        crust, core, cr_name="GPPVA_FSU2H", co_name="GM_GM1",
        method=method, n_bridge=40, core_compo_mode="strict_thermo_only",
        drop_composition=True,
    )

    assert info["i1"] == reference_info["i1"]
    assert info["i2"] == reference_info["i2"]
    assert matched.n == legacy["nb"].size
    for key in ("nb", "Q1", "Q3", "Q4", "Q6", "Q7"):
        np.testing.assert_allclose(
            getattr(matched, key), legacy[key], rtol=8e-7, atol=2e-9,
            err_msg=f"{method}: native/legacy mismatch in {key}",
        )

    # Source pressure/energy knots outside the selected bridge are retained.
    i1, i2 = info["i1"], info["i2"]
    np.testing.assert_array_equal(matched.nb[:i1 + 1], crust.nb[:i1 + 1])
    np.testing.assert_array_equal(matched.nb[-(core.n - i2):], core.nb[i2:])
    assert matched.n == (i1 + 1) + 40 + (core.n - i2)

    assert np.all(np.diff(matched.nb) > 0.0)
    assert np.all(np.diff(matched.p) > 0.0)
    assert np.all(np.diff(matched.eps) > 0.0)
    assert np.all(matched.Q1 > 0.0)
    assert info["native_report"]["cs2_min"] > 0.0
    assert info["native_report"]["cs2_max"] < 1.0
    assert info["first_law_relative_residual"] < 5e-9
    np.testing.assert_allclose(matched.Q6, matched.Q7, rtol=0.0, atol=1e-14)

    # This benchmark has an input crust composition table. Its omission is a
    # deliberate, explicitly requested approximation of the native v1 API.
    assert info["composition_dropped"] is True
    assert matched.composition_file_present is False
    assert matched.composition_available is False
    assert np.all(matched.Iphase == 0)
    for field in (matched.Yn, matched.Yp, matched.Ye, matched.Ymu,
                  matched.Aav, matched.Zav, matched.Yh):
        assert np.all(field == 0.0)
    assert info["is_exact_maxwell"] is False
    assert info["density_jump_fm3"] == 0.0
    assert "consistency" in info and "overall" in info["consistency"]


def test_sharp_maxwell_method_is_rejected_by_smooth_matcher_until_m5_integration():
    crust = read_compose(CRUST)
    core = read_compose(CORE)
    with pytest.raises(NativeFeatureUnavailable, match="not a first-order transition"):
        match_crust_core(
            crust, core, method="maxwell", drop_composition=True
        )


def test_native_cli_keeps_sharp_maxwell_fail_closed_until_m5_integration(tmp_path):
    crust = read_compose(CRUST)
    core = read_compose(CORE)
    cdir, kdir, out = tmp_path / "crust", tmp_path / "core", tmp_path / "out"
    cdir.mkdir(); kdir.mkdir(); out.mkdir()
    crust.write_internal_directory(cdir)
    core.write_internal_directory(kdir)
    result = subprocess.run(
        [str(executable()), "match-core", str(cdir), str(kdir), str(out), "maxwell", "40", "true"],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode != 0
    assert "maxwell is the sharp coexistence method" in result.stdout


def test_native_matching_requires_explicit_composition_drop():
    crust = read_compose(CRUST)
    core = read_compose(CORE)
    with pytest.raises(ValueError, match="drop_composition=True"):
        match_crust_core(crust, core, method="quintic_P")
    with pytest.raises(ValueError, match="drop_composition=True"):
        match_crust_core(core, crust, method="quintic_P")


def test_native_matching_rejects_unported_literature_composition_mode():
    crust = read_compose(CRUST)
    core = read_compose(CORE)
    with pytest.raises(NativeFeatureUnavailable, match="strict_thermo_only"):
        match_crust_core(
            crust, core, core_compo_mode="literature_beta_eq", drop_composition=True
        )


def test_native_matching_rejects_invalid_controls():
    crust = read_compose(CRUST)
    core = read_compose(CORE)
    with pytest.raises(ValueError, match="method"):
        match_crust_core(crust, core, method="maxwell_jump", drop_composition=True)
    with pytest.raises(ValueError, match="n_bridge"):
        match_crust_core(crust, core, n_bridge=7, drop_composition=True)
    with pytest.raises(ValueError, match="n_bridge"):
        match_crust_core(crust, core, n_bridge=True, drop_composition=True)


@pytest.mark.parametrize("mutation,match", [
    ("entropy", "zero entropy per baryon Q2"),
    ("q6q7", "Q6 and Q7"),
    ("beta", "beta-equilibrated"),
])
def test_native_matching_rejects_nonmatching_thermodynamic_state(mutation, match):
    crust = read_compose(CRUST)
    core = read_compose(CORE)
    if mutation == "entropy":
        crust = crust.copy()
        crust.Q2[0] = 1e-3
    elif mutation == "q6q7":
        crust = crust.copy()
        crust.Q6[0] += 1e-3
    else:
        crust = crust.copy()
        crust.beta_equilibrated = False
    with pytest.raises((NativeExecutionError, ValueError), match=match):
        match_crust_core(
            crust, core, method="quintic_P", drop_composition=True
        )
