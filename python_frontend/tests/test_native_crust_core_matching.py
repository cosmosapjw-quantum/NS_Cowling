from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from compose_eos import (
    EOSData,
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
METHODS = ("cubic_F", "cubic_P", "quintic_F", "quintic_P")


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
    assert matched.composition_file_present is True  # phase-only sidecar
    assert matched.composition_available is False
    i1, i2 = info["i1"], info["i2"]
    np.testing.assert_array_equal(matched.Iphase[:i1+1], crust.Iphase[:i1+1])
    np.testing.assert_array_equal(matched.Iphase[-(core.n-i2):], core.Iphase[i2:])
    assert np.all(matched.Iphase[i1+1:-(core.n-i2)] == info["native_report"]["bridge_phase_id"])
    for field in (matched.Yn, matched.Yp, matched.Ye, matched.Ymu,
                  matched.Aav, matched.Zav, matched.Yh):
        assert np.all(field == 0.0)
    assert info["is_exact_maxwell"] is False
    assert info["density_jump_fm3"] == 0.0
    assert "consistency" in info and "overall" in info["consistency"]


def test_native_css_is_single_speed_with_crust_side_first_order_jump():
    """CSS is one constant-c_s^2 bridge: continuous at core, Maxwell-jump at crust."""
    _require_reference_inputs()
    crust = read_compose(CRUST)
    core = read_compose(CORE)
    matched, info = match_crust_core(
        crust, core, cr_name="GPPVA_FSU2H", co_name="GM_GM1",
        method="css", n_bridge=40, core_compo_mode="strict_thermo_only",
        drop_composition=True,
    )

    report = info["native_report"]
    assert report["construction_kind"] == "css_single_speed_with_crust_density_jump"
    assert report["is_exact_maxwell"] is True
    assert report["seidov_applicable"] is True
    assert info["is_exact_maxwell"] is True
    assert info["phase_transition_order"] == "first_order"
    assert info["density_jump_fm3"] > 0.0

    # The CSS bridge has one and only one sound-speed value.
    cs2 = report["css_single_sound_speed"]
    assert 0.0 < cs2 <= 1.0
    assert report["cs2_left"] == pytest.approx(cs2, rel=1e-11, abs=1e-13)
    assert report["cs2_right"] == pytest.approx(cs2, rel=1e-11, abs=1e-13)
    assert report["cs2_min"] == pytest.approx(cs2, rel=1e-11, abs=1e-13)
    assert report["cs2_max"] == pytest.approx(cs2, rel=1e-11, abs=1e-13)
    assert info["dcs2_left"] is None
    assert info["dcs2_left_semantics"] == "not_applicable_across_first_order_density_jump"

    # The only sharp interface is crust -> bridge: P and mu are continuous,
    # while n and epsilon jump according to Delta epsilon = mu_t Delta n.
    i_low = report["css_transition_low_index"]
    i_high = report["css_transition_high_index"]
    assert i_high == i_low + 1
    assert matched.nb[i_high] > matched.nb[i_low]
    assert matched.p[i_high] == pytest.approx(matched.p[i_low], rel=1e-11, abs=1e-12)
    assert matched.mub[i_high] == pytest.approx(matched.mub[i_low], rel=1e-11, abs=1e-10)
    delta_n = matched.nb[i_high] - matched.nb[i_low]
    delta_eps = matched.eps[i_high] - matched.eps[i_low]
    assert delta_eps == pytest.approx(matched.mub[i_low] * delta_n, rel=2e-10, abs=1e-10)
    assert report["css_transition_epsilon_low"] == pytest.approx(matched.eps[i_low], rel=2e-10)
    assert report["css_transition_epsilon_high"] == pytest.approx(matched.eps[i_high], rel=2e-10)
    assert report["css_transition_delta_eps"] == pytest.approx(delta_eps, rel=2e-10)
    assert info["density_jump_fm3"] == pytest.approx(delta_n, rel=1e-12)
    assert info["seidov_applicable"] is True
    assert info["seidov_ratio"] == pytest.approx(report["seidov_ratio"], rel=1e-12)

    # The bridge joins continuously onto the selected core endpoint.
    i_core = report["css_core_bridge_index"]
    source_core_index = report["i2"]
    assert matched.nb[i_core] == pytest.approx(core.nb[source_core_index], rel=1e-12)
    assert matched.p[i_core] == pytest.approx(core.p[source_core_index], rel=2e-10)
    assert matched.eps[i_core] == pytest.approx(core.eps[source_core_index], rel=2e-8)
    mu_core_euler = (
        (core.Q6[source_core_index] + 1.0) * core.mn
        + core.Q1[source_core_index]
    )
    assert matched.mub[i_core] == pytest.approx(mu_core_euler, rel=2e-10)
    assert info["first_law_relative_residual"] < 5e-9


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


# ---------------------------------------------------------------------------
# §3.5 acceptance tests #9 and #10: root-solver failure / no-root case and
# a deliberately invalid case with nonpositive density jump. The native
# Fortran solve_single_css_match rejects both; these tests prove the rejection
# surfaces to the caller instead of returning a nominal success or a clipped
# solution.  (Tests #1-#8 and #11-#12 are already covered by
# test_native_css_is_single_speed_with_crust_side_first_order_jump above and
# test_apr_css_maxwell_audit.py; these two close the remaining gap.)
# ---------------------------------------------------------------------------


def _css_barotrope(*, rest_mass: float, stiffness: float, n_min: float,
                   n_max: float, count: int = 81) -> EOSData:
    """Exact cold barotrope P = K n^2, eps = m n + K n^2, mu = m + 2 K n.

    This is the canonical cold polytrope of index 2 used by the CSS tests; its
    thermodynamic identities hold to machine precision and cs^2 = dP/deps is
    monotone, so it is a well-behaved crust or core input for the single-speed
    CSS matcher.
    """
    nb = np.linspace(n_min, n_max, count)
    pressure = stiffness * nb ** 2
    epsilon = rest_mass * nb + pressure
    mu = (epsilon + pressure) / nb
    zeros = np.zeros_like(nb)
    q7 = epsilon / (nb * 939.0) - 1.0
    return EOSData(
        nb=nb, mn=939.0, mp=939.0,
        Q1=pressure / nb, Q2=zeros.copy(), Q3=mu / 939.0 - 1.0,
        Q4=zeros.copy(), Q5=zeros.copy(), Q6=q7.copy(), Q7=q7.copy(),
        q_extra=zeros.copy(), Iphase=np.zeros(nb.size, dtype=int),
        Yn=zeros.copy(), Yp=zeros.copy(), Ye=zeros.copy(), Ymu=zeros.copy(),
        Aav=zeros.copy(), Zav=zeros.copy(), Yh=zeros.copy(),
        Il=1, iT=1, iYq=1, inb_start=1, temperature_mev=0.0,
        charge_fraction_grid=0.0, q_extra_count=0,
        composition_file_present=False, composition_available=False,
        beta_equilibrated=True,
    )


def test_css_single_speed_rejects_endpoints_with_no_causal_solution():
    """§3.5 #9: when the crust/core endpoints admit no causal single-speed CSS
    solution, the native matcher must fail closed rather than clip a noncausal
    solution into the allowed interval or return a nominal success. The exact
    rejection point depends on where in the pipeline the causality violation is
    detected (endpoint sound-speed check, no-sign-change bracket, or noncausal
    solved root), so this test asserts failure with a causality-related message,
    not a specific substring. The key contract is that the matcher must NOT
    silently succeed on an acausal endpoint pair."""
    # A causal polytropic crust joined to an explicitly acausal cubic core
    # (P = K n^3 with stiffness so high that c_s^2 = dP/deps exceeds 1 at the
    # high-density end). The CSS matcher requires 0 < cs^2 <= 1 at the selected
    # core endpoint, so an acausal core endpoint is rejected.
    crust = _css_barotrope(rest_mass=939.0, stiffness=1.0, n_min=0.05, n_max=0.60)
    # acausal cubic core: cs^2 = 3 K n^2 / (m + 3 K n^2) -> 1 as n grows; with
    # stiffness large, cs^2 exceeds 1 well inside the grid.
    nb = np.linspace(0.30, 1.50, 81)
    stiffness = 4000.0
    rest_mass = 939.0
    pressure = stiffness * nb ** 3
    epsilon = rest_mass * nb + 0.5 * stiffness * nb ** 3
    mu = (epsilon + pressure) / nb
    zeros = np.zeros_like(nb)
    q7 = epsilon / (nb * 939.0) - 1.0
    core = EOSData(
        nb=nb, mn=939.0, mp=939.0,
        Q1=pressure / nb, Q2=zeros.copy(), Q3=mu / 939.0 - 1.0,
        Q4=zeros.copy(), Q5=zeros.copy(), Q6=q7.copy(), Q7=q7.copy(),
        q_extra=zeros.copy(), Iphase=np.zeros(nb.size, dtype=int),
        Yn=zeros.copy(), Yp=zeros.copy(), Ye=zeros.copy(), Ymu=zeros.copy(),
        Aav=zeros.copy(), Zav=zeros.copy(), Yh=zeros.copy(),
        Il=1, iT=1, iYq=1, inb_start=1, temperature_mev=0.0,
        charge_fraction_grid=0.0, q_extra_count=0,
        composition_file_present=False, composition_available=False,
        beta_equilibrated=True,
    )
    with pytest.raises((NativeExecutionError, ValueError)) as excinfo:
        match_crust_core(
            crust, core, cr_name="synth_crust", co_name="synth_core",
            method="css", n_bridge=40, core_compo_mode="strict_thermo_only",
            drop_composition=True,
        )
    # The rejection must cite a causality / sound-speed / no-causal-solution
    # failure; it must not be a silent success or an unrelated validation error.
    msg = str(excinfo.value).lower()
    assert ("causal" in msg or "sound speed" in msg or "no.*solution" in msg
            or "endpoint" in msg or "outside (0,1)" in msg or "outside (0, 1)" in msg), (
        f"rejection must cite causality/endpoint failure; got: {excinfo.value!r}")


def test_css_exhaustive_search_finds_positive_jump_in_overlapping_sources(tmp_path):
    """Density-domain overlap does not forbid another certified CSS window."""
    from compose_eos import write_compose
    from compose_eos.native import run_native
    crust = _css_barotrope(rest_mass=939.0, stiffness=1.0, n_min=0.05, n_max=0.50)
    core = _css_barotrope(rest_mass=939.0, stiffness=1.2, n_min=0.45, n_max=1.20)
    matched, info = match_crust_core(
        crust, core, method="css", n_bridge=40, drop_composition=True,
    )
    report = info["native_report"]
    assert report["window_candidates_rejected"] > 0
    assert report["window_pairs_examined"] > 1
    assert report["i1"] < crust.n - 1
    assert core.nb[report["i2"]] > crust.nb[report["i1"]]
    low, high = report["css_transition_low_index"], report["css_transition_high_index"]
    assert high == low + 1
    assert info["density_jump_fm3"] > 0
    assert matched.nb[high] > matched.nb[low]
    assert matched.eps[high] > matched.eps[low]
    assert 0 < report["css_single_sound_speed"] <= 1
    assert matched.p[high] == pytest.approx(matched.p[low], rel=1e-11, abs=1e-12)
    assert matched.mub[high] == pytest.approx(matched.mub[low], rel=1e-11, abs=1e-10)
    target = tmp_path / "overlap_css"
    write_compose(target, matched)
    assert run_native("verify", target, tmp_path / "verification.json")["overall"] == "PASS"


def test_css_rejects_when_every_endpoint_pair_has_invalid_density_order():
    """Every core knot is below every crust knot: no ordered window exists."""
    crust = _css_barotrope(rest_mass=939.0, stiffness=1.0, n_min=0.05, n_max=0.50)
    core = _css_barotrope(rest_mass=939.0, stiffness=1.2, n_min=0.01, n_max=0.04)
    assert core.nb[-1] < crust.nb[0]
    with pytest.raises(NativeExecutionError, match="exhaustive native endpoint search"):
        match_crust_core(crust, core, method="css", n_bridge=40, drop_composition=True)
