"""EOS-06 native CSS parameterized high-density branch and join verification.

EOS-06 acceptance (data/ns_cowling_cloud_agent_dag.v1.json):
  1. per-task source-backed parameter manifest
  2. native CLI action
  3. joined crust+core+CSS table not sidecar-only
  4. Delta_epsilon and first-law jump verified
  5. sound speed and causal checks
  6. native read-back

The analytic CSS kernel (fortran_eos/src/compose_eos_css.f90) already existed
and passed CI unit tests (M5A / E390-05A: css_segment.log). The genuine
PARTIAL_NATIVE_KERNEL_JOIN_NOT_RUN remainder was that the kernel was NOT
exposed as a CLI action, NO three-segment crust+core+CSS join existed, NO
per-task source-backed parameter manifest existed, and NO native read-back
was performed (the run_css_task orchestrator route was an explicit
fail-closed stub). This node closes that gap with:
  - data/eos390_css_segments.v1.json (source-backed parameter manifest,
    criterion 1)
  - fortran_eos/src/compose_eos_css_join.f90 (the join subroutine) + the
    'css-join' CLI action in compose_eos_fortran.f90 (criterion 2)
  - the joined table = matched crust+core knots + CSS tail knots, written
    as a single CompOSE table (NOT a sidecar; criterion 3)
  - the first-law jump identity Delta_epsilon = mu_t * Delta_n_B verified
    by the kernel and re-asserted on the joined table (criterion 4)
  - assembled sound speed 0 < cs2 <= 1 across the whole join (criterion 5)
  - the joined table written to a CompOSE dir and read back losslessly
    (criterion 6)
  - Python API join_crust_core_css wiring the CLI to EOSData round-trip.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))
from compose_eos import (  # noqa: E402
    read_compose, match_crust_core, join_crust_core_css, write_compose, EOSData,
)

CRUST = ROOT / "eos_generator" / "eos_data_normalized" / "GPPVA_FSU2H"
CORE = ROOT / "eos_generator" / "eos_data_normalized" / "GM_GM1"
CSS_MANIFEST = ROOT / "data" / "eos390_css_segments.v1.json"


def _require_reference_inputs():
    required = (CRUST / "eos.nb", CRUST / "eos.thermo",
                CORE / "eos.nb", CORE / "eos.thermo")
    missing = [str(p) for p in required if not p.is_file()]
    if missing:
        pytest.skip("benchmark EOS inputs unavailable: " + ", ".join(missing))


@pytest.fixture(scope="module")
def matched_table():
    """Run match_crust_core (quintic_P) once to produce the matched crust+core
    table that the CSS join appends to. Shared across the criterion tests."""
    _require_reference_inputs()
    crust = read_compose(CRUST)
    core = read_compose(CORE)
    matched, info = match_crust_core(
        crust, core, method="quintic_P", n_bridge=40,
        core_compo_mode="strict_thermo_only", drop_composition=True,
    )
    assert info["consistency"]["overall"]["pass_"], (
        "matched crust+core base must pass consistency for the CSS join test")
    return matched


def _css_params_from_manifest():
    """Load the source-backed benchmark CSS segment from the manifest."""
    m = json.loads(CSS_MANIFEST.read_text(encoding="utf-8"))
    seg = m["segments"]["benchmark_conformal_first_order"]
    return {
        "p_transition": float(seg["p_transition"]),
        "epsilon_low": float(seg["epsilon_low"]),
        "nb_low": float(seg["nb_low"]),
        "delta_epsilon": float(seg["delta_epsilon"]),
        "cs2": float(seg["cs2"]),
        "n_css": int(seg["n_css"]),
    }


# ---- Criterion 1: per-task source-backed parameter manifest ----

def test_css_parameter_manifest_is_source_backed_and_well_formed():
    """EOS-06 criterion 1: a source-backed CSS parameter manifest exists,
    with a real paper reference (Alford-Han-Prakash 2013) and a per-task
    segment carrying all five CSS parameters (p_transition, epsilon_low,
    nb_low, delta_epsilon, cs2) plus n_css and interface_class."""
    assert CSS_MANIFEST.is_file(), "data/eos390_css_segments.v1.json missing"
    m = json.loads(CSS_MANIFEST.read_text(encoding="utf-8"))
    assert m["schema"] == "ns-cowling.eos390-css-segments"
    ref = m["reference"]
    assert "Alford" in ref["citation"] and "Prakash" in ref["citation"], (
        "manifest reference is not Alford-Han-Prakash (source not backed)")
    assert ref["arxiv"] == "1302.4732"
    seg = m["segments"]["benchmark_conformal_first_order"]
    assert seg["source_reference"] and "Alford" in seg["source_reference"]
    for key in ("p_transition", "epsilon_low", "nb_low", "delta_epsilon", "cs2", "n_css"):
        assert key in seg, f"segment missing {key}"
    assert seg["interface_class"] == "first_order"
    assert 0 < float(seg["cs2"]) <= 1
    # the manifest must document the blocked models (EOS-03 source-blocked)
    assert "blocked_segments" in m
    assert set(m["blocked_segments"]["models"]) >= {
        "GPPVA_FSU2", "GPPVA_FSU2H", "GPPVA_NL3wrL55", "GPPVA_TM1e",
        "SPG_M3", "SPG_M5"}


# ---- Criterion 2: native CLI action ----

def test_css_join_cli_action_is_available():
    """EOS-06 criterion 2: the 'css-join' CLI action is exposed in the
    native Fortran executable (not just the analytic kernel)."""
    from compose_eos.native import executable
    import subprocess
    exe = executable()
    # the no-arg usage string lists the available actions
    r = subprocess.run([str(exe)], capture_output=True, text=True)
    usage = r.stdout + r.stderr
    assert "css-join" in usage, (
        "'css-join' not in the CLI usage -- the native CLI action is not exposed")


def test_css_join_python_api_invokes_native_action(matched_table):
    """EOS-06 criterion 2 (Python wiring): join_crust_core_css calls the
    native css-join CLI action and returns a joined EOSData + report."""
    params = _css_params_from_manifest()
    joined, info = join_crust_core_css(matched_table, **params)
    assert isinstance(joined, EOSData)
    assert info["operation"] == "css-join"
    assert info["n_total"] == info["n_matched"] + info["n_css"]
    assert info["n_matched"] == matched_table.n


# ---- Criterion 3: joined crust+core+CSS table not sidecar-only ----

def test_joined_table_is_three_segment_and_not_sidecar(matched_table):
    """EOS-06 criterion 3: the joined table is the matched crust+core knots
    followed by the CSS tail knots, as a SINGLE EOSData (not a sidecar file).
    The CSS tail knots have nb strictly greater than the matched top, and the
    joined grid spans the full density range (crust low -> CSS high)."""
    params = _css_params_from_manifest()
    joined, info = join_crust_core_css(matched_table, **params)
    # three-segment: matched (crust+core) knots + CSS tail knots
    assert info["n_css"] >= 4
    assert joined.n == info["n_total"]
    # the CSS tail knots extend past the matched core's top density
    matched_top_nb = float(matched_table.nb[-1])
    css_tail_nbs = [float(joined.nb[i]) for i in range(info["n_matched"], joined.n)]
    assert all(n > matched_top_nb for n in css_tail_nbs), (
        "CSS tail does not extend past the matched core top -- join failed")
    # the joined grid is strictly increasing (criterion 3 + the join invariant)
    assert np.all(np.diff(np.asarray(joined.nb, dtype=float)) > 0), (
        "joined nb grid is not strictly increasing")
    # NOT a sidecar: the joined EOSData carries the CSS rows inline; there is
    # no separate CSS-only table. The composition is absent (strict-thermo-only)
    # but the thermo rows (nb, Q1..Q7) include the CSS tail in the main array.
    assert not joined.composition_file_present
    assert not joined.composition_available


# ---- Criterion 4: Delta_epsilon and first-law jump verified ----

def test_css_transition_jump_identity_verified(matched_table):
    """EOS-06 criterion 4: the CSS transition jump identity
    Delta_epsilon = mu_t * Delta_n_B is verified at machine precision (the
    Fortran kernel enforces it at initialization, residual < 64*eps*scale).
    The first-law jump residual is reported and must be ~0."""
    params = _css_params_from_manifest()
    joined, info = join_crust_core_css(matched_table, **params)
    assert info["jump_verified"] is True, (
        "first-law jump identity not verified by the kernel")
    assert info["interface_class"] == 1, (
        "first-order interface class not set (Delta_epsilon > 0)")
    # machine-precision residual (the kernel uses 64*eps*scale tolerance)
    assert info["first_law_jump_residual"] < 1e-10, (
        f"first-law jump residual {info['first_law_jump_residual']} too large")
    # the high-phase state is the jump: epsilon_high = epsilon_low + delta_epsilon
    assert abs(info["epsilon_high"] - (info["epsilon_low"] + info["delta_epsilon"])) < 1e-6, (
        "epsilon_high != epsilon_low + delta_epsilon (jump not applied)")
    # nb_high > nb_low (the first-order jump increases density)
    assert info["nb_high"] > info["nb_low"]


def test_joined_table_first_law_holds(matched_table):
    """EOS-06 criterion 4 (first-law on the joined table): the joined table
    satisfies mu_B = F + P/nb at every row (the join did not break closure)."""
    params = _css_params_from_manifest()
    joined, info = join_crust_core_css(matched_table, **params)
    mn = float(joined.mn)
    nb = np.asarray(joined.nb, dtype=float)
    q1 = np.asarray(joined.Q1, dtype=float)
    q3 = np.asarray(joined.Q3, dtype=float)
    q6 = np.asarray(joined.Q6, dtype=float)
    P = q1 * nb               # code convention: P = q1*nb
    mu_B = (q3 + 1.0) * mn
    F = (q6 + 1.0) * mn
    # mu_B = F + P/nb  (first law)
    residual = np.abs(mu_B - (F + P / nb))
    assert np.all(residual < 1e-3 * np.maximum(1.0, np.abs(mu_B))), (
        f"first law violated on the joined table; max residual {residual.max()}")


# ---- Criterion 5: sound speed and causal checks ----

def test_joined_table_is_causal(matched_table):
    """EOS-06 criterion 5: the assembled sound speed on the joined table is
    causal: 0 < cs2 <= 1 across the whole join (crust low-density through the
    CSS high-density branch). The native report's cs2_min/cs2_max bound the
    joined table, not just the CSS branch."""
    params = _css_params_from_manifest()
    joined, info = join_crust_core_css(matched_table, **params)
    assert info["causal"] is True, "joined table flagged non-causal"
    assert info["cs2_min"] > 0.0, (
        f"nonpositive min sound speed {info['cs2_min']} on the joined table")
    assert info["cs2_max"] <= 1.0, (
        f"superluminal max sound speed {info['cs2_max']} on the joined table")
    # the CSS branch is constant-cs2 = 1/3 by construction; the joined cs2_max
    # is dominated by the core's higher-density sound speed, all subluminal.
    assert info["cs2_max"] < 1.0


def test_joined_table_passes_native_consistency_gate(matched_table):
    """EOS-06 criterion 5 (consistency): the joined crust+core+CSS table
    passes the native C1-C11 consistency gate (hard checks C1-C6/C10/C11).
    This re-uses the EOS-04 consistency guarantee on the three-segment table."""
    params = _css_params_from_manifest()
    joined, info = join_crust_core_css(matched_table, **params)
    cons = info["consistency"]
    assert cons["overall"]["pass_"] is True, (
        f"joined CSS table fails consistency: "
        f"failed_checks={cons['overall'].get('failed_checks')}")
    assert cons["C6"]["pass_"] is True   # dP/dn_B > 0
    assert cons["C10"]["pass_"] is True  # causal sound speed
    assert cons["C11"]["pass_"] is True  # dmu_B/dn_B > 0


# ---- Criterion 6: native read-back ----

def test_joined_table_native_read_back(tmp_path, matched_table):
    """EOS-06 criterion 6: the joined crust+core+CSS table is written to a
    CompOSE directory (eos.nb, eos.t, eos.thermo, eos.yq -- a full CompOSE
    table, not a sidecar) and read back losslessly: row count, density range,
    and the CSS high-density tail survive the round-trip."""
    params = _css_params_from_manifest()
    joined, info = join_crust_core_css(matched_table, **params)
    # write the joined table to a CompOSE directory
    compose_dir = tmp_path / "joined_compose"
    write_compose(compose_dir, joined)
    # the CompOSE dir is a full table, not a sidecar
    for f in ("eos.nb", "eos.t", "eos.thermo", "eos.yq"):
        assert (compose_dir / f).is_file(), f"{f} not written (not a full CompOSE table)"
    # read it back
    readback = read_compose(compose_dir)
    # lossless round-trip: row count + density range + CSS tail
    assert readback.n == joined.n, (
        f"read-back row count {readback.n} != joined {joined.n}")
    np.testing.assert_allclose(np.asarray(readback.nb, dtype=float),
                               np.asarray(joined.nb, dtype=float),
                               rtol=1e-12, atol=0)
    # the CSS high-density tail survived (max nb ~ the CSS branch top)
    assert float(readback.nb[-1]) > 1.3, (
        f"CSS tail lost on read-back; max nb {float(readback.nb[-1])}")
    # the first-law identity holds on the read-back table
    mn = float(readback.mn)
    nb = np.asarray(readback.nb, dtype=float)
    q1 = np.asarray(readback.Q1, dtype=float)
    q3 = np.asarray(readback.Q3, dtype=float)
    q6 = np.asarray(readback.Q6, dtype=float)
    P = q1 * nb
    mu_B = (q3 + 1.0) * mn
    F = (q6 + 1.0) * mn
    residual = np.abs(mu_B - (F + P / nb))
    assert np.all(residual < 1e-3 * np.maximum(1.0, np.abs(mu_B))), (
        "first law violated on the read-back table")
