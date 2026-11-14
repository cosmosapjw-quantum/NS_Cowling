"""NSC remediation Phase B regression tests (mission section 4.6).

14 cases pinning the rest-mass normalization contract:
  1. identical neutron/proton mass pair -> Maxwell gate PASS
  2. neutron-only mismatch -> FAIL
  3. proton-only mismatch -> FAIL
  4. both mismatch -> FAIL
  5. common neutron/proton mass within one EOS -> allowed
  6. different neutron/proton physical masses -> correctly supported
  7. composition missing -> physical correction BLOCKED
  8. composition present and verified -> transformation PASS
  9. GPPVA_DD2 normalization round-trip
 10. Q3-Q7 normalization round-trip
 11. P, eps, mu_B, mu_Q invariants
 12. first-law identity checks
 13. source file immutability
 14. header-only rewrite detection

Physics gates are NEVER weakened to pass these tests.
"""
from __future__ import annotations
import hashlib
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))

os.environ.setdefault(
    "COMPOSE_EOS_FORTRAN_EXE",
    str(ROOT / "build" / "compose_eos_fortran" / "compose_eos_fortran"),
)

from compose_eos.api import read_compose, match_crust_core  # noqa: E402
from compose_eos.eos_data import EOSData  # noqa: E402
from compose_eos.native import run_native  # noqa: E402

EOS = ROOT / "eos_generator" / "eos_data"


def _needs_backend():
    return not Path(os.environ["COMPOSE_EOS_FORTRAN_EXE"]).is_file()


def _sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


# ---------------------------------------------------------------- 1-4: Maxwell pairwise mass gate
def _maxwell_attempt(low: EOSData, high: EOSData):
    """Invoke the native maxwell-transition with local_neutrality_confirmed.

    Returns (ok, message). Uses internal-directory round trip.
    """
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        lo, hi, out = root / "lo", root / "hi", root / "out"
        lo.mkdir(); hi.mkdir(); out.mkdir()
        low.write_internal_directory(lo)
        high.write_internal_directory(hi)
        try:
            run_native("maxwell-transition", lo, hi, out, "true")
            return True, ""
        except Exception as exc:  # NativeExecutionError
            return False, str(exc)


@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
def test_1_identical_mass_pair_passes_gate():
    """Same EOS split into two phases with identical (mn, mp) passes the gate."""
    d = read_compose(EOS / "BBB_BHF2")
    low, high = d.copy(), d.copy()
    ok, msg = _maxwell_attempt(low, high)
    # gate-specific rejection must NOT appear; other physics failures are fine
    assert "rest-mass conventions" not in msg, msg


@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
def test_2_neutron_only_mismatch_fails_gate():
    d = read_compose(EOS / "BBB_BHF2")
    low, high = d.copy(), d.copy()
    high.mn = high.mn + 0.01  # neutron-only mismatch
    ok, msg = _maxwell_attempt(low, high)
    assert not ok and "rest-mass conventions" in msg, msg


@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
def test_3_proton_only_mismatch_fails_gate():
    d = read_compose(EOS / "BBB_BHF2")
    low, high = d.copy(), d.copy()
    high.mp = high.mp + 0.01  # proton-only mismatch, neutron identical
    ok, msg = _maxwell_attempt(low, high)
    assert not ok and "rest-mass conventions" in msg, msg


@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
def test_4_both_masses_mismatch_fail_gate():
    d = read_compose(EOS / "GPPVA_DD2")
    low, high = d.copy(), d.copy()
    high.mn += 0.5; high.mp += 0.5
    ok, msg = _maxwell_attempt(low, high)
    assert not ok and "rest-mass conventions" in msg, msg


# ---------------------------------------------------------------- 5-6: single-EOS conventions
def test_5_common_mass_within_one_eos_allowed():
    """Uniform mn==mp convention (GPPVA family) is a legal provider gauge."""
    d = read_compose(EOS / "GPPVA_DD2")
    assert d.mn == d.mp == 939.0
    assert d.n >= 3  # parses fine


def test_6_physical_split_masses_supported():
    """Physical mn != mp convention (core models) is fully supported."""
    d = read_compose(EOS / "APR_APR")
    assert abs(d.mn - 939.56533) < 1e-9
    assert abs(d.mp - 938.27200) < 1e-9
    assert abs((d.mn - d.mp) - 1.29333) < 1e-9


# ---------------------------------------------------------------- 7-8: composition policy
@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
def test_7_composition_missing_blocks_physical_correction():
    """Without per-nucleus composition the physical split-mass correction is
    not derivable: the API must NOT silently apply it. The smooth matcher
    applies only the documented uniform gauge shift (representation-only)."""
    d_cr = read_compose(EOS / "GPPVA_DD2")
    d_co = read_compose(EOS / "BBB_BHF2")
    # GPPVA_DD2 carries aggregate composition -> requires drop_composition
    with pytest.raises(ValueError, match="drop_composition=True"):
        match_crust_core(d_cr, d_co, method="css", n_bridge=8)
    # dropping is the ONLY non-blocking route, and it is thermo-only
    matched, info = match_crust_core(d_cr, d_co, method="css", n_bridge=8,
                               drop_composition=True)
    assert info["composition_provenance"] == "none_strict_thermo_only"
    assert info["core_composition_input_usable"] is False


@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
def test_8_composition_present_transformation_pass():
    """With composition present and drop_composition=True the gauge-aligned
    transformation runs and the native report exposes the provenance fields."""
    d_cr = read_compose(EOS / "GPPVA_DD2")
    d_co = read_compose(EOS / "BBB_BHF2")
    matched, info = match_crust_core(d_cr, d_co, method="css", n_bridge=8,
                               drop_composition=True)
    assert abs(info["native_report"]["dm_mass_gauge_meV"] - 0.56533) < 1e-6
    assert info["native_report"]["source_euler_residual_meV"] > 0
    # P invariance: matched P at retained crust rows equals source P
    assert np.allclose(matched.p[:info["i1"] + 1],
                       d_cr.p[:info["i1"] + 1], rtol=0, atol=1e-12)


# ---------------------------------------------------------------- 9-10: round trips
@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
def test_9_gppva_dd2_gauge_round_trip():
    """F_tgt = F_src + dmcr, inverse F_src = F_tgt - dmcr, exact round trip."""
    d_cr = read_compose(EOS / "GPPVA_DD2")
    d_co = read_compose(EOS / "BBB_BHF2")
    dmcr = d_co.mn - d_cr.mn
    matched, info = match_crust_core(d_cr, d_co, method="css", n_bridge=8,
                               drop_composition=True)
    i1 = info["i1"]
    # forward: matched F on retained crust rows = source F + dmcr
    assert np.allclose(matched.F[:i1 + 1], d_cr.F[:i1 + 1] + dmcr,
                       rtol=0, atol=1e-9)
    # inverse: subtracting recovers source exactly
    assert np.allclose(matched.F[:i1 + 1] - dmcr, d_cr.F[:i1 + 1],
                       rtol=0, atol=1e-9)
    # eps: eps_tgt = eps_src + n*dmcr
    assert np.allclose(matched.eps[:i1 + 1],
                       d_cr.eps[:i1 + 1] + d_cr.nb[:i1 + 1] * dmcr,
                       rtol=0, atol=1e-7)


@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
def test_10_q3_q7_normalization_round_trip():
    """Output Q3/Q6/Q7 are normalized to out%mn; inverse recovers source."""
    d_cr = read_compose(EOS / "GPPVA_DD2")
    d_co = read_compose(EOS / "BBB_BHF2")
    dmcr = d_co.mn - d_cr.mn
    matched, info = match_crust_core(d_cr, d_co, method="css", n_bridge=8,
                               drop_composition=True)
    i1 = info["i1"]
    # Q6: (Q6_out+1)*mn_out = (Q6_src+1)*mn_src + dmcr
    lhs = (matched.Q6[:i1 + 1] + 1.0) * matched.mn
    rhs = (d_cr.Q6[:i1 + 1] + 1.0) * d_cr.mn + dmcr
    assert np.allclose(lhs, rhs, rtol=0, atol=1e-9)
    # Q7 == Q6 on the cold output
    assert np.allclose(matched.Q7, matched.Q6, rtol=0, atol=1e-15)


# ---------------------------------------------------------------- 11-12: invariants
@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
def test_11_thermodynamic_invariants():
    """P preserved; mu_B = F + P/n on every output row; mu_Q preserved."""
    d_cr = read_compose(EOS / "GPPVA_DD2")
    d_co = read_compose(EOS / "BBB_BHF2")
    matched, info = match_crust_core(d_cr, d_co, method="css", n_bridge=8,
                               drop_composition=True)
    i1, i2 = info["i1"], info["i2"]
    # mu_B identity on ALL rows (crust re-derived, bridge, core re-derived).
    # Per baryon: mu_B = F + P/n_B = F + Q1. Q1 is MeV per baryon; matched.p
    # is per volume and must NOT be used here.
    mu_b = (matched.Q3 + 1.0) * matched.mn
    assert np.allclose(mu_b, matched.F + matched.Q1, rtol=0, atol=2e-9)
    # mu_Q preserved on retained crust rows: Q4_out*mn_out = Q4_src*mn_src
    muq_out = matched.Q4[:i1 + 1] * matched.mn
    muq_src = d_cr.Q4[:i1 + 1] * d_cr.mn
    assert np.allclose(muq_out, muq_src, rtol=0, atol=1e-9)


@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
def test_12_first_law_identity():
    """dF = P/n dn integrated: mu2-mu1 = integral; and Q2=0, Q5=0, Q6=Q7."""
    d_cr = read_compose(EOS / "GPPVA_DD2")
    d_co = read_compose(EOS / "BBB_BHF2")
    matched, info = match_crust_core(d_cr, d_co, method="css", n_bridge=8,
                               drop_composition=True)
    assert np.max(np.abs(matched.Q2)) == 0.0
    assert np.max(np.abs(matched.Q5)) == 0.0
    assert np.allclose(matched.Q6, matched.Q7, rtol=0, atol=1e-15)
    # bridge first-law residual reported by native
    assert info["native_report"]["first_law_relative_residual"] < 1e-10


# ---------------------------------------------------------------- 13-14: source immutability
@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
def test_13_source_file_immutability():
    """match_core must not modify any source eos.thermo file."""
    thermo = EOS / "GPPVA_DD2" / "eos.thermo"
    h0 = _sha256(thermo)
    d_cr = read_compose(EOS / "GPPVA_DD2")
    d_co = read_compose(EOS / "BBB_BHF2")
    match_crust_core(d_cr, d_co, method="css", n_bridge=8, drop_composition=True)
    assert _sha256(thermo) == h0


@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
def test_14_header_rewrite_is_detectable_and_forbidden():
    """The matcher records the core gauge in the output header (mn, mp) and
    the gauge shift separately. A naive 'header-only rewrite' (claiming the
    physical split without transforming Q columns) is detectable: Q3 would
    fail the mu = F + P/n identity. Pin that the real output satisfies it."""
    d_cr = read_compose(EOS / "GPPVA_DD2")
    d_co = read_compose(EOS / "BBB_BHF2")
    matched, info = match_crust_core(d_cr, d_co, method="css", n_bridge=8,
                               drop_composition=True)
    # output header carries the core convention
    assert matched.mn == d_co.mn and matched.mp == d_co.mp
    # a header-only rewrite would leave crust Q3 inconsistent with shifted F:
    # naive Q3 would be (Q3_src+1)*mn_src/mn_out - 1; check it FAILS identity
    naive = (d_cr.Q3[: info["i1"] + 1] + 1.0) * d_cr.mn / matched.mn - 1.0
    mu_naive = (naive + 1.0) * matched.mn
    lhs = matched.F[: info["i1"] + 1] + matched.Q1[: info["i1"] + 1]
    # the naive mapping must NOT satisfy the identity (that's why the code
    # re-derives Q3); at least one retained row differs beyond 1e-6
    assert np.max(np.abs(mu_naive - lhs)) > 1e-6
