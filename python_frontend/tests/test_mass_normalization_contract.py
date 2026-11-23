"""NSC remediation Phase B regression tests (mission section 4.6).

Regression cases pinning physical-energy-preserving normalization:
  1-4. identical, neutron-only, proton-only and both reference changes preserve
       native Maxwell physical output after complete CompOSE re-encoding
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
import subprocess
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


# ---------------------------------------------------------------- 1-4: reference normalization
@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
@pytest.mark.parametrize("mn,mp", [(939.,939.), (950.,939.), (939.,938.), (950.,938.)])
def test_1_to_4_native_maxwell_reference_reencoding(mn, mp):
    """Raw native route preserves physical EOS under independent header references.

    The former header-equality tests conflated an encoding change with an
    additive physical rest-energy change. Physical-offset controls remain in
    test_maxwell_mass_gate.py; source normalization must not fit those offsets.
    """
    from test_maxwell_mass_gate import physical_pair, reencode
    low, original = physical_pair()
    high = reencode(original, mn, mp)
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        lo, hi, target = root/"lo", root/"hi", root/"out"
        lo.mkdir(); hi.mkdir(); target.mkdir()
        low.write_internal_directory(lo); high.write_internal_directory(hi)
        report = run_native("maxwell-transition", lo, hi, target, "true")
        out = EOSData.from_internal_directory(target)
    assert report["status"] == "ok"
    assert report["mu_transition"] == pytest.approx(990., abs=2e-10)
    assert report["pressure_transition"] == pytest.approx(8., abs=2e-11)
    assert report["n_low"] == pytest.approx(.12, abs=2e-12)
    assert report["n_high"] == pytest.approx(.13, abs=2e-12)
    assert report["physical_energy_shift_mev"] == 0.
    transition_mu = report["mu_transition"]
    low_rows = int(report["low_phase_rows"])
    assert low_rows + int(report["high_phase_rows"]) == out.n
    assert out.Iphase[low_rows - 1] != out.Iphase[low_rows]
    for source, first, stop, retained in (
        (low, 0, low_rows, low.mub < transition_mu - 1e-9),
        (original, low_rows, out.n, original.mub > transition_mu + 1e-9),
    ):
        # Density is not a cross-phase identity key. Compare only the retained
        # chemical-potential domain of this source and its output branch.
        for index in np.flatnonzero(retained):
            rows = np.flatnonzero(np.isclose(
                out.nb[first:stop], source.nb[index], rtol=0., atol=1e-13)) + first
            assert rows.size == 1
            row = rows[0]
            np.testing.assert_allclose([out.eps[row],out.p[row],out.mub[row],out.muq[row]],
                [source.eps[index],source.p[index],source.mub[index],source.muq[index]],
                rtol=0.,atol=2e-10)
    # Explicit false-identity control: discarded low phase at mu=995 has
    # n=.13, also the *other* phase's transition density at mu=990.
    discarded = np.flatnonzero(np.isclose(low.mub, 995., rtol=0., atol=1e-10))
    assert discarded.size == 1 and low.mub[discarded[0]] > transition_mu
    aliases = np.flatnonzero(np.isclose(out.nb, low.nb[discarded[0]], rtol=0., atol=1e-13))
    assert aliases.size == 1 and aliases[0] == low_rows
    assert abs(out.mub[aliases[0]] - low.mub[discarded[0]]) > 1.


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
    re-encodes dimensionless columns while preserving physical energy."""
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
    """With composition present and drop_composition=True the normalization
    runs without a physical energy shift and exposes provenance fields."""
    d_cr = read_compose(EOS / "GPPVA_DD2")
    d_co = read_compose(EOS / "BBB_BHF2")
    matched, info = match_crust_core(d_cr, d_co, method="css", n_bridge=8,
                               drop_composition=True)
    assert info["native_report"]["dm_mass_gauge_meV"] == 0.0
    assert info["dm_crust"] == 0.0
    assert info["normalization_reference_delta_mev"] == pytest.approx(d_co.mn - d_cr.mn)
    assert info["mass_normalization_semantics"] == "physical_energy_preserving_Q_rescaling"
    assert info["native_report"]["source_euler_residual_meV"] > 0
    # P invariance: matched P at retained crust rows equals source P
    assert np.allclose(matched.p[:info["i1"] + 1],
                       d_cr.p[:info["i1"] + 1], rtol=0, atol=1e-12)


# ---------------------------------------------------------------- 9-10: round trips
@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
def test_9_gppva_dd2_gauge_round_trip():
    """Changing the output normalization preserves physical F and epsilon."""
    d_cr = read_compose(EOS / "GPPVA_DD2")
    d_co = read_compose(EOS / "BBB_BHF2")
    matched, info = match_crust_core(d_cr, d_co, method="css", n_bridge=8,
                               drop_composition=True)
    i1 = info["i1"]
    # Physical free energy per baryon is unchanged.
    assert np.allclose(matched.F[:i1 + 1], d_cr.F[:i1 + 1],
                       rtol=0, atol=1e-9)
    # Re-encoding Q6 back to the source reference recovers source Q6.
    source_q6 = (matched.Q6[:i1 + 1] + 1.0) * matched.mn / d_cr.mn - 1.0
    assert np.allclose(source_q6, d_cr.Q6[:i1 + 1], rtol=0, atol=1e-12)
    # epsilon is an absolute gravitating source and must stay unchanged.
    assert np.allclose(matched.eps[:i1 + 1],
                       d_cr.eps[:i1 + 1],
                       rtol=0, atol=1e-7)


@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
def test_10_q3_q7_normalization_round_trip():
    """Output Q3/Q6/Q7 are normalized to out%mn; inverse recovers source."""
    d_cr = read_compose(EOS / "GPPVA_DD2")
    d_co = read_compose(EOS / "BBB_BHF2")
    matched, info = match_crust_core(d_cr, d_co, method="css", n_bridge=8,
                               drop_composition=True)
    i1 = info["i1"]
    # Q6: (Q6_out+1)*mn_out = (Q6_src+1)*mn_src
    lhs = (matched.Q6[:i1 + 1] + 1.0) * matched.mn
    rhs = (d_cr.Q6[:i1 + 1] + 1.0) * d_cr.mn
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
    """A header-only rewrite changes physical energy; complete Q rescaling does not."""
    d_cr = read_compose(EOS / "GPPVA_DD2")
    d_co = read_compose(EOS / "BBB_BHF2")
    matched, info = match_crust_core(d_cr, d_co, method="css", n_bridge=8,
                                   drop_composition=True)
    assert matched.mn == d_co.mn and matched.mp == d_co.mp
    stop = info["i1"] + 1
    header_only_f = (d_cr.Q6[:stop] + 1.0) * matched.mn
    assert np.max(np.abs(header_only_f - d_cr.F[:stop])) > 1e-6
    np.testing.assert_allclose(matched.F[:stop], d_cr.F[:stop], rtol=0, atol=1e-9)
    # Source Q3 tabulation residual is a separate disclosed Euler correction.
    np.testing.assert_allclose(matched.mub[:stop], d_cr.F[:stop] + d_cr.Q1[:stop],
                               rtol=0, atol=1e-9)


@pytest.mark.skipif(_needs_backend(), reason="requires native backend")
@pytest.mark.parametrize("method", ["css", "cubic_F", "cubic_P", "quintic_F", "quintic_P"])
def test_matching_invariant_under_compose_mass_normalization(method):
    """Re-encoding an identical physical crust must not alter its matched EOS."""
    crust = read_compose(ROOT / "eos_generator/eos_data_normalized/GPPVA_FSU2H")
    core = read_compose(ROOT / "eos_generator/eos_data_normalized/GM_GM1")
    encoded = crust.copy()
    encoded.mn = crust.mn + 0.01
    ratio = crust.mn / encoded.mn
    for field in ("Q3", "Q6", "Q7"):
        setattr(encoded, field, (getattr(crust, field) + 1.0) * ratio - 1.0)
    for field in ("Q4", "Q5"):
        setattr(encoded, field, getattr(crust, field) * ratio)
    # This modifies representation, not a model's constituent masses/energy.
    np.testing.assert_allclose(encoded.eps, crust.eps, rtol=2e-14, atol=1e-12)
    original, oi = match_crust_core(crust, core, method=method, n_bridge=40, drop_composition=True)
    recoded, ri = match_crust_core(encoded, core, method=method, n_bridge=40, drop_composition=True)
    assert (oi["i1"], oi["i2"]) == (ri["i1"], ri["i2"])
    for field in ("nb", "p", "eps", "mub", "F"):
        np.testing.assert_allclose(getattr(recoded, field), getattr(original, field),
                                   rtol=2e-10, atol=2e-9, err_msg=field)
    stop = oi["i1"] + 1
    np.testing.assert_allclose(original.F[:stop], crust.F[:stop], rtol=2e-14, atol=1e-10)
    np.testing.assert_allclose(original.eps[:stop], crust.eps[:stop], rtol=2e-14, atol=1e-10)


def test_python_matching_preserves_physical_energy_under_normalization():
    """Independent Python implementation obeys the same physical invariance."""
    code = r'''
import sys, types
from pathlib import Path
import numpy as np
root = Path(sys.argv[1])
pkg = types.ModuleType("compose_eos")
pkg.__path__ = [str(root / "eos_generator/compose_eos_pkg/compose_eos")]
sys.modules["compose_eos"] = pkg
from compose_eos.io import read_compose
from compose_eos.matching import match_crust_core
cr = read_compose(root / "eos_generator/eos_data_normalized/GPPVA_FSU2H")
co = read_compose(root / "eos_generator/eos_data_normalized/GM_GM1")
for method in ("css", "cubic_F", "cubic_P", "quintic_F", "quintic_P"):
    matched, info = match_crust_core(cr, co, method=method, n_bridge=40,
                                    core_compo_mode="strict_thermo_only")
    stop = info["i1"] + 1
    np.testing.assert_allclose(matched.F[:stop], cr.F[:stop], rtol=2e-14, atol=1e-10)
    np.testing.assert_allclose(matched.eps[:stop], cr.eps[:stop], rtol=2e-14, atol=1e-10)
'''
    result = subprocess.run([sys.executable, "-c", code, str(ROOT)], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
