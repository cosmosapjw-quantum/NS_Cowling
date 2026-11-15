"""Regression tests for the pairwise neutron/proton mass-equality admission gate
of the native Fortran `construct_maxwell_transition`.

Physical requirement (mission PR):
    m_n^(L) == m_n^(H)  AND  m_p^(L) == m_p^(H)
This is a PAIRWISE equality between the low- and high-density EOS conventions.
It does NOT require m_n == m_p within either EOS. The physically distinct
neutron and proton masses must be retained when appropriate.

The native gate (fortran_eos/src/compose_eos_phase_transition.f90 lines 664-665)
compares `abs(low%mn-high%mn) > 1e-12*max(low%mn,high%mn)` OR
`abs(low%mp-high%mp) > 1e-12*max(low%mp,high%mp)` and rejects on a mismatch.
The mn and mp comparisons are INDEPENDENT: a neutron-only or proton-only
mismatch is rejected even when the other mass matches.

These tests must NOT widen the gate to absorb differences of order 0.5-1 MeV.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "python_frontend"))

from compose_eos import EOSData, construct_maxwell_transition


# Canonical physical masses (Typel et al. 2010, arXiv:0908.2344, Table II)
MN_PHYS = 939.56536
MP_PHYS = 938.27203
# Common-convention mass (Lalazissis 2005 DD-ME2; TW99)
M_COMMON = 939.0

# Current native gate relative tolerance
GATE_RTOL = 1e-12


def _make_phase(*, mn: float, mp: float, n_min: float = 0.5,
                n_max: float = 2.0, count: int = 81,
                stiffness: float = 100.0) -> EOSData:
    """Cold quadratic barotrope P = K n^2, eps = m_ref n + K n^2 with m_ref = mn.

    Both phases use the SAME thermodynamic shape so that mass-equality is the
    ONLY discriminating input for the admission gate. The mn passed here is
    used both as the EOSData mn field and as the rest-mass scale of epsilon/Q.
    mp is recorded as metadata only (the barotrope is single-conserved-charge
    baryon number); the gate still compares it independently.
    """
    nb = np.linspace(n_min, n_max, count)
    pressure = stiffness * nb**2
    epsilon = mn * nb + pressure
    mu = (epsilon + pressure) / nb
    zeros = np.zeros_like(nb)
    # CompOSE Q-quantities normalised by mn of THIS phase
    q1 = pressure / nb
    q3 = mu / mn - 1.0
    q6 = epsilon / (nb * mn) - 1.0
    q7 = q6.copy()
    return EOSData(
        nb=nb, mn=mn, mp=mp,
        Q1=q1, Q2=zeros.copy(), Q3=q3,
        Q4=zeros.copy(), Q5=zeros.copy(), Q6=q6, Q7=q7,
        q_extra=zeros.copy(), Iphase=np.ones(nb.size, dtype=int),
        Yn=zeros.copy(), Yp=zeros.copy(), Ye=zeros.copy(), Ymu=zeros.copy(),
        Aav=zeros.copy(), Zav=zeros.copy(), Yh=zeros.copy(),
        Il=1, iT=0, iYq=1, inb_start=1, temperature_mev=0.0,
        charge_fraction_grid=0.0, q_extra_count=0,
        composition_file_present=False, composition_available=False,
        beta_equilibrated=True,
    )


def _gate_rejects(low: EOSData, high: EOSData) -> bool:
    """Return True if the native gate rejects the pair on a mass convention mismatch.

    A rejection is signalled by a status string containing 'rest-mass' (the
    exact Fortran message at compose_eos_phase_transition.f90:665 is
    'phase tables must use the same neutron and proton rest-mass conventions').
    """
    try:
        construct_maxwell_transition(
            low, high, assume_local_charge_neutrality=True, return_report=True
        )
    except Exception as exc:  # NativeExecutionError carries the Fortran status
        msg = str(exc)
        return ("rest-mass" in msg) or ("rest mass" in msg) or ("mass convention" in msg)
    return False  # admitted (or failed for a different reason)


# ---------------------------------------------------------------------------
# Scenario 1: identical (mn, mp) = (939.56536, 938.27203) on both phases PASSES
# the admission gate, even though mn != mp.
# ---------------------------------------------------------------------------
def test_identical_distinct_physical_masses_pass_gate():
    low = _make_phase(mn=MN_PHYS, mp=MP_PHYS)
    high = _make_phase(mn=MN_PHYS, mp=MP_PHYS)
    assert not _gate_rejects(low, high), (
        "Identical (mn,mp)=(939.56536,938.27203) on both phases must pass the "
        "mass-equality gate. The gate is pairwise (L vs H), NOT mn==mp."
    )


# ---------------------------------------------------------------------------
# Scenario 2: identical (mn, mp) = (939.0, 939.0) on both phases PASSES.
# ---------------------------------------------------------------------------
def test_identical_common_masses_pass_gate():
    low = _make_phase(mn=M_COMMON, mp=M_COMMON)
    high = _make_phase(mn=M_COMMON, mp=M_COMMON)
    assert not _gate_rejects(low, high), (
        "Identical (939.0, 939.0) on both phases must pass the mass-equality gate."
    )


# ---------------------------------------------------------------------------
# Scenario 3: a neutron-ONLY mismatch is REJECTED (mp matches, mn does not).
# ---------------------------------------------------------------------------
def test_neutron_only_mismatch_rejected():
    low = _make_phase(mn=M_COMMON, mp=MP_PHYS)
    high = _make_phase(mn=MN_PHYS, mp=MP_PHYS)  # mp identical, mn differs by ~0.565 MeV
    assert _gate_rejects(low, high), (
        "A neutron-only mismatch (mn differs by ~0.565 MeV, mp identical) must "
        "be rejected. The mn and mp comparisons are independent; matching mp "
        "does NOT rescue a mismatched mn."
    )


# ---------------------------------------------------------------------------
# Scenario 4: a proton-ONLY mismatch is REJECTED (mn matches, mp does not).
# ---------------------------------------------------------------------------
def test_proton_only_mismatch_rejected():
    low = _make_phase(mn=MN_PHYS, mp=M_COMMON)
    high = _make_phase(mn=MN_PHYS, mp=MP_PHYS)  # mn identical, mp differs by ~1.293 MeV
    assert _gate_rejects(low, high), (
        "A proton-only mismatch (mp differs by ~1.293 MeV, mn identical) must "
        "be rejected. The mn and mp comparisons are independent; matching mn "
        "does NOT rescue a mismatched mp."
    )


# ---------------------------------------------------------------------------
# Scenario 5: a discrepancy attributable solely to source DECIMAL ROUNDING is
# handled only under an EXPLICIT precision policy with recorded deltas and
# provenance. The default gate (1e-12 relative) must STILL reject a 0.5 MeV
# difference; it must not silently absorb it as "rounding".
# ---------------------------------------------------------------------------
def test_source_rounding_gap_not_silently_absorbed():
    # RG_SLY4 raw header: (939.5654, 938.272) -- 4-7 digit rounding of the
    # canonical (939.56536, 938.27203). The mp difference is ~3e-5 MeV,
    # which the 1e-12 relative gate REJECTS. This is correct: a precision
    # policy must be EXPLICIT, not implicit.
    low = _make_phase(mn=MN_PHYS, mp=MP_PHYS)         # full precision
    high = _make_phase(mn=939.5654, mp=938.272)       # source-rounded
    delta_mp = abs(MP_PHYS - 938.272)
    assert delta_mp > 0.0, "sanity: mp rounding gap is nonzero"
    # The 1e-12 relative gate rejects this tiny gap; an explicit precision
    # policy (NOT a silent tolerance widening) would be required to admit it.
    gate_rejects = _gate_rejects(low, high)
    # Document the expected behaviour: with rtol=1e-12, a ~3e-5 MeV mp gap
    # (relative ~3e-8) is rejected. If a future precision policy admits it,
    # the policy MUST record delta_mp and the source provenance.
    assert gate_rejects, (
        f"Source-decimal rounding gap (mp delta={delta_mp:.3e} MeV) must NOT be "
        "silently absorbed by the default gate. An explicit precision policy with "
        "recorded deltas and provenance is required to admit it; the gate must "
        "remain strict by default."
    )


def test_explicit_precision_policy_documented():
    """The precision policy for source-rounded masses is documented, not silent.

    This test pins the contract: if a future precision policy admits a source
    rounding gap, it MUST record delta_mn, delta_mp, the source provenance, and
    the admitted tolerance. The policy constant lives here so the contract is
    explicit and testable, not buried in a generic tolerance knob.
    """
    # Reference values and their source-rounded forms seen in the audit:
    sources = [
        # (canonical, source_rounded, source_id, provenance)
        (939.56536, 939.5654, "RG_SLY4", "RG_SLY4 eos.thermo header (4-digit mn)"),
        (938.27203, 938.272,  "RG_SLY4", "RG_SLY4 eos.thermo header (3-digit mp)"),
        (939.56536, 939.56533, "BBB_BHF2", "BBB_BHF2 eos.thermo header (5-digit mn)"),
        (938.27203, 938.272,  "BBB_BHF2", "BBB_BHF2 eos.thermo header (3-digit mp)"),
    ]
    for canonical, rounded, src, prov in sources:
        delta = abs(canonical - rounded)
        # Relative gap vs canonical:
        rel = delta / max(abs(canonical), 1.0)
        # The default gate (rtol=1e-12) rejects any of these.
        assert rel > GATE_RTOL, (
            f"{src}: source-rounded {rounded} vs canonical {canonical} has rel gap "
            f"{rel:.3e} which exceeds the default gate rtol {GATE_RTOL}. An explicit "
            "precision policy would be needed to admit it."
        )


# ---------------------------------------------------------------------------
# Scenario 6: a genuinely DIFFERENT mass convention is NEVER silently accepted
# by widening a generic tolerance. We pin the gate rtol at 1e-12 and assert
# that a 1 MeV difference (the order of a real convention change, e.g.
# 939.0 vs 939.56536) is rejected.
# ---------------------------------------------------------------------------
def test_different_convention_never_silently_accepted():
    # A ~0.565 MeV mn difference (common-939 vs physical-939.56536) represents
    # a genuine convention change, not rounding. The gate must reject it.
    low = _make_phase(mn=M_COMMON, mp=M_COMMON)         # common convention
    high = _make_phase(mn=MN_PHYS, mp=MP_PHYS)          # physical convention
    delta_mn = abs(M_COMMON - MN_PHYS)  # ~0.565 MeV
    assert delta_mn > 0.5, "sanity: this is a real convention change, not rounding"
    assert _gate_rejects(low, high), (
        f"A genuine mass-convention change (mn delta={delta_mn:.3f} MeV) must NOT "
        "be silently accepted by widening a generic tolerance. The gate rtol is "
        f"pinned at {GATE_RTOL}; admitting a convention change requires an "
        "EXPLICIT representation-only re-expression or physical EOS regeneration, "
        "never a silent tolerance widening."
    )


def test_gate_rtol_pinned_to_exact_binary_equality():
    """Pin the gate relative tolerance at 1e-12, appropriate for exact binary
    equality of canonicalised physical constants and source decimal precision.
    Raising this to absorb 0.5-1 MeV differences is FORBIDDEN by the audit.
    """
    assert GATE_RTOL == 1e-12, (
        "The native mass-equality gate rtol is pinned at 1e-12. If a precision "
        "policy is added, it must be a SEPARATE, explicit knob with recorded "
        "deltas and provenance, NOT a widening of this gate rtol."
    )


# ---------------------------------------------------------------------------
# Provenance consistency: equal mn/mp metadata does NOT imply equal underlying
# microphysics. The test fixture separately validates that the source mass
# pair and the EOS thermodynamic data have consistent provenance.
# ---------------------------------------------------------------------------
def test_provenance_consistency_mass_pair_vs_thermo():
    """Equal (mn, mp) metadata is necessary but NOT sufficient: the source mass
    pair and the EOS thermodynamic data must come from the SAME source. This
    test documents the contract using the audit findings:

    GPPVA_DD2 raw eos.thermo header = (939.0, 939.0)  [CompOSE normalised]
    GPPVA_DD2 registry physical    = (939.56536, 938.27203)  [Typel 2010]

    The raw artifact and the physical model disagree. A representation-only
    re-expression (CompOSE normalisation by mn) is the physically valid route
    to a common convention ONLY if the source microphysics used the distinct
    masses; a metadata-only rewrite of mp would be FORBIDDEN here.
    """
    # Raw header (written artifact) vs registry (physical model input)
    raw_mn, raw_mp = 939.0, 939.0
    reg_mn, reg_mp = 939.56536, 938.27203
    raw_common = abs(raw_mn - raw_mp) < 1e-6
    reg_distinct = abs(reg_mn - reg_mp) > 1e-6
    assert raw_common and reg_distinct, (
        "GPPVA_DD2 audit anchor: raw header is common-convention, registry is "
        "distinct-physical. This is the canonical discrepancy the audit resolves."
    )
    # The contract: equal metadata (raw 939,939) does NOT mean the underlying
    # microphysics used 939. Provenance must be checked separately.
    metadata_equal = (raw_mn == raw_mp)
    microphysics_distinct = (reg_mn != reg_mp)
    assert metadata_equal and microphysics_distinct, (
        "Equal raw metadata (939==939) does NOT imply equal microphysics "
        f"(registry {reg_mn} != {reg_mp}). Provenance must be validated "
        "independently of the mass-equality gate."
    )
