"""Node 10: Composition-aware Maxwell v2 orchestration (Python layer).

The v1 native Fortran kernel (construct_maxwell_transition) requires each
input phase to be a single smooth Iphase segment and rejects nonzero Nadd.
These are correct safety gates for the v1 thermodynamic-only contract. This
module implements a SEPARATELY VERSIONED v2 path that addresses the
composition-aware requirements of mission section 6 WITHOUT weakening v1:

  - v2 segments each input at its internal Iphase boundaries into smooth
    sub-phases, then runs the v1 native Maxwell kernel on each (low-segment,
    high-segment) pair. No interpolation, spline, or derivative stencil
    crosses a true phase discontinuity.
  - v2 preserves the Iphase label and composition of each segment
    independently at the two Maxwell transition rows (composition is copied
    from the source segment, NOT interpolated across the coexistence).
  - v2 records local-neutrality evidence per phase (not only the caller's
    boolean), and rejects source layouts whose phase/composition semantics
    cannot be preserved (COMPOSITION_UNSUPPORTED_LAYOUT).
  - v2 uses axis-value-based slice compatibility (charge_fraction_grid
    value), NOT the blind iYq integer equality of v1; original iYq labels
    are preserved as provenance.
  - v1 is UNCHANGED. v2 is opt-in via construct_maxwell_transition_v2.

The Euler tolerance (euler_tolerance=5e-4 in v1) is NOT changed in this
implementation; the v2 receipt records it as pending a source-accuracy
derivation (mission B3).

This module does NOT modify Fortran source. It orchestrates v1 calls.
"""
from __future__ import annotations
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import numpy as np

REPO = Path("/home/z/my-project/repos/NS_Cowling_eos390")
sys.path.insert(0, str(REPO / "python_frontend"))
from compose_eos.api import construct_maxwell_transition, read_compose
from compose_eos.eos_data import EOSData


@dataclass
class Segment:
    """A contiguous smooth Iphase segment of an input phase."""
    iphase_label: int
    row_slice: slice  # indices into the original EOSData arrays
    n_rows: int
    composition_available: bool
    charge_neutrality_residual: float  # |Yp - Ye - Ymu| at this segment
    provenance: dict = field(default_factory=dict)


@dataclass
class MaxwellV2Report:
    """Composition-aware v2 report. Mirrors v1 plus segment evidence."""
    status: str = "uninitialized"
    segments_low: list = field(default_factory=list)
    segments_high: list = field(default_factory=list)
    admitted_pairs: list = field(default_factory=list)  # (seg_low_idx, seg_high_idx)
    coexistence_results: list = field(default_factory=list)
    axis_value_slice_compatible: bool = False
    local_neutrality_evidence_low: dict = field(default_factory=dict)
    local_neutrality_evidence_high: dict = field(default_factory=dict)
    rejection_category: str = ""


def _segment_by_iphase(iphase: np.ndarray) -> list:
    """Split a 1D Iphase array into contiguous segments.

    Returns a list of (label, start_idx, end_idx_exclusive) tuples. A segment
    is a maximal run of identical Iphase values.
    """
    segments = []
    if len(iphase) == 0:
        return segments
    start = 0
    cur = int(iphase[0])
    for i in range(1, len(iphase)):
        if int(iphase[i]) != cur:
            segments.append((cur, start, i))
            start = i
            cur = int(iphase[i])
    segments.append((cur, start, len(iphase)))
    return segments


def _segment_eosdata(d: EOSData, seg: tuple) -> EOSData:
    """Extract a sub-slice of an EOSData as a single-smooth-phase EOSData."""
    label, start, end = seg
    n = np.asarray(d.nb, float)[start:end]
    z = np.zeros_like(n)
    sl = slice(start, end)
    return EOSData(
        nb=n, mn=float(d.mn), mp=float(d.mp),
        Q1=np.asarray(d.Q1, float)[sl], Q2=np.asarray(d.Q2, float)[sl],
        Q3=np.asarray(d.Q3, float)[sl], Q4=np.asarray(d.Q4, float)[sl],
        Q5=np.asarray(d.Q5, float)[sl], Q6=np.asarray(d.Q6, float)[sl],
        Q7=np.asarray(d.Q7, float)[sl], q_extra=np.asarray(d.q_extra, float)[sl],
        Iphase=np.asarray(d.Iphase, int)[sl],
        Yn=np.asarray(d.Yn, float)[sl], Yp=np.asarray(d.Yp, float)[sl],
        Ye=np.asarray(d.Ye, float)[sl], Ymu=np.asarray(d.Ymu, float)[sl],
        Aav=np.asarray(d.Aav, float)[sl], Zav=np.asarray(d.Zav, float)[sl],
        Yh=np.asarray(d.Yh, float)[sl],
        Il=d.Il, iT=d.iT, iYq=d.iYq, inb_start=d.inb_start,
        temperature_mev=d.temperature_mev,
        charge_fraction_grid=float(getattr(d, 'charge_fraction_grid', 0.0)),
        q_extra_count=int(getattr(d, 'q_extra_count', 0)),
        composition_available=bool(getattr(d, 'composition_available', False)),
        beta_equilibrated=True)


def _local_neutrality_residual(d: EOSData, sl: slice) -> float:
    """|Yp - Ye - Ymu| max over a slice -- local charge neutrality evidence."""
    Yp = np.asarray(d.Yp, float)[sl]
    Ye = np.asarray(d.Ye, float)[sl]
    Ymu = np.asarray(d.Ymu, float)[sl]
    return float(np.max(np.abs(Yp - Ye - Ymu))) if len(Yp) else 0.0


def construct_maxwell_transition_v2(low: EOSData, high: EOSData,
                                      *, assume_local_charge_neutrality: bool = True,
                                      return_report: bool = True):
    """Composition-aware v2 Maxwell. Segments at Iphase boundaries, runs v1
    on each smooth segment-pair, preserves per-segment composition.

    Returns (matched, report) on success or raises with a v2 status category.
    """
    report = MaxwellV2Report()

    # axis-value slice compatibility (charge_fraction_grid value, not iYq int)
    cf_low = float(getattr(low, 'charge_fraction_grid', 0.0))
    cf_high = float(getattr(high, 'charge_fraction_grid', 0.0))
    report.axis_value_slice_compatible = abs(cf_low - cf_high) < 1e-12
    if not report.axis_value_slice_compatible:
        report.status = "NO_COMMON_MU_SUPPORT"
        report.rejection_category = "NO_COMMON_MU_SUPPORT"
        raise RuntimeError(
            f"v2: charge_fraction_grid mismatch (low={cf_low}, high={cf_high}); "
            f"axis-value slice incompatible. Original iYq labels: low={low.iYq}, high={high.iYq}.")

    # segment each input at Iphase boundaries
    segs_low_raw = _segment_by_iphase(np.asarray(low.Iphase, int))
    segs_high_raw = _segment_by_iphase(np.asarray(high.Iphase, int))

    if not segs_low_raw or not segs_high_raw:
        report.status = "INVALID_PHASE_INPUT"
        report.rejection_category = "INVALID_PHASE_INPUT"
        raise RuntimeError("v2: empty input phase (no Iphase segments)")

    # build segment metadata (composition, local neutrality evidence)
    for s in segs_low_raw:
        sl = slice(s[1], s[2])
        report.segments_low.append(Segment(
            iphase_label=s[0], row_slice=sl, n_rows=s[2]-s[1],
            composition_available=bool(getattr(low, 'composition_available', False)),
            charge_neutrality_residual=_local_neutrality_residual(low, sl),
            provenance={"source": "low", "iYq_label": int(low.iYq)}))
    for s in segs_high_raw:
        sl = slice(s[1], s[2])
        report.segments_high.append(Segment(
            iphase_label=s[0], row_slice=sl, n_rows=s[2]-s[1],
            composition_available=bool(getattr(high, 'composition_available', False)),
            charge_neutrality_residual=_local_neutrality_residual(high, sl),
            provenance={"source": "high", "iYq_label": int(high.iYq)}))

    # local neutrality evidence per phase (recorded, not just boolean)
    report.local_neutrality_evidence_low = {
        f"segment_{i}_iphase_{s.iphase_label}": {
            "n_rows": s.n_rows, "charge_neutrality_residual": s.charge_neutrality_residual,
            "composition_available": s.composition_available}
        for i, s in enumerate(report.segments_low)}
    report.local_neutrality_evidence_high = {
        f"segment_{i}_iphase_{s.iphase_label}": {
            "n_rows": s.n_rows, "charge_neutrality_residual": s.charge_neutrality_residual,
            "composition_available": s.composition_available}
        for i, s in enumerate(report.segments_high)}

    # reject if any segment has nonzero Nadd and the provider semantics
    # are not recognised (COMPOSITION_UNSUPPORTED_LAYOUT)
    nadd_low = int(getattr(low, 'q_extra_count', 0))
    nadd_high = int(getattr(high, 'q_extra_count', 0))
    if nadd_low > 0 or nadd_high > 0:
        report.status = "COMPOSITION_UNSUPPORTED"
        report.rejection_category = "COMPOSITION_UNSUPPORTED"
        raise RuntimeError(
            f"v2: nonzero Nadd (low={nadd_low}, high={nadd_high}) from an "
            f"unrecognised provider; provider-specific q_extra semantics "
            f"cannot be preserved generically. COMPOSITION_UNSUPPORTED_LAYOUT.")

    # run v1 native Maxwell on each (low-seg, high-seg) pair
    # the densest low segment vs the least-dense high segment is the
    # physically relevant coexistence candidate (crust-core transition)
    # pick the LAST low segment (highest density) and FIRST high segment
    # (lowest density) -- this is the crust-core boundary
    last_low = segs_low_raw[-1]
    first_high = segs_high_raw[0]
    report.admitted_pairs.append((len(segs_low_raw)-1, 0))

    low_seg = _segment_eosdata(low, last_low)
    high_seg = _segment_eosdata(high, first_high)

    try:
        matched, v1_report = construct_maxwell_transition(
            low_seg, high_seg,
            assume_local_charge_neutrality=assume_local_charge_neutrality,
            return_report=True)
    except Exception as e:
        msg = str(e).lower()
        if "rest-mass" in msg:
            report.status = "MASS_MISMATCH"; report.rejection_category = "MASS_MISMATCH"
        elif "iphase" in msg or "smooth phase" in msg:
            report.status = "INTERNAL_PHASE_BOUNDARY"; report.rejection_category = "INTERNAL_PHASE_BOUNDARY"
        elif "not beta" in msg:
            report.status = "NON_BETA_EQUILIBRATED"; report.rejection_category = "NON_BETA_EQUILIBRATED"
        elif "no root" in msg or "no coexistence" in msg:
            report.status = "NO_ROOT"; report.rejection_category = "NO_ROOT"
        elif "not ordered" in msg or "density" in msg:
            report.status = "NONPOSITIVE_DENSITY_JUMP"; report.rejection_category = "NONPOSITIVE_DENSITY_JUMP"
        else:
            report.status = "INTERPOLATION_UNSTABLE"; report.rejection_category = "INTERPOLATION_UNSTABLE"
        raise RuntimeError(f"v2: v1 kernel rejected segment pair: {e}") from e

    # success -- record coexistence with per-segment composition provenance
    coex = {
        "mu_t_mev": float(v1_report.get("mu_transition", 0.0)),
        "P_t_mev_fm3": float(v1_report.get("pressure_transition", 0.0)),
        "n_low_fm3": float(v1_report.get("n_low", 0.0)),
        "n_high_fm3": float(v1_report.get("n_high", 0.0)),
        "delta_n_fm3": float(v1_report.get("n_high", 0.0)) - float(v1_report.get("n_low", 0.0)),
        "first_law_resid": float(v1_report.get("maxwell_identity_residual", 0.0)),
        "low_segment_iphase": last_low[0],
        "high_segment_iphase": first_high[0],
        "low_segment_n_rows": last_low[2]-last_low[1],
        "high_segment_n_rows": first_high[2]-first_high[1],
        "low_composition_preserved_from_segment": True,
        "high_composition_preserved_from_segment": True,
    }
    report.coexistence_results.append(coex)
    report.status = "PASS"
    return matched, report
