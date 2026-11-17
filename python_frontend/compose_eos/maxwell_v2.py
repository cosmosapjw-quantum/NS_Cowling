"""Segmented thermo-only Maxwell construction with lossless retained tails.

Only the last smooth low-density segment and first smooth high-density segment
are submitted to the native Maxwell kernel. Remaining outer-crust and high-core
rows are retained verbatim in physical units, with separately recorded inherited
phase boundaries. Those boundaries are not silently interpolated or certified
as additional Maxwell transitions. Composition must be explicitly discarded.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any
import numpy as np

from compose_eos.api import construct_maxwell_transition, _thermo_only_view
from compose_eos.eos_data import EOSData
from compose_eos.native import NativeExecutionError

_ROW_FIELDS = (
    "nb", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "q_extra", "Iphase",
    "Yn", "Yp", "Ye", "Ymu", "Aav", "Zav", "Yh", "row_status",
    "Yn_present", "Yp_present", "Ye_present", "Ymu_present", "Abar_present",
    "Zbar_present", "Yh_present", "C3_applicable", "C4_applicable",
)


@dataclass
class Segment:
    """One contiguous source Iphase segment; missing composition is not zero."""
    iphase_label: int
    row_slice: slice
    n_rows: int
    composition_available: bool
    charge_neutrality_residual: float | None
    provenance: dict = field(default_factory=dict)


@dataclass
class MaxwellV2Report:
    """Native interface evidence, source retention, and uncertified boundaries."""
    status: str = "uninitialized"
    segments_low: list = field(default_factory=list)
    segments_high: list = field(default_factory=list)
    attempted_pairs: list = field(default_factory=list)
    admitted_pairs: list = field(default_factory=list)
    coexistence_results: list = field(default_factory=list)
    axis_value_slice_compatible: bool = False
    local_neutrality_evidence_low: dict = field(default_factory=dict)
    local_neutrality_evidence_high: dict = field(default_factory=dict)
    rejection_category: str = ""
    native_report: dict = field(default_factory=dict)
    transition_row_indices_zero_based: list = field(default_factory=list)
    output_phase_map: list = field(default_factory=list)
    retained_source_rows: dict = field(default_factory=dict)
    inherited_interfaces: list = field(default_factory=list)
    inherited_interfaces_status: str = "NONE"
    validation_scope: str = "native_boundary_segments_and_new_maxwell_interface"
    composition_dropped: bool = False
    source_thermodynamic_rows_preserved: bool = False


class MaxwellV2Error(NativeExecutionError):
    """A structured failure retaining the v2 evidence accumulated so far."""
    def __init__(self, message: str, *, report: MaxwellV2Report, category: str):
        super().__init__(message, category=category, operation="maxwell-transition-v2")
        self.report = report


def _fail(report, category, message):
    report.status = report.rejection_category = category
    raise MaxwellV2Error("v2: " + message, report=report, category=category)


def _segment_by_iphase(iphase: np.ndarray) -> list:
    if len(iphase) == 0:
        return []
    breaks = np.r_[0, np.flatnonzero(np.diff(iphase) != 0)+1, len(iphase)]
    return [(int(iphase[a]), int(a), int(b)) for a, b in zip(breaks[:-1], breaks[1:])]


def _segment_eosdata(d: EOSData, seg: tuple) -> EOSData:
    """Slice rows without losing support seals, masks, axes or beta evidence."""
    _, start, end = seg
    out = d.copy()
    for name in _ROW_FIELDS:
        setattr(out, name, np.asarray(getattr(d, name))[start:end].copy())
    for name in ("q_full", "q_extra_full"):
        arr = np.asarray(getattr(d, name))
        if arr.size:
            setattr(out, name, arr[..., start:end].copy())
    out.inb_start = d.inb_start + start
    return out


def _local_neutrality_residual(d: EOSData, sl: slice) -> float | None:
    """Evaluate total charge only where the complete nuclear/lepton inventory exists."""
    applicable = np.asarray(d.C4_applicable, bool)[sl]
    if not d.composition_available or not applicable.size or not np.all(applicable):
        return None
    yp, ye, ym = (np.asarray(getattr(d, name), float)[sl] for name in ("Yp", "Ye", "Ymu"))
    a, z, yh = (np.asarray(getattr(d, name), float)[sl] for name in ("Aav", "Zav", "Yh"))
    if np.any(a <= 0):
        return None
    return float(np.max(np.abs(yp+(z/a)*yh-ye-ym)))


def _physical_reencode(data: EOSData, mn: float, mp: float) -> EOSData:
    """Change the serialization reference, preserving each physical column."""
    out = data.copy()
    if data.mn != mn:
        for name in ("Q3", "Q6", "Q7"):
            setattr(out, name, (getattr(data, name)+1.)*data.mn/mn-1.)
        for name in ("Q4", "Q5"):
            setattr(out, name, getattr(data, name)*data.mn/mn)
    out.mn, out.mp = mn, mp
    return out


def _retain_tails(low, high, matched, raw_low, raw_high, report):
    prefix = _segment_eosdata(low, (0, 0, raw_low[-1][1]))
    suffix = _segment_eosdata(high, (0, raw_high[0][2], high.n))
    prefix = _physical_reencode(prefix, matched.mn, matched.mp)
    suffix = _physical_reencode(suffix, matched.mn, matched.mp)
    nlow_segments = len(raw_low)
    report.output_phase_map = []
    for side, segments, offset in (("low", raw_low, 0), ("high", raw_high, nlow_segments)):
        for i, (label, start, end) in enumerate(segments):
            report.output_phase_map.append({
                "output_iphase": offset+i+1, "source": side,
                "source_iphase": label, "source_rows_zero_based": [start, end-1],
                "boundary_segment": (side == "low" and i == len(segments)-1) or (side == "high" and i == 0),
            })
    for i, (_, a, b) in enumerate(raw_low[:-1]):
        prefix.Iphase[a:b] = i+1
    for i, (_, a, b) in enumerate(raw_high[1:], start=1):
        suffix.Iphase[a-raw_high[0][2]:b-raw_high[0][2]] = nlow_segments+i+1
    middle = matched.copy()
    middle.Iphase = np.where(matched.Iphase == 1, nlow_segments, nlow_segments+1)
    out = matched.copy()
    for name in _ROW_FIELDS:
        setattr(out, name, np.concatenate([getattr(part, name) for part in (prefix, middle, suffix)]))
    out.inb_start = 1
    out.source_dir = None  # a composite is not a sealed source directory
    out.composition_file_present = True
    out.composition_available = False
    out.row_status[:] = 1  # phase-only, no abundance inventory
    out.beta_equilibrated = bool(low.beta_equilibrated and high.beta_equilibrated)
    if np.any(np.diff(out.nb) <= 0) or not np.all(np.isfinite(out.nb)):
        _fail(report, "INVALID_RETAINED_SOURCE_ORDER", "retained source tails do not form an ordered density table")
    first = prefix.n + int(report.native_report["low_phase_rows"]) - 1
    report.transition_row_indices_zero_based = [first, first+1]
    report.retained_source_rows = {"low_prefix": prefix.n, "high_suffix": suffix.n}
    boundaries = np.flatnonzero(np.diff(out.Iphase) != 0)
    report.inherited_interfaces = [{
        "row_indices_zero_based": [int(i), int(i+1)],
        "left_output_iphase": int(out.Iphase[i]), "right_output_iphase": int(out.Iphase[i+1]),
        "status": "PRESERVED_SOURCE_BOUNDARY_NO_CROSS_PHASE_INTERPOLATION",
    } for i in boundaries if i != first]
    report.inherited_interfaces_status = (
        "PRESERVED_SOURCE_BOUNDARIES_NO_CROSS_PHASE_INTERPOLATION" if report.inherited_interfaces else "NONE")
    report.source_thermodynamic_rows_preserved = True
    return out


def construct_maxwell_transition_v2(low: EOSData, high: EOSData, *,
                                      assume_local_charge_neutrality: bool = False,
                                      drop_composition: bool = False,
                                      return_report: bool = True):
    """Match boundary segments and retain all other supplied source rows.

    Composition is dropped only by explicit request. A successful report proves
    the new interface and its two native boundary segments; inherited source
    interfaces remain explicitly unvalidated. It does not certify a surface EOS
    or a stellar solver's support for all phase boundaries.
    """
    report = MaxwellV2Report()
    if not isinstance(low, EOSData) or not isinstance(high, EOSData):
        raise TypeError("low and high must both be EOSData instances")
    if not isinstance(drop_composition, bool) or not isinstance(return_report, bool):
        raise TypeError("drop_composition and return_report must be bool")
    if not isinstance(assume_local_charge_neutrality, bool):
        raise TypeError("assume_local_charge_neutrality must be bool")
    if not assume_local_charge_neutrality:
        _fail(report, "LOCAL_NEUTRALITY_NOT_CONFIRMED", "caller must confirm local charge neutrality")
    has_composition = any(bool(getattr(d, flag)) for d in (low, high)
                          for flag in ("composition_available", "composition_file_present"))
    if has_composition and not drop_composition:
        _fail(report, "COMPOSITION_UNSUPPORTED", "COMPOSITION_UNSUPPORTED: native kernel is thermo-only; explicit drop_composition=True is required")
    cf_low, cf_high = float(low.charge_fraction_grid), float(high.charge_fraction_grid)
    report.axis_value_slice_compatible = (low.Il == high.Il and abs(cf_low-cf_high) < 1e-12)
    if not report.axis_value_slice_compatible:
        _fail(report, "SLICE_CONVENTION_MISMATCH", f"charge_fraction_grid mismatch or Il mismatch (low={cf_low}/{low.Il}, high={cf_high}/{high.Il}); axis-value slice incompatible")
    raw_low = _segment_by_iphase(np.asarray(low.Iphase, int))
    raw_high = _segment_by_iphase(np.asarray(high.Iphase, int))
    if not raw_low or not raw_high:
        _fail(report, "INVALID_PHASE_INPUT", "empty input phase (no Iphase segments)")
    for side, data, segments in (("low", low, raw_low), ("high", high, raw_high)):
        target = getattr(report, "segments_"+side)
        for label, start, end in segments:
            sl = slice(start, end)
            target.append(Segment(label, sl, end-start, data.composition_available,
                _local_neutrality_residual(data, sl),
                {"source": side, "iYq_label": data.iYq, "source_dir": data.source_dir}))
        setattr(report, "local_neutrality_evidence_"+side, {
            f"segment_{i}_iphase_{s.iphase_label}": {
                "n_rows": s.n_rows, "charge_neutrality_residual": s.charge_neutrality_residual,
                "composition_available": s.composition_available,
                "basis": "complete_inventory" if s.charge_neutrality_residual is not None else "unavailable_caller_assumption_only",
            } for i, s in enumerate(target)})
    if low.q_extra_count > 0 or high.q_extra_count > 0:
        _fail(report, "COMPOSITION_UNSUPPORTED", "COMPOSITION_UNSUPPORTED_LAYOUT: nonzero Nadd requires explicit provider semantics")
    low_input = _thermo_only_view(low) if drop_composition else low
    high_input = _thermo_only_view(high) if drop_composition else high
    report.composition_dropped = has_composition and drop_composition
    selected = (len(raw_low)-1, 0)
    report.attempted_pairs.append(selected)
    try:
        matched, native = construct_maxwell_transition(
            _segment_eosdata(low_input, raw_low[-1]),
            _segment_eosdata(high_input, raw_high[0]),
            assume_local_charge_neutrality=True, return_report=True)
    except NativeExecutionError as exc:
        _fail(report, exc.category or "UNKNOWN", f"native boundary-segment construction rejected: {exc}")
    except (ValueError, TypeError) as exc:
        _fail(report, "INPUT_CONTRACT_FAILURE", f"source support/input contract rejected: {exc}")
    required = ("mu_transition", "pressure_transition", "n_low", "n_high",
                "maxwell_identity_residual", "low_phase_rows", "high_phase_rows")
    if any(name not in native for name in required):
        _fail(report, "MALFORMED_NATIVE_REPORT", "required coexistence fields missing from native report")
    if native.get("retained_thermodynamic_rows") != "physical_columns_preserved_without_reconstruction":
        _fail(report, "MALFORMED_NATIVE_REPORT", "native backend does not declare lossless retained-row semantics; rebuild the current backend")
    report.native_report = dict(native)
    out = _retain_tails(low_input, high_input, matched, raw_low, raw_high, report)
    report.coexistence_results.append({
        "mu_t_mev": float(native["mu_transition"]),
        "P_t_mev_fm3": float(native["pressure_transition"]),
        "n_low_fm3": float(native["n_low"]), "n_high_fm3": float(native["n_high"]),
        "delta_n_fm3": float(native["n_high"])-float(native["n_low"]),
        "first_law_resid": float(native["maxwell_identity_residual"]),
        "low_segment_iphase": raw_low[-1][0], "high_segment_iphase": raw_high[0][0],
        "low_segment_n_rows": raw_low[-1][2]-raw_low[-1][1],
        "high_segment_n_rows": raw_high[0][2]-raw_high[0][1],
        "low_composition_preserved_from_segment": False,
        "high_composition_preserved_from_segment": False,
        "composition_status": "NOT_AVAILABLE_THERMO_ONLY", "output_composition_available": False,
    })
    report.admitted_pairs.append(selected)
    report.status = "PASS"
    return (out, report) if return_report else out
