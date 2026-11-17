"""Adapter: native crust+CSS+core match -> phase-segmented TOV/Cowling input.

The native ``match-core css`` table contains exactly one first-order interface
(crust row ``css_transition_low_index`` -> CSS row ``css_transition_high_index``,
equal P and mu_B, positive density and energy jump) and one continuous
CSS->core join. Intake preserves the source Iphase labels verbatim. The
current two-provider solver rejects retained internal source boundaries before
constructing interpolants. Solver role labels 1/2 are assigned only on a
private copy after that check. MaxwellEOS
re-checks pressure and chemical-potential continuity and the latent-energy
identity itself.

Scope: the CLD inner crusts used by the census start at n_B >= 1e-3 fm^-3 and
carry no outer crust, so the TOV surface is the crust table's minimum
pressure and the reported radius is a truncation radius.  The tables carry no
composition; this adapter makes no abundance or elastic-crust claim.
"""
from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from compose_eos import EOSData, read_compose

from .maxwell import MaxwellEOS
from .css_barotrope import CSSCoreBarotrope

CSS_KIND = "css_single_speed_with_crust_density_jump"
PRESSURE_REL_GATE = 1.0e-10   # scripts/match_report_admission.PRESSURE_REL_GATE


class UnsupportedSourcePhaseBoundary(ValueError):
    """Intake retained a boundary for which the solver has no declared law."""

    category = "UNSUPPORTED_RETAINED_SOURCE_PHASE_BOUNDARY"

    def __init__(self, boundaries):
        self.boundaries = boundaries
        super().__init__(f"{self.category}: source boundary pairs {boundaries}; "
                         "phase labels do not specify an interpolation or interface law")


def phase_table_from_css_match(compose_dir: str | Path, match_report: Mapping[str, Any]) -> EOSData:
    """Validate CSS indices and return the table with its serialized labels."""
    if match_report.get("construction_kind") != CSS_KIND:
        raise ValueError("match report is not a css single-speed construction")
    data = read_compose(str(compose_dir))
    try:
        fields = [match_report[key] for key in
                  ("css_transition_low_index", "css_transition_high_index", "n")]
        if any(isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer))
               for value in fields):
            raise ValueError("css transition index fields must be exact integers")
        lo, hi, n_rows = map(int, fields)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("match report lacks css transition index fields") from exc
    if hi != lo + 1 or lo < 1 or hi >= data.n - 1 or n_rows != data.n:
        raise ValueError("css transition index fields are inconsistent with the table")
    # P is re-formed as n_B*Q1 on read-back, so the two jump rows can differ by
    # an ulp; use the declared admission pressure gate (1e-10 relative).
    p_scale = max(abs(float(data.p[lo])), abs(float(data.p[hi])))
    if (abs(float(data.p[hi]) - float(data.p[lo])) > PRESSURE_REL_GATE * p_scale
            or not (data.nb[hi] > data.nb[lo] and data.eps[hi] > data.eps[lo])):
        raise ValueError("rows at the css transition index are not a positive first-order jump")
    expected = match_report.get("output_phase_ids")
    valid = (isinstance(expected, list) and len(expected) == data.n
             and type(match_report.get("bridge_phase_id")) is int
             and all(type(tag) is int for tag in expected)
             and np.array_equal(data.Iphase, expected))
    if not valid:
        raise ValueError("CSS serialized phase layout disagrees with native receipt")
    return data


def css_match_maxwell_eos(compose_dir: str | Path, match_report: Mapping[str, Any]
                          ) -> tuple[MaxwellEOS, dict[str, Any]]:
    """Build the phase-segmented EOS and its scope metadata."""
    table = phase_table_from_css_match(compose_dir, match_report)
    join = match_report.get("css_core_bridge_index")
    if (isinstance(join, (bool, np.bool_)) or not isinstance(join, (int, np.integer))
            or not match_report["css_transition_high_index"] < join < table.n - 1):
        raise ValueError("invalid css core join index")
    if "bridge_phase_id" in match_report:
        tag = match_report["bridge_phase_id"]
        hi = match_report["css_transition_high_index"]
        if (type(tag) is not int or not np.all(table.Iphase[hi:join] == tag)
                or np.any(table.Iphase[:hi] == tag) or np.any(table.Iphase[join:] == tag)):
            raise ValueError("CSS bridge phase layout disagrees with retained source labels")
    lo = match_report["css_transition_low_index"]
    boundaries = [(i - 1, i) for start, end in ((0, lo + 1), (join, table.n))
                  for i in range(start + 1, end) if table.Iphase[i] != table.Iphase[i - 1]]
    if boundaries:
        raise UnsupportedSourcePhaseBoundary(boundaries)
    # Provider roles belong only to the solver's private view. The intake
    # table and CompOSE source retain every original label and boundary.
    source_phase_ids = table.Iphase.tolist()
    table = copy.deepcopy(table)
    table.Iphase = np.where(np.arange(table.n) <= lo, 1, 2).astype(int)
    high = table.Iphase == 2
    provider = CSSCoreBarotrope(
        table.p[high], table.eps[high], table.nb[high], table.mub[high],
        join_index=int(join - match_report["css_transition_high_index"]),
        sound_speed_sq=match_report.get("css_single_sound_speed", float("nan")))
    eos = MaxwellEOS(table, high_provider=provider)
    meta = {
        "construction_kind": CSS_KIND,
        "source_phase_ids": source_phase_ids,
        "source_phase_boundaries": boundaries,
        "crust_rows": int(np.count_nonzero(table.Iphase == 1)),
        "css_and_core_rows": int(np.count_nonzero(table.Iphase == 2)),
        "css_core_join": "continuous in n_B, P, eps, mu_B; one-sided sound-speed derivative",
        "css_provider": "exact analytic constant-sound-speed bridge",
        "core_provider": "separate PCHIP(eps-P versus ln P), first-law density",
        "css_core_join_pressure_mev_fm3": provider.join_pressure,
        "css_core_cs2_left": provider.sound_speed_sq,
        "css_core_cs2_right": provider.core.evaluate(provider.join_pressure).sound_speed_sq,
        "low_branch_first_law_anchor": "interface",
        "low_branch_surface_density_drift_relative": eos.low_branch_surface_density_drift_relative,
        "high_branch_end_density_drift_relative": eos.high_branch_end_density_drift_relative,
        "pressure_cutoff_mev_fm3": float(eos.low.pressure[0]),
        "perturbation_boundary": "Delta P=0 at finite cutoff pressure; artificial truncation boundary",
        "lapse_normalization": "vacuum Schwarzschild value at cutoff; complete exterior not modeled",
        "frequency_scope": "truncated fluid model; no complete-star observable or stability claim",
        "surface_definition": ("truncation at the crust table minimum pressure "
                               "(no outer crust in the CLD reconstruction)"),
        "composition": "none; no abundance or elastic-crust claim",
    }
    return eos, meta
