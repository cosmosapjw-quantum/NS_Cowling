"""Append an actual homogeneous core branch without rewriting a whole-crust EOS.

The caller binds model/source identity. This module checks the supplied native
anchor, cold slice, phase continuity and every continuous P(mu_E) interval.
It does not infer physical support from metadata or extrapolate a functional.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import replace

import numpy as np

from .crust_assembly import ROW_FIELDS, _cold_slice, rebase_header, select_rows
from .eos_data import EOSData
from .whole_crust import WholeCrustBarotrope


def _same(a, b, label):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    scale = np.maximum(np.maximum(abs(a), abs(b)), 1e-300)
    if a.shape != b.shape or np.any(abs(a-b) > 128*np.finfo(float).eps*scale):
        raise ValueError(label + " physical state mismatch")


def extend_whole_crust(whole: EOSData, contract: dict, core: EOSData,
                      source_anchor: EOSData):
    """Return a copied EOS/contract with actual core rows above the old endpoint.

    source_anchor is a separately identified native row at the original upper
    density. It must have the same physical state and present composition as
    the old endpoint. New rows must remain in that same source phase. Their
    source phase code is mapped to the existing whole-EOS phase code.
    """
    for data in (whole, core, source_anchor):
        _cold_slice(data)
    if source_anchor.n != 1:
        raise ValueError("source anchor must be exactly one native row")
    for data in (core, source_anchor):
        if data.Il != whole.Il or data.temperature_mev != whole.temperature_mev:
            raise ValueError("native core and whole EOS cold slice conventions differ")
        if np.any(data.Q5 != 0):
            raise ValueError("core extension requires the existing cold lepton convention")
    old_model = WholeCrustBarotrope.from_contract(contract)
    _same([whole.nb[-1], whole.p[-1], whole.eps[-1]],
          [old_model.upper_n[-1], old_model.upper_p[-1], old_model.upper_epsilon[-1]],
          "old table/contract")
    anchor_residuals = {}
    for name in ("nb", "p", "eps", "mub"):
        a, b = float(getattr(whole,name)[-1]), float(getattr(source_anchor,name)[0])
        _same(a,b,"source anchor "+name)
        anchor_residuals[name] = abs(a-b)
    for name in ("Yn","Yp","Ye","Ymu"):
        mask = name+"_present"
        if getattr(whole,mask)[-1] and getattr(source_anchor,mask)[0]:
            _same(getattr(whole,name)[-1],getattr(source_anchor,name)[0],
                  "source anchor composition "+name)
    indices = np.flatnonzero(core.nb > whole.nb[-1])
    if indices.size == 0:
        raise ValueError("actual core rows do not extend density support")
    selected = select_rows(core, indices)
    if np.any(selected.Iphase != source_anchor.Iphase[0]):
        raise ValueError("new native phase requires an explicit phase-boundary contract")
    # Header rebasing preserves physical energies; old row arrays are untouched.
    tail = rebase_header(selected, mn=whole.mn, mp=whole.mp)
    tail = replace(tail, Iphase=np.full(tail.n,whole.Iphase[-1],dtype=whole.Iphase.dtype))
    combined = replace(whole, source_dir=None,
        composition_file_present=whole.composition_file_present or tail.composition_file_present,
        composition_available=whole.composition_available or tail.composition_available,
        beta_equilibrated=whole.beta_equilibrated and tail.beta_equilibrated,
        **{name:np.concatenate((getattr(whole,name),getattr(tail,name))) for name in ROW_FIELDS})
    _cold_slice(combined)
    record = deepcopy(contract)
    for key, values in (("n",tail.nb),("p",tail.p),("epsilon",tail.eps)):
        record["upper"][key] = list(record["upper"][key])+values.tolist()
    model = WholeCrustBarotrope.from_contract(record)
    first = len(contract["upper"]["n"])-1
    checks = [dict(index=i,**model.upper_intervals[i].validate(),
                   contract=model.upper_intervals[i].contract())
              for i in range(first,len(model.upper_intervals))]
    # Exact old-array preservation is the required structural identity here.
    for name in ROW_FIELDS:
        if not np.array_equal(getattr(combined,name)[:whole.n],getattr(whole,name)):
            raise AssertionError("old row changed: "+name)
    return combined, record, dict(
        status="ACTUAL_CORE_TAIL_APPENDED_AND_CONTINUOUSLY_CHECKED",
        old_rows=whole.n, appended_rows=tail.n, total_rows=combined.n,
        native_source_row_indices=indices.tolist(),
        anchor_absolute_residuals=anchor_residuals,
        native_source_phase=int(source_anchor.Iphase[0]),
        retained_output_phase=int(whole.Iphase[-1]),
        old_row_fields_preserved=list(ROW_FIELDS),
        preserved_contract_sections=["surface","outer","lower_bridge"],
        old_density_domain_fm3=list(old_model.density_domain),
        density_domain_fm3=list(model.density_domain),
        pressure_domain_MeV_fm3=list(model.pressure_domain),
        extension_intervals=checks,
        extension_cs2_min=min(c["min_cs2_over_c2"] for c in checks),
        extension_cs2_max=max(c["max_cs2_over_c2"] for c in checks),
        energy_reference_shift="NONE",
        composition_scope="native knot fractions retained; continuous consumer is equilibrium thermo-only",
        support_scope="actual exported rows and explicit continuous representation; no high-density extrapolation")
