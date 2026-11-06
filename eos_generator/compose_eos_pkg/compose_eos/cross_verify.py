r"""Cross-verification and production EoS builder for CompOSE tables.

Provides automated comparison between a thermodynamically reconstructed EoS table
and a first-principles literature-regenerated EoS table across all three stellar
regimes (outer crust, inner crust / pasta, and homogeneous core).
"""
from __future__ import annotations

from typing import Dict

import numpy as np

from .consistency import check_consistency
from .eos_data import EOSData
from .literature_solvers import (
    DD_MODELS,
    NL_MODELS,
    SPG_GDFM_PARAMS,
    dd_bulk_thermo,
    gdfm_bulk_thermo,
    lepton_thermo,
    nl_bulk_thermo,
    outer_crust_pressure,
    ME,
    MMU,
)


def cross_verify_eos(rec: EOSData, lit: EOSData, dirname: str = "") -> Dict[str, object]:
    """Compare a reconstructed :class:`EOSData` against its literature counterpart."""
    rep_rec = check_consistency(rec, verbose=False)
    rep_lit = check_consistency(lit, verbose=False)

    outer = np.where((rec.Aav > 0) & (rec.Yn == 0))[0]
    inner = np.where((rec.Aav > 0) & (rec.Yn > 0))[0]
    core = np.where(rec.is_homogeneous)[0]

    def _med_rel(a: np.ndarray, b: np.ndarray) -> float:
        if a.size == 0:
            return 0.0
        den = np.where(np.abs(b) > 1e-30, np.abs(b), 1.0)
        return float(np.median(np.abs(a - b) / den))

    return dict(
        dirname=dirname,
        n_points=rec.n,
        rec_C2_max=rep_rec["C2"]["max_rel"],
        lit_C2_max=rep_lit["C2"]["max_rel"],
        rec_C3_max=max(rep_rec["C3"]["crust"], rep_rec["C3"]["core"]),
        rec_C4_max=max(rep_rec["C4"]["crust"], rep_rec["C4"]["core"]),
        rec_C5_maxmul=rep_rec["C5"]["maxmul"],
        outer_crust_med_p_rel=_med_rel(rec.p[outer], lit.p[outer]),
        inner_crust_med_p_rel=_med_rel(rec.p[inner], lit.p[inner]),
        core_med_p_rel=_med_rel(rec.p[core], lit.p[core]),
        core_med_eps_rel=_med_rel(rec.eps[core], lit.eps[core]),
    )


__all__ = ["cross_verify_eos"]
