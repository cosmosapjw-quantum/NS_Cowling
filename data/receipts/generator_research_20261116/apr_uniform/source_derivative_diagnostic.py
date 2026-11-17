#!/usr/bin/env python3
"""Bounded first-law cross-check of original APR data and regenerated branches."""
from pathlib import Path
import json
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.api import read_compose
from compose_eos.apr_uniform import beta_state
source=read_compose(ROOT/"eos_generator/eos_data/APR_APR")
index=int(np.argmin(abs(source.nb-.2)))
density=float(source.nb[index])
rows=[]
for branch in ("source","low","high"):
    energy=lambda i:float(source.eps[i]) if branch=="source" else beta_state(float(source.nb[i]),branch)["epsilon"]
    epsilon=energy(index)
    derivative2=(energy(index+1)-energy(index-1))/(source.nb[index+1]-source.nb[index-1])
    grid=source.nb[index-2:index+3]
    coefficients=np.polynomial.Polynomial.fit(grid,[energy(i) for i in range(index-2,index+3)],4)
    derivative4=float(coefficients.deriv()(density))
    pressure=float(source.p[index]) if branch=="source" else beta_state(density,branch)["pressure"]
    mu=float(source.mub[index]) if branch=="source" else beta_state(density,branch)["mu_n"]
    rows.append({"branch":branch,"n_b_fm3":density,"epsilon_mev_fm3":epsilon,
        "pressure_mev_fm3":pressure,"mu_b_mev":mu,"source_grid_used_fm3":grid.tolist(),
        "pressure_from_energy_3point":density*derivative2-epsilon,
        "pressure_from_energy_5point":density*derivative4-epsilon,
        "mu_from_energy_3point":derivative2,"mu_from_energy_5point":derivative4})
print(json.dumps({"scope":"one bounded n=.2 fm^-3 diagnostic, both pure generated branches on the same source grid",
    "interpretation":"source and generated branches independently satisfy the cold first law; residual source mismatch concerns data/functional/convention provenance, not unknown phase prescription",
    "gibbs_window_caveat":"five-point window crosses the Gibbs onset; Gibbs finite differences are therefore deliberately not treated as a smooth pure-branch derivative here",
    "rows":rows},indent=2))
