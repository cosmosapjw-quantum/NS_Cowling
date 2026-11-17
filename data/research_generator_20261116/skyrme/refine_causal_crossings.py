from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq
from compose_eos.skyrme_uniform import PARAMETERS,sound_speed_ratio
results=[]
for model in PARAMETERS:
    estimates=[]
    for step in (1e-3,5e-4,2.5e-4,1.25e-4):
        root=brentq(lambda n:sound_speed_ratio(model,n,relative_step=step)-1,1.,2.,xtol=1e-10)
        estimates.append({"relative_derivative_step":step,"causal_crossing_density_fm3":root})
    row={"model_id":model,"scope":"numerical cs2/c2=1 crossing of uniform npemu model branch",
         "crossings":estimates,"last_refinement_change_fm3":abs(estimates[-1]["causal_crossing_density_fm3"]-estimates[-2]["causal_crossing_density_fm3"]),
         "other_stability_channels":"spin and finite-wave-number stability unassessed"}
    results.append(row)
    print(json.dumps(row))
Path("build/skyrme_research/causal_crossings.json").write_text(json.dumps(results,indent=2)+"\n")
