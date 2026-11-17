from pathlib import Path
import os,sys,json,time,hashlib
from dataclasses import asdict
r=Path("/home/cosmosapjw/ns_cowling_wt/crust_refinement_20261117")
sys.path.insert(0,str(r/"python_frontend"))
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background,solve_mode
from compose_eos.stellar_cowling import CowlingSettings
o=r/"data/crust_refinement_20261117/stellar"
d=json.loads((o/"calculations/G_TM1_CANONICAL_N064__primary/BACKGROUND.json").read_text())
eos=SharpContractEnthalpyEOS.from_path(r/d["contract"])
bg=solve_sharp_background(eos,d["background"]["central_density_fm3"],**d["background"]["config"])
settings=CowlingSettings(rtol=3e-10,atol=3e-13,max_step_h=.0005,root_xtol_hz=1e-10)
t=time.time()
mode=solve_mode(bg,[176.,178.],settings,conversion="slow")
record=dict(settings=asdict(settings),mode=mode,wall_time_seconds=time.time()-t,
source_sha256=hashlib.sha256((r/"python_frontend/compose_eos/stellar_sharp.py").read_bytes()).hexdigest(),
scope="Bounded low-root control triggered by 6.93e-5 normalized join residual; same background/closure/surface cutoff")
(o/"LOW_ROOT_TIGHT_CONTROL.json").write_text(json.dumps(record,indent=2)+"\n")
print(json.dumps({k:v for k,v in mode.items() if k not in ["profile","interface_junctions"]}))
