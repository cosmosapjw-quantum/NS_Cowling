from pathlib import Path
import os,sys,json,time,hashlib
from dataclasses import asdict
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/"python_frontend"))
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS,solve_sharp_background,solve_mode
from compose_eos.stellar_cowling import CowlingSettings
OUT=ROOT/"data/crust_refinement_20261117/stellar"
d=json.loads((OUT/"calculations/G_TM1_CANONICAL_N032__primary/BACKGROUND.json").read_text())
eos=SharpContractEnthalpyEOS.from_path(ROOT/d["contract"])
bg=solve_sharp_background(eos,d["background"]["central_density_fm3"],**d["background"]["config"])
settings=CowlingSettings(rtol=3e-10,atol=3e-13,max_step_h=.0005,root_xtol_hz=1e-10)
results=[]
for conversion,bracket in [("slow",[176.,178.]),("rapid",[1903.,1905.])]:
    t=time.time();mode=solve_mode(bg,bracket,settings,conversion=conversion)
    record=dict(settings=asdict(settings),mode=mode,wall_time_seconds=time.time()-t,
        source_sha256=hashlib.sha256((ROOT/"python_frontend/compose_eos/stellar_sharp.py").read_bytes()).hexdigest(),
        scope="Only N32 roots outside 1e-6 join gate; same background/closure/surface cutoff")
    (OUT/("N32_"+conversion.upper()+"_TIGHT_CONTROL.json")).write_text(json.dumps(record,indent=2)+"\n")
    results.append({k:v for k,v in mode.items() if k not in ["profile","interface_junctions"]})
    print(json.dumps(results[-1]),flush=True)
