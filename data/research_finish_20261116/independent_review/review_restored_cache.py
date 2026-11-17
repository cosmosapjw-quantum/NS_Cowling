#!/usr/bin/env python3
"""Independent cache-poisoning and constrained-curvature counterexamples."""
from pathlib import Path
from dataclasses import FrozenInstanceError
import sys,json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.restored_spg_candidate import RestoredSPGCandidate,_local_minimum_check
from compose_eos.candidate_gdfm_cld import SurfaceParameters
module=ROOT/'python_frontend/compose_eos/restored_spg_candidate.py';sha=hashlib.sha256(module.read_bytes()).hexdigest()
m=RestoredSPGCandidate('SPG_M1');original=m.state(.002);frozen=False
try:m.surface=SurfaceParameters('SPG_M1',1.2,20.,'alternate','f'*64)
except FrozenInstanceError:frozen=True
original['energy_per_baryon']=-1234.;original['all_three_attempts'][0]['coordinates'][0]=-1234.
new=m.state(.002)
result={'source_sha256':sha,'configuration_frozen':bool(frozen),'returned_dict_poison_prevented':bool(new['energy_per_baryon']!=-1234.),'nested_dict_poison_prevented':bool(new['all_three_attempts'][0]['coordinates'][0]!=-1234.)}
Path(__file__).with_name('RESTORED_CACHE_GREEN.json').write_text(json.dumps(result,indent=2)+'\n');print(result)
lo=np.array([1e-7,.035,.001]);hi=np.array([1.,.25,.499]);z=np.array([1.,.1,.1])
g=lambda z:np.array([1-z[0],z[1]-.1,z[2]-.1])
try:
 answer=_local_minimum_check(g,z,lo,hi);boundary={'saddle_accepted':True,'reported_stability':answer}
except RuntimeError as exc:boundary={'saddle_accepted':False,'error':str(exc)}
boundary.update(source_sha256=sha,energy='-(t-1)^2/2+(ni-.1)^2/2+(yi-.1)^2/2',gradient_at_boundary=g(z).tolist(),inward_energy_change_at_dt_minus_0p001=-.000001/2)
out=Path(__file__).with_name('MARGINAL_BOUNDARY_CURVATURE_REVIEW.json')
if out.exists() and json.loads(out.read_text()).get('saddle_accepted'):
 out.rename(out.with_name('MARGINAL_BOUNDARY_CURVATURE_RED.json'))
out.write_text(json.dumps(boundary,indent=2)+'\n');print(boundary)
assert all(v for k,v in result.items() if k!='source_sha256')

assert not boundary["saddle_accepted"], "marginal-boundary inward saddle admitted"
