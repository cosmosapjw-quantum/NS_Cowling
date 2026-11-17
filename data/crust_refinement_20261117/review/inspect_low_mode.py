"""Read-only diagnostic of the sampled low-frequency mode, not an acceptance test."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.stellar_sharp import SharpContractEnthalpyEOS
from compose_eos.stellar_cowling import C_KM_S
b=ROOT/'data/crust_refinement_20261117/stellar/calculations/G_TM1_CANONICAL_N064__primary'
mode=json.loads((b/'slow_ROOT_00.json').read_text());bg=json.loads((b/'BACKGROUND.json').read_text())
eos=SharpContractEnthalpyEOS.from_path(ROOT/bg['contract'])
phi0=.5*np.log1p(-2*bg['background']['mass_geom_km']/bg['background']['radius_km'])
def point(q):
 s=eos.state_at_enthalpy(q['h'],q['side'])
 fL=C_KM_S/(2*np.pi)*np.sqrt(6*s['cs2_over_c2'])*np.exp(phi0-q['h'])/q['r_km']
 return dict(**q,n_fm3=s['n'],cs2_over_c2=s['cs2_over_c2'],local_l2_Lamb_frequency_hz=float(fL))
pairs=[dict(left=point(a),right=point(c)) for a,c in zip(mode['profile'][:-1],mode['profile'][1:]) if a['W']*c['W']<0]
sample=[point(q) for q in mode['profile']]
result=dict(scope='Physical diagnostics for one sampled provisional root; W-node count alone does not assign mode class',
 frequency_hz=mode['frequency_hz'],matching_residual=mode['matching_residual'],
 join_mismatch=mode['relative_join_mismatch'],sign_changes=pairs,
 sampled_min_cs2=min(q['cs2_over_c2'] for q in sample),
 source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [b/'slow_ROOT_00.json',b/'BACKGROUND.json',ROOT/bg['contract']]})
Path(__file__).with_name('LOW_ROOT_INSPECTION.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
