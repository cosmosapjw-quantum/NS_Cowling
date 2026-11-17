#!/usr/bin/env python3
"""Review the actually executed coarse public TF CLI calculation."""
from pathlib import Path
import sys,json,gzip,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.tf_composition import read_tf_composition
from compose_eos.tf_selection import select_geometry_candidates
base=ROOT/'data/research_finish_20261116/generation/actual_tf_cli';dispatch=json.loads((base/'DISPATCH.json').read_text())
out=base/'generated/GPPVA_FSU2H';receipt=json.loads((out/'RECEIPT.json').read_text())
data=read_tf_composition(out);native=data['thermodynamics']
assert dispatch['status']=='REQUESTED_ROWS_NATIVE_VERIFIED' and dispatch['returncode']==0
assert native.n==3 and np.allclose(native.nb,[.095,.0975,.1],rtol=0,atol=1e-15)
checks={n:hashlib.sha256((out/n).read_bytes()).hexdigest()==h for n,h in receipt['export']['output_sha256'].items()}
assert all(checks.values())
source={n:hashlib.sha256(Path(n).read_bytes()).hexdigest()==h for n,h in dispatch['source_sha256'].items()}
assert all(source.values())
F=[(3,'droplet'),(2,'droplet'),(1,'droplet'),(2,'hole'),(3,'hole'),(3,'uniform')]
nodes=[];attempt_count=0;failures=0
for points in [16,24]:
 records=json.loads(gzip.decompress((out/f'attempts_N{points:03d}.json.gz').read_bytes()))
 assert [r['n_b'] for r in records]==[.095,.0975,.1]
 for r in records:
  decision=select_geometry_candidates(r['attempts'],F,n_b=r['n_b'])
  selected=decision['selected'];assert selected is not None
  assert selected['model_id']=='GPPVA_FSU2H' and selected['points']==points
  attempt_count+=len(r['attempts']);failures+=sum(a.get('execution')=='FAILED' for a in r['attempts'])
  nodes.append(dict(points=points,n_B=r['n_b'],attempts=len(r['attempts']),
    selected_energy_per_baryon=selected['energy_per_baryon'],
    selected_status=selected['numerical_status'],KKT=selected['max_chemical_KKT_residual_MeV'],
    geometry=r['selected_geometry']))
# Native homogeneous fields are available only for explicitly uniform profiles.
assert native.Yn_present.tolist()==[False,True,True]
assert native.Yp_present.tolist()==[False,True,True]
assert native.Yh_present.tolist()==[False,True,True]
result=dict(scope='Actual end-to-end public TF CLI coarse smoke (N16/24), not N192/384 convergence or full-crust completion',
 native_rows=native.n,mesh_density_nodes=len(nodes),attempt_count=attempt_count,preserved_execution_failures=failures,
 phase_codes=native.Iphase.tolist(),native_free_masks=native.Yn_present.tolist(),
 baryon_integral_error=float(max(abs(data['fractions'][:,0]+data['fractions'][:,1]-1))),
 charge_integral_error=float(max(abs(data['fractions'][:,1]-data['fractions'][:,2]))),
 output_hash_matches=checks,dispatcher_source_hash_matches=source,nodes=nodes,
 runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
Path(__file__).with_name('ACTUAL_TF_CLI_REVIEW.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
