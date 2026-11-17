#!/usr/bin/env python3
"""Re-evaluate four frozen N384 seam profiles with the final band operator.

No nuclear density/radius minimization is performed. Only meson/electron
fields are recomputed for the stored neutron/proton profiles and radius.
All outputs remain in this fresh evidence directory.
"""
from pathlib import Path
import gzip,json,hashlib,platform
import numpy as np
import scipy
from compose_eos.spatial_tf_equilibrium import ReducedCellFunctional,RadialGrid

ROOT=Path(__file__).resolve().parents[3]
BASE=Path(__file__).resolve().parent
INPUT=ROOT/'data/research_finish_20261116/crust_assembly'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
paths=['python_frontend/compose_eos/spatial_tf_equilibrium.py',
       'python_frontend/compose_eos/spatial_tf.py',
       'python_frontend/compose_eos/gppva_rmf.py',
       'python_frontend/compose_eos/cold_leptons.py',
       str(Path(__file__).relative_to(ROOT))]
before={p:sha(ROOT/p) for p in paths}
if before['python_frontend/compose_eos/spatial_tf.py']!='2d1f1edc72b170b01194f2705e62905383a2d4292bd84cc58067bc2cc1067c02':
    raise RuntimeError('unexpected spatial TF operator identity')
if before['python_frontend/compose_eos/spatial_tf_equilibrium.py']!='8be3c6756f1b04a38cb59fdeeacb5c0aed5c98d8337c8ed14c64e7091501bdca':
    raise RuntimeError('unexpected screened electron operator identity')
old_evidence=INPUT/'final_operator_profile_verification.json'
old_evidence_sha=sha(old_evidence)
input_sha={str(old_evidence.relative_to(ROOT)):old_evidence_sha}
rows=[]
for model in ('GPPVA_DD2','GPPVA_TW','GPPVA_FSU2R','GPPVA_FSU2H'):
    path=INPUT/f'tf_interfaces/{model}/INTERFACE.json.gz'
    input_sha[str(path.relative_to(ROOT))]=sha(path)
    old=json.loads(gzip.decompress(path.read_bytes()))['tf_endpoint']
    grid=RadialGrid(old['points'],3);f=ReducedCellFunctional(model,grid,old['n_b'])
    old_nn=np.array(old['profile']['nn']);old_np=np.array(old['profile']['np'])
    x=np.r_[old_nn*grid.sqrtw/old['n_b'],old_np*grid.sqrtw/old['n_b'],np.log(old['radius_fm'])]
    nn,np_,ne,radius,muq,eres=f.unpack(x)
    s=f.evaluate(nn,np_,ne,radius)
    active=[nn>max(1e-10,old['n_b']*1e-7),np_>max(1e-10,old['n_b']*1e-7)]
    weights=np.r_[grid.weights[active[0]],grid.weights[active[1]]]
    chemical=np.r_[s['mu_n'][active[0]],(s['mu_p']+muq)[active[1]]]
    mu=float(weights@chemical/sum(weights));p=old['n_b']*mu-s['epsilon']
    kkt=[]
    for mask,chem in zip(active,[s['mu_n'],s['mu_p']+muq]):
        kkt.extend(abs(chem[mask]-mu));kkt.extend(np.maximum(mu-chem[~mask],0.))
    kkt.extend(abs(s['mu_e']-muq))
    r=dict(model=model,points=old['points'],original_interface_sha256=sha(path),
      n_b=old['n_b'],energy_difference_MeV_fm3=s['epsilon']-old['epsilon'],
      pressure_difference_MeV_fm3=p-old['pressure'],mu_difference_MeV=mu-old['mu_B'],
      nucleon_density_roundtrip_change_fm3=float(max(max(abs(nn-old_nn)),max(abs(np_-old_np)))),
      electron_density_max_change_fm3=float(max(abs(ne-np.array(old['profile']['ne'])))),
      electron_residual=eres,max_chemical_KKT_residual_MeV=float(max(kkt)),
      radial_residual_MeV=abs(s['d_epsilon_d_log_radius']/old['n_b']),
      baryon_constraint_error=abs(s['mean_baryon_density']-old['n_b']),
      charge_constraint_error=abs(s['mean_charge_density']),
      status='FROZEN_PROFILE_FINAL_BAND_OPERATOR_VERIFIED')
    if any(not np.isfinite(v) for v in r.values() if isinstance(v,(float,int))):
        raise RuntimeError('nonfinite reevaluation result')
    if (abs(r['energy_difference_MeV_fm3'])>1e-10 or abs(r['pressure_difference_MeV_fm3'])>1e-10
        or abs(r['mu_difference_MeV'])>1e-7 or max(kkt)>=.02 or r['radial_residual_MeV']>=2e-4
        or r['baryon_constraint_error']>1e-12 or r['charge_constraint_error']>1e-12
        or r['nucleon_density_roundtrip_change_fm3']>1e-15):
        raise RuntimeError(str(r))
    rows.append(r)
    print(model,r['energy_difference_MeV_fm3'],r['pressure_difference_MeV_fm3'],r['mu_difference_MeV'],flush=True)
after={p:sha(ROOT/p) for p in paths}
if before!=after:raise RuntimeError('code changed during verification')
if sha(old_evidence)!=old_evidence_sha:raise RuntimeError('previous sealed evidence changed')
root_receipt=INPUT/'tf_interfaces/RECEIPT.json'
input_sha[str(root_receipt.relative_to(ROOT))]=sha(root_receipt)
result={'claim':'Stored nuclear profiles/radii re-evaluated using final band meson operator and final electron operator; no nuclear minimization',
  'original_root_code_sha256':json.loads(root_receipt.read_text())['code_sha256'],
  'current_code_sha256':before,'source_unchanged_during_execution':True,
  'previous_sealed_evidence_unchanged':True,'input_sha256':input_sha,
  'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
  'rows':rows}
(BASE/'BAND_PROFILE_RECHECK.json').write_text(json.dumps(result,indent=2)+'\n')
snapshot={p:{'sha256':before[p],'content':(ROOT/p).read_text()} for p in paths}
with gzip.GzipFile(filename=str(BASE/'SOURCE_SNAPSHOT.json.gz'),mode='wb',mtime=0) as stream:
    stream.write((json.dumps(snapshot,indent=2)+'\n').encode())
print('COMPLETE',len(rows),'frozen profiles; prior sealed evidence unchanged',flush=True)
