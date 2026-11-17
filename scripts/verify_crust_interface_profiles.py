#!/usr/bin/env python3
"""Re-evaluate prior seam profiles with the final electron implementation."""
from pathlib import Path
import gzip,json,hashlib
import numpy as np
from compose_eos.spatial_tf_equilibrium import ReducedCellFunctional,RadialGrid
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'data/research_finish_20261116/crust_assembly'
rows=[]
for model in ('GPPVA_DD2','GPPVA_TW','GPPVA_FSU2R','GPPVA_FSU2H'):
    path=BASE/f'tf_interfaces/{model}/INTERFACE.json.gz'
    old=json.loads(gzip.decompress(path.read_bytes()))['tf_endpoint']
    grid=RadialGrid(old['points'],3);f=ReducedCellFunctional(model,grid,old['n_b'])
    nn=np.array(old['profile']['nn']);np_=np.array(old['profile']['np'])
    x=np.r_[nn*grid.sqrtw/old['n_b'],np_*grid.sqrtw/old['n_b'],np.log(old['radius_fm'])]
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
    r=dict(model=model,original_interface_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
      n_b=old['n_b'],energy_difference_MeV_fm3=s['epsilon']-old['epsilon'],
      pressure_difference_MeV_fm3=p-old['pressure'],mu_difference_MeV=mu-old['mu_B'],
      electron_density_max_change_fm3=float(max(abs(ne-np.array(old['profile']['ne'])))),
      electron_residual=eres,max_chemical_KKT_residual_MeV=float(max(kkt)),
      radial_residual_MeV=abs(s['d_epsilon_d_log_radius']/old['n_b']),
      baryon_constraint_error=abs(s['mean_baryon_density']-old['n_b']),
      charge_constraint_error=abs(s['mean_charge_density']),
      status='FROZEN_PROFILE_CURRENT_OPERATOR_VERIFIED')
    if abs(r['energy_difference_MeV_fm3'])>1e-10 or abs(r['mu_difference_MeV'])>1e-7 or max(kkt)>=.02 or r['radial_residual_MeV']>=2e-4:
        raise RuntimeError(str(r))
    rows.append(r);print(model,r['energy_difference_MeV_fm3'],r['mu_difference_MeV'],flush=True)
paths=['python_frontend/compose_eos/spatial_tf_equilibrium.py','python_frontend/compose_eos/spatial_tf.py','python_frontend/compose_eos/gppva_rmf.py','scripts/verify_crust_interface_profiles.py']
result={'claim':'same stored profiles evaluated with final conservative-flux electron operator; no new minimization',
  'original_root_code_sha256':json.loads((BASE/'tf_interfaces/RECEIPT.json').read_text())['code_sha256'],
  'current_code_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},'rows':rows}
(BASE/'final_operator_profile_verification.json').write_text(json.dumps(result,indent=2)+'\n')
