#!/usr/bin/env python3
"""Generate limited diagnostic SPG candidate CLD samples, preserving all starts."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.api import read_compose
from compose_eos.candidate_gdfm_cld import MASS, solve
from compose_eos.gdfm_uniform import MODEL_IDS, beta_state

STARTS=[[.5,7.,.15,.12],[.8,10.,.16,.2],[.25,6.,.14,.3]]
TARGETS=[.005,.015,.03,.045,.055]
# A diagnostic residual threshold, not a source-fidelity tolerance.
GRADIENT_LIMIT=1e-4


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    result={'status':'RUNNING','scope':'candidate spherical CLD; pure neutron gas and electrons, no outer crust; finite starts only',
            'source_calibrated':False, 'gradient_limit':GRADIENT_LIMIT,'models':{}}
    dest=args.output/'candidate_cld.json'
    dest.write_text(json.dumps(result,indent=2)+'\n')
    for model in MODEL_IDS:
        source_dir=ROOT/'eos_generator/eos_data'/model
        source=read_compose(source_dir)
        valid=np.flatnonzero((source.Iphase==7)&(source.nb<.06))
        rows=[]
        for target in TARGETS:
            idx=int(valid[np.argmin(np.abs(source.nb[valid]-target))]); nb=float(source.nb[idx])
            best,attempts=solve(model,nb,STARTS)
            t,r,ni,y=best['search_coordinates']; lo=max(.06,nb*1.001)
            bounds=[(1e-5,.99999),(2.,40.),(lo,.25),(.005,.45)]
            active=[j for j,(v,(a,b)) in enumerate(zip((t,r,ni,y),bounds)) if min(v-a,b-v)<1e-5*(b-a)]
            grad=max(map(abs,best['transformed_gradient']))
            margins=best['state']; violations=[]
            if not best['success']: violations.append('OPTIMIZER_FAILED')
            if active: violations.append('SEARCH_BOUND_ACTIVE')
            if grad>GRADIENT_LIMIT: violations.append('STATIONARITY_UNRESOLVED')
            if margins['muon_insertion_margin_MeV']<0: violations.append('MUON_INSERTION_UNSTABLE')
            if margins['gas_proton_insertion_margin_MeV']<0: violations.append('GAS_PROTON_INSERTION_UNSTABLE')
            homogeneous=beta_state(model,nb,include_muons=True)
            epsilon=best['energy']
            pressure=nb*margins['mu_n_gas']-epsilon if not violations else None
            f,r,ni,y=best['coordinates']; yh=f*ni/nb
            rows.append({'source_row_0based':idx,'n_B':nb,'status':'INTERIOR_STATIONARY_CANDIDATE' if not violations else 'HOLD',
                'violations':violations,'search_bound_active':active,'max_abs_transformed_gradient':grad,
                'epsilon_MeV_fm3':epsilon,'pressure_envelope_MeV_fm3':pressure,
                'homogeneous_epsilon_MeV_fm3':homogeneous['epsilon'],
                'delta_epsilon_vs_homogeneous_MeV_fm3':epsilon-homogeneous['epsilon'],
                'composition':{'cluster_baryon_fraction':yh,'free_neutron_fraction':(1-f)*margins['n_gas']/nb,
                    'electron_fraction':margins['n_e']/nb,'muon_fraction_assumed':0.,'cluster_proton_fraction':y},
                'source':{'pressure_MeV_fm3':float(nb*source.Q1[idx]),
                    'epsilon_MeV_fm3':float(nb*source.mn*(1+source.Q6[idx])),
                    'electron_fraction':float(source.Ye[idx])},'best':best,'attempts':attempts})
        result['models'][model]={'source_sha256':{name:hashlib.sha256((source_dir/name).read_bytes()).hexdigest()
            for name in ['eos.nb','eos.thermo','eos.compo']},'rows':rows}
        dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    result['status']='COMPLETED_DIAGNOSTIC_ONLY'
    dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')

if __name__=='__main__': main()
