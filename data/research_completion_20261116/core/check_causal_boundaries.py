#!/usr/bin/env python3
"""Density-resolved causal crossings of declared core functionals; no table gate."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'python_frontend'))
from compose_eos.core_regeneration import state_solver
records={}
for model,bracket in {'APR_APR':(.85,.91),'RG_KDE0v1':(1.39,1.44),'RG_SLY4':(1.24,1.29)}.items():
    solve,_=state_solver(model); trials=[]
    for rel in (2e-4,1e-4,5e-5):
        def cs2(n):
            h=n*rel
            dmu=(-solve(n+2*h)['mu_n']+8*solve(n+h)['mu_n']-8*solve(n-h)['mu_n']+solve(n-2*h)['mu_n'])/(12*h)
            return n*dmu/solve(n)['mu_n']
        n=brentq(lambda n:cs2(n)-1,*bracket,xtol=2e-12)
        trials.append({'relative_difference_step':rel,'n_fm3':n,'cs2_over_c2':cs2(n),
                       'cs2_below':cs2(n*.999),'cs2_above':cs2(n*1.001)})
    records[model]={'roots':trials,'root_spread_fm3':float(np.ptp([x['n_fm3'] for x in trials]))}
out={'meaning':'Cold equilibrated derivative c_s^2/c^2 = n_B dmu_B/dn_B / mu_B. Numerical local causal boundary, not a complete microscopic validity proof.',
     'method':'Fourth-order chemical-potential derivative and Brent root; three step sizes; no CompOSE values.',
     'models':records,'runner_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('CAUSAL_BOUNDARIES.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
print(json.dumps(records,indent=2))
