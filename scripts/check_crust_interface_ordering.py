#!/usr/bin/env python3
"""Two-branch local Gibbs ordering from actually reoptimized bracketing rows."""
from pathlib import Path
import gzip,json
import numpy as np
from scipy.optimize import brentq
from compose_eos.bsk22_extended_representation import ExtendedAnchoredBSk22
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'data/research_finish_20261116/crust_assembly'
d=json.loads((ROOT/'data/research_completion_20261116/outer_crust/generated/RECEIPT.json').read_text())['drip']
bsk=ExtendedAnchoredBSk22(outer_density=d['n_b'],outer_pressure=d['pressure'],outer_mu_b=d['g'],outer_epsilon=d['epsilon'])
records=[]
for kind,models in [('tf',[f'GPPVA_{m}' for m in ['DD2','TW','FSU2R','FSU2H']]),
                    ('restored_spg',[f'SPG_M{i}' for i in range(1,6)])]:
    for model in models:
        if kind=='tf':
            p=BASE/f'tf_interfaces/{model}'
            receipt=json.loads(gzip.decompress((p/'INTERFACE.json.gz').read_bytes()))
            endpoint=receipt['tf_endpoint'];nk='n_b';mk='mu_B';pk='pressure'
            rows=[r for r in json.loads(gzip.decompress((p/'attempts.json.gz').read_bytes()))
                  if r['points']==384 and r['numerical_status']=='STATIONARY']
        else:
            p=BASE/f'restored_spg/{model}'
            receipt=json.loads((p/'RECEIPT.json').read_text())
            endpoint=receipt['spg_endpoint'];nk='n_B';mk='mu_B';pk='pressure'
            rows=json.loads(gzip.decompress((p/'all_reoptimized_states.json.gz').read_bytes()))
        n=endpoint[nk]
        left=max((r for r in rows if r[nk]<n*(1-5e-5)),key=lambda r:r[nk])
        right=min((r for r in rows if r[nk]>n*(1+5e-5)),key=lambda r:r[nk])
        sides=[]
        for side,row in [('left',left),('right',right)]:
            bn=brentq(lambda q:bsk.state(q)['pressure']-row[pk],bsk.n0,.01,xtol=1e-14)
            b=bsk.state(bn)
            sides.append(dict(side=side,n_inner=row[nk],n_bsk=bn,pressure=row[pk],
                delta_mu_MeV=row[mk]-b['mu_b_euler']))
        slope=(right[pk]-left[pk])/(right[nk]-left[nk])
        bn=receipt['bsk_endpoint']['n_b']
        passed=sides[0]['delta_mu_MeV']>0>sides[1]['delta_mu_MeV'] and slope>0 and n>bn
        if not passed:raise RuntimeError(str((model,sides,slope)))
        records.append(dict(model=model,variant=kind,status='TWO_BRANCH_LOCAL_ORDERING_NUMERICALLY_CHECKED',
            sides=sides,inner_dP_dn_finite_difference=slope,
            bsk_dP_dn_at_interface=bsk.state(bn)['pressure_density_derivative'],
            analytic_d_delta_mu_dP_at_interface=1/n-1/bn,
            limitation='Selected branch pair only. Global lower convex envelope over every nucleus, geometry, seed, or functional is not proved.'))
        print(model,sides[0]['delta_mu_MeV'],sides[1]['delta_mu_MeV'],slope,flush=True)
(BASE/'local_branch_ordering.json').write_text(json.dumps({'records':records},indent=2)+'\n')
