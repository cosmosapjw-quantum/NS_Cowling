"""Read physical TF composition without inventing a free/cluster partition."""
from __future__ import annotations
import gzip,json
from pathlib import Path
import numpy as np
from .api import read_compose


def read_tf_composition(directory: str | Path) -> dict:
    """Return native thermo plus bound/unbound-inclusive TF fractions/profiles.

    Physical species conservation is checked from radial density integrals.
    CompOSE free-particle/group masks remain unavailable for nonuniform rows.
    The companion is bound to the table by density, energy, pressure and Ye.
    """
    root=Path(directory);metadata=json.loads((root/'composition.json').read_text())
    if metadata.get('composition_representation')!='spatial_TF_species_densities':
        raise ValueError('not a spatial TF composition representation')
    data=read_compose(root/'compose')
    states=json.loads(gzip.decompress((root/'selected_states.json.gz').read_bytes()))
    if len(states)!=data.n or len(metadata['rows'])!=data.n:
        raise ValueError('TF/table row count mismatch')
    totals=[];profiles=[]
    for i,(state,row) in enumerate(zip(states,metadata['rows'])):
        for source,actual in [('n_b',data.nb[i]),('epsilon',data.eps[i]),('pressure',data.p[i]),('Y_e',data.Ye[i])]:
            if not np.isclose(state[source],actual,rtol=3e-12,atol=2e-13):
                raise ValueError(f'TF companion/table identity mismatch in {source}')
        if not np.isclose(row['n_b'],state['n_b'],rtol=3e-12,atol=2e-13):
            raise ValueError('TF companion row density mismatch')
        p={k:np.asarray(v,float) for k,v in state['profile'].items()}
        size=state['points']
        if any(p[k].shape!=(size,) for k in ('x','weights','nn','np','ne')):
            raise ValueError('invalid TF radial profile shape')
        if any(not np.all(np.isfinite(p[k])) for k in ('x','weights','nn','np','ne')):
            raise ValueError('nonfinite TF profile')
        if any(np.any(p[k]<0) for k in ('nn','np','ne')):
            raise ValueError('negative TF particle density')
        if np.any(np.diff(p['x'])<=0) or p['x'][0]<=0 or p['x'][-1]>=1:
            raise ValueError('invalid TF radial coordinates')
        w=p['weights'];n=float(w@(p['nn']+p['np']));charge=float(w@(p['np']-p['ne']))
        if np.any(w<=0) or abs(sum(w)-1)>1e-12 or abs(n-state['n_b'])>2e-10 or abs(charge)>2e-10:
            raise ValueError('TF profile baryon/charge constraint mismatch')
        yn=float(w@p['nn'])/n;yp=float(w@p['np'])/n;ye=float(w@p['ne'])/n
        for name,value in [('Y_n_total',yn),('Y_p_total',yp),('Y_e',ye)]:
            if not np.isclose(row[name],value,rtol=2e-9,atol=2e-11):
                raise ValueError('TF companion fraction/profile mismatch')
        totals.append([yn,yp,ye,0.]);profiles.append(p)
    return {'thermodynamics':data,'representation':metadata['composition_representation'],
        'density_fm3':data.nb.copy(),'species':['n_total','p_total','e','mu'],
        'fractions':np.asarray(totals),'profiles':profiles,'metadata':metadata,
        'free_cluster_partition':'NOT_DEFINED_FOR_DIFFUSE_TF_PROFILE'}
