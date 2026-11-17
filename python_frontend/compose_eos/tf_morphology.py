"""Descriptive radial-density morphology, separate from physical admissibility.

A half-contrast occupancy is a declared shape diagnostic, not a sharp phase
interface or a criterion for rejecting a stationary variational candidate.
"""
from __future__ import annotations
import numpy as np

def _occupancy(values):
    values=np.asarray(values,float)
    contrast=float(np.ptp(values))
    if contrast<=0:return {'contrast':contrast,'intervals':[],'normalized_endpoints':None,
                            'interior_dense_interval':False,'multiple_dense_intervals':False}
    unit=(values-values.min())/contrast;mask=unit>.5
    starts=np.flatnonzero(mask&~np.r_[False,mask[:-1]])
    ends=np.flatnonzero(mask&~np.r_[mask[1:],False])
    return {'contrast':contrast,'normalized_endpoints':[float(unit[0]),float(unit[-1])],
            'intervals':[[int(a),int(b)] for a,b in zip(starts,ends)],
            'interior_dense_interval':bool(len(starts)==1 and not mask[0] and not mask[-1]),
            'multiple_dense_intervals':bool(len(starts)>1)}

def profile_morphology(row):
    """Identify strong shell-like cases using BOTH total baryons and protons."""
    p=row['profile'];protons=np.asarray(p['np'],float);neutrons=np.asarray(p['nn'],float)
    if (protons.ndim!=1 or len(protons)<2 or neutrons.shape!=protons.shape or
        not np.isfinite(protons).all() or not np.isfinite(neutrons).all()):
        raise ValueError('finite matching radial nucleon profiles required')
    baryons=protons+neutrons
    weights=np.asarray(p.get('weights',np.ones(len(baryons))),float)
    if weights.shape!=baryons.shape or not np.isfinite(weights).all() or np.any(weights<=0):
        raise ValueError('finite positive matching radial weights required')
    weights=weights/weights.sum()
    scale=float(weights@baryons)
    if scale<=0:raise ValueError('positive mean baryon density required')
    species={'nn':neutrons,'np':protons}
    if 'ne' in p:species['ne']=np.asarray(p['ne'],float)
    differences={}
    for name,values in species.items():
        if values.shape!=baryons.shape or not np.isfinite(values).all():
            raise ValueError('finite matching species profiles required')
        absolute=float(np.ptp(values));mean=float(weights@values)
        differences[name]={'absolute_fm3':absolute,
            'relative':absolute/max(abs(mean),1e-30),
            'uniform':bool(absolute<=max(1e-18,1e-12*scale) or
                           absolute<=1e-4*abs(mean))}
    if float(np.ptp(baryons))/scale<1e-4:
        uniform=all(item['uniform'] for item in differences.values())
        return {'classification':'UNIFORM' if uniform else 'NEAR_UNIFORM_COMPOSITION_PROFILE',
            'species_contrast':differences,
            'threshold':'Each species relative contrast <=1e-4 or absolute <=max(1e-18,1e-12*nB) fm^-3; total-baryon contrast <1e-4',
            'scope':'A descriptive numerical shape; near-uniform modulation is not evidence for a resolved new physical phase.'}
    b=_occupancy(baryons);q=_occupancy(protons);dimension=row['dimension']
    if dimension==1:
        label='PLANAR_REPEATED_OR_TRANSLATED' if b['multiple_dense_intervals'] or b['interior_dense_interval'] else 'PLANAR_SIMPLE'
    elif b['interior_dense_interval'] and q['interior_dense_interval']:
        label='NONCANONICAL_CYLINDRICAL_PROFILE' if dimension==2 else 'NONCANONICAL_SPHERICAL_PROFILE'
    elif b['multiple_dense_intervals'] and q['multiple_dense_intervals']:
        label='MULTIPLE_RADIAL_DENSE_INTERVALS'
    elif b['interior_dense_interval']!=q['interior_dense_interval'] or b['multiple_dense_intervals']!=q['multiple_dense_intervals']:
        label='BARYON_PROTON_SHAPE_DIFFERENCE'
    else:label='SIMPLE_RADIAL_PROFILE'
    return {'classification':label,'baryon_halfcontrast':b,'proton_halfcontrast':q,
        'scope':'Descriptive half-contrast occupancy only; retain all candidate energies and profiles.'}


# Local output codes describe sampled profiles, not source-provider phase IDs.
TF_PHASE_CODES={'uniform':3,'droplet':21,'rod':22,'slab':23,'tube':24,'bubble':25,
    'noncanonical_cylindrical_profile':26,'noncanonical_spherical_profile':27,
    'near_uniform_profile':28,'noncanonical_radial_profile':29}

def geometry_label(row):
    """Conservative export label; ambiguous profiles retain sparse composition."""
    result=profile_morphology(row);kind=result['classification'];dimension=int(row['dimension'])
    if kind=='UNIFORM':return 'uniform'
    if kind=='NEAR_UNIFORM_COMPOSITION_PROFILE':return 'near_uniform_profile'
    if dimension==1:return 'slab'
    if kind in ('NONCANONICAL_CYLINDRICAL_PROFILE','NONCANONICAL_SPHERICAL_PROFILE'):
        return kind.lower()
    if kind in ('MULTIPLE_RADIAL_DENSE_INTERVALS','BARYON_PROTON_SHAPE_DIFFERENCE'):
        return 'noncanonical_radial_profile'
    p=row['profile'];dense=(p['nn'][0]+p['np'][0])>(p['nn'][-1]+p['np'][-1])
    if dimension not in (2,3):raise ValueError('radial dimension must be1,2,3')
    return {2:'rod' if dense else 'tube',3:'droplet' if dense else 'bubble'}[dimension]
