"""Conservative initial profiles for finite screened-TF candidate searches."""
from __future__ import annotations
import copy
import math
import numpy as np

def weak_cluster_seed(uniform,dimension,amplitude,radius):
    """Return a clustered initial guess with unchanged mean species densities.

    Input is a stationary npe uniform-row representation on the solver's
    uniform radial grid. Positive amplitude gives a dense centre; negative
    amplitude gives a hollow centre. All densities and radius remain free in
    the subsequent solve. This constructs an initial guess, not a new EOS.
    """
    if dimension not in (1,2,3):raise ValueError('dimension must be1,2,3')
    if not math.isfinite(amplitude) or not 0<abs(amplitude)<1:
        raise ValueError('signed amplitude must have magnitude between0 and1')
    if not math.isfinite(radius) or radius<=0:raise ValueError('positive finite radius required')
    n=float(uniform['n_b']);yn=float(uniform['Y_n']);yp=float(uniform['Y_p'])
    if not all(math.isfinite(v) for v in [n,yn,yp]) or n<=0 or min(yn,yp)<0 or abs(yn+yp-1)>1e-8:
        raise ValueError('finite normalized neutron/proton means required')
    row=copy.deepcopy(uniform);row['dimension']=dimension;row['radius_fm']=float(radius)
    row['numerical_status']='INITIAL_GUESS_ONLY';row['optimizer_success']=False
    p=row['profile'];x=np.asarray(p['x'],float)
    if x.ndim!=1 or len(x)<2 or not np.allclose(x,(np.arange(len(x))+.5)/len(x),rtol=0,atol=1e-13):
        raise ValueError('uniform radial midpoint grid required')
    w=np.diff(np.linspace(0,1,len(x)+1)**dimension)
    shape=np.exp(-(x/.45)**2);shape-=w@shape
    p['weights']=w.tolist()
    p['np']=(yp*n*(1+amplitude*shape)).tolist()
    p['nn']=(yn*n*(1+.3*amplitude*shape)).tolist()
    p['ne']=np.full(len(x),yp*n).tolist()
    return row
