"""Initial-family recovery preserves the actual density/charge contract."""
import copy
import numpy as np
import pytest
from compose_eos.tf_search import weak_cluster_seed

def uniform_row():
    n=.095;yp=.02;points=16
    return {'n_b':n,'Y_p':yp,'Y_n':1-yp,'points':points,'dimension':3,'radius_fm':30.,
         'profile':{'x':((np.arange(points)+.5)/points).tolist(),
                    'nn':[n*(1-yp)]*points,'np':[n*yp]*points,'ne':[n*yp]*points}}

@pytest.mark.parametrize('dimension',[1,2,3])
@pytest.mark.parametrize('amplitude',[.03,.3,-.03,-.3])
def test_seed_preserves_means_and_declared_orientation(dimension,amplitude):
    old=uniform_row();saved=copy.deepcopy(old);row=weak_cluster_seed(old,dimension,amplitude,15.)
    p=row['profile'];w=np.array(p['weights']);nn=np.array(p['nn']);np_=np.array(p['np']);ne=np.array(p['ne'])
    assert abs(w@(nn+np_)-old['n_b'])<3e-17
    assert abs(w@(np_-ne))<1e-18
    assert (nn[0]+np_[0]-nn[-1]-np_[-1])*amplitude>0
    assert min(nn)>=0 and min(np_)>=0
    assert row['radius_fm']==15. and row['dimension']==dimension
    assert row['numerical_status']=='INITIAL_GUESS_ONLY' and not row['optimizer_success']
    assert old==saved

@pytest.mark.parametrize('dimension,amplitude,radius',[(0,.3,15.),(2,1.,15.),(2,float('nan'),15.),(2,.3,0.)])
def test_invalid_seed_parameters_fail_closed(dimension,amplitude,radius):
    with pytest.raises(ValueError):weak_cluster_seed(uniform_row(),dimension,amplitude,radius)
