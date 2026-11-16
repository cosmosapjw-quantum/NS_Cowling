"""Explicit surface prescriptions must produce independently re-equilibrated cells."""
from dataclasses import FrozenInstanceError
import hashlib,json
import numpy as np
import pytest
from compose_eos import candidate_gdfm_cld as cld
from compose_eos.spg_crust_generator import neutron_drip,sequence

def test_surface_parameters_are_immutable_and_model_specific():
    p=cld.SurfaceParameters('SPG_M1',1.03,12.,'explicit_test','a'*64)
    with pytest.raises(FrozenInstanceError):p.b=20.
    with pytest.raises(ValueError,match='model'):
        cld.reduced_cell('SPG_M2',.02,[.5,.15,.1],surface=p)

def test_default_surface_cache_tracks_content_identity(tmp_path,monkeypatch):
    path=tmp_path/'data/research_generator_20261116/candidate_surface/candidate_surface_fit.json'
    path.parent.mkdir(parents=True)
    data={'models':{'SPG_M1':{'all_mass':{'sigma0_MeV_fm2':1.1,'b':13.}}}}
    path.write_text(json.dumps(data));monkeypatch.setattr(cld,'project_root',lambda:tmp_path)
    first=cld.default_surface_parameters('SPG_M1')
    data['models']['SPG_M1']['all_mass']['b']=17.;path.write_text(json.dumps(data))
    second=cld.default_surface_parameters('SPG_M1')
    assert first.b==13. and second.b==17.
    assert first.source_sha256!=second.source_sha256
    assert second.source_sha256==hashlib.sha256(path.read_bytes()).hexdigest()

def test_explicit_default_preserves_functional_and_alternative_gradient():
    p=cld.default_surface_parameters('SPG_M3')
    z=np.array([.3,.14,.1,.02])
    e,g,_=cld.reduced_two_phase_cell('SPG_M3',.03,z)
    ee,gg,_=cld.reduced_two_phase_cell('SPG_M3',.03,z,surface=p)
    assert e==ee;np.testing.assert_array_equal(g,gg)
    alt=cld.SurfaceParameters('SPG_M3',1.08,10.,'explicit_test','b'*64)
    _,grad,_=cld.reduced_two_phase_cell('SPG_M3',.03,z,surface=alt)
    numeric=[]
    for j,h in enumerate([1e-5,1e-6,1e-6,1e-5]):
        a=z.copy();b=z.copy();a[j]+=h;b[j]-=h
        numeric.append((cld.reduced_two_phase_cell('SPG_M3',.03,a,surface=alt)[0]-cld.reduced_two_phase_cell('SPG_M3',.03,b,surface=alt)[0])/(2*h))
    np.testing.assert_allclose(grad,numeric,rtol=0,atol=2e-6)

def test_surface_change_reoptimizes_composition_and_first_law():
    p=cld.SurfaceParameters('SPG_M2',1.04,14.,'explicit_test','c'*64)
    q=cld.SurfaceParameters('SPG_M2',1.06,7.,'alternative_test','d'*64)
    n=.02;h=2e-6
    a,_=cld.solve_two_phase('SPG_M2',n,surface=p)
    b,_=cld.solve_two_phase('SPG_M2',n,surface=q)
    assert a['status']==b['status']=='STATIONARY_CANDIDATE'
    assert abs(a['composition']['Ye']-b['composition']['Ye'])>1e-5
    assert a['surface_parameters']['source_sha256']=='c'*64
    lo,_=cld.solve_two_phase('SPG_M2',n-h,surface=q)
    hi,_=cld.solve_two_phase('SPG_M2',n+h,surface=q)
    np.testing.assert_allclose((hi['energy']-lo['energy'])/(2*h),b['mu_B_MeV'],rtol=0,atol=2e-5)
    assert neutron_drip('SPG_M2',surface=p)['surface_parameters']['source_sha256']=='c'*64

@pytest.mark.parametrize('sigma,b',[(0,5),(float('nan'),5),(1,-16),(1,float('inf'))])
def test_invalid_explicit_surface_rejected(sigma,b):
    with pytest.raises(ValueError):cld.SurfaceParameters('SPG_M1',sigma,b,'test','e'*64)


def test_native_export_preserves_physical_support_findings(monkeypatch,tmp_path):
    # A source-preserving table must survive a declared acausal segment, while
    # conservation errors remain implementation failures.
    import importlib
    from pathlib import Path
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[2]/'scripts'))
    runner=importlib.import_module('reoptimize_spg_surfaces')
    states=sequence('SPG_M1',[.01,.02,.03])
    data=runner.to_eos_data(states)
    monkeypatch.setattr(runner,'write_compose',lambda *args:None)
    monkeypatch.setattr(runner,'read_compose',lambda *args:data)
    checks={f'C{i}':{'pass_':i<=5} for i in range(1,12)}
    monkeypatch.setattr(runner,'check_consistency',lambda *args,**kwargs:checks)
    monkeypatch.setattr(runner,'physical_domain_report',lambda *args:{'segments':'retained'})
    result=runner.export_native(tmp_path,states)
    assert not result['native_checks']['C10']['pass_']
    assert result['physical_support']=={'segments':'retained'}
    checks['C3']['pass_']=False
    with pytest.raises(RuntimeError,match='constraint C3'):
        runner.export_native(tmp_path,states)
