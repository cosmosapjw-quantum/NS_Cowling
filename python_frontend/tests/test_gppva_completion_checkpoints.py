"""Scheduler checkpoint contracts; synthetic rows are not physics evidence."""
import gzip,json,sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import complete_gppva_geometry as runner

def fake(model,n,points,family,*args):
    d,seed=runner.FAMILIES[family]
    return [dict(model_id=model,n_b=n,points=points,dimension=d,seed=seed,requested_family=family,
        numerical_status='UNIFORM_STATIONARY',muon_threshold_crossed=False,epsilon=n*950.,energy_per_baryon=950.,
        baryon_density_contrast=0.,baryon_constraint_error=0.,charge_constraint_error=0.,scalar_pde_residual=0.,vector_pde_residual=0.,electron_pde_residual=0.,
        profile={'nn':[.048]*points,'np':[.002]*points,'ne':[.002]*points,'weights':[1/points]*points,'x':[(i+.5)/points for i in range(points)]})]

def test_complete_node_checkpoint_has_all_families_and_identity(monkeypatch,tmp_path):
    monkeypatch.setattr(runner,'solve_attempt',fake)
    runner.solve_node('TEST',.05,8,{g:[] for g in runner.FAMILIES},None,[],tmp_path,{'synthetic':'unit_contract'})
    path=next(tmp_path.rglob('*.json.gz'));saved=json.loads(gzip.decompress(path.read_bytes()))
    assert len(saved['attempts'])==6
    assert saved['status']=='LOWEST_EVALUATED_STATIONARY_CANDIDATE'
    assert saved['source_sha256']=={'synthetic':'unit_contract'}
    assert not list(tmp_path.rglob('*.tmp.gz'))

def test_completed_family_survives_next_attempt_crash(monkeypatch,tmp_path):
    def crash_later(model,n,points,family,*args):
        if family=='rod':raise RuntimeError('simulated process-side interruption')
        return fake(model,n,points,family,*args)
    monkeypatch.setattr(runner,'solve_attempt',crash_later)
    with pytest.raises(RuntimeError,match='interruption'):
        runner.solve_node('TEST',.05,8,{g:[] for g in runner.FAMILIES},None,[],tmp_path,{'synthetic':'unit_contract'})
    path=next(tmp_path.rglob('*.json.gz'));saved=json.loads(gzip.decompress(path.read_bytes()))
    assert saved['status']=='RUNNING'
    assert saved['completed_families']==['droplet']
    assert len(saved['attempts'])==1
