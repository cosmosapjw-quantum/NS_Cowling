"""The public runner retains a failed cold start and executes weak-profile recovery."""
from pathlib import Path
import importlib.util
import gzip,json
import numpy as np


def test_public_tf_runner_recovers_failed_family_and_keeps_failure(tmp_path,monkeypatch):
    path=Path(__file__).resolve().parents[2]/'scripts/build_gppva_tf.py'
    spec=importlib.util.spec_from_file_location('tf_recovery_builder',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    def solve(model,n,**kw):
        if (kw['dimension'],kw['seed'])==(2,'droplet') and kw['initial'] is None:
            raise RuntimeError('vector field line search on cold trial')
        points=kw['points'];w=np.diff(np.linspace(0,1,points+1)**kw['dimension'])
        return dict(model_id=model,n_b=n,epsilon=n*940.,energy_per_baryon=940.,
            points=points,dimension=kw['dimension'],seed=kw['seed'],
            numerical_status='UNIFORM_STATIONARY',muon_threshold_crossed=False,
            baryon_constraint_error=0.,charge_constraint_error=0.,
            scalar_pde_residual=0.,vector_pde_residual=0.,electron_pde_residual=0.,
            baryon_density_contrast=0.,Y_n=.98,Y_p=.02,radius_fm=30.,
            profile={'x':((np.arange(points)+.5)/points).tolist(),'weights':w.tolist(),
                     'nn':[.98*n]*points,'np':[.02*n]*points,'ne':[.02*n]*points})
    monkeypatch.setattr(module,'solve_cell',solve)
    monkeypatch.setattr(module,'export',lambda model,states,out:
                        {'status':'EXPORTED_NATIVE_READBACK_VERIFIED','row_count':len(states)})
    result=module.run_model('GPPVA_DD2',[.0625],[8],tmp_path)
    assert result['rows']==1,'cold failure was not followed by the available weak-profile search'
    raw=json.loads(gzip.decompress((tmp_path/'GPPVA_DD2/attempts_N008.json.gz').read_bytes()))
    attempts=raw[0]['attempts']
    assert any(a.get('exception')=='RuntimeError' for a in attempts)
    assert any(a.get('initial_kind')=='WEAK_CLUSTER_RECOVERY' for a in attempts)
