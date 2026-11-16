"""A finite-density HOLD must not become an interpolatable native interval."""
from pathlib import Path
import importlib.util


def test_unresolved_density_splits_native_exports(tmp_path, monkeypatch):
    path=Path(__file__).resolve().parents[2]/"scripts/build_gppva_tf.py"
    spec=importlib.util.spec_from_file_location("tf_gap_builder",path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    def solve(model,n,**kw):
        if n==.04 and (kw['dimension'],kw['seed'])==(2,'droplet'):
            raise RuntimeError('unresolved electron field')
        return dict(n_b=n,epsilon=n*940.,energy_per_baryon=940.,
            points=kw['points'],dimension=kw['dimension'],seed=kw['seed'],
            numerical_status='STATIONARY',muon_threshold_crossed=False,
            baryon_constraint_error=0.,charge_constraint_error=0.,
            scalar_pde_residual=0.,vector_pde_residual=0.,electron_pde_residual=0.)
    exported=[]
    def export(model,states,out):
        exported.append([s['n_b'] for s in states])
        return {'status':'EXPORTED_NATIVE_READBACK_VERIFIED','row_count':len(states)}
    monkeypatch.setattr(module,'solve_cell',solve)
    monkeypatch.setattr(module,'physical_geometry',lambda row:'droplet')
    monkeypatch.setattr(module,'export',export)
    result=module.run_model('GPPVA_DD2',[.01,.02,.03,.04,.05,.06,.07],[16],tmp_path)
    assert exported==[[.01,.02,.03],[.05,.06,.07]], 'HOLD density was bridged by one native table'
    assert result['rows']==6 and result['requested_rows']==7
