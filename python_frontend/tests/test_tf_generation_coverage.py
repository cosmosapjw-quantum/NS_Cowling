"""Regression: a failed geometry must not disappear from a successful export."""
from pathlib import Path
import importlib.util

def test_main_tf_generator_retains_unresolved_geometry(tmp_path, monkeypatch):
    path=Path(__file__).resolve().parents[2]/"scripts/build_gppva_tf.py"
    spec=importlib.util.spec_from_file_location("tf_coverage_builder",path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    def solve(model,n,**kw):
        if (kw["dimension"],kw["seed"])==(2,"droplet"):
            raise RuntimeError("screened electron line search")
        return dict(n_b=n,epsilon=n*940.,energy_per_baryon=940.,
                    points=kw["points"],dimension=kw["dimension"],seed=kw["seed"],
                    numerical_status="STATIONARY",muon_threshold_crossed=False,
                    baryon_constraint_error=0.,charge_constraint_error=0.,
                    scalar_pde_residual=0.,vector_pde_residual=0.,electron_pde_residual=0.)
    exported=[]
    monkeypatch.setattr(module,"solve_cell",solve)
    monkeypatch.setattr(module,"physical_geometry",lambda row:"droplet")
    monkeypatch.setattr(module,"export",lambda model,states,out:
                        exported.extend(states) or {"status":"TEST_EXPORT","row_count":len(states)})
    result=module.run_model("GPPVA_DD2",[.0625],[16],tmp_path)
    assert exported==[], "an unresolved rod family was silently omitted"
    assert result["rows"]==0
