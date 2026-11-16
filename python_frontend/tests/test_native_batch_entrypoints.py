"""Public batch entry points carry native evidence and never run Python physics."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]

def load_batch():
    spec=importlib.util.spec_from_file_location('native_core_batch',ROOT/'scripts/generate_core_eos.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def test_native_batch_import_has_no_python_physics_dependency():
    # The package frontend may import array containers; forbid physical solvers.
    code="""import runpy,sys
runpy.run_path('scripts/generate_core_eos.py',run_name='import_only')
blocked=('compose_eos.core_regeneration','compose_eos.apr_uniform','compose_eos.apr_gibbs',
         'compose_eos.bbb_uniform','compose_eos.gm1_uniform','compose_eos.skyrme_uniform',
         'compose_eos.gdfm_uniform','compose_eos.spatial_tf','compose_eos.candidate_gdfm_cld')
assert not any(name in sys.modules for name in blocked), sorted(set(blocked)&set(sys.modules))
"""
    subprocess.run([sys.executable,'-c',code],cwd=ROOT,check=True)

def test_historical_physics_requires_explicit_oracle_request():
    for script in ('reoptimize_spg_surfaces.py','build_gppva_tf.py'):
        result=subprocess.run([sys.executable,str(ROOT/'scripts'/script),'unused'],
                              text=True,capture_output=True,cwd=ROOT)
        assert result.returncode!=0
        assert 'historical Python research oracle' in result.stderr
        assert 'compose-eos generate' in result.stderr

def test_native_batch_failure_preserves_dispatch(monkeypatch,tmp_path):
    batch=load_batch()
    calls=[]
    def fail(model,directory,**options):
        calls.append(model);directory.mkdir()
        (directory/'DISPATCH.json').write_text('{"status":"CHILD_FAILED","returncode":7}')
        return 7
    monkeypatch.setattr(batch,'run_generation',fail)
    assert batch.main([str(tmp_path/'batch'),'--models','GM_GM1','--points','3'])==7
    result=json.loads((tmp_path/'batch/SUMMARY.json').read_text())
    assert result['status']=='INCOMPLETE'
    assert calls==['GM_GM1']
    assert (tmp_path/'batch/GM_GM1/DISPATCH.json').read_text()=='{"status":"CHILD_FAILED","returncode":7}'

def test_native_scalar_receipt_rejects_nonfinite(tmp_path):
    batch=load_batch()
    path=tmp_path/'state.csv';path.write_text('n_b,epsilon\n0.1,nan\n')
    import pytest
    with pytest.raises(ValueError,match='nonfinite'):
        batch._native_rows(path)
