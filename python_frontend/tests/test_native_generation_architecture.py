"""Native production architecture tests; legacy Python physics is an oracle only."""
import json
import subprocess
import sys
from pathlib import Path
from compose_eos import generation_dispatch as dispatch

ROOT=Path(__file__).resolve().parents[2]

def test_public_generation_command_is_native_for_all_families(tmp_path):
    requests=[
      dispatch.generation_request('GM_GM1',points=3),
      dispatch.generation_request('SPG_M1',points=17),
      dispatch.generation_request('GPPVA_DD2',densities=[.03,.04,.05],resolutions=[16,24])]
    for request in requests:
        command=dispatch._command(ROOT,tmp_path/'generated',request)
        assert command[0] != sys.executable
        assert Path(command[0]).name == 'run-generator.sh'
        assert Path(command[1]).suffix == '.nml'
        assert command[2] == str(tmp_path/'generated')

def test_catalog_does_not_load_python_physical_solvers():
    code="from compose_eos.generation_dispatch import generator_catalog; import sys,json; c=generator_catalog(); print(json.dumps({'n':len(c),'loaded':[m for m in sys.modules if m.startswith('compose_eos.')]}))"
    p=subprocess.run([sys.executable,'-c',code],capture_output=True,text=True,check=True)
    result=json.loads(p.stdout)
    assert result['n']==18
    forbidden={'compose_eos.core_regeneration','compose_eos.apr_uniform','compose_eos.apr_gibbs','compose_eos.gdfm_uniform','compose_eos.spatial_tf','compose_eos.candidate_gdfm_cld'}
    assert not forbidden.intersection(result['loaded'])
