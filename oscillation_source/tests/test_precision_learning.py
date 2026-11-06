"""Selected-workspace checks. They never select the worked solution implicitly."""
from pathlib import Path
import os,subprocess
from tools.precision_study import ROOT, _fpm, run_one, load_pair

def selected_workspace():
    value=os.environ.get('NS_WORKSPACE','learner')
    if value not in {'learner','solutions','calculation'}:raise ValueError('Unknown workspace')
    return value

def test_selected_cell_thermodynamics():
    p=subprocess.run([_fpm(),'test','cell_checks','--profile','release'],
          cwd=ROOT/selected_workspace(),text=True,capture_output=True,timeout=120)
    assert p.returncode==0,p.stdout+'\n'+p.stderr
    assert '12 checks passed' in p.stdout

def test_selected_precision_pair(tmp_path):
    for alpha,name in [(0.,'base'),(.1,'variant')]:
        run_one(tmp_path/name,{'alpha':alpha,'rtol':1e-11,'mode_rtol':1e-11,'step_cap':.05},selected_workspace())
    q=load_pair(tmp_path/'base',tmp_path/'variant')
    assert q['identity_consistent']
    assert q['core_function_unchanged']
    # Deliberately no assertion on a preferred sign: this checks the pair, not its full numerical allowance.
    assert q['reference']['f_hz']>0 and q['variant']['p1_hz']>0
