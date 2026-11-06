"""Numerical evidence handling and direct-cell solver checks, not learner grading."""
from pathlib import Path
import json,os
import pytest
from tools.precision_study import run_one, load_pair, compare_sequence, StudyError, ROOT

@pytest.fixture(scope='module')
def pair(tmp_path_factory):
    out=tmp_path_factory.mktemp('direct-cell')
    a=run_one(out/'a0',{'alpha':0.0,'rtol':1e-10,'mode_rtol':1e-10,'step_cap':.05})
    b=run_one(out/'a1',{'alpha':.1,'rtol':1e-10,'mode_rtol':1e-10,'step_cap':.05})
    return out,a,b

def test_direct_pair_closes_boundary_and_identity(pair):
    d,a,b=pair
    q=load_pair(d/'a0',d/'a1')
    assert q['identity_consistent']
    assert q['delta_f_hz']==pytest.approx(-.01044847414,abs=2e-7)
    assert q['delta_p1_hz']==pytest.approx(-.46974177113,abs=2e-6)
    for r in [a,b]:
        assert abs(r['f_residual'])<1e-8 and abs(r['p1_residual'])<1e-8
        assert abs(r['f_surface_residual'])<1e-14
        assert abs(r['p1_surface_residual'])<1e-14
        assert abs(r['mass_kg']/r['target_mass_kg']-1)<1e-12

def test_identical_run_is_reused_without_rerun(pair):
    d,a,b=pair
    p=d/'a0/execution.json';before=p.stat().st_mtime_ns
    r=run_one(d/'a0',{'alpha':0.0,'rtol':1e-10,'mode_rtol':1e-10,'step_cap':.05})
    assert p.stat().st_mtime_ns==before
    assert r['f_hz']==a['f_hz']

def test_changed_input_cannot_reuse_output(pair):
    d,a,b=pair
    with pytest.raises(StudyError):run_one(d/'a0',{'alpha':.2,'rtol':1e-10,'mode_rtol':1e-10,'step_cap':.05})

def test_failure_never_becomes_complete(tmp_path):
    with pytest.raises(StudyError):run_one(tmp_path/'bad',{'p_cut_pa':1e5})
    rec=json.loads((tmp_path/'bad/execution.json').read_text())
    assert rec['status']=='failed' and rec['returncode']!=0

def test_pair_requires_same_controls(pair,tmp_path):
    import shutil
    d,a,b=pair;shutil.copytree(d/'a1',tmp_path/'altered')
    p=tmp_path/'altered/result.json';r=json.loads(p.read_text());r['target_mass_kg']*=1.001;p.write_text(json.dumps(r))
    with pytest.raises(StudyError):load_pair(d/'a0',tmp_path/'altered')

def test_false_identity_cannot_resolve():
    rows=[dict(delta_f_hz=x,delta_p1_hz=2*x,identity_consistent=False,core_function_unchanged=True) for x in [.1,.101,.1011]]
    assert compare_sequence(rows)['status']=='identity_unresolved'

def test_diverging_changes_cannot_resolve():
    rows=[dict(delta_f_hz=x,delta_p1_hz=2*x,identity_consistent=True,core_function_unchanged=True) for x in [.1,.2,1.]]
    assert compare_sequence(rows)['status']=='convergence_unresolved'

def test_exact_null_has_absolute_error_handling():
    rows=[dict(delta_f_hz=0.,delta_p1_hz=0.,identity_consistent=True,core_function_unchanged=True)]*3
    q=compare_sequence(rows)
    assert q['modes']['f']['verdict']=='below_resolution'
    assert q['modes']['p1']['verdict']=='below_resolution'

def test_matching_background_must_agree(pair,tmp_path):
    import shutil
    from tools.precision_study import sha
    d,a,b=pair; dest=tmp_path/'wrong-background';shutil.copytree(d/'a1',dest)
    p=dest/'result.json';r=json.loads(p.read_text())
    r['f_match_radius_mismatch_m']=100.0
    p.write_text(json.dumps(r))
    rec=dest/'execution.json';e=json.loads(rec.read_text())
    e['outputs']['result.json']=sha(p);rec.write_text(json.dumps(e))
    with pytest.raises(StudyError,match='background'):load_pair(d/'a0',dest)

def test_unfinished_learner_has_no_solution_fallback(tmp_path):
    with pytest.raises(StudyError):run_one(tmp_path/'learner',{'alpha':0.0},workspace='learner')
    e=json.loads((tmp_path/'learner/execution.json').read_text())
    assert e['workspace']=='learner' and e['status']=='failed'
    assert 'not_implemented' in ((tmp_path/'learner/stdout.log').read_text()+(tmp_path/'learner/stderr.log').read_text())

def test_manifest_does_not_admit_arbitrary_targets(tmp_path):
    from tools.runtime import copy_sources,validate_packages,DependencyViolation
    copy_sources(ROOT,tmp_path/'source','solutions')
    p=tmp_path/'source/solutions/fpm.toml'
    p.write_text(p.read_text().replace('main = "cell_solve.f90"','main = "other.f90"'))
    with pytest.raises(DependencyViolation):validate_packages(tmp_path/'source','solutions')

def test_changed_mass_across_levels_is_not_refinement():
    rows=[]
    for mass in [1.,1.,1.1]:
        rows.append(dict(delta_f_hz=.1,delta_p1_hz=.2,identity_consistent=True,
            core_function_unchanged=True,reference={'target_mass_kg':mass},variant={'target_mass_kg':mass}))
    assert compare_sequence(rows)['status']=='controls_failed'
