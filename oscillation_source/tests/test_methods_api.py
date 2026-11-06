"""The provided-method path must work independently of exercise completion."""
from pathlib import Path
import importlib, json, math, shutil
import pytest
from tools import precision_study as old
ROOT=Path(__file__).resolve().parents[1]

def api():
    assert (ROOT/'tools/stellar_methods.py').is_file(), 'No callable method API is provided'
    return importlib.import_module('tools.stellar_methods')

def test_numerical_library_is_available_without_solution_dependency():
    for name in ['methods','calculation']:
        p=ROOT/name/'fpm.toml';assert p.is_file(), f'{name} library is missing'
        assert 'solutions' not in p.read_text() and 'learner' not in p.read_text()
    assert (ROOT/'calculation/src/physics/ns_cowling.f90').is_file()

def test_nonfinite_difference_is_never_convergent():
    rows=[dict(delta_f_hz=v,delta_p1_hz=1.,identity_consistent=True,core_function_unchanged=True) for v in [.1,.101,float('nan')]]
    assert old.compare_sequence(rows)['status']=='invalid_numerical_data'

def test_unknown_method_and_nonfinite_controls_are_rejected(tmp_path):
    m=api()
    with pytest.raises(ValueError):m.solve_pair(tmp_path/'a',level='magic')
    with pytest.raises(ValueError):m.solve_star(tmp_path/'b',alpha=float('nan'))
    assert not (tmp_path/'a').exists() and not (tmp_path/'b').exists()

def test_repeated_settings_are_not_a_refinement():
    m=api();r={'f_hz':2.,'p1_hz':4.,'rtol':1e-9,'mode_rtol':1e-9,'step_cap':.1,'f_bracket_hz':1e-12,'p1_bracket_hz':1e-12}
    rows=[{'reference':r.copy(),'variant':r.copy()} for i in range(4)]
    assert m.absolute_convergence(rows)['status']=='invalid_refinement_schedule'

def test_cancellation_cannot_hide_divergent_absolute_frequencies():
    m=api();rows=[]
    for i,v in enumerate([100.,101.,104.,116.]):
        base={'f_hz':v,'p1_hz':2*v,'rtol':10**(-7-2*i),'mode_rtol':10**(-7-2*i),'step_cap':.2/2**i,'f_bracket_hz':1e-11,'p1_bracket_hz':1e-11}
        rows.append({'reference':base,'variant':{**base,'f_hz':v-.01,'p1_hz':2*v-.4}})
    assert m.absolute_convergence(rows)['status']=='absolute_convergence_unresolved'

def test_plan_can_be_created_without_scientific_runs(tmp_path):
    m=api();r=m.run_study(tmp_path/'study',max_new_runs=0)
    assert r['status']=='pending' and len(r['pending_runs'])==21
    assert r['physical_claim_status']=='pending_additional_checks'
    assert not list((tmp_path/'study').glob('*/result.json'))

@pytest.fixture(scope='module')
def computed_pair(tmp_path_factory):
    import os
    known=os.environ.get('NS_VALIDATED_STUDY')
    if known:
        p=Path(known)
        return p/'L2_a0.0',p/'L2_a0.1'
    p=tmp_path_factory.mktemp('callable-method-pair')
    api().solve_pair(p,level='fine')
    return p/'reference',p/'variant'

def test_ready_method_produces_two_identified_modes(computed_pair):
    a,b=computed_pair;q=old.load_pair(a,b)
    assert q['identity_consistent'] and q['core_function_unchanged']
    assert q['delta_f_hz']==pytest.approx(-.01044847414,abs=2e-7)
    assert q['delta_p1_hz']==pytest.approx(-.46974177113,abs=2e-6)

def test_output_profile_tampering_is_detected(computed_pair,tmp_path):
    a,b=computed_pair;dest=tmp_path/'tampered';shutil.copytree(a,dest)
    p=dest/'result_background.csv';p.write_text(p.read_text()+'\n')
    with pytest.raises(old.StudyError,match='Output hash mismatch'):old._verified_result(dest)

def test_input_tampering_is_detected(computed_pair,tmp_path):
    a,b=computed_pair;dest=tmp_path/'tampered';shutil.copytree(a,dest)
    p=dest/'input.nml';p.write_text(p.read_text().replace('alpha=0.0','alpha=0.1'))
    with pytest.raises(old.StudyError,match='Input hash mismatch'):old._verified_result(dest)

def test_nan_result_rejected_even_when_file_checksum_is_updated(computed_pair,tmp_path):
    a,b=computed_pair;dest=tmp_path/'nan';shutil.copytree(a,dest)
    p=dest/'result.json';r=json.loads(p.read_text());r['p1_residual']=float('nan');p.write_text(json.dumps(r))
    rec=dest/'execution.json';e=json.loads(rec.read_text());e['outputs']['result.json']=old.sha(p);rec.write_text(json.dumps(e))
    with pytest.raises(old.StudyError,match='Nonfinite'):old._verified_result(dest)

def test_original_failed_sequence_is_still_rejected():
    rows=[dict(delta_f_hz=a,delta_p1_hz=b,identity_consistent=True,core_function_unchanged=True)
          for a,b in zip([-.0318299,-.559565,-.144152],[-.674666,.596514,3.531938])]
    q=old.compare_sequence(rows)
    assert q['modes']['f']['verdict']=='below_resolution'
    assert q['modes']['p1']['verdict']=='convergence_unresolved'

def test_library_sources_do_not_call_solution_tree():
    for base in ['methods','calculation']:
        for p in (ROOT/base).rglob('*.f90'):
            assert 'solutions/' not in p.read_text().lower()

def test_new_study_requires_all_evidence(tmp_path):
    m=api();p=tmp_path/'study';r=m.run_study(p,max_new_runs=1)
    assert len(r['completed_runs'])==1 and len(r['pending_runs'])==20
    assert r['physical_claim_status']=='pending_additional_checks'
    record=p/'L0_a0.0/execution.json';stamp=record.stat().st_mtime_ns
    q=m.run_study(p,max_new_runs=0)
    assert q['new_runs_this_call']==0 and record.stat().st_mtime_ns==stamp

def test_interrupted_attempt_requires_explicit_retry_and_is_preserved(tmp_path):
    m=api();p=tmp_path/'study';m.run_study(p,max_new_runs=0)
    dest=p/'L0_a0.0';dest.mkdir();(dest/'execution.json').write_text('{"status":"partial"}')
    (dest/'partial.txt').write_text('original interrupted evidence')
    with pytest.raises(old.StudyError,match='Incomplete'):m.run_study(p,max_new_runs=1)
    r=m.run_study(p,max_new_runs=1,retry_failed=True)
    retained=list((p/'failed_attempts').glob('*/partial.txt'))
    assert len(retained)==1 and retained[0].read_text()=='original interrupted evidence'
    assert len(r['completed_runs'])==1

def test_unauthorized_folders_are_not_implicit_workspaces():
    with pytest.raises(old.StudyError):old.source_digest('../solutions')

def test_partial_supplementary_run_does_not_block_pending_report(computed_pair,tmp_path):
    """Duplicate a valid pair as a report fixture, not as a convergence experiment."""
    m=api();out=tmp_path/'report';out.mkdir()
    recipe={'workspace':'calculation','groups':{g:[] for g in m.REQUIRED_GROUPS}}
    for i in range(4):
        for alpha,source in zip((0.0,0.1),computed_pair):
            name=f'L{i}_a{alpha}';shutil.copytree(source,out/name)
            rec=json.loads((source/'execution.json').read_text())
            recipe['groups']['main'].append((name,rec['parameters']))
    for alpha in (0.0,0.1):
        name=f'rk4_a{alpha}';p=out/name;p.mkdir()
        (p/'execution.json').write_text('{"status":"running"}')
        (p/'original_partial.txt').write_text('preserve interrupted computation')
        recipe['groups']['method'].append((name,{}))
    result=m._report(out,recipe)
    assert result['status']=='pending'
    assert len(result['pending_runs'])==2
    assert result['physical_claim_status']=='pending_additional_checks'
    assert (out/'rk4_a0.1/original_partial.txt').read_text()=='preserve interrupted computation'
