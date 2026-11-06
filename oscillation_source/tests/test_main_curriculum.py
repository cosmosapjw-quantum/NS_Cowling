"""The callable path is part of the course, not a deferred implementation task."""
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]

def test_eos_and_boundary_value_methods_are_in_their_main_chapters():
    eos=(ROOT/'notes/chapters/ch09.md').read_text()
    bvp=(ROOT/'notes/chapters/ch05.md').read_text()
    assert 'table_pressure_coordinate' in eos and 'cell_state' in eos
    assert 'det(' in bvp and 'ns_stellar_methods' in bvp
    plan=(ROOT/'curriculum/course_plan.md').read_text()
    assert '# 선택 심화' not in plan
    assert 'run_study' in plan and '16회' in plan

def test_provided_assignments_select_calculation_by_default():
    for name in ('A18','A20','A21'):
        p=subprocess.run([sys.executable,'tools/run_assignment.py',name,'--plan'],cwd=ROOT,
                         capture_output=True,text=True,check=True)
        r=json.loads(p.stdout)
        assert r['workspace']=='calculation'
        assert r['missing_files']==[]

def test_physics_exercise_is_not_silently_replaced():
    p=subprocess.run([sys.executable,'tools/run_assignment.py','A09','--plan'],cwd=ROOT,
                     capture_output=True,text=True,check=True)
    r=json.loads(p.stdout)
    assert r['workspace']=='learner'
    assert any('learner/src/modes/ns_cowling.f90'==s for s in r['source_paths'])

def test_eos_task_precedes_controlled_comparison_and_all_paths_exist():
    rows=json.loads((ROOT/'curriculum/assignments.json').read_text())['assignments']
    ids=[r['id'] for r in rows]
    assert ids.index('A20')<ids.index('A19')<ids.index('A21')<ids.index('P0')
    assert next(r for r in rows if r['id']=='A20')['stage']=='ch09'
    for r in rows:
        assert (ROOT/r['cue_source']).is_file()
        for p in r['source_paths']:assert (ROOT/p).is_file(),(r['id'],p)
