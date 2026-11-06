"""Table evaluation uses supplied methods, without learner implementations."""
import importlib
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]

def api():
    assert (ROOT/'tools/inspect_table.py').is_file(), 'No callable provided-table inspection exists'
    return importlib.import_module('tools.inspect_table')

def test_provided_table_accepts_consistent_input(tmp_path):
    r=api().inspect_table(ROOT/'fixtures/tables/hermite_consistent.csv',tmp_path/'table')
    assert r['status']=='ok'
    assert r['rows']==3 and r['admitted_intervals']==2 and r['rejected_intervals']==0
    assert r['model_ready'] is True

def test_provided_table_rejects_unadmitted_real_input(tmp_path):
    r=api().inspect_table(ROOT/'data/normalized/compose_134_si.csv',tmp_path/'table')
    assert r['rows']==197 and r['admitted_intervals']==124 and r['rejected_intervals']==72
    assert r['status']=='matching_inconsistent' and r['model_ready'] is False

def test_table_inspection_preserves_existing_output(tmp_path):
    p=tmp_path/'keep';p.mkdir();(p/'mine').write_text('mine')
    with pytest.raises(FileExistsError):api().inspect_table(ROOT/'fixtures/tables/hermite_consistent.csv',p)
    assert (p/'mine').read_text()=='mine'
