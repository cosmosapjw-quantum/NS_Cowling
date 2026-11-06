"""Export keeps supplied methods and explicitly selected physical functions."""
from pathlib import Path
import importlib
import pytest
ROOT=Path(__file__).resolve().parents[1]

def module():
    assert (ROOT/'tools/make_calculation.py').is_file(), 'No explicit physics-to-methods connection exists'
    return importlib.import_module('tools.make_calculation')

def test_existing_output_is_not_overwritten(tmp_path):
    m=module();p=tmp_path/'keep';p.mkdir();(p/'mine').write_text('mine')
    with pytest.raises(FileExistsError):m.create(p)
    assert (p/'mine').read_text()=='mine'

def test_missing_physics_is_not_completed_from_examples(tmp_path):
    m=module();p=tmp_path/'physical';p.mkdir()
    with pytest.raises(FileNotFoundError):m.create(tmp_path/'new',p)
    assert not (tmp_path/'new').exists()

def test_student_physics_is_copied_exactly(tmp_path):
    m=module();src=tmp_path/'physics';src.mkdir()
    for name in m.PHYSICS:
        source=next((ROOT/'learner/src').rglob(name))
        (src/name).write_bytes(source.read_bytes()+b'\n! MY PHYSICS MARKER\n')
    out=tmp_path/'my calculation';m.create(out,src)
    assert not (out/'solutions').exists()
    for name in m.PHYSICS:assert (out/'calculation/src/physics'/name).read_bytes()==(src/name).read_bytes()
    assert (out/'methods/fpm.toml').is_file()
