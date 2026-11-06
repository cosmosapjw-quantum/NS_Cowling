from pathlib import Path
import os
import shutil
import pytest
from tests.support.runner import CaseRunner, load_contract

ROOT = Path(__file__).resolve().parents[1]

@pytest.fixture
def root():
    return ROOT

@pytest.fixture
def copy_root(tmp_path):
    out = tmp_path / 'source copy with spaces'
    shutil.copytree(ROOT, out, ignore=shutil.ignore_patterns('.git', '.runtime', 'build', '__pycache__', '.pytest_cache', 'evidence'))
    return out

@pytest.fixture
def runner(root, tmp_path):
    return CaseRunner(root, tmp_path / 'runs', workspace=os.environ.get('NS_WORKSPACE', 'solutions'))

@pytest.fixture
def run_case(runner):
    return runner.run_case

@pytest.fixture
def contract(root):
    return lambda fixture_id: load_contract(root, fixture_id)
